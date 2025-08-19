## ================= HY_vs_HO 전용 GSEA 재런(파워↑) =================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(clusterProfiler); library(msigdbr); library(enrichplot); library(ggplot2)
})

## ---- 설정 ----
mast_rds   <- "/data/project/diabetes_LYH/tanya/rds/MAST/HY_vs_HO_MAST.rds"  # 필요시 경로 확인
out_root   <- "GSEA_results_HYvsHO_boost"
dir.create(out_root, FALSE, TRUE)

remove_mt   <- TRUE
remove_ribo <- FALSE   # RPL/RPS/MRPL/MRPS 제거
remove_hsp  <- FALSE   # HSP/HSPA/HSPB 등 제거
minGS       <- 10
maxGS       <- 1000   # 너무 큰 세트는 제외
padj_show   <- 0.25   # HYvsHO는 GSEA 관행(탐색)대로 0.25도 보고

## ---- 유틸 ----
needs_symbol_map <- function(genes) mean(startsWith(genes, "ENSG")) > 0.8
map_ensg_to_symbol <- function(genes){
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db 설치 필요: BiocManager::install('org.Hs.eg.db')")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ids <- sub("\\..*$","", genes)
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype="ENSEMBL", columns="SYMBOL")
  tibble(ENSEMBL=ids) |>
    left_join(distinct(as_tibble(m), ENSEMBL, .keep_all=TRUE), by="ENSEMBL") |>
    mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL)) |>
    pull(SYMBOL)
}

drop_genes <- function(df, gene_col="gene"){
  patts <- c()
  if (isTRUE(remove_mt))   patts <- c(patts, "^(MT-|MTRNR|mt-|mt-Rnr)")
  if (isTRUE(remove_ribo)) patts <- c(patts, "^(RPL|RPS|MRPL|MRPS)\\b")
  if (isTRUE(remove_hsp))  patts <- c(patts, "^HSP")
  if (length(patts)==0) return(df)
  keep <- !grepl(paste(patts, collapse="|"), df[[gene_col]], ignore.case=TRUE)
  message("  [filter] removed genes: ", sum(!keep))
  df[keep, , drop=FALSE]
}

get_msig <- function(collection){
  switch(collection,
    "H"            = msigdbr(species="Homo sapiens", collection="H"),
    "C2:REACTOME"  = msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"),
    "C5:BP"        = msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"),
    stop("unknown collection")
  )
}

make_rank <- function(df, method=c("fc_x_p","sign_p","fc_only")){
  method <- match.arg(method)
  fc_col <- if ("avg_log2FC" %in% names(df)) "avg_log2FC" else "avg_logFC"
  p_col  <- if ("p_val"     %in% names(df)) "p_val"     else "p_val_adj"
  df <- mutate(df, puse=pmax(.data[[p_col]], .Machine$double.xmin))
  score <- switch(method,
    "fc_x_p" = df[[fc_col]] * -log10(df$puse),
    "sign_p" = sign(df[[fc_col]]) * -log10(df$puse),
    "fc_only"= df[[fc_col]]
  )
  v <- score; names(v) <- df$gene
  sort(v, decreasing=TRUE)
}

run_one <- function(collection, method){
  message(sprintf("\n== [%s | %s] ==", collection, method))
  # 데이터 로딩/정리
  df <- readRDS(mast_rds)
  if (!"gene" %in% names(df)) df <- tibble::rownames_to_column(df, "gene")
  # ID 매핑
  if (needs_symbol_map(df$gene)) df$gene <- map_ensg_to_symbol(df$gene)
  # 필터링
  df <- drop_genes(df, "gene") |> distinct(gene, .keep_all=TRUE) |> filter(!is.na(gene))
  # 랭킹
  geneList <- make_rank(df, method=method)
  message("  [rank] genes=", length(geneList), " | top head=", paste(head(names(geneList),3), collapse=", "))
  # gene sets
  msig <- get_msig(collection)
  t2g  <- select(msig, gs_name, gene_symbol)
  # GSEA (fgsea 엔진; multilevel 사용)
  gsea <- GSEA(geneList=geneList, TERM2GENE=t2g,
               minGSSize=minGS, maxGSSize=maxGS,
               pvalueCutoff=1, pAdjustMethod="BH",
               seed=TRUE, by="fgsea", verbose=FALSE)
  res  <- as.data.frame(gsea)
  sig05 <- sum(res$p.adjust < 0.05, na.rm=TRUE)
  sig25 <- sum(res$p.adjust < 0.25, na.rm=TRUE)
  message(sprintf("  [GSEA] terms=%d | FDR<0.05: %d | FDR<0.25: %d | minFDR=%.3g",
                  nrow(res), sig05, sig25, min(res$p.adjust, na.rm=TRUE)))
  pref <- file.path(out_root, sprintf("HY_vs_HO__%s__%s__", gsub(":","",collection), method))
  write_tsv(arrange(res, p.adjust, desc(abs(NES))), paste0(pref, "GSEA_result.tsv"))

  # 탐색용 diverging barplot (padj<0.25 우선)
  if (nrow(res)>0){
    base <- filter(res, !is.na(p.adjust))
    pool <- if (any(base$p.adjust < padj_show)) filter(base, p.adjust < padj_show) else base
    n_each <- max(1, floor(20/2))
    up   <- pool |> filter(NES>=0) |> arrange(p.adjust, desc(NES)) |> slice_head(n=n_each)
    down <- pool |> filter(NES<0)  |> arrange(p.adjust, NES)       |> slice_head(n=n_each)
    sel  <- bind_rows(down, up); if (nrow(sel)>0){
      sel$label <- sel$Description
      sel$label <- gsub("^(GOBP_|GO:BP_|GO_BP_|REACTOME_|R-HSA-\\d+_|R-HSA-)","", sel$label)
      sel$label <- stringr::str_to_sentence(gsub("_"," ", sel$label))
      sel$dir   <- ifelse(sel$NES>=0, "Up","Down")
      lim <- max(abs(sel$NES), na.rm=TRUE)
      p <- ggplot(sel, aes(x=reorder(label, NES), y=NES, fill=dir)) +
        geom_col(width=0.8) + coord_flip(clip="off") +
        geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
        scale_y_continuous(limits=c(-lim, lim), expand=expansion(mult=c(0.05,0.08))) +
        scale_fill_manual(values=c(Down="#4575b4", Up="#d73027")) +
        labs(x=NULL,y="NES",
             title=sprintf("HY_vs_HO · %s", ifelse(collection=="H","Hallmark",collection)),
             subtitle=sprintf("Top 20 (%s)", if (any(base$p.adjust<padj_show)) paste0("padj<",padj_show) else "best available")) +
        theme_bw(base_size=11) +
        theme(legend.position="top", axis.text.y=element_text(size=9, lineheight=0.95),
              panel.grid.major.y=element_blank(), plot.margin=margin(5.5,14,5.5,5.5))
      ggsave(paste0(pref, "barplot_diverging.pdf"), p, width=8, height=max(6, 0.45*nrow(sel)))
    }
  }
  invisible(res)
}

## ---- 실행 순서 ----
# 1) Hallmark 먼저(검정 개수↓) → 세 가지 랭킹 시도
res_H_fcxp  <- run_one("H", "fc_x_p")
res_H_signp <- run_one("H", "sign_p")
res_H_fconly<- run_one("H", "fc_only")

# 2) Reactome / GO:BP도 동일하게 시도(원한다면 주석 해제)
res_R_fcxp  <- run_one("C2:REACTOME", "fc_x_p")
res_G_fcxp  <- run_one("C5:BP",       "fc_x_p")

message("\n== done ==\nOutput: ", normalizePath(out_root))
