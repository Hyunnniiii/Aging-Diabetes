## =========================================
## MAST DEG RDS → GSEA (Reactome / GO:BP 전용)
##  - 모든 플롯: PDF 저장
##  - padj<0.25가 0개면 fgsea 재시도(강화 파라미터) + FC-랭크 대안도 저장
##  - Diverging barplot(빨강=Down, 파랑=Up) 저장
##  - 결과는 한 폴더(out_dir)에 모음
## =========================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(clusterProfiler); library(enrichplot); library(msigdbr)
  library(ggplot2)
})

## === 설정 ===
deg_dir <- "/data/project/diabetes_LYH/tanya/rds/MAST"   # *_MAST.rds 위치
out_dir <- "GSEA_results_Reactome_GO"                    # 결과 폴더
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

remove_mt <- TRUE   # MT-/MTRNR/mt-/mt-Rnr 유전자 제거 여부

## === 헬퍼들 ===
infer_contrast <- function(path){
  x <- toupper(basename(path))
  m <- stringr::str_match(x, "(HY|HO|DO)[-_ ]?VS[-_ ]?(HY|HO|DO)")
  if (!is.na(m[1,1])) paste0(m[1,2], "_vs_", m[1,3]) else tools::file_path_sans_ext(basename(path))
}

# ENSG → SYMBOL (데이터가 이미 심볼이면 스킵)
needs_symbol_map <- function(genes) mean(startsWith(genes, "ENSG")) > 0.8
map_ensg_to_symbol <- function(genes){
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Ensembl ID 감지됨. BiocManager::install('org.Hs.eg.db') 후 재시작하세요.")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ids <- sub("\\..*$","", genes)
  m <- AnnotationDbi::select(org.Hs.eg.db, keys=ids, keytype="ENSEMBL", columns="SYMBOL")
  dplyr::tibble(ENSEMBL=ids) |>
    dplyr::left_join(dplyr::distinct(tibble::as_tibble(m), ENSEMBL, .keep_all=TRUE), by="ENSEMBL") |>
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL)) |>
    dplyr::pull(SYMBOL)
}

# MT 유전자 제거 (패턴: MT-|MTRNR|mt-|mt-Rnr, 대소문자 무시)
drop_mt_genes <- function(df, gene_col = "gene", verbose = TRUE) {
  stopifnot(gene_col %in% names(df))
  patt <- "^(MT-|MTRNR|mt-|mt-Rnr)"
  keep <- !grepl(patt, df[[gene_col]], ignore.case = TRUE)
  if (verbose) message("  [filter] MT genes removed: ", sum(!keep))
  df[keep, , drop = FALSE]
}

# Diverging barplot (빨강=Down, 파랑=Up)
make_barplot_diverging <- function(resdf, contrast, collection,
                                   topn = 20, padj_cut = 0.25, wrap = 45) {
  stopifnot(all(c("NES","p.adjust") %in% names(resdf)))
  name_col <- if ("Description" %in% names(resdf)) "Description" else if ("ID" %in% names(resdf)) "ID" else "ID"

  base <- dplyr::filter(resdf, !is.na(p.adjust))
  sig  <- dplyr::filter(base, p.adjust < padj_cut)
  pool <- if (nrow(sig) > 0) sig else base

  n_each <- max(1, floor(topn/2))
  up   <- pool %>% dplyr::filter(NES >= 0) %>% dplyr::arrange(p.adjust, dplyr::desc(NES)) %>% dplyr::slice_head(n = n_each)
  down <- pool %>% dplyr::filter(NES < 0)  %>% dplyr::arrange(p.adjust, NES)             %>% dplyr::slice_head(n = n_each)
  sel  <- dplyr::bind_rows(down, up)
  if (nrow(sel) == 0) return(NULL)

  lab <- sel[[name_col]]
  lab <- gsub("^(GOBP_|GO_BP_|REACTOME_|R-HSA-)", "", lab)
  lab <- gsub("_", " ", lab)
  lab <- stringr::str_to_sentence(lab)
  lab <- stringr::str_wrap(lab, width = wrap)

  sel$label <- lab
  sel$dir   <- ifelse(sel$NES >= 0, "Up", "Down")
  lim <- max(abs(sel$NES), na.rm = TRUE)

  plt <- ggplot2::ggplot(sel, ggplot2::aes(x = reorder(label, NES), y = NES, fill = dir)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::scale_y_continuous(limits = c(-lim, lim)) +
    ggplot2::scale_fill_manual(values = c("Down" = "#d73027", "Up" = "#4575b4")) +
    ggplot2::labs(
      x = NULL, y = "NES",
      title = paste0(contrast, " · ", collection),
      subtitle = paste0("Top ", nrow(sel), " (padj<", padj_cut, " 우선) · 좌우 대칭 스케일")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  list(plot = plt, n = nrow(sel))
}

# 기본 랭크: sign(FC) * -log10(p_raw or adj)
build_ranks <- function(mast_rds){
  message("  [rank] loading: ", mast_rds)
  df <- readRDS(mast_rds)
  if (!"gene" %in% names(df)) df <- tibble::rownames_to_column(df, "gene")
  stopifnot(any(c("avg_log2FC","avg_logFC") %in% names(df)))
  fc_col <- if ("avg_log2FC" %in% names(df)) "avg_log2FC" else "avg_logFC"
  p_col  <- if ("p_val"     %in% names(df)) "p_val" else if ("p_val_adj" %in% names(df)) "p_val_adj" else NA
  stopifnot(!is.na(p_col))

  if (needs_symbol_map(df$gene)) df$gene <- map_ensg_to_symbol(df$gene)
  if (isTRUE(remove_mt)) df <- drop_mt_genes(df, gene_col = "gene", verbose = TRUE)

  df <- df |>
    dplyr::mutate(p_for_rank = pmax(.data[[p_col]], .Machine$double.xmin),
                  score      = sign(.data[[fc_col]]) * -log10(p_for_rank)) |>
    dplyr::arrange(dplyr::desc(abs(score))) |>
    dplyr::distinct(gene, .keep_all=TRUE) |>
    dplyr::filter(!is.na(gene), is.finite(score))

  message("  [rank] genes in rank = ", nrow(df))
  v <- df$score; names(v) <- df$gene
  sort(v, decreasing = TRUE)
}

# 대안 랭크: FC만 (p가 평평할 때 유용)
build_fc <- function(mast_rds){
  df <- readRDS(mast_rds)
  if (!"gene" %in% names(df)) df <- tibble::rownames_to_column(df, "gene")
  if (needs_symbol_map(df$gene)) df$gene <- map_ensg_to_symbol(df$gene)
  if (isTRUE(remove_mt)) df <- drop_mt_genes(df, gene_col = "gene", verbose = TRUE)
  fc_col <- if ("avg_log2FC" %in% names(df)) "avg_log2FC" else "avg_logFC"
  v <- df[[fc_col]]; names(v) <- df$gene
  sort(v, decreasing = TRUE)
}

# 최신 msigdbr 호출(경고 제거)
get_msig <- function(collection){
  if (collection=="C2:REACTOME") return(msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"))
  if (collection=="C5:BP")       return(msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"))
  stop("unknown collection: ", collection)
}

# enrich 커브 폴백
make_fallback_enrich <- function(term, t2g, geneList, title, outfile){
  gs <- t2g |> dplyr::filter(gs_name == term) |> pull(gene_symbol)
  plt <- fgsea::plotEnrichment(gs, sort(geneList, decreasing=TRUE)) + ggtitle(title)
  ggsave(outfile, plt, width=7, height=5, device="pdf")
}

## === 메인: 한 비교 처리 (Reactome/GO:BP만) ===
gsea_and_plots_one <- function(mast_rds,
                               collections = c("C2:REACTOME","C5:BP"),
                               showCategory = 20){
  contrast <- infer_contrast(mast_rds)
  message("\n========== [", contrast, "] ==========")

  geneList   <- build_ranks(mast_rds)
  geneListFC <- build_fc(mast_rds)  # 대안 랭크
  readr::write_tsv(tibble(gene=names(geneList), score=as.numeric(geneList)),
                   file.path(out_dir, paste0(contrast, "__rank.tsv")))

  for (collection in collections) {
    message("  [GSEA] collection = ", collection)
    pref <- file.path(out_dir, paste0(contrast, "__", gsub(":","",collection), "__"))
    msig <- get_msig(collection)
    t2g  <- dplyr::select(msig, gs_name, gene_symbol)

    # clusterProfiler::GSEA (fgsea 엔진)
    gseaRes <- tryCatch({
      GSEA(geneList = geneList, TERM2GENE = t2g,
           minGSSize = 10, maxGSSize = 2000,
           pvalueCutoff = 1, pAdjustMethod = "BH",
           seed = TRUE, by = "fgsea", verbose = FALSE)
    }, error=function(e){ message("  [ERR] GSEA: ", e$message); NULL })
    if (is.null(gseaRes)) next

    resdf <- as.data.frame(gseaRes)
    sig_n <- sum(resdf$p.adjust < 0.25, na.rm=TRUE)
    message(sprintf("  [GSEA] terms=%d (padj<0.25: %d)", nrow(resdf), sig_n))
    readr::write_tsv(dplyr::arrange(resdf, p.adjust, dplyr::desc(abs(NES))),
                     paste0(pref, "GSEA_result.tsv"))

    # padj<0.25가 0개면: fgsea(강화) + FC랭크도 같이 저장
    if (sig_n == 0) {
      message("  [RERUN] fgsea relaxed + FC-rank")
      pw <- split(msig$gene_symbol, msig$gs_name)
      fg  <- fgsea::fgsea(pathways=pw, stats=geneList,   minSize=10, maxSize=2000, nperm=100000)
      fg  <- fg[order(fg$padj, -abs(fg$NES)), ]
      fgfc<- fgsea::fgsea(pathways=pw, stats=geneListFC, minSize=10, maxSize=2000, nperm=100000)
      fgfc<- fgfc[order(fgfc$padj, -abs(fgfc$NES)), ]
      readr::write_tsv(fg,   paste0(pref, "FGSEA_relaxed.tsv"))
      readr::write_tsv(fgfc, paste0(pref, "FGSEA_relaxed_FCrank.tsv"))
    }

    # 플롯: enrichment(Top+/Top−), diverging barplot, dotplot(Up/Down), cnetplot
    if (nrow(resdf) > 0) {
      gseaRes@result$.sign <- ifelse(gseaRes@result$NES >= 0, "Up", "Down")
      top_up   <- resdf |> dplyr::filter(NES>0) |> arrange(p.adjust, dplyr::desc(NES)) |> dplyr::slice_head(n=1) |> pull(ID)
      top_down <- resdf |> dplyr::filter(NES<0) |> arrange(p.adjust, NES)             |> dplyr::slice_head(n=1) |> pull(ID)

      if (length(top_up)) {
        ok <- tryCatch({
          p1 <- gseaplot2(gseaRes, geneSetID = top_up, pvalue_table = FALSE,
                          title = paste0(collection, " · Top+ : ", top_up))
          ggsave(paste0(pref, "TopPos_enrichment.pdf"), p1, width=7, height=5, device="pdf"); TRUE
        }, error=function(e){ message("  [warn] gseaplot2(+): ", e$message); FALSE })
        if (!ok) make_fallback_enrich(top_up, t2g, geneList,
                                      paste0(collection, " · Top+ : ", top_up),
                                      paste0(pref, "TopPos_enrichment.pdf"))
      }
      if (length(top_down)) {
        ok <- tryCatch({
          p2 <- gseaplot2(gseaRes, geneSetID = top_down, pvalue_table = FALSE,
                          title = paste0(collection, " · Top- : ", top_down))
          ggsave(paste0(pref, "TopNeg_enrichment.pdf"), p2, width=7, height=5, device="pdf"); TRUE
        }, error=function(e){ message("  [warn] gseaplot2(-): ", e$message); FALSE })
        if (!ok) make_fallback_enrich(top_down, t2g, geneList,
                                      paste0(collection, " · Top- : ", top_down),
                                      paste0(pref, "TopNeg_enrichment.pdf"))
      }

      # (새) Diverging barplot
      db <- tryCatch(
        make_barplot_diverging(resdf, contrast, collection, topn = showCategory, padj_cut = 0.25, wrap = 40),
        error = function(e){ message("  [warn] diverging barplot: ", e$message); NULL }
      )
      if (!is.null(db)) {
        ggsave(paste0(pref, "barplot_diverging.pdf"),
               db$plot, width = 8.5, height = max(4, 0.35 * db$n), device = "pdf")
      }

      # Up/Down dotplot (데이터 없으면 자동 스킵)
      for (sg in c("Up","Down")) {
        dp <- tryCatch({
            sub <- gseaRes; sub@result <- subset(gseaRes@result, .sign == sg)
            if (nrow(sub@result) == 0) NULL else
              enrichplot::dotplot(sub, showCategory = min(showCategory, nrow(sub@result)))
          },
          error=function(e){ message("  [warn] dotplot(", sg, "): ", e$message); NULL }
        )
        if (!is.null(dp)) {
          ggsave(paste0(pref, "dotplot_", sg, ".pdf"),
                 dp, width = 9, height = 7, device = "pdf")
        }
      }

      # cnetplot (항목 있을 때만)
      cn <- NULL
      try(cn <- cnetplot(gseaRes, showCategory = min(10, nrow(resdf)),
                         foldChange = geneList, circular = FALSE), silent=TRUE)
      if (!is.null(cn)) ggsave(paste0(pref, "cnetplot.pdf"), cn, width=9, height=7, device="pdf")
    }
    message("  [save] ", pref, "*.pdf / *.tsv")
  }
}

## === 실행 루프 (Reactome/GO:BP만) ===
rds_files <- list.files(deg_dir, pattern = "_MAST\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("RDS가 없습니다: ", deg_dir)
message("대상 파일: ", length(rds_files), "개 → ", paste(basename(rds_files), collapse=", "))

for (f in rds_files) {
  tryCatch(
    gsea_and_plots_one(f, collections = c("C2:REACTOME","C5:BP"), showCategory = 20),
    error = function(e) message("\n[SKIP] ", basename(f), " : ", e$message)
  )
}

message("\n=== 완료 ===\n출력 폴더: ", normalizePath(out_dir))
