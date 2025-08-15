suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr); library(ggplot2)
})

# ==== 설정 ====
rds_path <- "/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds"
ct_col   <- "celltype"      # 메타데이터의 셀타입 컬럼명
fc_cut   <- 0.58
padj_cut <- 0.05
min_cells_per_ct <- 50       # 셀타입 최소 세포 수
min_pb_samples   <- 6        # pseudobulk 최소 샘플 수
out_pdf  <- "celltypeDEG_allContrasts_Run1.pdf"

# ==== Load ====
seu <- readRDS(rds_path)
counts <- tryCatch(GetAssayData(seu, assay="RNA", layer="counts"),
                   error=function(e) GetAssayData(seu, assay="RNA", slot="counts"))
meta <- seu@meta.data
stopifnot(ncol(counts) == nrow(meta), all(colnames(counts) == rownames(meta)))

# ==== Run1만 필터 ====
meta_f <- meta %>% filter(orig.ident == "Run1")
counts_f <- counts[, rownames(meta_f)]
meta_f$batch     <- factor(meta_f$orig.ident, levels = c("Run1","Run2_1","Run2_2"))
meta_f$group     <- factor(meta_f$group, levels = c("Diabetes_Old","Healthy_Old","Healthy_Young"))
meta_f$sample_id <- interaction(meta_f$patient, meta_f$batch, drop = TRUE)

# ==== 결과 누적 ====
all_res <- list()

# ==== cell type loop ====
for (ct in sort(unique(meta_f[[ct_col]]))) {
  sub_cells <- rownames(meta_f)[meta_f[[ct_col]] == ct]
  if (length(sub_cells) < min_cells_per_ct) next

  sub_meta   <- meta_f[sub_cells, , drop=FALSE]
  sub_counts <- counts_f[, sub_cells, drop=FALSE]

  # pseudobulk (patient×batch)
  idx <- split(seq_len(ncol(sub_counts)), sub_meta$sample_id, drop=TRUE)
  if (length(idx) < min_pb_samples) next
  pb_counts <- do.call(cbind, lapply(idx, function(i) Matrix::rowSums(sub_counts[, i, drop=FALSE])))
  colnames(pb_counts) <- names(idx)
  pb_counts <- pb_counts[!grepl("^MT-", rownames(pb_counts), ignore.case=TRUE), , drop=FALSE]

  samples_info <- sub_meta %>%
    group_by(sample_id) %>%
    summarise(patient=first(patient), group=first(group), batch=first(batch), .groups="drop") %>%
    distinct(sample_id, .keep_all=TRUE) %>%
    as.data.frame()

  keep_ids <- intersect(samples_info$sample_id, colnames(pb_counts))
  if (length(keep_ids) < 2) next
  samples_info <- samples_info[match(keep_ids, samples_info$sample_id), , drop=FALSE]
  pb_counts    <- pb_counts[, keep_ids, drop=FALSE]
  rownames(samples_info) <- samples_info$sample_id

  # edgeR → voom → limma (Run1은 배치 한 레벨 → ~0+group)
  y <- DGEList(counts=pb_counts, samples=samples_info)
  y$samples$group   <- droplevels(y$samples$group)
  keep <- filterByExpr(y, group=y$samples$group, min.count=10)
  y <- y[keep, , keep.lib.sizes=FALSE]
  if (nrow(y) < 50) next
  y <- calcNormFactors(y)

  design <- model.matrix(~ 0 + group, data=y$samples)
  v <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)

  # 대비 생성
  have <- colnames(design)
  contr <- list()
  if (all(c("groupHealthy_Old","groupDiabetes_Old") %in% have)) 
    contr$HO_vs_DO <- "groupHealthy_Old - groupDiabetes_Old"
  if (all(c("groupHealthy_Young","groupDiabetes_Old") %in% have)) 
    contr$HY_vs_DO <- "groupHealthy_Young - groupDiabetes_Old"
  if (all(c("groupHealthy_Old","groupHealthy_Young") %in% have)) 
    contr$HO_vs_HY <- "groupHealthy_Old - groupHealthy_Young"
  if (!length(contr)) next

  cont <- makeContrasts(contrasts = unlist(contr), levels = colnames(coef(fit)))
  fit2 <- eBayes(contrasts.fit(fit, cont))

  for (nm in colnames(cont)) {
    tt <- topTable(fit2, coef=nm, n=Inf, sort.by="P")
    tt$Gene <- rownames(tt)
    tt$cell_type <- ct
    tt$contrast  <- nm
    all_res[[length(all_res)+1]] <- tt
  }
}

# ==== 합치고 색 지정 ====
df_all <- bind_rows(all_res)
df_all$color <- "grey80"
df_all$color[df_all$logFC >  fc_cut & df_all$adj.P.Val < padj_cut] <- "red"
df_all$color[df_all$logFC < -fc_cut & df_all$adj.P.Val < padj_cut] <- "blue"

# ==== 라벨 없음 버전 ====
cairo_pdf(out_pdf, width = 10, height = 6, onefile = TRUE)
for (nm in unique(df_all$contrast)) {
  df_sub <- df_all %>% filter(contrast == nm)

    # (옵션) 라벨 달 유의 DEG만 현재 대비로 필터
  label_df_sub <- df_sub %>% dplyr::filter(color %in% c("red","blue"))

  p <- ggplot(df_sub, aes(x=cell_type, y=logFC, color=color)) +
    geom_jitter(width=0.3, alpha=0.6, size=1) +
    scale_color_identity() +
    theme_minimal(base_size=12) +
    labs(title=paste0("[Run1] DEG across cell types: ", nm),
         x="Cell type", y="log2 Fold Change") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  print(p)

  # 라벨 있음 버전
  q <- ggplot(df_sub, aes(x=cell_type, y=logFC, color=color)) +
    geom_jitter(width=0.3, alpha=0.6, size=1) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      data = label_df_sub, aes(label = Gene),
      size = 2, max.overlaps = Inf, min.segment.length = 0
    ) +
    theme_minimal(base_size=12) +
    labs(title=paste("[Run1] DEG across cell types (labeled):", nm),
         x="Cell type", y="log2 Fold Change") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  print(q)
}
dev.off()
message("Saved: ", out_pdf)
