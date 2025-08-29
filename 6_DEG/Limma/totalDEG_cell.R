# ============================================================
## Date | 2025-08-29 | Yoon JiHyun
## Anchor 기반 integration 결과 파일로 이후 분석 진행 (Harmony 추가로 처리한 파일은 다루지 않음)

## Description | 모든 샘플에서 celltype별로 환자×배치 pseudobulk 생성
## edgeR(TMM) - limma-voom으로 HY-HO-DO 대비 DEG를 계산, celltype별 logFC 분포를 plot으로 저장
# ============================================================


suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr); library(ggplot2)
})

# ==== User 설정 ====
rds_path <- "/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds"
fc_cut  <- 0.58
padj_cut <- 0.05

# ==== Load ====
seu <- readRDS(rds_path)
counts <- tryCatch(GetAssayData(seu, assay="RNA", layer="counts"),
                   error=function(e) GetAssayData(seu, assay="RNA", slot="counts"))
meta <- seu@meta.data
stopifnot(ncol(counts) == nrow(meta), all(colnames(counts) == rownames(meta)))

# ==== Common factors ====
batch_levels <- c("Run1","Run2_1","Run2_2")
meta$batch     <- factor(meta$orig.ident, levels = batch_levels)
meta$group     <- factor(meta$group, levels = c("Diabetes_Old","Healthy_Old","Healthy_Young"))
meta$sample_id <- interaction(meta$patient, meta$batch, drop = TRUE)

# ==== Cell type 컬럼 직접 지정 ====
ct_col <- "celltype"  # 메타데이터에서 실제 cell type 컬럼명

# ==== 결과 저장용 ====
all_res <- list()

# ==== cell type loop ====
for (ct in sort(unique(meta[[ct_col]]))) {
  sub_cells <- rownames(meta)[meta[[ct_col]] == ct]
  if (length(sub_cells) < 50) next

  sub_meta <- meta[sub_cells, , drop=FALSE]
  sub_counts <- counts[, sub_cells, drop=FALSE]

  # pseudobulk
  idx <- split(seq_len(ncol(sub_counts)), sub_meta$sample_id, drop=TRUE)
  if (length(idx) < 6) next
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
  pb_counts <- pb_counts[, keep_ids, drop=FALSE]
  rownames(samples_info) <- samples_info$sample_id

  # edgeR → voom → limma
  y <- DGEList(counts=pb_counts, samples=samples_info)
  y$samples$group   <- factor(y$samples$group,   levels=levels(meta$group))
  y$samples$batch   <- factor(y$samples$batch,   levels=levels(meta$batch))
  y$samples$patient <- factor(y$samples$patient)

  keep <- filterByExpr(y, group=y$samples$group, min.count=10)
  y <- y[keep, , keep.lib.sizes=FALSE]
  if (nrow(y) < 50) next
  y <- calcNormFactors(y)

  design <- model.matrix(~ 0 + group + batch, data=y$samples)
  v <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)

  have <- colnames(design)
  contr <- list()
  if (all(c("groupHealthy_Old","groupDiabetes_Old") %in% have)) 
    contr$HO_vs_DO <- "groupHealthy_Old - groupDiabetes_Old"
  if (all(c("groupHealthy_Young","groupDiabetes_Old") %in% have)) 
    contr$HY_vs_DO <- "groupHealthy_Young - groupDiabetes_Old"
  if (all(c("groupHealthy_Old","groupHealthy_Young") %in% have)) 
    contr$HO_vs_HY <- "groupHealthy_Old - groupHealthy_Young"
  if (!length(contr)) next

  cont <- makeContrasts(contrasts=unlist(contr), levels=colnames(coef(fit)))
  fit2 <- eBayes(contrasts.fit(fit, cont))

  for (nm in colnames(cont)) {
    tt <- topTable(fit2, coef=nm, n=Inf, sort.by="P")
    tt$Gene <- rownames(tt)
    tt$cell_type <- ct
    tt$contrast <- nm
    all_res[[length(all_res)+1]] <- tt
  }
}

# ==== 결과 합치기 ====
df_all <- bind_rows(all_res)

# ==== 색상 지정 ====
df_all$color <- "grey80"
df_all$color[df_all$logFC >  fc_cut & df_all$adj.P.Val < padj_cut] <- "red"
df_all$color[df_all$logFC < -fc_cut & df_all$adj.P.Val < padj_cut] <- "blue"


# ==== PDF에 3개 대비 모두 출력 ====
cairo_pdf("celltypeDEG_allContrasts.pdf", width=10, height=6, onefile=TRUE)

for (nm in unique(df_all$contrast)) {
  df_sub <- df_all %>% dplyr::filter(contrast == nm)

  # (옵션) 라벨 달 유의 DEG만 현재 대비로 필터
  label_df_sub <- df_sub %>% dplyr::filter(color %in% c("red","blue"))

  # 라벨 없음 버전
  p <- ggplot(df_sub, aes(x=cell_type, y=logFC, color=color)) +
    geom_jitter(width=0.3, alpha=0.6, size=1) +
    scale_color_identity() +
    theme_minimal(base_size=12) +
    labs(title=paste("DEG across cell types:", nm),
         x="Cell type", y="log2 Fold Change") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  print(p)

  # 라벨 있음 버전 (원하면만 사용)
  q <- ggplot(df_sub, aes(x=cell_type, y=logFC, color=color)) +
    geom_jitter(width=0.3, alpha=0.6, size=1) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      data = label_df_sub, aes(label = Gene),
      size = 2, max.overlaps = Inf, min.segment.length = 0
    ) +
    theme_minimal(base_size=12) +
    labs(title=paste("DEG across cell types (labeled):", nm),
         x="Cell type", y="log2 Fold Change") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  print(q)
}

dev.off()