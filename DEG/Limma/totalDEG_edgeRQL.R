# ============================================================
# Pseudobulk (patient × batch) 공통 전처리  → 두 파이프라인 최소 출력
#   ① limma-voom: voom_mean_variance.pdf + Volcano_limma.pdf
#   ② edgeR QL GLM: BCV.pdf + Volcano_edger.pdf
#   DEG 테이블: DEG_HO_vs_DO.tsv, DEG_HY_vs_DO.tsv, DEG_HO_vs_HY.tsv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr)
  library(ggplot2); library(ggrepel)
})

# --------------------------
# 0) Paths & Output Folders
# --------------------------
rds_path <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT_with_markers.rds"
out_limma <- "deg_out_limma_voom_patientXbatch"
out_edger <- "deg_out_edger_ql_patientXbatch"
dir.create(out_limma, showWarnings = FALSE)
dir.create(out_edger, showWarnings = FALSE)

# --------------------------
# 1) Load Seurat & counts
# --------------------------
seu <- readRDS(rds_path)
counts <- tryCatch(
  GetAssayData(seu, assay = "RNA", layer = "counts"),   # Seurat v5
  error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts") # v4 fallback
)
stopifnot(inherits(counts, "dgCMatrix") || is.matrix(counts))
meta <- seu@meta.data
stopifnot(ncol(counts) == nrow(meta))

# ------------------------------------------
# 2) Define batch, group, sample_id (공통)
# ------------------------------------------
batch_levels <- c("Run1","Run2_1","Run2_2")
meta$batch  <- factor(meta$orig.ident, levels = batch_levels)
meta$group  <- factor(meta$group, levels = c("Diabetes_Old","Healthy_Old","Healthy_Young"))
meta$sample_id <- interaction(meta$patient, meta$batch, drop = TRUE)
stopifnot(all(colnames(counts) == rownames(meta)))

# ------------------------------------
# 3) Pseudobulk aggregation (공통)
# ------------------------------------
idx <- split(seq_len(ncol(counts)), meta$sample_id, drop = TRUE)
pb_counts <- do.call("cbind", lapply(idx, function(i) Matrix::rowSums(counts[, i, drop = FALSE])))

## ★ 추가: 미토콘드리아 유전자 제거 (대/소문자 무시)
rm_mt <- grepl("^MT-", rownames(pb_counts), ignore.case = TRUE)
pb_counts <- pb_counts[!rm_mt, , drop = FALSE]
colnames(pb_counts) <- names(idx)

# -----------------------------------------
# 4) samples_info for pseudobulk (공통)
# -----------------------------------------
samples_info <- meta %>%
  group_by(sample_id) %>%
  summarise(patient = first(patient),
            group   = first(group),
            batch   = first(batch), .groups = "drop") %>%
  as.data.frame()
pb_counts <- pb_counts[, samples_info$sample_id]
rownames(samples_info) <- samples_info$sample_id

# -----------------------------------------
# 5) edgeR container + filtering + TMM (공통)
# -----------------------------------------
y <- DGEList(counts = pb_counts, samples = samples_info)
y$samples$group   <- factor(y$samples$group, levels = levels(meta$group))
y$samples$batch   <- factor(y$samples$batch, levels = batch_levels)
y$samples$patient <- factor(y$samples$patient)

keep <- filterByExpr(y, group = y$samples$group, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# ------------------------------
# 6) Design (공통): ~ 0 + group + batch
# ------------------------------
design <- model.matrix(~ 0 + group + batch, data = y$samples)

# ------------------------------
# Helper: Volcano (DEG만 라벨)
# ------------------------------
plot_volcano <- function(tt, title, fc_cut = 0.58, padj_col = "padj") {
  # padj 통일 보장
  if (!(padj_col %in% names(tt))) {
    if ("adj.P.Val" %in% names(tt)) tt[[padj_col]] <- tt$adj.P.Val
    else if ("FDR" %in% names(tt))  tt[[padj_col]] <- tt$FDR
    else if ("P.Value" %in% names(tt)) tt[[padj_col]] <- p.adjust(tt$P.Value, "BH")
  }
  tt$Gene <- if ("Gene" %in% names(tt)) tt$Gene else rownames(tt)
  tt$Gene[is.na(tt$Gene) | tt$Gene == ""] <- rownames(tt)[is.na(tt$Gene) | tt$Gene == ""]

  tt$Regulation <- "None"
  tt$Regulation[tt$logFC >  fc_cut & tt[[padj_col]] < 0.05] <- "Up"
  tt$Regulation[tt$logFC < -fc_cut & tt[[padj_col]] < 0.05] <- "Down"

  tt$neglog10 <- -log10(pmax(tt[[padj_col]], .Machine$double.xmin))
  lab_df <- subset(tt, Regulation != "None")

  ggplot(tt, aes(logFC, neglog10)) +
    geom_point(data = subset(tt, Regulation == "None"),
               color = "grey70", alpha = 0.35, size = 1) +
    geom_point(data = subset(tt, Regulation == "Up"),
               color = "red", size = 1.4) +
    geom_point(data = subset(tt, Regulation == "Down"),
               color = "blue", size = 1.4) +
    geom_text_repel(data = lab_df,
                    aes(label = Gene),
                    size = 2.5, max.overlaps = Inf, min.segment.length = 0) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    labs(title = title, x = expression(log[2]*" FC"), y = expression(-log[10]*" adj.P")) +
    theme_minimal(base_size = 12)
}

# ============================================================
# ① limma-voom (진단: voom mean-variance 1장 + Volcano 3페이지)
# ============================================================
{
  # 모델 적합
  v <- limma::voomWithQualityWeights(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)

  cont <- makeContrasts(
    HO_vs_DO = groupHealthy_Old   - groupDiabetes_Old,
    HY_vs_DO = groupHealthy_Young - groupDiabetes_Old,
    HO_vs_HY = groupHealthy_Old   - groupHealthy_Young,
    levels = colnames(coef(fit))
  )
  fit2 <-limma::eBayes(limma::contrasts.fit(fit, cont), trend=TRUE, robust=TRUE)

  res_limma <- list(
    HO_vs_DO = topTable(fit2, coef = "HO_vs_DO", n = Inf, sort.by = "P"),
    HY_vs_DO = topTable(fit2, coef = "HY_vs_DO", n = Inf, sort.by = "P"),
    HO_vs_HY = topTable(fit2, coef = "HO_vs_HY", n = Inf, sort.by = "P")
  )
  for (nm in names(res_limma)) {
    tt <- res_limma[[nm]]
    tt$Gene <- rownames(tt); rownames(tt) <- NULL
    write.table(tt, file = file.path(out_limma, paste0("DEG_", nm, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Volcano (DEG 라벨)
  pdf(file.path(out_limma, "Volcano_limma.pdf"), width = 7, height = 6)
  print(plot_volcano(res_limma$HO_vs_DO, "Volcano (limma): HO vs DO"))
  print(plot_volcano(res_limma$HY_vs_DO, "Volcano (limma): HY vs DO"))
  print(plot_volcano(res_limma$HO_vs_HY, "Volcano (limma): HO vs HY"))
  dev.off()
}

# ---- (limma-voom) PCA: batch 제거(=group 보호) 후 샘플 구조 확인 ----
Xg <- model.matrix(~ 0 + group, data = y$samples)  # group 보호
E_clean <- removeBatchEffect(v$E, batch = y$samples$batch, design = Xg)

set.seed(42)
pca <- prcomp(t(E_clean), scale. = TRUE)
pca_df <- data.frame(pca$x[,1:2],
                     group = y$samples$group,
                     batch = y$samples$batch,
                     sample_id = rownames(y$samples))

p <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = batch, label = sample_id)) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 2.6, max.overlaps = 30) +
  labs(title = "PCA (limma-voom, removeBatchEffect | design=~0+group)",
       subtitle = "시각화 전용: 통계는 ~0+group+batch 결과 사용",
       x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2,1]),
       y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2,2])) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_limma, "PCA_limma.png"), p, width = 7.5, height = 6.2, dpi = 300)

# ==========================================
# ② edgeR QL GLM (진단: BCV 1장 + Volcano 3페이지)
# ==========================================
{
  # 1) 분산 추정 + BCV 진단
  y2 <- edgeR::estimateDisp(y, design, robust = TRUE)
  cairo_pdf(file.path(out_edger, "BCV.pdf"), width = 6.5, height = 5.5, onefile = FALSE)
  edgeR::plotBCV(y2)
  dev.off()

  # 2) QL GLM 적합 (★ DGEGLM 객체)
  fitq <- edgeR::glmQLFit(y2, design, robust = TRUE)

  # 3) 대비 만들기
  contr_edger <- limma::makeContrasts(
    HO_vs_DO = groupHealthy_Old   - groupDiabetes_Old,
    HY_vs_DO = groupHealthy_Young - groupDiabetes_Old,
    HO_vs_HY = groupHealthy_Old   - groupHealthy_Young,
    levels = design
  )

  # 4) 테스트 + 저장
  res_edger <- list()
  for (nm in colnames(contr_edger)) {
    qlf <- edgeR::glmQLFTest(fitq, contrast = contr_edger[, nm])
    tt  <- edgeR::topTags(qlf, n = Inf)$table
    tt$Gene <- rownames(tt); rownames(tt) <- NULL
    res_edger[[nm]] <- tt
    write.table(tt, file = file.path(out_edger, paste0("DEG_", nm, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # 5) Volcano (FDR 컬럼 사용)
  cairo_pdf(file.path(out_edger, "Volcano_edger.pdf"), width = 7, height = 6, onefile = TRUE)
  print(plot_volcano(res_edger$HO_vs_DO, "Volcano (edgeR QL): HO vs DO", padj_col = "FDR"))
  print(plot_volcano(res_edger$HY_vs_DO, "Volcano (edgeR QL): HY vs DO", padj_col = "FDR"))
  print(plot_volcano(res_edger$HO_vs_HY, "Volcano (edgeR QL): HO vs HY", padj_col = "FDR"))
  dev.off()
}
