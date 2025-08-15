# ===============================================
# Run2-1 & Run2-2 전용 DEG (현재 방식)
#  - 샘플 단위: patient × batch (Run2_1, Run2_2 유지)
#  - Design: ~ 0 + group + batch
#  - 출력: 최소(진단 + Volcano + DEG 테이블)
# ===============================================

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr)
  library(ggplot2); library(ggrepel)
})

# ---- 0) IO ----
rds_path  <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT_with_markers.rds"
out_limma <- "deg_out_Run2_limma_voom"
out_edger <- "deg_out_Run2_edger_ql"
dir.create(out_limma, showWarnings = FALSE)
dir.create(out_edger, showWarnings = FALSE)

# ---- 1) Load ----
seu <- readRDS(rds_path)
counts <- tryCatch(
  GetAssayData(seu, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
)
meta <- seu@meta.data
stopifnot(ncol(counts) == nrow(meta))

# ---- 2) Run2-1 / Run2-2만 선택 ----
meta$batch <- meta$orig.ident
keep_r2 <- meta$batch %in% c("Run2_1","Run2_2")
meta2   <- meta[keep_r2, , drop = FALSE]
counts2 <- counts[, keep_r2, drop = FALSE]

# 그룹 레벨(해석 편의: DO, HO, HY 순)
meta2$group <- factor(meta2$group,
                      levels = c("Diabetes_Old","Healthy_Old","Healthy_Young"))
meta2$batch <- factor(meta2$batch, levels = c("Run2_1","Run2_2"))

# ---- 3) sample_id = patient × batch → pseudobulk ----
meta2$sample_id <- interaction(meta2$patient, meta2$batch, drop = TRUE)

idx <- split(seq_len(ncol(counts2)), meta2$sample_id, drop = TRUE)
pb_counts <- do.call("cbind", lapply(idx, function(i) Matrix::rowSums(counts2[, i, drop = FALSE])))

# (옵션) 미토 유전자 제외하려면 TRUE
remove_mt <- FALSE
if (remove_mt) {
  rm_mt <- grepl("^MT-", rownames(pb_counts), ignore.case = TRUE)
  pb_counts <- pb_counts[!rm_mt, , drop = FALSE]
}

# 샘플 메타
samples_info <- meta2 |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(patient = dplyr::first(patient),
                   group   = dplyr::first(group),
                   batch   = dplyr::first(batch), .groups = "drop") |>
  as.data.frame()
pb_counts <- pb_counts[, samples_info$sample_id]
rownames(samples_info) <- samples_info$sample_id

# ---- 4) DGEList + 필터 + TMM ----
y <- DGEList(counts = pb_counts, samples = samples_info)
y$samples$group <- droplevels(y$samples$group)
y$samples$batch <- droplevels(y$samples$batch)

keep <- filterByExpr(y, group = y$samples$group, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# ---- 5) Design: ~ 0 + group + batch ----
design <- model.matrix(~ 0 + group + batch, data = y$samples)

# ---- helper: Volcano (DEG만 라벨) ----
plot_volcano <- function(tt, title, padj_col = c("adj.P.Val","FDR"), fc_cut = 0.58) {
  padj_col <- padj_col[padj_col %in% names(tt)][1]
  if (is.na(padj_col)) { padj_col <- "padj"; tt[[padj_col]] <- p.adjust(tt$P.Value, "BH") }
  if (!("Gene" %in% names(tt))) tt$Gene <- rownames(tt)
  tt$Regulation <- "None"
  tt$Regulation[tt$logFC >  fc_cut & tt[[padj_col]] < 0.05] <- "Up"
  tt$Regulation[tt$logFC < -fc_cut & tt[[padj_col]] < 0.05] <- "Down"
  tt$neglog10 <- -log10(pmax(tt[[padj_col]], .Machine$double.xmin))
  lab_df <- subset(tt, Regulation != "None")

  ggplot(tt, aes(logFC, neglog10)) +
    geom_point(data = subset(tt, Regulation == "None"), color = "grey70", alpha = 0.35, size = 1) +
    geom_point(data = subset(tt, Regulation == "Up"),   color = "red",  size = 1.4) +
    geom_point(data = subset(tt, Regulation == "Down"), color = "blue", size = 1.4) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = Gene),
                             size = 2.5, max.overlaps = Inf, min.segment.length = 0) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    labs(title = title, x = expression(log[2]*" FC"), y = expression(-log[10]*" adj.P")) +
    theme_minimal(base_size = 12)
}

# =========================
# A) limma-voom (Run2_1 + Run2_2)
# =========================
{
  # 진단 1: voom mean-variance
  pdf(file.path(out_limma, "voom_mean_variance.pdf"), width = 6.2, height = 5.2)
  invisible(voom(y, design, plot = TRUE))
  dev.off()

  # 중복 환자 있으면 block=patient
  rep_counts <- table(y$samples$patient)
  use_block <- any(rep_counts >= 2)

  if (use_block) {
    v0 <- voom(y, design, plot = FALSE)
    dupcor <- duplicateCorrelation(v0, design = design, block = y$samples$patient)
    v <- voom(y, design, plot = FALSE, block = y$samples$patient, correlation = dupcor$consensus)
    fit <- lmFit(v, design, block = y$samples$patient, correlation = dupcor$consensus)
  } else {
    v <- voom(y, design, plot = FALSE)
    fit <- lmFit(v, design)
  }

  cont <- makeContrasts(
    HO_vs_DO = groupHealthy_Old   - groupDiabetes_Old,
    HY_vs_DO = groupHealthy_Young - groupDiabetes_Old,
    HO_vs_HY = groupHealthy_Old   - groupHealthy_Young,
    levels = colnames(coef(fit))
  )
  fit2 <- eBayes(contrasts.fit(fit, cont))

  res_limma <- list(
    HO_vs_DO = topTable(fit2, coef = "HO_vs_DO", n = Inf, sort.by = "P"),
    HY_vs_DO = topTable(fit2, coef = "HY_vs_DO", n = Inf, sort.by = "P"),
    HO_vs_HY = topTable(fit2, coef = "HO_vs_HY", n = Inf, sort.by = "P")
  )
  for (nm in names(res_limma)) {
    tt <- res_limma[[nm]]; tt$Gene <- rownames(tt); rownames(tt) <- NULL
    write.table(tt, file = file.path(out_limma, paste0("DEG_", nm, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # 진단 2: SA plot
  pdf(file.path(out_limma, "SA_plot.pdf"), width = 6.2, height = 5.2)
  plotSA(fit2, main = "limma-voom: SA plot (Run2_1+Run2_2)")
  dev.off()

  # Volcano(3페이지)
  pdf(file.path(out_limma, "Volcano_limma_Run2.pdf"), width = 7, height = 6)
  print(plot_volcano(res_limma$HO_vs_DO, "Volcano (limma): HO vs DO", padj_col = "adj.P.Val"))
  print(plot_volcano(res_limma$HY_vs_DO, "Volcano (limma): HY vs DO", padj_col = "adj.P.Val"))
  print(plot_volcano(res_limma$HO_vs_HY, "Volcano (limma): HO vs HY", padj_col = "adj.P.Val"))
  dev.off()
}

# =========================
# B) edgeR QL GLM (Run2_1 + Run2_2)
# =========================
{
  y2 <- estimateDisp(y, design, robust = TRUE)

  # 진단: BCV
  pdf(file.path(out_edger, "BCV.pdf"), width = 6.2, height = 5.2)
  plotBCV(y2)
  dev.off()

  fitq <- glmQLFit(y2, design, robust = TRUE)
  contr <- makeContrasts(
    HO_vs_DO = groupHealthy_Old   - groupDiabetes_Old,
    HY_vs_DO = groupHealthy_Young - groupDiabetes_Old,
    HO_vs_HY = groupHealthy_Old   - groupHealthy_Young,
    levels = design
  )

  res_edger <- list()
  for (nm in colnames(contr)) {
    qlf <- glmQLFTest(fitq, contrast = contr[, nm])
    tt  <- topTags(qlf, n = Inf)$table
    tt$Gene <- rownames(tt); rownames(tt) <- NULL
    res_edger[[nm]] <- tt
    write.table(tt, file = file.path(out_edger, paste0("DEG_", nm, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # Volcano(3페이지)
  pdf(file.path(out_edger, "Volcano_edger_Run2.pdf"), width = 7, height = 6)
  print(plot_volcano(res_edger$HO_vs_DO, "Volcano (edgeR QL): HO vs DO", padj_col = "FDR"))
  print(plot_volcano(res_edger$HY_vs_DO, "Volcano (edgeR QL): HY vs DO", padj_col = "FDR"))
  print(plot_volcano(res_edger$HO_vs_HY, "Volcano (edgeR QL): HO vs HY", padj_col = "FDR"))
  dev.off()
}

message("Done. Outputs:\n",
        "- ", out_limma, ": voom_mean_variance.pdf, SA_plot.pdf, Volcano_limma_Run2.pdf, 3×DEG_*.tsv\n",
        "- ", out_edger, ": BCV.pdf, Volcano_edger_Run2.pdf, 3×DEG_*.tsv")
