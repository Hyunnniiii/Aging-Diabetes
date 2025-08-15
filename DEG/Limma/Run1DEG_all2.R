# ===============================================
# Run1 전용 pseudobulk(DEG): 현재 방식 (limma-voom + edgeR-QL)
#  - 샘플 단위: patient × batch (Run1만 → 사실상 patient)
#  - Design: ~ 0 + group  (HY vs DO)
#  - 출력:
#     limma: voom_mean_variance.pdf, SA_plot.pdf, Volcano_limma_Run1.pdf, DEG_DO_vs_HY.tsv
#     edgeR: BCV.pdf, Volcano_edger_Run1.pdf,       DEG_DO_vs_HY.tsv
# ===============================================

suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr)
  library(ggplot2); library(ggrepel); library(tibble)
})

# ---- 0) IO ----
rds_path  <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT_with_markers.rds"
out_limma <- "deg_out_Run1_limma_voom"
out_edger <- "deg_out_Run1_edger_ql"
dir.create(out_limma, showWarnings = FALSE)
dir.create(out_edger, showWarnings = FALSE)

# ---- 1) Seurat & counts ----
seu <- readRDS(rds_path)
counts <- tryCatch(
  GetAssayData(seu, assay = "RNA", layer = "counts"),   # Seurat v5
  error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
)
stopifnot(inherits(counts, "dgCMatrix") || is.matrix(counts))

meta <- seu@meta.data
meta$batch <- meta$orig.ident
stopifnot(ncol(counts) == nrow(meta))

# ---- 2) Run1만 선택 + HY/DO만 남기기 ----
keep_run1 <- meta$batch == "Run1"
meta1  <- meta[keep_run1, , drop = FALSE]
counts1 <- counts[, keep_run1, drop = FALSE]

# group 레벨 정리 (HY/DO만)
meta1$group <- droplevels(factor(meta1$group,
                                 levels = c("Healthy_Young","Diabetes_Old")))
stopifnot(all(levels(meta1$group) %in% c("Healthy_Young","Diabetes_Old")))

# ---- 3) sample_id = patient × batch (Run1만이라 patient와 동일) ----
meta1$sample_id <- interaction(meta1$patient, meta1$batch, drop = TRUE)

# pseudobulk (환자×배치별 합산)
idx <- split(seq_len(ncol(counts1)), meta1$sample_id, drop = TRUE)
pb_counts <- do.call("cbind", lapply(idx, function(i) Matrix::rowSums(counts1[, i, drop = FALSE])))

# (옵션) 미토 유전자 제외하려면 TRUE로
remove_mt <- FALSE
if (remove_mt) {
  rm_mt <- grepl("^MT-", rownames(pb_counts), ignore.case = TRUE)
  pb_counts <- pb_counts[!rm_mt, , drop = FALSE]
}

# 샘플 메타
samples_info <- meta1 |>
  group_by(sample_id) |>
  summarise(patient = first(patient),
            group   = first(group), .groups = "drop") |>
  as.data.frame()
pb_counts <- pb_counts[, samples_info$sample_id]
rownames(samples_info) <- samples_info$sample_id

# ---- 4) DGEList + 필터 + TMM ----
y <- DGEList(counts = pb_counts, samples = samples_info)
y$samples$group <- droplevels(y$samples$group)  # HY / DO
keep <- filterByExpr(y, group = y$samples$group, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# 디자인: ~ 0 + group  (배치는 한 레벨뿐)
design <- model.matrix(~ 0 + group, data = y$samples)

# ---- 공통 Volcano 함수 (DEG만 라벨) ----
plot_volcano <- function(tt, title, padj_col = c("adj.P.Val","FDR"), fc_cut = 0.58) {
  padj_col <- padj_col[padj_col %in% names(tt)][1]
  if (is.na(padj_col)) {
    padj_col <- "padj"
    tt[[padj_col]] <- if ("P.Value" %in% names(tt)) p.adjust(tt$P.Value, "BH") else NA_real_
  }
  tt$Gene <- if ("Gene" %in% names(tt)) tt$Gene else rownames(tt)
  tt$Regulation <- "None"
  tt$Regulation[tt$logFC >  fc_cut & tt[[padj_col]] < 0.05] <- "Up"
  tt$Regulation[tt$logFC < -fc_cut & tt[[padj_col]] < 0.05] <- "Down"
  tt$neglog10 <- -log10(pmax(tt[[padj_col]], .Machine$double.xmin))
  lab_df <- subset(tt, Regulation != "None")

  ggplot(tt, aes(logFC, neglog10)) +
    geom_point(data = subset(tt, Regulation == "None"), color = "grey70", alpha = 0.35, size = 1) +
    geom_point(data = subset(tt, Regulation == "Up"),   color = "red",  size = 1.4) +
    geom_point(data = subset(tt, Regulation == "Down"), color = "blue", size = 1.4) +
    geom_text_repel(data = lab_df, aes(label = Gene), size = 2.5, max.overlaps = Inf) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    labs(title = title, x = expression(log[2]*" FC"), y = expression(-log[10]*" adj.P")) +
    theme_minimal(base_size = 12)
}

# =========================
# A) limma-voom (Run1)
# =========================
{
  # 진단 1: voom mean-variance
  pdf(file.path(out_limma, "voom_mean_variance.pdf"), width = 6.2, height = 5.2)
  invisible(voom(y, design, plot = TRUE))
  dev.off()

  v   <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)

  # 대비: DO_vs_HY (DO - HY)
  cont <- makeContrasts(DO_vs_HY = groupDiabetes_Old - groupHealthy_Young,
                        levels = colnames(coef(fit)))
  fit2 <- eBayes(contrasts.fit(fit, cont))

  res_limma <- topTable(fit2, coef = "DO_vs_HY", n = Inf, sort.by = "P")
  res_limma$Gene <- rownames(res_limma); rownames(res_limma) <- NULL
  write.table(res_limma, file = file.path(out_limma, "DEG_DO_vs_HY.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # 진단 2: SA plot
  pdf(file.path(out_limma, "SA_plot.pdf"), width = 6.2, height = 5.2)
  plotSA(fit2, main = "limma-voom: SA plot (Run1)")
  dev.off()

  # Volcano
  pdf(file.path(out_limma, "Volcano_limma_Run1.pdf"), width = 7, height = 6)
  print(plot_volcano(res_limma, "Volcano (limma-voom, Run1): DO vs HY", padj_col = "adj.P.Val"))
  dev.off()
}

# =========================
# B) edgeR QL GLM (Run1)
# =========================
{
  y2 <- estimateDisp(y, design, robust = TRUE)

  # 진단: BCV
  pdf(file.path(out_edger, "BCV.pdf"), width = 6.2, height = 5.2)
  plotBCV(y2)
  dev.off()

  fitq <- glmQLFit(y2, design, robust = TRUE)
  contr <- makeContrasts(DO_vs_HY = groupDiabetes_Old - groupHealthy_Young,
                         levels = design)
  qlf <- glmQLFTest(fitq, contrast = contr[, "DO_vs_HY"])
  res_edger <- topTags(qlf, n = Inf)$table
  res_edger$Gene <- rownames(res_edger); rownames(res_edger) <- NULL
  write.table(res_edger, file = file.path(out_edger, "DEG_DO_vs_HY.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Volcano
  pdf(file.path(out_edger, "Volcano_edger_Run1.pdf"), width = 7, height = 6)
  print(plot_volcano(res_edger, "Volcano (edgeR QL, Run1): DO vs HY", padj_col = "FDR"))
  dev.off()
}

message("Done. Outputs:\n",
        "- ", out_limma, ": voom_mean_variance.pdf, SA_plot.pdf, Volcano_limma_Run1.pdf, DEG_DO_vs_HY.tsv\n",
        "- ", out_edger, ": BCV.pdf, Volcano_edger_Run1.pdf, DEG_DO_vs_HY.tsv")

