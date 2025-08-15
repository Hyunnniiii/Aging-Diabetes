suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr); library(ggplot2)
})

# ===== User 옵션 =====
rds_path <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT_with_markers.rds"
out_dir  <- "deg_out_noBatch_patientMerge_ALL3"
fc_cut   <- 0.58     # |log2FC|
padj_cut <- 0.05     # FDR
remove_mt <- FALSE   # TRUE면 미토 유전자 제거
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===== Load =====
seu <- readRDS(rds_path)
counts <- tryCatch(
  GetAssayData(seu, assay = "RNA", layer = "counts"),   # Seurat v5
  error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
)
meta <- seu@meta.data
stopifnot(ncol(counts) == nrow(meta), all(colnames(counts) == rownames(meta)))

# ===== 세 그룹만 사용 (DO/HO/HY) =====
meta2   <- subset(meta, group %in% c("Diabetes_Old","Healthy_Old","Healthy_Young"))
meta2$group <- factor(meta2$group, levels = c("Diabetes_Old","Healthy_Old","Healthy_Young"))
counts2 <- counts[, rownames(meta2), drop = FALSE]

# ===== 환자 단위로 모든 런 합치기 (patient-merge) =====
meta2$sample_id_patient <- factor(meta2$patient)
idx_pat <- split(seq_len(ncol(counts2)), meta2$sample_id_patient, drop = TRUE)
pb_counts <- do.call(cbind, lapply(idx_pat, function(i) Matrix::rowSums(counts2[, i, drop = FALSE])))
colnames(pb_counts) <- names(idx_pat)

if (remove_mt) {
  rm_mt <- grepl("^MT-", rownames(pb_counts), ignore.case = TRUE)
  pb_counts <- pb_counts[!rm_mt, , drop = FALSE]
}

samples_info <- meta2 %>%
  group_by(sample_id_patient) %>%
  summarise(patient = first(patient),
            group   = first(group), .groups = "drop") %>%
  as.data.frame()
rownames(samples_info) <- samples_info$sample_id_patient
pb_counts <- pb_counts[, rownames(samples_info), drop = FALSE]

# ===== edgeR → voom(QW) → limma (배치항 없음) =====
y <- edgeR::DGEList(pb_counts, samples = samples_info)
keep <- edgeR::filterByExpr(y, group = samples_info$group, min.count = 10)
y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])

design <- model.matrix(~ 0 + group, data = y$samples)
v   <- limma::voomWithQualityWeights(y, design, plot = FALSE)
fit <- limma::lmFit(v, design)

# 3개 대비
cont <- limma::makeContrasts(
  HO_vs_DO = groupHealthy_Old   - groupDiabetes_Old,
  HY_vs_DO = groupHealthy_Young - groupDiabetes_Old,
  HO_vs_HY = groupHealthy_Old   - groupHealthy_Young,
  levels = colnames(coef(fit))
)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)

# ===== 결과 저장 =====
save_tt <- function(obj, coef_name, out_dir) {
  tt <- limma::topTable(obj, coef = coef_name, n = Inf, sort.by = "P")
  tt$Gene <- rownames(tt); rownames(tt) <- NULL
  fn <- file.path(out_dir, paste0("DEG_", coef_name, "_noBatch_patientMerge.tsv"))
  write.table(tt, file = fn, sep = "\t", quote = FALSE, row.names = FALSE)
  tt
}
tt_list <- list(
  HO_vs_DO = save_tt(fit2, "HO_vs_DO", out_dir),
  HY_vs_DO = save_tt(fit2, "HY_vs_DO", out_dir),
  HO_vs_HY = save_tt(fit2, "HO_vs_HY", out_dir)
)

# ===== Volcano helper (라벨 없이 색만) =====
plot_volcano <- function(tt, title, fc_cut = 0.58, padj_cut = 0.05, padj_col = "adj.P.Val") {
  if (!(padj_col %in% names(tt))) {
    if ("adj.P.Val" %in% names(tt)) padj_col <- "adj.P.Val"
    else if ("FDR" %in% names(tt))  padj_col <- "FDR"
    else if ("P.Value" %in% names(tt)) { tt$adj.P.Val <- p.adjust(tt$P.Value, "BH"); padj_col <- "adj.P.Val" }
    else stop("No p-values in table.")
  }
  df <- tt
  df$Reg <- "None"
  df$Reg[df$logFC >  fc_cut & df[[padj_col]] < padj_cut] <- "Up"
  df$Reg[df$logFC < -fc_cut & df[[padj_col]] < padj_cut] <- "Down"
  df$neglog10 <- -log10(pmax(df[[padj_col]], .Machine$double.xmin))

  ggplot(df, aes(logFC, neglog10)) +
    geom_point(data = subset(df, Reg == "None"), color = "grey75", alpha = 0.45, size = 1) +
    geom_point(data = subset(df, Reg == "Up"),   color = "red",   alpha = 0.9,  size = 1.4) +
    geom_point(data = subset(df, Reg == "Down"), color = "blue",  alpha = 0.9,  size = 1.4) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2) +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2) +
    labs(title = title, x = expression(log[2]*" FC"), y = expression(-log[10]*" FDR")) +
    theme_minimal(base_size = 12)
}

# ===== Volcano PDF (3 페이지) =====
cairo_pdf(file.path(out_dir, "Volcano_noBatch_patientMerge_ALL3.pdf"), width = 7, height = 6, onefile = TRUE)
print(plot_volcano(tt_list$HO_vs_DO, "Volcano: HO vs DO (no batch, patient-merged)",
                   fc_cut = fc_cut, padj_cut = padj_cut, padj_col = "adj.P.Val"))
print(plot_volcano(tt_list$HY_vs_DO, "Volcano: HY vs DO (no batch, patient-merged)",
                   fc_cut = fc_cut, padj_cut = padj_cut, padj_col = "adj.P.Val"))
print(plot_volcano(tt_list$HO_vs_HY, "Volcano: HO vs HY (no batch, patient-merged)",
                   fc_cut = fc_cut, padj_cut = padj_cut, padj_col = "adj.P.Val"))
dev.off()

# ===== 간단 요약 =====
summ_line <- function(tt, name) {
  paste0(name,
         " | genes=", nrow(tt),
         " | Up=",   sum(tt$adj.P.Val < padj_cut & tt$logFC >  fc_cut, na.rm = TRUE),
         " | Down=", sum(tt$adj.P.Val < padj_cut & tt$logFC < -fc_cut, na.rm = TRUE))
}
cat(
  "\n[no-batch patient-merge]\n",
  summ_line(tt_list$HO_vs_DO, "HO_vs_DO"), "\n",
  summ_line(tt_list$HY_vs_DO, "HY_vs_DO"), "\n",
  summ_line(tt_list$HO_vs_HY, "HO_vs_HY"), "\n", sep = ""
)
message("Saved to: ", out_dir)
