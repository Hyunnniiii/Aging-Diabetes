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

# ===== 진단(환자 합본 버전: group만) =====
## ===== 진단(환자-합본: batch 무시, group만) =====
# 컷오프(화산도/DEG 요약과 동일하게)
fc_cut  <- 0.58   # |log2FC| >= 0.58 (≈1.5배)
padj_cut <- 0.05  # FDR < 0.05

# 1) 환자 기준 pseudobulk 만들기 (batch 합치기)
meta$patient <- factor(meta$patient)
idx_pat <- split(seq_len(ncol(counts)), meta$patient, drop = TRUE)
pb_counts_pat <- do.call("cbind", lapply(idx_pat, function(i) Matrix::rowSums(counts[, i, drop = FALSE])))

# 미토 유전자 제거(일관성 유지)
rm_mt_pat <- grepl("^MT-", rownames(pb_counts_pat), ignore.case = TRUE)
pb_counts_pat <- pb_counts_pat[!rm_mt_pat, , drop = FALSE]

# 2) 환자 수준 메타데이터
samples_info_pat <- meta %>%
  group_by(patient) %>%
  summarise(group = first(group), .groups = "drop") %>%
  as.data.frame()
# 열 정렬 일치
pb_counts_pat <- pb_counts_pat[, as.character(samples_info_pat$patient)]
rownames(samples_info_pat) <- as.character(samples_info_pat$patient)

cat("\n[patient-merged] Group counts:\n")
print(table(samples_info_pat$group))

# 3) edgeR 컨테이너 + 필터/정규화
y0 <- edgeR::DGEList(counts = pb_counts_pat, samples = samples_info_pat)
y0$samples$group <- factor(y0$samples$group, levels = levels(meta$group))

keep <- edgeR::filterByExpr(y0, group = y0$samples$group, min.count = 10)
y_diag  <- y0[keep, , keep.lib.sizes = FALSE]
y_diag  <- edgeR::calcNormFactors(y_diag)

cat("Genes before/after filter:", nrow(y0), "/", nrow(y_diag), "\n")
print(y_diag$samples[, c("lib.size","norm.factors","group")])

# 4) design(~0+group) 구성 및 full-rank 확인
design_only <- model.matrix(~ 0 + group, data = y_diag$samples)
cat("Design uses batch?  FALSE  | is.fullrank:", limma::is.fullrank(design_only), "\n")

# 5) 대비 자동 구성 + DEG 개수 집계
mk_contr <- function(cols) {
  lst <- list()
  if (all(c("groupHealthy_Old","groupDiabetes_Old")   %in% cols)) lst$HO_vs_DO <- "groupHealthy_Old - groupDiabetes_Old"
  if (all(c("groupHealthy_Young","groupDiabetes_Old") %in% cols)) lst$HY_vs_DO <- "groupHealthy_Young - groupDiabetes_Old"
  if (all(c("groupHealthy_Old","groupHealthy_Young")  %in% cols)) lst$HO_vs_HY <- "groupHealthy_Old - groupHealthy_Young"
  lst
}

count_sig <- function(design, yobj, fc_cut, padj_cut) {
  if (!limma::is.fullrank(design)) return(NULL)
  v   <- limma::voom(yobj, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  contr <- mk_contr(colnames(design))
  if (!length(contr)) return(NULL)
  cont <- limma::makeContrasts(contrasts = unlist(contr), levels = colnames(coef(fit)))
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  do.call(rbind, lapply(colnames(cont), function(nm) {
    tt <- limma::topTable(fit2, coef = nm, n = Inf, sort.by = "P")
    data.frame(
      contrast    = nm,
      genes_tested= nrow(tt),
      n_sig       = sum(tt$adj.P.Val < padj_cut & abs(tt$logFC) > fc_cut, na.rm = TRUE)
    )
  }))
}

res_only <- count_sig(design_only, y_diag, fc_cut, padj_cut)
cat("\n[patient-merged] 유의 DEG 개수(|logFC|>", fc_cut, ", FDR <", padj_cut, ")\n", sep = "")
print(res_only)
cat("==== 진단 끝 ====\n")
