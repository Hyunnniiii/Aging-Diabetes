suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(limma)
  library(Matrix); library(dplyr)
  library(ggplot2); library(ggrepel); library(tibble)
})

if (!exists("fc_cut",   inherits = TRUE)) fc_cut   <- 0.58
if (!exists("padj_cut", inherits = TRUE)) padj_cut <- 0.05

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

# ======= 여기에 붙여넣기 =======
# ---- pseudobulk 만들고 나서 바로 아래 블록 추가 ----

# 0) cbind로 만든 컬럼명 붙이기 + samples_info 생성/정렬(매우 중요)
colnames(pb_counts) <- names(idx)   # ← 없으면 subscript 오류납니다

samples_info <- meta2 %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(patient = dplyr::first(patient),
                   group   = dplyr::first(group),
                   batch   = dplyr::first(batch), .groups = "drop") %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  as.data.frame()

keep_ids <- intersect(samples_info$sample_id, colnames(pb_counts))
if (length(keep_ids) < 2) stop("pseudobulk 샘플 교집합이 너무 적습니다.")
samples_info <- samples_info[match(keep_ids, samples_info$sample_id), , drop = FALSE]
pb_counts    <- pb_counts[, keep_ids, drop = FALSE]
rownames(samples_info) <- samples_info$sample_id

# 1) 그룹×배치 분포(Confounding 확인)
cat("\n[Run2] Group × Batch 분포:\n")
print(table(samples_info$group, samples_info$batch))

# 2) 필터/정규화 후 라이브러리 정보
y0 <- edgeR::DGEList(counts = pb_counts, samples = samples_info)
keep <- edgeR::filterByExpr(y0, group = samples_info$group, min.count = 10)
y  <- y0[keep, , keep.lib.sizes = FALSE]
y  <- edgeR::calcNormFactors(y)
cat("Genes before/after filter:", nrow(y0), "/", nrow(y), "\n")
print(y$samples[, c("lib.size","norm.factors","group","batch")])

# 3) 디자인 랭크(설계가 추정 가능한지)
use_batch <- length(unique(samples_info$batch)) > 1
design_b  <- if (use_batch) model.matrix(~ 0 + group + batch, data = samples_info) else model.matrix(~ 0 + group, data = samples_info)
cat("Design uses batch? ", use_batch, " | is.fullrank:", limma::is.fullrank(design_b), "\n")

# 4) 대비 목록 만들기(존재하는 그룹만)
mk_contr <- function(cols) {
  lst <- list()
  if (all(c("groupHealthy_Old","groupDiabetes_Old")   %in% cols)) lst$HO_vs_DO <- "groupHealthy_Old - groupDiabetes_Old"
  if (all(c("groupHealthy_Young","groupDiabetes_Old") %in% cols)) lst$HY_vs_DO <- "groupHealthy_Young - groupDiabetes_Old"
  if (all(c("groupHealthy_Old","groupHealthy_Young")  %in% cols)) lst$HO_vs_HY <- "groupHealthy_Old - groupHealthy_Young"
  lst
}

# 5) 배치 포함 vs 제외 DEG 수 비교(Confounding 진단)
count_sig <- function(design) {
  if (!limma::is.fullrank(design)) return(NULL)
  v   <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  contr <- mk_contr(colnames(design))
  if (!length(contr)) return(NULL)
  cont <- limma::makeContrasts(contrasts = unlist(contr), levels = colnames(coef(fit)))
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  out <- lapply(colnames(cont), function(nm) {
    tt <- limma::topTable(fit2, coef = nm, n = Inf, sort.by = "P")
    data.frame(contrast = nm,
               genes_tested = nrow(tt),
               n_sig = sum(tt$adj.P.Val < padj_cut & abs(tt$logFC) > fc_cut, na.rm = TRUE))
  })
  do.call(rbind, out)
}

res_with_batch    <- if (use_batch) count_sig(model.matrix(~ 0 + group + batch, data = samples_info)) else NULL
res_without_batch <- count_sig(model.matrix(~ 0 + group,        data = samples_info))

cat("\n[Run2] 유의 DEG 개수 비교(|logFC|>", fc_cut, ", FDR <", padj_cut, ")\n", sep = "")
if (!is.null(res_with_batch)) {
  cat(" with batch:\n"); print(res_with_batch)
} else {
  cat(" with batch: (배치 레벨 1개 → 동일)\n")
}
cat(" without batch:\n"); print(res_without_batch)
cat("==== 진단 끝 ====\n")
