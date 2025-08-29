# ============================================================
## Date | 2025-08-29 | Yoon JiHyun
## Harmony로 추가 보정한 결과 파일로 이후 분석 진행

## Description | HY·HO·DO 세 그룹 간 차이를 MAST로 DEG 분석 
## 각 결과를 TSV/RDS로 저장하며 유의 유전자 수를 콘솔에 요약
# ============================================================


suppressPackageStartupMessages({library(Seurat)
library(MAST)
library(Matrix)
library(dplyr)
library(readr)})

seu <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')
DefaultAssay(seu) <- 'RNA'

# HY/HO/DO가 들어있는 컬럼명 지정
grp_col <- 'group'
stopifnot(grp_col %in% colnames(seu@meta.data))
Idents(seu) <- factor(seu@meta.data[[grp_col]],
                      levels = c("Healthy_Young","Healthy_Old","Diabetes_Old"))

# ---- 1) 헬퍼: 한 비교를 MAST로 수행 ----
run_mast <- function(obj, ident1, ident2,
                     out_prefix,
                     min_pct = 0.10,       # 각 그룹에서 최소 10% 셀에서 발현된 유전자만 테스트 (연산량/잡음 감소)
                     logfc_thresh = 0,     # 사전 필터링 제거 (MAST에 맡김)
                     latent = c("nFeature_RNA","percent.mt","Run"),
                     assay = "RNA", slot = "data") {

  message(sprintf("[MAST] %s vs %s  | latent: %s",
                  ident1, ident2, paste(latent, collapse = ", ")))

  # 존재하는 latent만 사용
  latent <- intersect(latent, colnames(obj@meta.data))

  # 실제 비교에 포함될 셀 수 확인
  n1 <- sum(Idents(obj) == ident1)
  n2 <- sum(Idents(obj) == ident2)
  message(sprintf("Cells: %s=%d, %s=%d", ident1, n1, ident2, n2))

  # FindMarkers with MAST
  res <- FindMarkers(
    object = obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    assay = assay,
    slot = slot,
    min.pct = min_pct,
    logfc.threshold = logfc_thresh,
    latent.vars = latent,
    verbose = FALSE
  )

  # 결과 정리: gene 열 추가 및 정렬
  res <- res %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(abs(ifelse("avg_log2FC" %in% names(.), avg_log2FC,
                                       ifelse("avg_logFC" %in% names(.), avg_logFC, 0)))))

  # 저장
  out_tsv <- paste0(out_prefix, "_MAST.tsv")
  readr::write_tsv(res, out_tsv)
  message(sprintf("Saved: %s  (Sig FDR<0.05: %d)", out_tsv, sum(res$p_val_adj < 0.05, na.rm = TRUE)))

  return(res)
}

# ---- 2) 실행: 3가지 대비 ----
dir.create("MAST_DEG", showWarnings = FALSE)

res_HY_DO <- run_mast(seu, "Healthy_Young", "Diabetes_Old",
                      out_prefix = "MAST_DEG/HY_vs_DO",
                      latent = c("nFeature_RNA","percent.mt","Run"))

res_HY_HO <- run_mast(seu, "Healthy_Young", "Healthy_Old",
                      out_prefix = "MAST_DEG/HY_vs_HO",
                      latent = c("nFeature_RNA","percent.mt","Run"))

res_HO_DO <- run_mast(seu, "Healthy_Old", "Diabetes_Old",
                      out_prefix = "MAST_DEG/HO_vs_DO",
                      latent = c("nFeature_RNA","percent.mt","Run"))

# ---- 3) 요약 출력 ----
summ <- data.frame(
  contrast = c("HY vs DO", "HY vs HO", "HO vs DO"),
  n_sig_FDR_0.05 = c(sum(res_HY_DO$p_val_adj < 0.05, na.rm = TRUE),
                     sum(res_HY_HO$p_val_adj < 0.05, na.rm = TRUE),
                     sum(res_HO_DO$p_val_adj < 0.05, na.rm = TRUE))
)
print(summ)

# 1) 개별 RDS 백업 (재사용 쉬움)
saveRDS(res_HY_DO, "MAST_DEG/HY_vs_DO_MAST.rds")
saveRDS(res_HY_HO, "MAST_DEG/HY_vs_HO_MAST.rds")
saveRDS(res_HO_DO, "MAST_DEG/HO_vs_DO_MAST.rds")
