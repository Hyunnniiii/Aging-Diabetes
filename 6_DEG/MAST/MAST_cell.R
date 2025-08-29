# ============================================================
## Date | 2025-08-29 | Yoon JiHyun
## Harmony로 추가 보정한 결과 파일로 이후 분석 진행

## Description | 각 celltype별로 HY·HO·DO 간 MAST DEG 수행, 각 결과를 TSV·RDS로 저장
# ============================================================


suppressPackageStartupMessages({
  library(Seurat); library(MAST); library(Matrix); library(dplyr); library(readr)
})

seu <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')
DefaultAssay(seu) <- 'RNA'

grp_col <- 'group'
Idents(seu) <- factor(seu[[grp_col]],
                      levels = c("Healthy_Young","Healthy_Old","Diabetes_Old"))

# --- celltype 한 개에 대해 한 대비 돌리기 ---
run_mast_celltype <- function(obj, celltype_value,
                              ident1, ident2,
                              latent = c("nFeature_RNA","percent.mt","Run"),
                              min_pct = 0.10, logfc_thresh = 0,
                              assay = "RNA", layer = "data") {

sub <- subset(obj, subset = celltype == !!celltype_value)
if (ncol(sub) == 0) { message("[SKIP] no cells in ", celltype_value); return(NULL) }

  Idents(sub) <- factor(sub[[grp_col]],
                        levels = c("Healthy_Young","Healthy_Old","Diabetes_Old"))

  n1 <- sum(Idents(sub) == ident1); n2 <- sum(Idents(sub) == ident2)
  if (n1 == 0 || n2 == 0) {
    message(sprintf("[SKIP] %s | %s vs %s (cells: %d vs %d)", celltype_value, ident1, ident2, n1, n2))
    return(NULL)
  }

  latent <- intersect(latent, colnames(sub@meta.data))
  message(sprintf("[MAST][%s] %s vs %s | cells %d/%d | latent: %s",
                  celltype_value, ident1, ident2, n1, n2, paste(latent, collapse=", ")))

  res <- FindMarkers(
    object = sub,
    ident.1 = ident1, ident.2 = ident2,
    test.use = "MAST",
    assay = assay, layer = layer,     # Seurat v5 권장
    min.pct = min_pct, logfc.threshold = logfc_thresh,
    latent.vars = latent, verbose = FALSE
  ) %>% tibble::rownames_to_column("gene")

  # 저장
  base <- sprintf("%s__%s_vs_%s",
                gsub("[^A-Za-z0-9._-]+","_", celltype_value), ident1, ident2)
  out_dir <- "/home/tanya0721/PROJECT/DEG/MAST"
  out_tsv <- file.path(out_dir, paste0(base, "_MAST.tsv"))
  out_rds <- file.path(out_dir, paste0(base, "_MAST.rds"))
  write_tsv(res, out_tsv)
  saveRDS(res, out_rds)
  message("Saved: ", basename(out_tsv), " / ", basename(out_rds))}

  cts <- sort(unique(seu$celltype))

for (ct in cts) {
  run_mast_celltype(seu, ct, "Healthy_Young", "Diabetes_Old")
  run_mast_celltype(seu, ct, "Healthy_Young", "Healthy_Old")
  run_mast_celltype(seu, ct, "Healthy_Old",  "Diabetes_Old")
}