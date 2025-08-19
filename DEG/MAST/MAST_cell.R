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


# 경로/파일명 바꿔도 됨
SCRIPT=/home/tanya0721/PROJECT/DEG/MAST/MAST_cell.R
LOG=/home/tanya0721/PROJECT/DEG/MAST/MAST_cell_$(date +%Y%m%d_%H%M).log

# 30코어(0-29)로 고정 + BLAS 스레드 1개로 제한 + DEG 환경에서 실행
nohup env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
taskset -c 0-29 mamba run -n DEG Rscript "$SCRIPT" > "$LOG" 2>&1 &

echo "PID: $!  |  log: $LOG"

#######
tail -f "$LOG"     # 로그 실시간 보기 (Ctrl+C로 종료)
ps -p $! -o pid,pcpu,pmem,cmd   # 방금 띄운 PID로 상태 확인


LOG=/home/tanya0721/PROJECT/DEG/MAST/MAST_cell_$(date +%Y%m%d_%H%M).log
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
taskset -c 0-29 mamba run -n DEG Rscript "$SCRIPT" > "$LOG" 2>&1 < /dev/null &
echo "PID: $!  |  log: $LOG"
sleep 1
ps -o pid,ppid,stat,pcpu,pmem,etime,cmd -p $!    # PPID가 1이면 세션과 분리 OK
tail -n 30 -f "$LOG"