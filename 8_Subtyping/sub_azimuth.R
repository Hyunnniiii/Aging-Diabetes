## =========================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | T/NK/Mono subset의 RNA counts + meta.data 추출해 Azimuth reference 파일 생성 
## 웹에서 Azimuth 돌린 결과 파일 매핑 후 UMAP으로 확인
## =========================================


# --- v5 reference 생성 ---
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# --- 원본 v5 객체 (RNA counts의 원천) ---
orig_path <- "/data/project/diabetes_LYH/tanya/rds/sct_annotated_1.rds"
orig <- readRDS(orig_path)
stopifnot("RNA" %in% Assays(orig))

rna_counts_all <- tryCatch({ GetAssayData(orig, assay = "RNA", layer = "counts") })
stopifnot(inherits(rna_counts_all, "dgCMatrix"))

# --- 라인리지별 서브 통합 객체에서 바코드/메타만 사용 ---
lineages    <- c("T", "NK", "Mono")
input_root  <- "Subtyping"
output_root <- "AzimuthReferences"
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

for (lin in lineages) {
  message(sprintf("\n[%s] ---", lin))
  sub_rds <- file.path(input_root, lin, sprintf("subset_%s_SCTintegrated.rds", lin))
  stopifnot(file.exists(sub_rds))
  sub <- readRDS(sub_rds)

  # 이 라인리지의 셀 목록
  cells <- colnames(sub)
  stopifnot(length(cells) > 0)

  # 원본과 교집합(혹시 이름 불일치 대비)
  cells_use <- cells[cells %in% colnames(orig)]
  if (length(cells_use) < length(cells)) {
    warning(sprintf("[%s] %d cells not found in original; using %d cells.",
                    lin, length(setdiff(cells, cells_use)), length(cells_use)))
  }
  stopifnot(length(cells_use) > 0)

  # 원본에서 RNA counts 슬라이스(희소행렬 유지)
  counts <- rna_counts_all[, cells_use, drop = FALSE]

  # 메타데이터는 sub에서 동일 바코드만 정렬하여 사용
  meta <- sub@meta.data[cells_use, , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(meta)))

  # payload 저장 (이건 v5/v4 무관: 단순 리스트)
  saveRDS(
    list(counts = counts, meta = meta),
    file = file.path(output_root, sprintf("v4_payload_%s.rds", lin))
  )
  message(sprintf("[%s] wrote payload: %s",
                  lin, normalizePath(file.path(output_root, sprintf("v4_payload_%s.rds", lin)))))
}
message("\n[v5] Done. Now run step #2 in a Seurat v4 environment.")


####### v4 객체 만들기 #######
suppressPackageStartupMessages({
  library(Seurat)        # <-- 반드시 v4가 로드된 환경
  library(Matrix)
})

lineages    <- c("T", "NK", "Mono")
input_root  <- "AzimuthReferences"
output_root <- "AzimuthReferences"

for (lin in lineages) {
  payload_path <- file.path(input_root, sprintf("v4_payload_%s.rds", lin))
  stopifnot(file.exists(payload_path))
  payload <- readRDS(payload_path)

  # 안전 점검
  stopifnot(inherits(payload$counts, "dgCMatrix"))
  stopifnot(is.data.frame(payload$meta))
  stopifnot(identical(colnames(payload$counts), rownames(payload$meta)))

  # v4 Seurat 객체 생성 (Assay 클래스)
  seu <- CreateSeuratObject(
    counts   = payload$counts,
    meta.data= payload$meta,
    project  = paste0("AziRef_", lin)
  )

  # 저장
  out <- file.path(output_root, sprintf("v4_seurat_aziref_%s.rds", lin))
  saveRDS(seu, out)
  message(sprintf("[%s] saved v4 Seurat: %s", lin, normalizePath(out)))
}

message("\n[v4] All done.")


##### web에서 Azimuth 돌린 후 ####
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

lineages <- c("T","NK","Mono")
sub_root <- "Subtyping"                           # subset_{LIN}_SCTintegrated.rds 위치
azi_root <- file.path(sub_root, "Azimuth_sub")    # Azimuth_sub/azimuth_pred_<LIN>.tsv
out_root <- "AzimuthPlots"; dir.create(out_root, FALSE, TRUE)

strip <- function(x) {                            # 접미사 -1 / _1 제거
  x <- as.character(x)
  x <- sub("-\\d+$","",x)
  x <- sub("_[0-9]+$","",x)
  x
}

for (lin in lineages) {
  message(sprintf("\n[%s] ---", lin))

  # 1) 라인리지별 Seurat subset 로드
  obj <- readRDS(file.path(sub_root, lin, sprintf("subset_%s_SCTintegrated.rds", lin)))

  # 2) 해당 라인리지 TSV 로드
  anno_path <- file.path(azi_root, sprintf("azimuth_pred_%s.tsv", lin))
  stopifnot(file.exists(anno_path))
  anno <- read.delim(anno_path, check.names = FALSE)

  # 3) 바코드 매칭 (우선 exact, 안 맞으면 접미사 제거 매칭)
  m <- match(colnames(obj), anno$cell)                           # ← 바코드 컬럼 = "cell"
  if (any(is.na(m))) m <- match(strip(colnames(obj)), strip(anno$cell))

  matched <- sum(!is.na(m))
  message(sprintf("[%s] matched %d / %d", lin, matched, ncol(obj)))
  stopifnot(matched > 0)

  # 4) 라벨 붙이고 UMAP 저장
  obj$azimuth_l2 <- anno$predicted.celltype.l2[m]                # ← 라벨 컬럼 = "predicted.celltype.l2"
  Idents(obj) <- "azimuth_l2"

  p <- DimPlot(obj, reduction="umap", group.by="azimuth_l2",
               label=TRUE, repel=TRUE) + coord_fixed()
  ggsave(file.path(out_root, sprintf("Azimuth_l2_%s.pdf", lin)),
         p, width=8, height=8, units="in")
}
message("\nDone.")
