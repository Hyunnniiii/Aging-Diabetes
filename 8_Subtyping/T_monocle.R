## =========================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | 원본 rds의 RNA raw counts로 T subset을 재구성한 뒤 monocle3로 UMAP·클러스터링·그래프를 학습
## NaiveT 스코어/최대클러스터 기반 루트를 자동 지정해 pseudotime을 계산, 결과(CDS·UMAP·pseudotime·그룹별 박스플롯)를 저장
## =========================================


# Subtyping/Monocle3_T/run_monocle3_T_no_wrappers.R
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(monocle3)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
ts  <- function() format(Sys.time(), "%H:%M:%S")
msg <- function(...) cat(sprintf("[%s] %s\n", ts(), sprintf(...)))

## ===== 경로/설정 =====
PARENT_RDS <- "/data/project/diabetes_LYH/tanya/rds/sct_annotated_1.rds"
T_RDS      <- "/data/project/diabetes_LYH/tanya/rds/subset_T_SCTintegrated.rds"

out_dir <- file.path("Subtyping","Monocle3_T"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

## ===== 헬퍼 =====
# 부모의 RNA raw counts를 서브셋 객체 RNA assay에 주입해서 확실한 count 행/열이름 확보
rebuild_rna_from_parent <- function(subobj, parent, assay="RNA") {
  cells <- colnames(subobj)
  if (!all(cells %in% colnames(parent))) {
    stop("Some cells in subset are not present in parent.")
  }
  # Seurat v5 Assay5 counts layer 사용
  if (inherits(parent[[assay]], "Assay5")) {
    stopifnot("counts" %in% Layers(parent[[assay]]))
    cnt <- GetAssayData(parent, assay=assay, layer="counts")[, cells, drop = FALSE]
  } else {
    cnt <- GetAssayData(parent, assay=assay, slot ="counts")[, cells, drop = FALSE]
  }
  subobj[[assay]] <- CreateAssayObject(counts = cnt)
  DefaultAssay(subobj) <- assay
  subobj
}

# Seurat → monocle3 cell_data_set 변환
to_cds <- function(seu, assay="RNA"){
  # counts 추출
  counts <- if (inherits(seu[[assay]], "Assay5")) {
    GetAssayData(seu, assay=assay, layer="counts")
  } else {
    GetAssayData(seu, assay=assay, slot ="counts")
  }
  # gene/cell 메타 구성
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Counts matrix must have row and column names.")
  }
  gene_metadata <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )
  cell_metadata <- seu@meta.data
  # CDS 생성
  new_cell_data_set(counts,
                    cell_metadata = cell_metadata,
                    gene_metadata = gene_metadata)
}

## ===== 로드 & counts 보정 =====
msg("Loading T subset...")
sobj <- readRDS(T_RDS)

msg("Loading parent (for raw RNA counts)...")
parent <- readRDS(PARENT_RDS)

msg("Rebuilding RNA assay from parent raw counts...")
sobj <- rebuild_rna_from_parent(sobj, parent, assay="RNA")
DefaultAssay(sobj) <- "RNA"

## ===== naïve T marker 스코어(루트셀 힌트) =====
naive_markers <- list(c("SELL","CCR7","TCF7","IL7R","LTB"))
present <- intersect(unique(unlist(naive_markers)), rownames(sobj[["RNA"]]))
if (length(present) == 0) {
  msg("⚠️ Naive markers not found in RNA features. Will use fallback for root selection.")
}
sobj <- tryCatch({
  AddModuleScore(sobj, features = naive_markers, name = "NaiveT", assay = "RNA")
}, error = function(e) { msg("AddModuleScore failed: %s", e$message); sobj })

# 보기 좋게 group/cluster 컬럼 파악
grp_col  <- if ("group" %in% colnames(sobj@meta.data)) "group" else NULL
clus_col <- if ("seurat_clusters" %in% colnames(sobj@meta.data)) "seurat_clusters" else NULL

## ===== Seurat → Monocle3 =====
msg("Converting Seurat → cell_data_set ...")
cds <- to_cds(sobj, assay="RNA")

# gene_short_name 보장(위에서 이미 넣었지만 안전망)
if (is.null(rowData(cds)$gene_short_name)) {
  rowData(cds)$gene_short_name <- rownames(cds)
}

# 컬러링용 메타 유지
if (!is.null(grp_col))  colData(cds)$group  <- colData(cds)[[grp_col]]
if (!is.null(clus_col)) colData(cds)$seurat_clusters <- colData(cds)[[clus_col]]

## ===== Monocle3 파이프라인 =====
msg("preprocess_cds ...")
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")

msg("reduce_dimension (UMAP) ...")
cds <- reduce_dimension(cds, reduction_method = "UMAP")

msg("cluster_cells ...")
cds <- cluster_cells(cds, reduction_method = "UMAP")

msg("learn_graph ...")
cds <- learn_graph(cds)

## ===== root cell 자동 지정 =====
msg("Selecting root cells ...")
root_cells <- NULL
nmcol <- grep("^NaiveT", colnames(colData(cds)), value = TRUE)
if (length(nmcol) > 0) {
  sc  <- colData(cds)[[nmcol[1]]]
  thr <- stats::quantile(sc, 0.98, na.rm = TRUE)
  root_cells <- names(sc)[which(sc >= thr)]
}

## ===== root cell 자동 지정 (수정된 폴백) =====
# 2) 폴백: 가장 큰 클러스터의 중심(UMAP) 근처 셀 일부
if (length(root_cells) < 10) {
  msg("NaiveT-based root insufficient → fallback to largest cluster medoid.")
  # clusters()는 벡터를 반환. reduction_method 명시!
  clus <- monocle3::clusters(cds, reduction_method = "UMAP")  # factor(named by cell ids)
  if (is.null(names(clus))) names(clus) <- colnames(cds)      # 안전장치

  big_c <- names(sort(table(clus), decreasing = TRUE))[1]
  sel   <- names(clus)[clus == big_c]

  umap <- reducedDims(cds)$UMAP[sel, , drop = FALSE]
  center <- colMeans(umap)
  d <- sqrt(rowSums((umap - matrix(center, nrow(umap), ncol(umap), byrow = TRUE))^2))

  keep_n <- min(200, nrow(umap))      # 최대 200개 루트셀
  root_cells <- sel[order(d)][seq_len(keep_n)]
}

msg("order_cells ... (root n=%d)", length(root_cells))
cds <- order_cells(cds, root_cells = root_cells)

## ===== 저장 & 플롯 =====
msg("Saving CDS and plots ...")
saveRDS(cds, file.path(out_dir, "T_monocle3_cds.rds"))

# UMAP + graph + cluster
p1 <- plot_cells(
  cds,
  color_cells_by = if (!is.null(clus_col)) "seurat_clusters" else "cluster",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
) + ggtitle("T cells: UMAP + principal graph (clusters)")
ggsave(file.path(out_dir, "T_umap_graph_clusters.png"), p1, width = 8, height = 6, dpi = 300)

# Pseudotime
p2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
) + ggtitle("T cells: pseudotime")
ggsave(file.path(out_dir, "T_pseudotime.png"), p2, width = 8, height = 6, dpi = 300)

# Group(있으면)
if (!is.null(grp_col)) {
  p3 <- plot_cells(
    cds,
    color_cells_by = "group",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3
  ) + ggtitle("T cells: group")
  ggsave(file.path(out_dir, "T_umap_group.png"), p3, width = 8, height = 6, dpi = 300)
}

msg("Done → %s", normalizePath(out_dir))



pt <- monocle3::pseudotime(cds)
grp <- colData(cds)$group
pdf('Tmonocle_boxplot.pdf')
boxplot(pt ~ grp, main="Pseudotime by group")
dev.off()

br <- monocle3::principal_graph(cds)$UMAP   # 그래프
# 각 셀의 가장 가까운 edge/leaf를 받아오는 헬퍼
bp <- monocle3::branch_nodes(cds)
leaf_cells <- monocle3::choose_cells(cds)  # 또는 monocle3::partitions/ clusters 사용
# 간단 대안: clusters(cds, "UMAP")로 큰 클러스터별 그룹 비율 비교
tab <- prop.table(table(clusters(cds, "UMAP"), colData(cds)$group), 1)
print(round(tab,3))