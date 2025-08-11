library(Seurat)
library(ggplot2)
library(plyr)

obj <- readRDS('SCT&Harmony_mkrs.rds')

# 
new_cluster_names <- c(
  "CD14 Mono",       # 0
  "CD4 TCM",         # 1
  "CD4 TCM",         # 2
  "CD8 TEM",         # 3
  "NK",              # 4
  "naive B",         # 5
  "CD8 naive",       # 6
  "CD14 Mono",       # 7
  "CD14 Mono",       # 8
  "NK",              # 9
  "naive B",         # 10
  "CD8 TEM",         # 11
  "intermediate B",  # 12
  "intermediate B",  # 13
  "CD14 Mono",       # 14
  "Treg",            # 15
  "CD16 Mono",       # 16
  "MAIT",            # 17
  "CD4 TCM",         # 18
  "CD4 TCM",         # 19
  "CD14 Mono",       # 20
  "CD8 TEM",         # 21
  "pDC",             # 22
  "plasmablast",     # 23
  "cDC2",            # 24
  "remove"             # 25
)

obj$celltype <- mapvalues(
  x    = as.character(obj$cluster_r1.0), 
  from = as.character(0:25), 
  to   = new_cluster_names
)

# 1) celltype 벡터 가져오기
ct <- obj$celltype

# 2) "remove"가 아닌 셀만 남길 논리 벡터 생성
keep_cells <- ct != "remove"

# 3) [,] 인덱싱으로 객체 서브셋
obj_filtered <- obj[, keep_cells]

# 4. UMAP 등 시각화 (celltype 기준)
pdf("SH_final_annotation1.pdf", width=8, height=6)
DimPlot(obj_filtered, reduction = 'umap_r1.0', group.by = "celltype", 
        label = TRUE, repel = TRUE) +
theme(
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  ) +
  coord_fixed(ratio = 1)
dev.off()

saveRDS(obj_filtered, file = "SCT&Harmony_annot.rds")


# heatmap, dotplot으로 canonical marker 발현 확인
library(Seurat)
library(ggplot2)
library(patchwork)

# 1) 객체 불러오기
seurat_obj <- readRDS("SCT&Harmony_annot.rds")

# 4) 마커 리스트 준비 (각 cell type마다 ≥2개)
marker_list <- list(
  `intermediate Mono`   = c("FCGR3A","CD14"),
  `CD4 TCM`             = c("LEF1","IL7R","CCR7"),
  `NK`                  = c("NCAM1","NKG7","GNLY"),
  `classical Mono`      = c("CD14","LYZ"),
  `naive CD4`           = c("CCR7","SELL","TCF7"),
  `CD8 TEM`             = c("GZMK","S100A4"),
  `Intermediate B`      = c("CD27","TNFRSF13B"),
  `naive CD8`           = c("SELL","LEF1","TCF7"),
  `naive B`             = c("TCL1A","IGHD"),
  `γδ-T`                = c("TRDC","TRGC1"),
  `Memory B`            = c("CD27","IGHG1"),
  `non-classical Mono`  = c("FCGR3A","MS4A7"),
  `Treg`                = c("FOXP3","IL2RA","CTLA4"),
  `pDC`                 = c("LILRA4","CLEC4C","IL3RA"),
  `Plasma cell`         = c("SDC1","MZB1","XBP1"),
  `cDC`                 = c("CLEC9A","XCR1","IRF8")
)
all_markers <- unique(unlist(marker_list))

# 3) RNA assay normalize
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", verbose = TRUE)

# 4) 메타데이터 셋업
Idents(seurat_obj)       <- seurat_obj$celltype    # 'celltype' 컬럼이 annotation
DefaultAssay(seurat_obj) <- "RNA"

# 5) DotPlot
dp <- DotPlot(seurat_obj, assay = "RNA", features = all_markers, 
              group.by = 'celltype') +
  RotatedAxis() +
  labs(title = "DotPlot of canonical markers") +
  theme(plot.title = element_text(hjust = 0.5))

# 6) Heatmap 준비
seurat_obj <- ScaleData(seurat_obj, assay="RNA", features = all_markers, verbose = TRUE)
ht <- DoHeatmap(
  seurat_obj,
  assay    = "RNA",
  features = all_markers,
  group.by = "celltype",
  size     = 3
) +
  scale_fill_gradientn(colors = c("navy","white","firebrick")) +
  labs(title = "Heatmap of canonical markers") +
  theme(plot.title = element_text(hjust = 0.5))

# 7) PDF 출력
pdf("finalcluster_Dot_and_Heatmap.pdf", width = 12, height = 10)
dp + ht + plot_layout(ncol = 1, heights = c(1,2))
dev.off()


# data상 celltype별 marker gene expression 확인
library()