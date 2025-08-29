# ============================================================
## Date | 2025-08-28 | Yoon JiHyun
## Anchor 기반 integration 결과 파일로 이후 분석 진행

## Description | 최종 annotation 결과 매핑, cell 수 적은 cluster는 삭제, 새로운 rds 파일로 저장
## 결과 umap으로 확인 후 canonical marker expression을 dotplot & heatmap으로 최종 celltype별 확인
# ============================================================


seurat_obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/integrated_SCT_with_markers.rds')

# 
new_cluster_names <- c(
  "intermediate Mono",           # 0
  "CD4 TCM",           # 1
  "NK",          # 2
  "classical Mono",# 3
  "naive CD4",           # 4
  "CD8 TEM",           # 5
  "Intermediate B",# 6
  "naive CD8",           # 7
  "naive B",     # 8
  "gamma_delta T",           # 9
  "Memory B",    # 10
  "non-classical Mono",# 11
  "Treg",           # 12 
  "pDC",         # 13
  "Plasma cell", # 14
  "cDC",         # 15
  "Remove"       # 16
)

seurat_obj$celltype <- plyr::mapvalues(
  seurat_obj$seurat_clusters,
  from = as.character(0:16),
  to   = new_cluster_names
)

# 3. "Remove" celltype 제외
seurat_obj_filtered <- subset(seurat_obj, subset = celltype != "Remove")

# 4. UMAP 등 시각화 (celltype 기준)
pdf("final_annotation1.pdf", width=8, height=6)
DimPlot(seurat_obj_filtered, group.by = "celltype", label = TRUE, repel = TRUE)
dev.off()

saveRDS(seurat_obj_filtered, file = "sct_annotated_1.rds")


# heatmap, dotplot으로 canonical marker 발현 확인
library(Seurat)
library(ggplot2)
library(patchwork)

# 1) 객체 불러오기
seurat_obj <- readRDS("sct_annotated_final.rds")

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
pdf("CanonicalMarkers_Dot_and_Heatmap.pdf", width = 12, height = 10)
dp + ht + plot_layout(ncol = 1, heights = c(1,2))
dev.off()


