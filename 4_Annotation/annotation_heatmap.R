# ============================================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | annotation된 celltype별 marker gene expression 정규화 후 heatmap으로 확인
# ============================================================


library(Seurat)
library(dplyr)
library(ggplot2)

# 1) 객체 불러오기
seurat_obj <- readRDS("sct_annotated_1.rds")

# 2) 클러스터별 top10 마커 유전자 추출
markers_df <- seurat_obj@misc$markers
top10 <- markers_df %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  pull(gene) %>% 
  unique()

# 3) 관심 유전자만 subset
sub <- subset(seurat_obj, features = top10)

# 4) RNA assay 로 전환 후 정규화·스케일
DefaultAssay(sub) <- "RNA"

sub <- NormalizeData(
  sub,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = 1e4
)
sub <- ScaleData(
  sub,
  assay   = "RNA",
  verbose = FALSE
)

# 5) Idents() 를 클러스터로 설정
Idents(sub) <- sub$celltype

# 6) heatmap 그리기
pdf('cluster_heatmap.pdf', width=10, height=25)
DoHeatmap(
  object   = sub,
  features = top10,
  assay    = "RNA",
  slot     = "scale.data",
  group.by = "celltype"
) +
  scale_fill_gradientn(colors = c("navy","white","firebrick")) +
  theme(axis.text.y = element_text(size = 6))

dev.off()