# ============================================================
## Date | 2025-08-28 | Yoon JiHyun
## Anchor 기반 integration 결과 파일로 이후 분석 진행

## Description | 각 cluster별 marker gene 추출 후 top10, top20 gene heatmap으로 발현 확인
## marker 정보까지 포함해서 새로운 rds 파일로 저장
# ============================================================


library(Seurat)
library(ggplot2)
library(dplyr)

obj <- readRDS("integrated_SCT_markergene.rds")
obj <- PrepSCTFindMarkers(object =obj, assay='SCT')

# marker 찾기
markers <- FindAllMarkers(
  object         = obj,
  assay          = "SCT",
  slot           = "data",
  only.pos       = TRUE,
  min.pct        = 0.25,
  logfc.threshold = 0.2
)

# 유의미한 것만 필터
markers <- markers %>% 
  filter(p_val_adj < 0.05)

# 클러스터별 Top20 뽑기 (avg_log2FC 기준)
top20 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

# metadata 에 marker 결과 저장 (@misc 슬롯 사용)
obj@misc$markers <- markers
obj@misc$top20   <- top20

saveRDS(obj, file = "integrated_SCT_with_markers.rds")



######### top20 gene Heatmap ########
library(Seurat)
library(cowplot)  # grid 정렬에 유용
library(ggplot2)

# 1) 통합 객체 불러오기
obj <- readRDS("integrated_SCT_with_markers.rds")

# 2) Top20 gene 리스트 추출
features <- unique(obj@misc$top20$gene)

# 3) 전체 클러스터 Top20 히트맵
p_all <- DoHeatmap(
  object   = obj,
  features = features,
  assay     = "SCT",
  slot      = "data",
  size      = 3  # 폰트 크기
) + 
  scale_fill_gradientn(
    colours = c("navy", "white", "firebrick"),
    limits  = c(-2, 2)
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

# 4) PDF로 저장
pdf("cluster_top20_heatmap.pdf", width = 10, height = 30)
print(p_all)
dev.off()


######### top10 gene Heatmap ########
# 1) 이미 저장된 top20에서 클러스터별로 top10만 추려서 features 선정
top10 <- obj@misc$top20 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

features <- unique(top10$gene)

# 2) 히트맵 (화면 출력)
t10 <- DoHeatmap(
  object   = obj,
  features = features,
  assay    = "SCT",
  slot     = "scale.data",
  size     = 3
) +
  scale_fill_gradientn(
    colours = c("navy", "white", "firebrick"),
    limits = c(-2, 2)
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

pdf("cluster_top10_heatmap.pdf", width = 20, height = 30)
print(t10)
dev.off()
library(dplyr)
top20 <- read.csv("top20_markers_by_cluster.csv")