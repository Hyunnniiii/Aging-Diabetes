library(Seurat)
library(ggplot2)
library(dplyr)

# 1) Harmony 처리된 오브젝트 불러오기
obj <- readRDS("SCT&Harmony.rds")

# 2) SCT assay 의 counts 만 RNA assay 로 복사
rna_assay <- CreateAssayObject(
  counts = GetAssayData(obj, assay = "SCT", layer = "counts")
)

# 3) 방금 만든 RNA assay 의 data 슬롯을 SCT assay 의 data 로 채워주기
rna_assay <- SetAssayData(
  object = rna_assay,
  layer   = "data",
  new.data = GetAssayData(obj, assay = "SCT", layer = "data")
)

# 4) 완성된 RNA assay 를 obj 에 붙여넣기
obj@assays[["RNA"]] <- rna_assay

# 확인: 이제 RNA assay 에 counts/data 가 들어있어야 함
print(obj[["RNA"]])

# 5) PrepSCTFindMarkers → FindAllMarkers
DefaultAssay(obj) <- "SCT"

obj <- PrepSCTFindMarkers(
  object = obj,
  assay  = "SCT"
)

markers <- FindAllMarkers(
  object          = obj,
  assay           = "SCT",
  slot            = "data",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.2
)

head(markers)

# 3) 유의미한 것만 필터
markers <- markers %>% 
  filter(p_val_adj < 0.05)

# 4) 클러스터별 Top20 뽑기 (avg_log2FC 기준)
top20 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()

# 6) metadata 에 marker 결과 저장 (@misc 슬롯 사용)
obj@misc$markers <- markers
obj@misc$top20   <- top20

saveRDS(obj, file = "SCT&Harmony_mkrs.rds")




library(Seurat)
library(cowplot)  # grid 정렬에 유용
library(ggplot2)

# 1) 통합 객체 불러오기
obj <- readRDS("SCT&Harmony_mkrs.rds")


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



library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# 1) 객체 불러오기
obj <- readRDS("SCT&Harmony_annot.rds")

Idents(obj) <- obj$cluster_r1.0

# 2) 이미 저장된 top20에서 클러스터별로 top10만 추려서 features 선정
top10 <- obj@misc$top20_r1.0 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()
features <- unique(top10$gene)

DefaultAssay(obj) <- "SCT"

# 3) ScaleData 로 상대 비교용 스케일링
obj <- ScaleData(
  object   = obj,
  assay    = "SCT",
  features = features,
  verbose  = FALSE
)


# 3) 히트맵 (화면 출력)
t10 <- DoHeatmap(
  object   = obj,
  features = features,
  assay    = "SCT",
  group.by = "cluster_r1.0",
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

pdf("sh_cluster_top10_heatmap.pdf", width = 15, height = 20)
print(t10)
dev.off()



library(dplyr)

# 1) 클러스터별로 gene 열을 콤마로 결합한 데이터프레임 생성
cluster_markers <- top20 %>%
  group_by(cluster) %>%
  summarise(
    genes = paste(gene, collapse = ","),
    .groups = "drop"
  )

# 결과 확인
print(cluster_markers)
# A tibble: N x 2
#   cluster genes
#   <fct>   <chr>
# 1 0       GENE1,GENE2,GENE3,...
# 2 1       GENE4,GENE5,GENE6,...
# ...

# 2) CSV 파일로 저장 (클러스터 번호, 유전자 리스트 컬럼)
write.csv(
  cluster_markers,
  file = "sh_cluster_top20_markers.csv",
  row.names = FALSE,
  quote = FALSE
)

cluster_markers <- top20 %>%
  group_by(cluster) %>%
  summarise(
    genes = paste(gene, collapse = "\t"),   # 쉼표 → 세미콜론
    .groups = "drop"
  )

# 3) 탭 구분 텍스트(.txt)로 저장하고 싶다면:
write.table(
  cluster_markers,
  file = "cluster_top20_markers.txt",
  sep = ",",
  row.names = FALSE,
  quote = FALSE
)


genes_of_interest <- c("IGHD","IGHM","TCL1A","FCER2","CXCR4",
  "MME","CD24","CD38","GNG7",
  "CD27","TNFRSF13B","CD80","CD86","GPR183")

pdf('B cell marker.pdf', width=20, height=20)
VlnPlot(
  object = obj,
  features = genes_of_interest,
  assay = "SCT",
  slot = "data",
  group.by = "cluster_r1.0",
  pt.size = 0
) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()