# cluster 세분화
obj <- readRDS('SCT&Harmony_mkrs.rds')

# 2) 해상도 올려서 클러스터링
obj <- FindClusters(
  obj,
  graph.name ='integrated_snn',
  resolution = 1.0,   # ← 여기 값을 조절 (기본 0.5 → 1.0, 1.2 등)
  algorithm  = 1
)

# 3) UMAP 다시 그리기
pdf('1.0Clusters.pdf')
DimPlot(
  obj,
  reduction = "umap",
  label     = TRUE,
  pt.size   = 0.5
) +
  ggtitle("Seurat_clusters_1.0")
dev.off()


# 원하는 gene list 발현 확인
library(Seurat)

obj <- readRDS('SCT&Harmony_mkrs.rds')

# 1) 확인하고 싶은 유전자 리스트
genes_to_check <- c(
  "IL1B","EREG","NLRP3","LRMDA","VCAN",
  "PLXDC2","PID1","SLC8A1","ACSL1","CXCL8",
  "AZIN1-AS1","NAMPT","RBM47","MCTP1","RNF144B"
)

# 2) obj 에 들어 있는 전체 유전자 이름 벡터 추출
all_genes <- rownames(obj)

# 3) 존재하는 유전자와 빠진 유전자 분리
present_genes <- genes_to_check[genes_to_check %in% all_genes]
missing_genes <- setdiff(genes_to_check, present_genes)

# 4) 결과 출력
cat("### 존재하는 유전자 (", length(present_genes), "개) ###\n")
print(present_genes)

cat("\n### 찾을 수 없는 유전자 (", length(missing_genes), "개) ###\n")
print(missing_genes)

pdf("monocyte marker.pdf", width = 8, height = 5)
DotPlot(
  object    = obj,
  assay     = "SCT",
  features  = genes_to_check,
  group.by  = "cluster_r1.0"   # 또는 "seurat_clusters"
) +
  scale_size_area() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  )
dev.off()



library(dplyr)

# 1) 클러스터별로 gene 열을 묶어서 데이터프레임 생성
cluster_markers_r1.0 <- top20_r1.0 %>%
  group_by(cluster) %>%
  summarise(
    genes = paste(gene, collapse = ","),  # 내부 구분자: 콤마
    .groups = "drop"
  )

# 2) CSV 파일로 저장
write.csv(
  cluster_markers_r1.0,
  file      = "cluster_r1.0_top20_markers.csv",
  row.names = FALSE,
  quote     = FALSE
)

# 3) 탭-구분 텍스트(.txt)로 저장 (내부 구분자: \t)
cluster_markers_r1.0 <- cluster_markers_r1.0 %>%
  mutate(genes = gsub(",", "\t", genes))

write.table(
  cluster_markers_r1.0,
  file      = "cluster_r1.0_top20_markers.txt",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE,
  col.names = TRUE
)





monocyte_markers <- list(
  Classical    = c("CD14","S100A8","S100A9","LYZ","CCR2","FCN1"),
  Intermediate = c("HLA-DRA","CD74","AIF1"),
  NonClassical = c("CX3CR1","MS4A7","TNFSF10","LILRB2"),
  Inflammatory = c("IL1B","NLRP3","CXCL8","S100A12","TNF")
)

# 예: DotPlot으로 비교
pdf('monocyte.pdf')
DotPlot(obj, features = monocyte_markers, assay = "SCT", group.by = "cluster_r1.0") +
  scale_size_area() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()


Mono_markers <- c(
  "CD14","CCR2","CCR5", "SELL",
  "CD68","ITGAX",
  "FCGR3A","CX3CR1", "HLA-DRA"
)
pdf('mono.pdf')
DotPlot(obj, features = Mono_markers, assay = "SCT", group.by = "cluster_r1.0") +
  scale_size_area() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()


markers_B <- list(
  Naive        = c("IGHD","IGHM","TCL1A","FCER2","CXCR4"),
  Intermediate = c("MME","CD24","CD38","GNG7"),
  Memory       = c("CD27","TNFRSF13B","CD80","CD86","GPR183")
)

# 4) DotPlot으로 비교
pdf("Bcell_markers.pdf", width=10, height=7)
DotPlot(
  obj,
  features = markers_B,
  group.by = 'cluster_r1.0',
  assay    = "SCT",
  dot.scale = 8
) +
  scale_size_area() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  )
dev.off()