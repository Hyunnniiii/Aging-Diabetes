library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(tidyverse)

# 1. 불러오기
seu <- readRDS("integrated_SCT_with_markers.rds")

# 2. raw count 추출
counts <- seu@assays$RNA@layers$counts
rownames(counts) <- rownames(seu@assays$RNA)


# 3. SCE 변환
sce <- SingleCellExperiment(list(logcounts = seu@assays$SCT@data))

# 4. Reference 준비
ref <- celldex::HumanPrimaryCellAtlasData()

# 5. SingleR 실행
singler.res <- SingleR(test = sce, ref = ref, labels = ref$label.main, 
                      de.method = 'wilcox')

singler.res %>% head()


# singler.res@pruned.labels 또는 singler.res$pruned.labels (DataFrame일 경우)
seu$SingleR <- singler.res$pruned.labels

# PDF 파일로 출력 시작
pdf("umap_SingleR.pdf", width=8, height=6)  # 원하는 파일명, 크기로 지정

# 원하는 플롯 실행
DimPlot(seu, group.by = "SingleR", label = TRUE) + ggtitle("SingleR annotation")

# 저장 끝내기
dev.off()

# singleR 결과 파일로 저장
# saveRDS(seu, "seurat_with_singler.rds") -> celltype 등 annotation 포함 Seurat 객체
# saveRDS(singler.res, "singler.res.rds") -> SingleR 분석 전체 결과


library(Seurat)
library(pheatmap)

# 1. 불러오기
seu <- readRDS("seurat_with_singler.rds")
singler.res <- readRDS("singler.res.rds")   # 저장했을 경우만!

# 2. cluster id와 score matrix 추출
cluster_ids <- seu$seurat_clusters
scores <- singler.res$scores

# cluster별로 평균 구하기
mean_scores <- aggregate(scores, by = list(cluster = cluster_ids), FUN = mean)
rownames(mean_scores) <- mean_scores$cluster
mean_scores <- mean_scores[,-1]   # cluster column 제거 (matrix만 남기기)
mean_scores <- t(mean_scores)     # 열(=cell type)이 더 많은 경우 보기 좋게 전치도 가능


# PDF 파일로 출력 시작
pdf("singleR_Cluster_Heatmap.pdf", width=8, height=6)

# heatmap 그리기
pheatmap(mean_scores,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 7,
         fontsize_col = 8,
         border_color = "grey60",
         main = "SingleR Cluster-wise Score Heatmap"
)
dev.off()