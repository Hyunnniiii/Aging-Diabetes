# ============================================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | Batch 간 차이가 심하게 나타나는 cluster0의 샘플/그룹/배치별 셀 수를 집계, QC 통계 요약 및 Boxplot
## 상위 마커 추출해 FeaturePlot으로 확인, 그 안에서 percent.mt=0 비율과 배치별 분포·셀 수 계산·출력
# ============================================================


###### cluster0 QC 확인
library(Seurat)
seu <- readRDS('sct_annotated_1.rds')

# cluster 0만 subset
cluster0 <- subset(seu, idents = "0")
# 샘플별 cell 개수
table(cluster0$orig.ident)
# 그룹별(celltype/condition 등) cell 개수
table(cluster0$group)
table(cluster0$condition)

# cluster별 nCount_RNA, nFeature_RNA, percent.mt의 요약 통계
library(dplyr)

qc_summary <- seu@meta.data %>%
  group_by(cluster = seu$seurat_clusters) %>%
  summarise(
    cell_n = n(),
    nCount_RNA_mean = mean(nCount_RNA),
    nCount_RNA_median = median(nCount_RNA),
    nFeature_RNA_mean = mean(nFeature_RNA),
    nFeature_RNA_median = median(nFeature_RNA),
    percent.mt_mean = mean(percent.mt),
    percent.mt_median = median(percent.mt)
  )

print(qc_summary)



###### Boxplot
library(ggplot2)
# nCount_RNA boxplot by cluster
pdf('cluster0_qcplot.pdf')
ggplot(seu@meta.data, aes(x = seurat_clusters, y = nCount_RNA)) +
  geom_boxplot() +
  xlab("Cluster") + ylab("nCount_RNA") + theme_bw()

# percent.mt boxplot by cluster
ggplot(seu@meta.data, aes(x = seurat_clusters, y = percent.mt)) +
  geom_boxplot() +
  xlab("Cluster") + ylab("percent.mt") + theme_bw()
dev.off()

# 이미 Seurat에서 markers를 계산한 경우
markers_0 <- subset(seu@misc$markers, cluster == "0")

# avg_log2FC가 큰 순서대로 TOP 10 gene 추출
top10_marker_0 <- markers_0[order(-markers_0$avg_log2FC), ][1:10, ]
print(top10_marker_0)

pdf('cluster0_marker.pdf', width=7, height=7)  # 반드시 width=height
genes <- c("RNF144B", "PID1", 'LRMDA', 'TNFAIP2', 'RBM47',
           'CSF3R', 'SLC11A1', 'WDFY3', 'SLC8A1', 'PLXDC2') # top10 marker
for(g in genes) {
  print(FeaturePlot(seu, features = g) + coord_fixed())
}
dev.off()



###### cluster0 cell 수 확인
# meta.data 로딩
md <- seu@meta.data

# cluster 0 전체 cell 수
total_c0 <- sum(md$cluster == "0")

# cluster 0 중 percent.mt == 0인 cell 수
zero_mt_c0 <- sum(md$cluster == "0" & md$percent.mt == 0)

# 비율 계산
ratio_zero_mt <- zero_mt_c0 / total_c0 * 100

# 결과 출력
cat("Cluster 0 전체 cell 수:", total_c0, "\n")
cat("percent.mt == 0인 cell 수:", zero_mt_c0, "\n")
cat("비율:", round(ratio_zero_mt, 2), "%\n")



###### batch별 percent.mt 비율 확인
library(ggplot2)

# cluster 0 cell subset
cluster0_cells <- subset(seu, idents = "0")

# percent.mt 분포를 batch별로 boxplot으로 시각화
pdf('cluster0_mt_by_batch.pdf')
ggplot(cluster0_cells@meta.data, aes(x = orig.ident, y = percent.mt)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  xlab("Batch (orig.ident)") +
  ylab("Percent Mitochondrial Genes") +
  ggtitle("Cluster 0 Percent.mt by Batch")
dev.off()



###### cluster 0 cell 중 batch별 cell 수 구하기
cluster0_batch_table <- table(md$orig.ident[md$cluster == "0"])

# 결과를 보기 좋게 데이터프레임으로 변환
cluster0_batch_df <- data.frame(
  Batch = names(cluster0_batch_table),
  Cluster0_Cells = as.numeric(cluster0_batch_table)
)

print(cluster0_batch_df)