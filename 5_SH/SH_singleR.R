# ============================================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | Harmony까지 적용한 data에 대한 singleR 분석 및 umap으로 결과 확인
# ============================================================


library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(tidyverse)

# 1. 불러오기
seu <- readRDS("SCT&Harmony.rds")

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
options(ggrepel.max.overlaps = Inf)

# PDF 파일로 출력 시작
pdf("SH_umap_SingleR.pdf", width=8, height=6)  # 원하는 파일명, 크기로 지정

# 원하는 플롯 실행
DimPlot(seu, group.by = "SingleR", label = TRUE, repel = TRUE) + ggtitle("SingleR annotation")

# 저장 끝내기
dev.off()

# singleR 결과 파일로 저장
saveRDS(seu, file = "SCT&Harmony.rds")