# ============================================================
## Date | 2025-08-28 | Yoon JiHyun
## Run1 / Run2-1 / Run2-2 sample 존재, 코드에서는 순서대로 Run1/2/3으로 지정
## batch 보정 결과 비교를 위해 anchor 기반 integration과 Harmony 모두 시도
       
## Description | QC 완료된 data에 SCTransform, FindIntegrationAnchors, IntegrateData 적용
## PCA, clustering 후 umap으로 분포 확인
## 결과 새로운 rds 파일로 저장
# ============================================================


library(Seurat)
library(dplyr)

# 1) QC 완료된 객체 불러오기
qc_files <- c("QCRun1.rds", "QCRun2.rds", "QCRun3.rds")
projects <- c("Run1", "Run2", "Run3")

seu_list <- vector("list", length(qc_files))
names(seu_list) <- projects
for (i in seq_along(qc_files)) {
  seu_list[[i]] <- readRDS(qc_files[i])
}

# 2) SCTransform 정규화
seu_sct <- lapply(seu_list, function(x) {
  SCTransform(x, verbose = TRUE)
})

# 3) 통합용 feature 뽑기 & 준비
features_int <- SelectIntegrationFeatures(object.list = seu_sct, nfeatures = 3000)
seu_sct <- PrepSCTIntegration(object.list     = seu_sct,
                              anchor.features  = features_int,
                              verbose          = TRUE)

# 4) Anchor 찾기
anchors <- FindIntegrationAnchors(object.list            = seu_sct,
                                  normalization.method   = "SCT",
                                  anchor.features        = features_int,
                                  verbose                = TRUE)

# 5) 데이터 통합 (SCT)
integrated <- IntegrateData(anchorset            = anchors,
                            normalization.method = "SCT",
                            verbose              = TRUE)

# 6) 차원 축소 & 클러스터링
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = TRUE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = TRUE)
integrated <- FindClusters(integrated, resolution = 0.5, verbose = TRUE)


# 7) 결과 확인
library(ggplot2)
pdf("Integrated_SCT_UMAP.pdf", width = 8, height = 6)
print(
  DimPlot(integrated, reduction = "umap", group.by = "orig.ident") +
    ggtitle("Integrated UMAP (SCT)")
)
dev.off()
# 필요하면 저장

saveRDS(integrated, "integrated_SCT.rds")


library(Seurat)
library(ggplot2)

# 1) 통합된 SCT 결과 불러오기
integrated <- readRDS("integrated_SCT.rds")

# 2) PDF 장치 열기 (가로 12, 세로 6 인치)
pdf("SCT_split_UMAP.pdf", width = 12, height = 6)

# 3) orig.ident 기준으로 3개 패싯 UMAP
DimPlot(
  object    = integrated,
  reduction = "umap",
  group.by  = "orig.ident",
  split.by  = "group",   # Diabetes_Old / Healthy_Old / Healthy_Young
  ncol      = 3,
  label     = FALSE
) + ggtitle("orig.ident")

# 4) patient 기준으로 3개 패싯 UMAP
DimPlot(
  object    = integrated,
  reduction = "umap",
  group.by  = "patient", # DO1…HY5
  split.by  = "group",
  ncol      = 3,
  label     = FALSE
) + ggtitle("patient")

# 5) 장치 닫기 → 파일 저장 완료
dev.off()
