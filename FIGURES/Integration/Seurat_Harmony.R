library(Seurat)
library(harmony)
library(ggplot2)

# 1) QC 완료된 객체 불러오기
qc_files <- c("QCRun1.rds", "QCRun2.rds", "QCRun3.rds")
projects <- c("Run1", "Run2", "Run3")
seu_list <- setNames(lapply(qc_files, readRDS), projects)

# 2) SCTransform (각 샘플별 정규화)
seu_sct_list <- lapply(seu_list, SCTransform, verbose = TRUE)

# 3) 병합(Merge) — 하나의 Seurat 객체로 합치기
seu_merged <- merge(
  x = seu_sct_list[[1]],           # Run1 객체
  y = seu_sct_list[-1],            # Run2, Run3 객체 리스트
  add.cell.ids = projects,         # Run1, Run2, Run3
  project = "MergedHarmony"        # (옵션) merged object 에 붙일 project 이름
)

# 3.1) 기본 assay를 SCT로 맞추기
DefaultAssay(seu_merged) <- "SCT"

# 3.2) 발현 변동 유전자(HVG) 2,000개 뽑기
seu_merged <- FindVariableFeatures(
  seu_merged,
  selection.method = "vst",
  nfeatures        = 3000,
  verbose          = TRUE
)

# 4) PCA
seu_merged <- RunPCA(seu_merged, verbose = TRUE)

# 5) Harmony batch correction
#    orig.ident 컬럼(Run1/Run2/Run3)에 따라 보정
seu_merged <- RunHarmony(
  object        = seu_merged,
  group.by.vars = "orig.ident",
  assay.use     = "SCT",
  verbose       = TRUE
)

# 6) UMAP/클러스터링 (Harmony embedding 사용)
seu_merged <- RunUMAP(seu_merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu_merged <- FindNeighbors(seu_merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu_merged <- FindClusters(seu_merged, resolution = 0.5, verbose = TRUE)

# 7) UMAP 결과 PDF로 저장
pdf("Harmony_Integration_UMAP.pdf", width = 8, height = 6)
print(
  DimPlot(seu_merged, reduction = "umap", group.by = "orig.ident") +
    ggtitle("Harmony-corrected UMAP")
)
dev.off()

saveRDS(seu_merged, file = "Integrateded_Harmony.rds")


library(Seurat)
library(ggplot2)

# 통합된 Harmony 결과 불러오기
seu_merged <- readRDS("Integrateded_Harmony.rds")

# PDF 장치 열기 (가로 12, 세로 6)
pdf("Harmony_split_UMAP.pdf", width = 12, height = 6)

# 1) orig.ident 기준으로 3개 패싯 UMAP
DimPlot(
  object    = seu_merged,
  reduction = "umap",
  group.by  = "orig.ident",  # 색상 기준: Run1/Run2/Run3
  split.by  = "group",       # 패싯 기준: Diabetes_Old / Healthy_Old / Healthy_Young
  ncol      = 3,
  label     = FALSE
) + ggtitle("orig.ident (Harmony)")

# 2) patient 기준으로 3개 패싯 UMAP
DimPlot(
  object    = seu_merged,
  reduction = "umap",
  group.by  = "patient",     # 색상 기준: DO1…HY5
  split.by  = "group",
  ncol      = 3,
  label     = FALSE
) + ggtitle("patient (Harmony)")

# 장치 닫기 — 파일 저장 완료
dev.off()
