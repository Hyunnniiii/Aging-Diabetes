library(Seurat)
library(scmap)
library(SingleCellExperiment)
library(Matrix)

# Reference counts matrix 
load('scmap/PBMC_ref.RData')  

# celltype 정보 추출 (sce_ref가 이미 메모리에 있다고 가정)
celltypes <- as.character(colData(sce_ref)$celltype)
names(celltypes) <- rownames(colData(sce_ref))   # cell barcode (colData rownames)

# cell barcode 순서 맞추기 (중요!)
celltypes_matched <- celltypes[colnames(sm)]

sce_ref_new <- SingleCellExperiment(
  assays = list(counts = sm),    # sm이 진짜 raw counts matrix!
  colData = DataFrame(cell_type1 = celltypes_matched)
)

# celltype 분포 확인
table(colData(sce_ref_new)$cell_type1)

# counts 값 점검
sum(assay(sce_ref_new, "counts"))
Matrix::nnzero(assay(sce_ref_new, "counts"))


# 2. query 파일 불러오기
query_path <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT.rds"
query <- readRDS(query_path)

# SingleCellExperiment 변환
sce_query <- as.SingleCellExperiment(query, assay='data')

# 3. feature_symbol 필수 컬럼 지정
rowData(sce_ref)$feature_symbol <- rownames(sce_ref)
rowData(sce_query)$feature_symbol <- rownames(sce_query)

# 4. celltype label 컬럼 명확히 할당(필요시)
# colData(sce_ref)$celltype <- colData(sce_ref)$celltype.l2

# 5. reference feature selection & indexing
sce_ref <- selectFeatures(sce_ref, suppress_plot = TRUE)
sce_ref <- indexCluster(sce_ref, cluster_col = "celltype")

# 6. query 예측/label transfer
scmap_results <- scmapCluster(
  projection = sce_query,
  index_list = list(ref = metadata(sce_ref)$scmap_cluster_index)
)

# 7. 예측 결과 query(Seurat) meta에 저장
query$scmap_predicted <- scmap_results$combined_labs

# 8. UMAP 플롯 (PDF 저장)
pdf('umap_scmap.pdf', width=8, height=6)
print(
  DimPlot(query, group.by = "scmap_predicted", reduction = "umap", label = TRUE) +
    ggtitle("scmap annotation")
)
dev.off()