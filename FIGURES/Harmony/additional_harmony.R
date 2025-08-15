library(harmony)
library(Seurat)

# 1) 이미 통합까지 완료된 오브젝트(예: integrated)
seu <- readRDS('/data/project/diabetes_LYH/tanya/rds/integrated_SCT.rds')

# 2) PCA 다시 수행 (integrated assay 기준)
DefaultAssay(seu) <- "integrated"
seu <- RunPCA(seu, verbose = TRUE)

# 3) Harmony 적용 (batch 변수: orig.ident)
seu <- RunHarmony(
  object         = seu,
  group.by.vars  = "orig.ident",    # 배치 변수명
  reduction.use  = "pca",           # PCA embedding 사용
  dims           = 1:30
)

# 4) Harmony embedding 기반 UMAP, clustering
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

saveRDS(seu, 'SCT&Harmony.rds')

seu <- readRDS("SCT&Harmony.rds")
Assays(seu)  # "RNA" "SCT" "integrated"
Reductions(seu)  # "pca" "harmony" "umap"
head(colnames(seu@meta.data))  # 확인


# 5) 결과 시각화
pdf('after_harmony.pdf')
DimPlot(seu, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Harmony UAMP by Batch")

DimPlot(seu, reduction = "umap", label = TRUE) +
  ggtitle("Harmony Clustering")
dev.off()

pdf('after_harmony_by_group.pdf', width =20)
DimPlot(
  seu, reduction = "umap",     
  split.by = "orig.ident", label = TRUE,             
  pt.size = 0.7              
) + ggtitle("HARMONY SPLIT UMAP by Batch") + coord_fixed()
dev.off()


## 과보정 확인
# 이미 Harmony 이후 객체가 seu 에 들어 있으므로
harmony_obj      <- seu                # Harmony 적용 후
preharm_obj_path <- "/data/project/diabetes_LYH/tanya/rds/integrated_SCT.rds"  # Harmony 직전
preharm_obj      <- readRDS(preharm_obj_path)  # Harmony 적용 전

# 시각적 점검
pdf('UMAP_check.pdf')
DimPlot(preharm_obj, group.by = "orig.ident", reduction = "umap") + ggtitle("Pre-Harmony: by batch")
DimPlot(harmony_obj,  group.by = "orig.ident", reduction = "umap") + ggtitle("Post-Harmony: by batch")

DimPlot(preharm_obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE) + ggtitle("Pre-Harmony: by cluster")
DimPlot(harmony_obj,  group.by = "seurat_clusters", reduction = "umap", label = TRUE) + ggtitle("Post-Harmony: by cluster")
dev.off()



############## 정량적 확인: lisi
## 0. 패키지 로드 ###########
library(lisi)
library(ggplot2)
library(dplyr)


## 1. compute_lisi() 수정 ##
calc_lisi <- function(seu, label, redn = "umap") {
  emb    <- seu@reductions[[redn]]@cell.embeddings[, 1:2]
  meta   <- seu@meta.data[, label, drop = FALSE]
  lisi_df <- compute_lisi(emb, meta, label_col = label)
  return(lisi_df[[label]])
}

## 2. 지표 계산 ############
batch_var <- "orig.ident"
cell_var  <- "seurat_clusters"

# iLISI (배치 혼합도)
ilisi_pre  <- calc_lisi(preharm_obj, batch_var)
ilisi_post <- calc_lisi(harmony_obj, batch_var)

# cLISI (클러스터/세포타입 보존도)
clisi_pre  <- calc_lisi(preharm_obj, cell_var)
clisi_post <- calc_lisi(harmony_obj, cell_var)


## 3. 데이터프레임 생성 ###
# iLISI 전후
df_ilisi <- data.frame(
  Score  = c(ilisi_pre, ilisi_post),
  Metric = "iLISI",
  Status = rep(c("Pre","Post"), each = length(ilisi_pre))
)
# cLISI 전후
df_clisi <- data.frame(
  Score  = c(clisi_pre, clisi_post),
  Metric = "cLISI",
  Status = rep(c("Pre","Post"), each = length(clisi_pre))
)

# 합치기
df_lisi <- bind_rows(df_ilisi, df_clisi)

## 4. 박스플롯 그리기 ######
pdf('iLISI&cLISI.pdf')
ggplot(df_lisi, aes(x = Status, y = Score, fill = Status)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_x_discrete(limits = c("Pre", "Post")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(title = "iLISI & cLISI (Pre vs Post Harmony)",
       x = NULL, y = "LISI Score")
dev.off()


############## 정량적 확인: kBET, Silhouette
library(kBET)   
library(cluster)   

## 1. 파라미터 설정 ###################
batch_var <- "orig.ident"       # 배치 라벨
cell_var  <- "seurat_clusters"  # 클러스터/세포타입 라벨
dims_kbet <- 1:10               # kBET 에 사용할 embedding 차원
test_size <- 500                # kBET 서브샘플 크기

## 2. kBET 계산 함수 ##################
run_kbet <- function(seu, batch_col, redn, dims, test_size=500) {
  # 1) embedding matrix 추출
  mat   <- Embeddings(seu, reduction = redn)[, dims]
  batch <- seu@meta.data[[batch_col]]
  # 2) kBET 실행 (do.pca=FALSE: 이미 embedding 된 상태)
  set.seed(123)
  res <- kBET(mat, batch, do.pca = FALSE)
  # 3) acceptance rate 반환
  return(mean(res$summary$kBET.pvalue > 0.05))
}

# Pre-Harmony (PCA embedding 사용)
if (!"pca" %in% Reductions(preharm_obj)) {
  preharm_obj <- RunPCA(preharm_obj, dims = dims_kbet, verbose = FALSE)
}
kbet_pre  <- run_kbet(preharm_obj, batch_var, redn = "pca",   dims = dims_kbet)

# Post-Harmony (harmony embedding 사용)
kbet_post <- run_kbet(harmony_obj, batch_var, redn = "harmony", dims = dims_kbet)

cat(sprintf("kBET acceptance rate:\n  Pre-Harmony : %.3f\n  Post-Harmony: %.3f\n",
            kbet_pre, kbet_post))

########################################
## 3. Silhouette width 계산 함수 ######
########################################
run_silhouette <- function(seu, label_col, redn, dims) {
  # 1) embedding + 거리행렬 계산
  emb   <- Embeddings(seu, reduction = redn)[, dims]
  distm <- dist(emb)
  # 2) 클러스터 벡터
  clusters <- as.numeric(factor(seu@meta.data[[label_col]]))
  # 3) silhouette 계산 후 평균 폭 반환
  sil <- silhouette(clusters, distm)
  return(summary(sil)$avg.width)
}

# Pre-Harmony silhouette (PCA embedding 사용)
sil_pre  <- run_silhouette(preharm_obj, cell_var, redn = "pca",     dims = dims_kbet)
# Post-Harmony silhouette (Harmony embedding 사용)
sil_post <- run_silhouette(harmony_obj,  cell_var, redn = "harmony", dims = dims_kbet)

cat(sprintf("Silhouette width (avg):\n  Pre-Harmony : %.3f\n  Post-Harmony: %.3f\n",
            sil_pre, sil_post))



# umap resolution 1.0 추가
# 1) 객체 로드
obj <- readRDS("SCT&Harmony_mkrs.rds")

# 2) 원본 UMAP과 클러스터 백업
#  – 원본 UMAP 좌표를 "umap_r0.5" 에 저장
obj[["umap_r0.5"]]      <- obj@reductions$umap_r0.5
#  – 원본 클러스터 아이덴티티(해상도0.5) 를 "cluster_r0.5" 에 저장
obj$cluster_r0.5        <- Idents(obj)

# 3) 새로 해상도 1.0 클러스터링 (기존 그래프 사용)
obj <- FindClusters(
  obj,
  graph.name = "integrated_snn",
  resolution = 1.0,
  algorithm  = 1
)
#  – 해상도 1.0 클러스터 아이덴티티를 "cluster_r1.0" 에 저장
obj$cluster_r1.0        <- Idents(obj)

# 4) 새 UMAP 계산 (원본 umap 건드리지 않고 "umap_r1.0" 에 저장)
obj <- RunUMAP(
  obj,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "umap_r1.0"
)

# 5) 결과 확인
#   – 원본 UMAP: obj@reductions$umap_r0.5
#   – 새 UMAP   : obj@reductions$umap_r1.0
#   – 원본 클러스터: obj$cluster_r0.5
#   – 새 클러스터  : obj$cluster_r1.0

# 6) 필요시 저장
saveRDS(obj, file = "SCT&Harmony_mkrs.rds")
