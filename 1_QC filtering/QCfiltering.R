# ============================================================
## Date | 2025-08-28 | Yoon JiHyun
## Run1 / Run2-1 / Run2-2 sample 존재, 코드에서는 순서대로 Run1/2/3으로 지정
       
## Description | QC-Run1은 plot으로 분포 확인하고 임의로 cut 지정, Run2, Run3은 MAD로 cut off.
## 이후 plot으로 분포 확인하고 새로운 rds 파일로 저장
# ============================================================


library(Seurat)
library(tools)
library(dplyr)
library(ggplot2)
library(patchwork)  

files    <- c("newRun1.rds", "newRun2.rds", "newRun3.rds")
projects <- c("Run1",         "Run2",         "Run3")

# fixed params for Run1
fixed_params <- list(
  min_features = 500,  max_features = 3000,
  min_counts   = 500,  max_counts   = 7000,
  max_mito     = 7
)

# 우리가 보여줄 세 가지 feature
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
res_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f       <- files[i]
  project <- projects[i]
  
  # 1) 원본 RDS 읽기
  obj <- readRDS(f)
  
  # 2) Seurat 객체 생성
  seu <- CreateSeuratObject(
    counts    = obj$counts,
    meta.data = obj$meta,
    project   = project
  )
  
  # 3) 필터링 전 셀 수
  n_before <- ncol(seu)
  
  # 4) percent.mt 재계산
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # 5) QC 기준 결정
  if (project == "Run1") {
    p <- fixed_params
  } else {
    feat <- seu$nFeature_RNA; cnts <- seu$nCount_RNA; mito <- seu$percent.mt
    med_f <- median(feat); mad_f <- mad(feat)
    med_c <- median(cnts); mad_c <- mad(cnts)
    med_m <- median(mito); mad_m <- mad(mito)
    p <- list(
      min_features = med_f - 3*mad_f,
      max_features = med_f + 3*mad_f,
      min_counts   = med_c - 3*mad_c,
      max_counts   = med_c + 3*mad_c,
      max_mito     = med_m + 3*mad_m
    )
  }
  
  # 6) QC 필터링
  seu_filtered <- subset(
    seu,
    subset = nFeature_RNA >  p$min_features  &
             nFeature_RNA <  p$max_features  &
             nCount_RNA   >  p$min_counts    &
             nCount_RNA   <  p$max_counts    &
             percent.mt   <  p$max_mito
  )
  
  # 7) 필터링 후 셀 수
  n_after <- ncol(seu_filtered)
  
  # 8) 결과 요약 저장
  res_list[[i]] <- data.frame(
    Run          = project,
    BeforeCells  = n_before,
    AfterCells   = n_after,
    RetainedFrac = round(n_after / n_before * 100, 2)
  )

  # 9) QC 완료된 객체를 디스크에 저장
  saveRDS(seu_filtered,
          file = paste0("QCRun", i, ".rds"))}
  
  # ——— 여기서부터 플롯 추가 ———
  
  # a) 원본 seu에 passedQC flag 달기
  kept_cells <- colnames(seu_filtered)
  seu$passedQC <- ifelse(colnames(seu) %in% kept_cells, "Kept", "Removed")
 

  # 1) VlnPlot을 리스트로 생성
  p_after <- VlnPlot(seu, features = features,
                       ncol = 3, combine = FALSE)

  # 2) 각 패널에 컷오프 라인 추가
  # features[1] == "nFeature_RNA"
p_after[[1]] <- p_after[[1]] +
  geom_hline(yintercept = p$min_features, color = "red", linetype = "dashed") +
  geom_hline(yintercept = p$max_features, color = "red", linetype = "dashed")

  # features[2] == "nCount_RNA"
p_after[[2]] <- p_after[[2]] +
  geom_hline(yintercept = p$min_counts,   color = "red", linetype = "dashed") +
  geom_hline(yintercept = p$max_counts,   color = "red", linetype = "dashed")

  # features[3] == "percent.mt"
p_after[[3]] <- p_after[[3]] +
  geom_hline(yintercept = p$max_mito,     color = "red", linetype = "dashed")

# 3) 다시 패치워크로 결합 + 제목
p_after <- wrap_plots(p_after, ncol = 3) +
  plot_annotation(title = paste(project, "- QC Cutoffs"))
  
# d) Scatter plots: 제거된 셀 빨강
  df <- seu@meta.data %>% mutate(passedQC = factor(passedQC, levels = c("Kept","Removed")))
  p_sc1 <- ggplot(df, aes(x = nCount_RNA, y = percent.mt, color = passedQC)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Kept"="black","Removed"="red")) +
    ggtitle("nCount_RNA vs percent.mt")
  p_sc2 <- ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA, color = passedQC)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Kept"="black","Removed"="red")) +
    ggtitle("nCount_RNA vs nFeature_RNA")
  
# e) PDF로 저장
  pdf(paste0("QC_", project, ".pdf"), width = 14, height = 10)
  print((p_before | p_after) / (p_sc1 | p_sc2))
  dev.off()}





# ——— QC plot 생성 루프 끝난 다음에 붙여주세요 ———

# 각 Run별 분포 시각화
for (i in seq_along(files)) {
  project <- projects[i]
  obj     <- readRDS(files[i])
  
  # Seurat 객체 생성 & mito 재계산
  seu <- CreateSeuratObject(counts = obj$counts, meta.data = obj$meta, project = project)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # QC 컷오프 기준 다시 계산
  if (project == "Run1") {
    p <- fixed_params
  } else {
    feat  <- seu$nFeature_RNA
    cnts  <- seu$nCount_RNA
    mito  <- seu$percent.mt
    med_f <- median(feat); mad_f <- mad(feat)
    med_c <- median(cnts); mad_c <- mad(cnts)
    med_m <- median(mito); mad_m <- mad(mito)
    p <- list(
      min_features = med_f - 3*mad_f,
      max_features = med_f + 3*mad_f,
      min_counts   = med_c - 3*mad_c,
      max_counts   = med_c + 3*mad_c,
      max_mito     = med_m + 3*mad_m
    )
  }
  
  # 메타데이터 꺼내기
  df <- seu@meta.data
  
  # 히스토그램 + 컷오프 라인
  p_dist_feat <- ggplot(df, aes(nFeature_RNA)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "black") +
    geom_vline(xintercept = p$min_features, color = "red", linetype = "dashed") +
    geom_vline(xintercept = p$max_features, color = "red", linetype = "dashed") +
    ggtitle(paste(project, "- nFeature_RNA")) + theme_minimal()
  
  p_dist_count <- ggplot(df, aes(nCount_RNA)) +
    geom_histogram(bins = 100, fill = "lightgreen", color = "black") +
    geom_vline(xintercept = p$min_counts, color = "red", linetype = "dashed") +
    geom_vline(xintercept = p$max_counts, color = "red", linetype = "dashed") +
    ggtitle(paste(project, "- nCount_RNA")) + theme_minimal()
  
  p_dist_mito <- ggplot(df, aes(percent.mt)) +
    geom_histogram(bins = 100, fill = "lightpink", color = "black") +
    geom_vline(xintercept = p$max_mito, color = "red", linetype = "dashed") +
    ggtitle(paste(project, "- percent.mt")) + theme_minimal()
  
  # PDF로 저장
  pdf(paste0("Dist_", project, ".pdf"), width = 14, height = 5)
  print(p_dist_feat | p_dist_count | p_dist_mito)
  dev.off()
}

 