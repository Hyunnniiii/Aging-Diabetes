## =========================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | 전체 파일에서 T/NK/Mono를 각각 분리해 배치별 SCTransform → SCT-anchor 통합 → PCA/UMAP/클러스터링(res=0.2–1.0) 
## UMAP과 RDS 저장, 해당 라인리지 UMAP 위에 기존 celltype 주석을 확인용으로 겹쳐 저장
## =========================================


# ===== Seurat v5: Lineage-wise subtyping (T / NK / Mono run separately) =====
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(future)
})

# 안전: 병렬 off + 글로벌 사이즈 넉넉히 (필요시 조정)
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

set.seed(1234)

# --- 입력 RDS ---
obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/sct_annotated_1.rds')
message("Seurat: ", packageVersion("Seurat"), " | SeuratObject: ", packageVersion("SeuratObject"))
stopifnot("RNA" %in% Assays(obj))
stopifnot("counts" %in% Layers(obj[["RNA"]]))

# --- supertype 없으면 만들어두기(있으면 이 블록은 그냥 지나감) ---
if (!"supertype" %in% colnames(obj@meta.data)) {
  t_lvls    <- c("CD4 TCM","naive CD4","CD8 TEM","naive CD8", "Treg")
  nk_lvls   <- c("NK")
  mono_lvls <- c("CD14 Mono","CD16 Mono","intermediate Mono")
  obj$supertype <- NA_character_
  obj$supertype[obj$celltype %in% t_lvls]    <- "T"
  obj$supertype[obj$celltype %in% nk_lvls]   <- "NK"
  obj$supertype[obj$celltype %in% mono_lvls] <- "Mono"
  obj$supertype <- factor(obj$supertype, levels = c("T","NK","Mono"))
}

stopifnot("supertype" %in% colnames(obj@meta.data))
stopifnot("orig.ident" %in% colnames(obj@meta.data))  # 배치 열 필요

# --- 공통 설정 ---
out_root <- "Subtyping"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
res_default <- 0.4

# 라인리지별 TCR/Ig HVG 제거 여부(보통 T만 TRUE)
remove_receptors_by_group <- c("T" = FALSE, "NK" = FALSE, "Mono" = FALSE)

# === 함수: 한 라인리지만 처리 ===
process_lineage <- function(obj, lineage, remove_receptors = FALSE) {
  msg <- function(...) message(sprintf("[%s] %s", lineage, sprintf(...)))

  out_dir <- file.path(out_root, lineage)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 1) 서브셋
  Idents(obj) <- "supertype"
  sub <- subset(obj, idents = lineage)
  stopifnot(ncol(sub) > 0)
  msg("cells: %d", ncol(sub))

  # RNA counts 레이어가 서브셋에 없으면 부모에서 복사(드문 케이스 대비)
  if (!"counts" %in% Layers(sub[["RNA"]])) {
    cnt_parent <- GetAssayData(obj, assay = "RNA", layer = "counts")[, colnames(sub), drop = FALSE]
    sub <- SetAssayData(sub, assay = "RNA", layer = "counts", new.data = cnt_parent)
  }

  # 2) 배치별 Split & SCTransform (새 어세이명으로 충돌 방지)
  sub.list <- SplitObject(sub, split.by = "orig.ident")
  glm_ok <- requireNamespace("glmGamPoi", quietly = TRUE)

  sub.list <- lapply(sub.list, function(s1) {
    DefaultAssay(s1) <- "RNA"
    SCTransform(
      s1,
      method              = if (glm_ok) "glmGamPoi" else "poisson",
      variable.features.n = 5000,        # 교집합 늘리기
      new.assay.name      = "SCTsub",    # 기존 SCT와 분리
      verbose             = TRUE
    )
  })

  # 3) 앵커 피처: 실제 scale.data 교집합만 사용
  scale_genes_list   <- lapply(sub.list, function(x) rownames(GetAssayData(x, assay = "SCTsub", layer = "scale.data")))
  common_scale_genes <- Reduce(intersect, scale_genes_list)
  msg("common scale.genes: %d", length(common_scale_genes))

  features <- SelectIntegrationFeatures(sub.list, nfeatures = 3000)
  features <- intersect(features, common_scale_genes)
  if (remove_receptors) {
    features <- setdiff(features, grep("^TR[ABDG]|^IG[HKL]", features, value = TRUE))
  }
  msg("anchor features after filters: %d", length(features))
  if (length(features) < 1200) {
    msg("WARNING: anchor features 적음(%d). variable.features.n을 더 키우는 걸 고려하세요.", length(features))
  }

  # 4) Prep → Anchors → Integrate (SCT 통합)
  sub.list <- PrepSCTIntegration(sub.list, anchor.features = features, verbose = TRUE)
  anchors  <- FindIntegrationAnchors(
    sub.list,
    normalization.method = "SCT",
    anchor.features      = features,
    dims                 = 1:30,
    verbose              = TRUE
  )
  sub.int <- suppressWarnings(
    IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)
  )

  # 5) 차원축소/클러스터 + 저장/플롯
  DefaultAssay(sub.int) <- "integrated"
  sub.int <- RunPCA(sub.int, npcs = 50, verbose = TRUE)
  sub.int <- RunUMAP(sub.int, dims = 1:30, verbose = TRUE)
  sub.int <- FindNeighbors(sub.int, dims = 1:30, verbose = TRUE)

  for (res in resolutions) {
    sub.int <- FindClusters(sub.int, resolution = res, verbose = TRUE)
    p <- DimPlot(sub.int, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.2) +
         ggtitle(sprintf("%s — Clusters (res=%.1f)", lineage, res))
    ggsave(file.path(out_dir, sprintf("UMAP_%s_clusters_res%.1f.pdf", lineage, res)),
           p, width = 7, height = 6, units = "in")
  }

  Idents(sub.int) <- sub.int[[sprintf("integrated_snn_res.%.1f", 0.4)]][,1]

  ggsave(file.path(out_dir, sprintf("UMAP_%s_by_batch.pdf", lineage)),
         DimPlot(sub.int, reduction = "umap", group.by = "orig.ident", pt.size = 0.2) +
           ggtitle(sprintf("%s — Batch (orig.ident)", lineage)),
         width = 7, height = 6, units = "in")

  # 필요하면 patient/group 오버레이도 추가 가능
  if ("patient" %in% colnames(sub.int@meta.data)) {
    ggsave(file.path(out_dir, sprintf("UMAP_%s_by_patient.pdf", lineage)),
           DimPlot(sub.int, reduction = "umap", group.by = "patient", pt.size = 0.2) +
             ggtitle(sprintf("%s — patient", lineage)),
           width = 7, height = 6, units = "in")
  }

  saveRDS(sub.int, file = file.path(out_dir, sprintf("subset_%s_SCTintegrated.rds", lineage)))
  
  msg("done → %s", normalizePath(out_dir))
  invisible(sub.int)
}

# === 실제 실행: T / NK / Mono 각각 따로 ===
lineages <- c("T","NK","Mono")
for (lin in lineages) {
  process_lineage(
    obj,
    lineage = lin,
    remove_receptors = isTRUE(remove_receptors_by_group[[lin]])
  )
}

message("All done. Outputs in: ", normalizePath(out_root))


## Umap 상에서 원래 annotation 확인
## ===== Check original annotations on lineage UMAPs =====
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})

lineages <- c("T","NK","Mono")
for (lin in lineages) {
  rds_path <- file.path("/data/project/diabetes_LYH/tanya/rds", sprintf("subset_%s_SCTintegrated.rds", lin))
  stopifnot(file.exists(rds_path))
  sub.int <- readRDS(rds_path)

  um <- Embeddings(sub.int, "umap"); xr <- range(um[,1], finite=TRUE); yr <- range(um[,2], finite=TRUE)
  fix_axes <- function(p) p + coord_fixed(xlim=xr, ylim=yr) + theme(aspect.ratio=1)

  # (선택) 클러스터 아이덴티티를 0.4로 맞춰두기
  idcol <- sprintf("integrated_snn_res.%.1f", 0.4)
  if (idcol %in% colnames(sub.int@meta.data)) {
    Idents(sub.int) <- sub.int[[idcol]][,1]
  }

  stopifnot("celltype" %in% colnames(sub.int@meta.data))
  out_dir <- file.path("Subtyping", lin, "AnnoCheck")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ## (1) 기본: celltype 라벨 UMAP
  p_main <- DimPlot(sub.int, reduction = "umap", group.by = "celltype", pt.size = 0.2) +
            ggtitle(sprintf("%s — celltype (labeled)", lin))
  p_main <- LabelClusters(p_main, id = "celltype", repel = TRUE)
  p_main <- fix_axes(p_main)

  ## (2) group별 facet + 라벨
  p_group <- NULL
  if ("group" %in% colnames(sub.int@meta.data)) {
    plist <- DimPlot(sub.int, reduction = "umap", group.by = "celltype",
                     split.by = "group", pt.size = 0.2, combine = FALSE)
    plist <- lapply(plist, function(p) fix_axes(LabelClusters(p, id = "celltype", repel = TRUE)))
    ncol_fac <- min(3, length(plist))
    p_group <- wrap_plots(plist, ncol = ncol_fac) +
               plot_annotation(title = sprintf("%s — celltype faceted by group", lin))
  }

  ## (3) orig.ident별 facet + 라벨
  p_orig <- NULL
  if ("orig.ident" %in% colnames(sub.int@meta.data)) {
    plist2 <- DimPlot(sub.int, reduction = "umap", group.by = "celltype",
                      split.by = "orig.ident", pt.size = 0.2, combine = FALSE)
    plist2 <- lapply(plist2, function(p) fix_axes(LabelClusters(p, id = "celltype", repel = TRUE)))
    ncol_fac2 <- min(3, length(plist2))
    p_orig <- wrap_plots(plist2, ncol = ncol_fac2) +
              plot_annotation(title = sprintf("%s — celltype faceted by orig.ident", lin))
  }

  ## (4) 한 PDF에 순서대로 저장 (페이지 1: 기본, 2: group, 3: orig.ident)
  pdf_path <- file.path(out_dir, sprintf("UMAP_%s_celltype_ALL.pdf", lin))
  pdf(pdf_path, width = 10, height = 8)
  for (pg in list(p_main, p_group, p_orig)) if (!is.null(pg)) print(pg)
  dev.off()
}


