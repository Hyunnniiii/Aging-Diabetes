suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(MAST)
})

## ===== 설정 =====
PARENT_RDS <- "/data/project/diabetes_LYH/tanya/rds/sct_annotated_1.rds"
lineages   <- c("T","NK","Mono")
GROUP_COL  <- "group"
grp_lvls   <- c("Healthy_Young","Healthy_Old","Diabetes_Old")
min_cells  <- 30
latent_candidates <- c("nCount_RNA","nFeature_RNA","percent.mt","orig.ident","S.Score","G2M.Score")

contrasts <- list(
  HO_vs_HY = c("Healthy_Old", "Healthy_Young"),
  DO_vs_HY = c("Diabetes_Old","Healthy_Young"),
  DO_vs_HO = c("Diabetes_Old","Healthy_Old")
)

sub_rds_path <- function(lin) file.path("/data/project/diabetes_LYH/tanya/rds", sprintf("subset_%s_SCTintegrated.rds", lin))
out_root <- file.path("Subtyping","MAST_DEG")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## ===== 유틸 =====
# 부모 counts로 RNA assay 재구성 (cells 이름 매칭 필수)
rebuild_rna_from_parent <- function(subobj, parent, assay="RNA") {
  cells <- colnames(subobj)
  stopifnot(all(cells %in% colnames(parent)))
  stopifnot("counts" %in% Layers(parent[[assay]]))
  cnt <- GetAssayData(parent, assay=assay, layer="counts")[, cells, drop = FALSE]
  subobj[[assay]] <- CreateAssayObject(counts = cnt)  # 클래스는 환경에 따라 Assay/Assay5
  DefaultAssay(subobj) <- assay
  subobj
}

# Assay 클래스에 맞춰 정규화 + (필요 시) JoinLayers
prep_rna_for_mast <- function(subobj, assay="RNA") {
  # 정규화('data' 생성 또는 data slot 채움)
  DefaultAssay(subobj) <- assay
  # Assay5 여부 판단
  is_a5 <- inherits(subobj[[assay]], "Assay5")
  if (is_a5) {
    if (!("data" %in% Layers(subobj[[assay]]))) {
      subobj <- NormalizeData(subobj)  # 'data' layer 생성
    }
    # 레이어 조인(Assay5에서만)
    subobj[[assay]] <- JoinLayers(subobj[[assay]])
  } else {
    # 구버전 Assay: layer 개념 없음, NormalizeData가 data slot 채움
    if (!"data" %in% slotNames(subobj[[assay]])) {
      subobj <- NormalizeData(subobj)
    }
  }
  list(obj=subobj, is_a5=is_a5)
}

# FindMarkers(MAST) 안전 호출: Assay5면 layer="data", 아니면 layer 인자 생략
run_mast_pair <- function(obj, g1, g2, latent_vars, assay="RNA", is_a5=FALSE) {
  if (is_a5) {
    fm <- FindMarkers(
      obj, ident.1=g1, ident.2=g2,
      test.use="MAST", assay=assay, layer="data",
      logfc.threshold=0, min.pct=0.1,
      latent.vars=latent_vars
    )
  } else {
    fm <- FindMarkers(
      obj, ident.1=g1, ident.2=g2,
      test.use="MAST", assay=assay,            # layer 인자 없이
      logfc.threshold=0, min.pct=0.1,
      latent.vars=latent_vars
    )
  }
  tibble::rownames_to_column(fm, "gene") |>
    mutate(contrast = paste0(g1, "_vs_", g2))
}

## ===== 메인 =====
message("Loading parent...")
parent <- readRDS(PARENT_RDS)

for (lin in lineages) {
  message(sprintf("\n[%s] ===", lin))
  sobj <- readRDS(sub_rds_path(lin))

  # RNA 재구성(부모 counts) → 정규화(+필요시 JoinLayers)
  sobj <- rebuild_rna_from_parent(sobj, parent, assay="RNA")
  prep  <- prep_rna_for_mast(sobj, assay="RNA")
  sobj  <- prep$obj
  is_a5 <- prep$is_a5

  # 그룹 팩터/Idents
  stopifnot(GROUP_COL %in% colnames(sobj@meta.data))
  sobj[[GROUP_COL]] <- factor(as.character(sobj[[GROUP_COL]][,1]), levels = grp_lvls)
  Idents(sobj) <- sobj[[GROUP_COL]][,1]

  # latent vars (있는 것만)
  latent_vars <- intersect(latent_candidates, colnames(sobj@meta.data))

  # 셀 수 체크
  tab <- table(Idents(sobj))
  present <- names(tab)[tab >= min_cells]
  run_list <- Filter(function(pp) all(pp %in% present), contrasts)
  if (!length(run_list)) {
    message(sprintf("[%s] 대비 없음(각 그룹 ≥ %d cells 필요).", lin, min_cells))
    next
  }

  # 대비별 실행
  res_list <- lapply(names(run_list), function(nm){
    pp <- run_list[[nm]]
    run_mast_pair(sobj, pp[1], pp[2], latent_vars, assay="RNA", is_a5=is_a5)
  })
  res_all <- dplyr::bind_rows(res_list)

  # 저장
  out_lin <- file.path(out_root, lin); dir.create(out_lin, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(res_all, file.path(out_lin, sprintf("MAST_DEG_%s_allContrasts.csv", lin)))
  readr::write_csv(as.data.frame(tab) |>
                     dplyr::rename(group=Var1, cells=Freq) |>
                     dplyr::mutate(lineage=lin),
                   file.path(out_lin, sprintf("cells_per_group_%s.csv", lin)))

  # 요약
  summary <- res_all |>
    dplyr::group_by(contrast) |>
    dplyr::summarise(
      n_padj_0.05 = sum(p_val_adj < 0.05, na.rm = TRUE),
      n_padj_0.1  = sum(p_val_adj < 0.1,  na.rm = TRUE),
      .groups="drop"
    ) |>
    dplyr::mutate(lineage=lin)
  readr::write_csv(summary, file.path(out_lin, sprintf("summary_%s.csv", lin)))

  message(sprintf("[%s] done → %s", lin, normalizePath(out_lin)))
}

message("\nAll done → ", normalizePath(out_root))


