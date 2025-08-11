# 0) 패키지 로드
library(Seurat)
library(edgeR)
library(dplyr)
library(purrr)
library(tibble)

# 1) Seurat 객체 불러오기
seurat_obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')

# 2) 원시 counts (sparse dgCMatrix) 및 meta 준비
counts_raw <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
meta_all   <- seurat_obj@meta.data %>%
  mutate(
    cell     = rownames(.),
    patient  = patient,        # HY1–5, HO1–5, DO1–5
    group    = group,          # Healthy_Young / Healthy_Old / Diabetes_Old
    batch    = orig.ident,     # Run1, Run2_1, Run2_2
    celltype = celltype
  ) %>%
  select(cell, patient, group, batch, celltype)

# 3) Run1 세포만 필터링
meta_sel   <- meta_all %>% filter(batch == "Run1")
cells_run1 <- meta_sel$cell

# 4) Run1 counts 만 base matrix 로 변환
counts_all <- as.matrix(counts_raw[, cells_run1, drop = FALSE])

# 5) pseudobulk 생성 함수 (batch 변수 제거, group 만 사용)
make_pseudobulk <- function(ct) {
  sub_meta     <- meta_sel %>% filter(celltype == ct)
  samples_meta <- sub_meta %>% distinct(patient, group)
  if (nrow(samples_meta) < 2) return(NULL)  # 최소 2명 이상 필요

  # Run1 matrix 에서 바로 subset → matrix 보장
  mat_ct    <- counts_all[, sub_meta$cell, drop = FALSE]
  pb_mat    <- rowsum(t(mat_ct), group = sub_meta$patient)
  pb_counts <- t(pb_mat)

  pb_meta <- samples_meta %>%
    rename(sample = patient) %>%
    mutate(group = factor(group, levels = c("Healthy_Young", "Diabetes_Old"))) %>%
    arrange(sample)

  # 순서 동기화
  keep_samps <- intersect(colnames(pb_counts), pb_meta$sample)
  pb_counts  <- pb_counts[, keep_samps, drop = FALSE]
  pb_meta    <- pb_meta %>%
                  filter(sample %in% keep_samps) %>%
                  slice(match(keep_samps, sample))

  list(counts = pb_counts, meta = pb_meta)
}

# 6) DE 분석 함수 (design ~ group)
analyze_celltype <- function(ct) {
  pb <- make_pseudobulk(ct)
  if (is.null(pb)) return(NULL)

  y <- DGEList(counts = pb$counts, samples = pb$meta)

  # 필터링 & TMM 정규화
  keep <- filterByExpr(y, group = y$samples$group)
  y    <- y[keep, , keep.lib.sizes = FALSE] %>% calcNormFactors()

  # model.matrix: ~ group
  design <- model.matrix(~ group, data = y$samples)

  # dispersion & GLM fitting
  y   <- estimateDisp(y, design)
  fit <- glmFit(y, design)

  # global test
  res_global <- topTags(glmLRT(fit), n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(celltype  = ct, comparison = "global")

  # HY vs DO contrast
  contr  <- makeContrasts(DO_vs_HY = groupDiabetes_Old, levels = design)
  res_pw <- topTags(glmLRT(fit, contrast = contr), n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(celltype  = ct, comparison = "DO_vs_HY")

  bind_rows(res_global, res_pw)
}

# 7) 전체 cell type 분석 실행
celltypes   <- unique(meta_sel$celltype)
deg_results <- map(celltypes, analyze_celltype) %>%
               compact() %>%
               bind_rows()






# ── 0) 패키지 설치/로드 ──
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

# ── 1) volcano plot 함수 재정의 ──
make_volcano <- function(df, ct, fc_thresh = 0.58, fdr_thresh = 0.05) {
  df2 <- df %>%
    filter(comparison == "DO_vs_HY") %>%    # Run1은 DO_vs_HY 하나만
    mutate(
      regulation = case_when(
        logFC >  fc_thresh & FDR < fdr_thresh  ~ "Up",
        logFC < -fc_thresh & FDR < fdr_thresh  ~ "Down",
        TRUE                                   ~ "None"
      )
    )
  
  ggplot(df2, aes(x = logFC, y = -log10(FDR), color = regulation)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
      values = c(Up = "red", Down = "blue", None = "gray"),
      breaks = c("Up","Down")
    ) +
    geom_text_repel(
      data = df2 %>% filter(regulation %in% c("Up","Down")),
      aes(label = gene),
      max.overlaps   = Inf,
      box.padding    = 0.3,
      point.padding  = 0.2
    ) +
    labs(
      title = ct,
      x     = expression(log[2]~FC),
      y     = expression(-log[10]~FDR),
      color = "Regulation"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}

# ── 2) Run1 DO_vs_HY 결과만 뽑아서 plot 리스트 생성 ──
df_pw <- deg_results %>% filter(comparison == "DO_vs_HY")
plots  <- df_pw %>%
  group_split(celltype) %>%
  set_names(map_chr(., ~ unique(.x$celltype))) %>%
  map(~ make_volcano(.x, unique(.x$celltype)))

# ── 3) PDF로 저장 ──
pdf("volcano_Run1_DOvsHY_updown.pdf", width = 6, height = 6)
for(p in plots) print(p)
dev.off()

