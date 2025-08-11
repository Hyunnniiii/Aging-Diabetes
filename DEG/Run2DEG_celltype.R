# ── 0) 패키지 로드 ──
library(Seurat)
library(edgeR)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggrepel)

# ── 1) Seurat 객체 불러오기 ──
seurat_obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')

# ── 2) raw counts + metadata 준비 ──
counts_raw <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
meta_all   <- seurat_obj@meta.data %>%
  mutate(
    cell     = rownames(.),
    patient  = patient,         # HY1–5, HO1–5, DO1–5
    group    = group,           # Healthy_Young / Healthy_Old / Diabetes_Old
    batch    = orig.ident,      # Run1, Run2_1, Run2_2
    celltype = celltype
  ) %>%
  select(cell, patient, group, batch, celltype)

# ── 3) Run2_1 및 Run2_2 샘플만 필터링 ──
meta_sel   <- meta_all %>% filter(batch %in% c("Run2_1", "Run2_2"))
cells_run2 <- meta_sel$cell

# ── 4) 해당 세포들만 base matrix 로 변환 ──
counts_all <- as.matrix(counts_raw[, cells_run2, drop = FALSE])

# ── 5) pseudobulk 생성 함수 ──
make_pseudobulk <- function(ct) {
  sub_meta     <- meta_sel %>% filter(celltype == ct)
  samples_meta <- sub_meta %>% distinct(patient, group)
  if (nrow(samples_meta) < 2) return(NULL)

  mat_ct    <- counts_all[, sub_meta$cell, drop = FALSE]
  pb_mat    <- rowsum(t(mat_ct), group = sub_meta$patient)
  pb_counts <- t(pb_mat)

  pb_meta <- samples_meta %>%
    rename(sample = patient) %>%
    mutate(group = factor(group,
                          levels = c("Healthy_Young", "Healthy_Old", "Diabetes_Old"))) %>%
    arrange(sample)

  keep_samps <- intersect(colnames(pb_counts), pb_meta$sample)
  pb_counts  <- pb_counts[, keep_samps, drop = FALSE]
  pb_meta    <- pb_meta %>% filter(sample %in% keep_samps) %>%
                           slice(match(keep_samps, sample))

  list(counts = pb_counts, meta = pb_meta)
}

# ── 6) DE 분석 함수: ~ group + 3-way contrasts ──
analyze_celltype <- function(ct) {
  pb <- make_pseudobulk(ct)
  if (is.null(pb)) return(NULL)

  y <- DGEList(counts = pb$counts, samples = pb$meta)
  keep <- filterByExpr(y, group = y$samples$group)
  y    <- y[keep, , keep.lib.sizes = FALSE] %>% calcNormFactors()

  design <- model.matrix(~ group, data = y$samples)
  y      <- estimateDisp(y, design)
  fit    <- glmFit(y, design)

  # 글로벌 테스트
  res_global <- topTags(glmLRT(fit), n = Inf)$table %>%
    rownames_to_column("gene") %>%
    mutate(celltype  = ct, comparison = "global")

  # 세 가지 대비 정의
  contrs <- list(
    HO_vs_HY = makeContrasts(groupHealthy_Old, levels = design),
    DO_vs_HY = makeContrasts(groupDiabetes_Old, levels = design),
    DO_vs_HO = makeContrasts(groupDiabetes_Old - groupHealthy_Old, levels = design)
  )
  res_pw <- imap_dfr(contrs, ~ {
    tt <- glmLRT(fit, contrast = .x) %>% topTags(n = Inf)
    tt$table %>%
      rownames_to_column("gene") %>%
      mutate(celltype  = ct,
             comparison = .y)
  })

  bind_rows(res_global, res_pw)
}

# ── 7) 모든 cell type DE 실행 ──
celltypes   <- unique(meta_sel$celltype)
deg_results <- map(celltypes, analyze_celltype) %>% compact() %>% bind_rows()

# ── 8) volcano plot 함수 재정의 (Up/Down 모두 표시) ──
make_volcano <- function(df, ct, fc_thresh = 0.58, fdr_thresh = 0.05) {
  df2 <- df %>%
    filter(comparison != "global") %>%
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
      breaks = c("Up","Down"), labels = c("Up","Down")
    ) +
    geom_text_repel(
      data = df2 %>% filter(regulation %in% c("Up","Down")),
      aes(label = gene),
      max.overlaps   = Inf,
      box.padding    = 0.3,
      point.padding  = 0.2
    ) +
    facet_wrap(~ comparison, ncol = 1) +
    labs(
      title = ct,
      x     = expression(log[2]~FC),
      y     = expression(-log[10]~FDR),
      color = "Regulation"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 8)
    )
}

# ── 9) volcano plot 생성 & PDF 저장 ──
df_pw <- deg_results %>% filter(comparison != "global")
plots <- df_pw %>%
  group_split(celltype) %>%
  set_names(map_chr(., ~ unique(.x$celltype))) %>%
  map(~ make_volcano(.x, unique(.x$celltype)))

pdf("volcano_Run2_allContrasts_updown.pdf", width = 6, height = 8)
for (p in plots) print(p)
dev.off()

