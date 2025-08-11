# ── 0) 패키지 로드 ──
library(Seurat)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(tibble)

# ── 1) Seurat 객체 + raw counts/metadata 준비 ──
seurat_obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')
counts_raw <- GetAssayData(seurat_obj, assay="RNA", layer="counts")
meta_all   <- seurat_obj@meta.data %>%
  mutate(cell = rownames(.), batch = orig.ident) %>%
  select(cell, patient, group, batch)

# ── 2) Run2_1 & Run2_2 샘플만 필터링 ──
meta_sel   <- meta_all %>%
  filter(batch %in% c("Run2_1","Run2_2"))
cells_sel  <- meta_sel$cell
counts_sel <- as.matrix(counts_raw[, cells_sel, drop = FALSE])

# ── 3) 전체 pseudobulk (celltype 무시!) ──
pb_mat_all    <- rowsum(t(counts_sel), group = meta_sel$patient)
pb_counts_all <- t(pb_mat_all)

pb_meta_all <- meta_sel %>%
  distinct(patient, group) %>%
  rename(sample = patient) %>%
  # 여기서 꼭 3개 그룹을 levels 에 모두 넣어주세요
  mutate(group = factor(group,
                        levels = c("Healthy_Young","Healthy_Old","Diabetes_Old"))) %>%
  arrange(sample)

# ── 4) edgeR DE 분석 ──
y_all      <- DGEList(counts = pb_counts_all, samples = pb_meta_all)
keep_all   <- filterByExpr(y_all, group = y_all$samples$group)
y_all      <- y_all[keep_all, , keep.lib.sizes = FALSE] %>% calcNormFactors()
design_all <- model.matrix(~ group, data = y_all$samples)
y_all      <- estimateDisp(y_all, design_all)
fit_all    <- glmFit(y_all, design_all)

# 3가지 대비 정의 & 결과 합치기
contrs <- list(
  HO_vs_HY = makeContrasts(groupHealthy_Old,                     levels = design_all),
  DO_vs_HY = makeContrasts(groupDiabetes_Old,                    levels = design_all),
  DO_vs_HO = makeContrasts(groupDiabetes_Old - groupHealthy_Old, levels = design_all)
)

res_all <- imap_dfr(contrs, ~ {
  tt <- glmLRT(fit_all, contrast = .x) %>% topTags(n = Inf)
  tt$table %>%
    rownames_to_column("gene") %>%
    mutate(
      comparison = .y,
      regulation = case_when(
        logFC >  0.58 & FDR < 0.05 ~ "Up",
        logFC < -0.58 & FDR < 0.05 ~ "Down",
        TRUE                        ~ "None"
      ),
      negLogFDR = -log10(FDR)
    )
})

# ── 5) volcano plot 3페이지짜리 PDF 저장 ──
pdf("Run2_allDEG.pdf", width = 6, height = 8)
for (comp in unique(res_all$comparison)) {
  df <- filter(res_all, comparison == comp)
  print(
    ggplot(df, aes(x = logFC, y = negLogFDR, color = regulation)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c(Up="red", Down="blue", None="gray")) +
      geom_hline(yintercept = -log10(0.05), linetype="dashed") +
      geom_vline(xintercept = c(-0.58, 0.58),  linetype="dashed") +
      geom_text_repel(
        data = df %>% filter(regulation != "None"),
        aes(label = gene),
        max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2
      ) +
      labs(
        title = paste0("Run2 pseudobulk: ", comp),
        x     = expression(log[2]~FC),
        y     = expression(-log[10]~FDR),
        color = "Regulation"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 12)
      )
  )
}
dev.off()
