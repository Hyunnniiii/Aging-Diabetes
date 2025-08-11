# ── 0) 패키지 로드 ──
library(Seurat)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)

# ── 1) Seurat 객체 + raw counts/metadata 준비 ──
seurat_obj <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')
counts_raw <- GetAssayData(seurat_obj, assay="RNA", layer="counts")
meta_all   <- seurat_obj@meta.data %>%
  mutate(cell = rownames(.), batch = orig.ident) %>%
  select(cell, patient, group, batch)

# ── 2) (예시) Run1만 합칠 경우, Run1 샘플 필터링 ──
meta_sel   <- meta_all %>% filter(batch == "Run1")
cells_sel  <- meta_sel$cell
counts_sel <- as.matrix(counts_raw[, cells_sel, drop=FALSE])

# ── 3) 전체 pseudobulk (celltype 무시!) ──
#    환자별로 모든 세포 합산
pb_mat_all    <- rowsum(t(counts_sel), group = meta_sel$patient)
pb_counts_all <- t(pb_mat_all)   # genes × patients

#    sample 메타데이터
pb_meta_all <- meta_sel %>%
  distinct(patient, group) %>%
  rename(sample = patient) %>%
  mutate(group = factor(group, levels=c("Healthy_Young","Diabetes_Old"))) %>%
  arrange(sample)

# ── 4) edgeR DE 분석 ──
y_all <- DGEList(counts = pb_counts_all, samples = pb_meta_all)

# 4.1) 필터링 & 정규화
keep_all <- filterByExpr(y_all, group = y_all$samples$group)
y_all    <- y_all[keep_all, , keep.lib.sizes=FALSE] %>% calcNormFactors()

# 4.2) 디자인 매트릭스
design_all <- model.matrix(~ group, data = y_all$samples)

# 4.3) 분산 추정 & GLM 적합
y_all <- estimateDisp(y_all, design_all)
fit_all <- glmFit(y_all, design_all)

# 4.4) contrast 정의 & LRT
ctr <- makeContrasts(DO_vs_HY = groupDiabetes_Old, levels = design_all)
lrt_all <- glmLRT(fit_all, contrast = ctr)
res_all  <- topTags(lrt_all, n=Inf)$table %>%
  rownames_to_column("gene") %>%
  mutate(
    regulation = case_when(
      logFC >  0.58 & FDR < 0.05 ~ "Up",
      logFC < -0.58 & FDR < 0.05 ~ "Down",
      TRUE                         ~ "None"
    ),
    negLogFDR = -log10(FDR)
  )

# ── 5) volcano plot ──
pdf('Run1_allDEG.pdf')
ggplot(res_all, aes(x = logFC, y = negLogFDR, color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(Up="red", Down="blue", None="gray")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = c(-0.58, 0.58),  linetype="dashed") +
  labs(
    title = "Overall pseudobulk DO vs HY Volcano (Run1)",
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
dev.off()
