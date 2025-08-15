library(ggplot2)
library(patchwork)

# meta 데이터는 이미 준비되어 있다고 가정
# seu <- readRDS("sct_annotated_1.rds")
seu$batch <- factor(seu$orig.ident, levels = c("Run1","Run2_1","Run2_2"))
meta <- seu@meta.data

# 1) 각 plot 정의
p_group <- ggplot(meta, aes(x = factor(group,
                                       levels = c("Healthy_Young","Healthy_Old","Diabetes_Old")),
                            fill = celltype)) +
  geom_bar(position = "fill", width = 0.8) +
  coord_flip() +
  labs(title = "A. HY vs HO vs DO (Group)",
       x = NULL, y = "Proportion") +
  theme_minimal(base_size = 12)

p_patient <- ggplot(meta, aes(x = factor(patient,
                                         levels = c(
                                           "HY1","HY2","HY3","HY4","HY5",
                                           "HO1","HO2","HO3","HO4","HO5",
                                           "DO1","DO2","DO3","DO4","DO5")),
                              fill = celltype)) +
  geom_bar(position = "fill", width = 0.8) +
  coord_flip() +
  labs(title = "B. Patients",
       x = NULL, y = "Proportion") +
  theme_minimal(base_size = 12)

p_batch <- ggplot(meta, aes(x = batch, fill = celltype)) +
  geom_bar(position = "fill", width = 0.8) +
  coord_flip() +
  labs(title = "C. Run1 vs Run2_1 vs Run2_2 (Batch)",
       x = NULL, y = "Proportion") +
  theme_minimal(base_size = 12)

# 2) 패치워크로 결합
combined <- (p_group + p_patient + p_batch) +
  plot_layout(ncol = 1,
              heights = c(1, 2, 1),
              guides = "collect") &
  theme(legend.position = "right")

# 3) PDF로 저장
pdf("composition_combined.pdf", width = 14, height = 10)
print(combined)
dev.off()


# 2) UMAP split-by: 그룹 / 환자 / 배치
pdf("UMAP_splits.pdf", width = 12, height = 8)

# A. split by group
DimPlot(seu, reduction = "umap",
        group.by = "celltype", split.by = "group",
        label = FALSE, pt.size = 0.5) + NoLegend()

# B. split by patient
DimPlot(seu, reduction = "umap",
        group.by = "celltype", split.by = "patient",
        label = FALSE, pt.size = 0.5) + NoLegend()

# C. split by batch
DimPlot(seu, reduction = "umap",
        group.by = "celltype", split.by = "orig.ident",
        label = FALSE, pt.size = 0.5) + NoLegend()

dev.off()


########
library(Seurat)
library(ggplot2)
library(patchwork)

# ── 준비 ────────────────────────────────────────────────────────
seu <- readRDS("sct_annotated_1.rds")
seu$batch <- factor(seu$orig.ident, levels = c("Run1","Run2_1","Run2_2"))
Idents(seu) <- "celltype"

# ── A. Group (3개 패싯, 한 행에 3열) ──────────────────────────────
pdf("UMAP_by_group_fixed.pdf", width = 12, height = 4)
DimPlot(seu,
        reduction = "umap",
        split.by  = "group",
        ncol      = 3,        # 한 행에 3개
        pt.size   = 0.5) +
  coord_fixed() +           # 1:1 비율 유지
  theme_minimal(base_size = 12) +
  labs(title = "UMAP by Group",
       x = "UMAP 1", y = "UMAP 2")
dev.off()

# ── B. Patient (15개 패싯, 한 행에 5열 → 3행으로 자동 분산) ─────────
pdf("UMAP_by_patient_fixed.pdf", width = 15, height = 12)
DimPlot(seu,
        reduction = "umap",
        split.by  = "patient",
        ncol      = 5,        # 한 행에 5개씩
        pt.size   = 0.5) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(title = "UMAP by Patient",
       x = "UMAP 1", y = "UMAP 2")
dev.off()

# ── C. Batch (3개 패싯, 한 행에 3열) ───────────────────────────────
pdf("UMAP_by_batch_fixed.pdf", width = 12, height = 4)
DimPlot(seu,
        reduction = "umap",
        split.by  = "batch",
        ncol      = 3,
        pt.size   = 0.5) +
  coord_fixed() +
   theme_minimal(base_size = 12) +
  labs(title = "UMAP by Batch",
       x = "UMAP 1", y = "UMAP 2")
dev.off()




library(Seurat)

# seurat_obj: 이미 UMAP까지 계산돼 있는 객체라고 가정
DefaultAssay(seu) <- "integrated"   # 발현 행렬 사용
pdf('featureplot_0.pdf', width=8, height=8)
FeaturePlot(seu, features = c("TNFAIP2","RNF144B"), blend = TRUE)
dev.off()



# 1) 각 플롯에 aspect ratio 고정
p1 <- FeaturePlot(seu, features = "CD14")   + theme(aspect.ratio = 1)
p2 <- FeaturePlot(seu, features = "FCGR3A")   + theme(aspect.ratio = 1)
p3 <- FeaturePlot(seu, features = c("CD14","FCGR3A"), blend = TRUE) +
      theme(aspect.ratio = 1)

# 2) PDF 장치 열 때 적절한 크기 지정 (여기서는 3패널 + 범례 공간 고려)
pdf("featureplots_fixed.pdf", width = 12, height = 4)

# 3) patchwork 로 붙이기: 각 패널 동일 너비, 마지막에 범례용 스페이서
(p1 + p2 + p3 + plot_spacer()) +
  plot_layout(ncol = 4,
              widths = c(1, 1, 1, 0.3)) & 
  theme(legend.position = "right")

dev.off()