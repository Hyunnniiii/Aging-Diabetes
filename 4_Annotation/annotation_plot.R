# ============================================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | singleR 결과 umap으로 확인, 특정 cluster만 확인, 특정 gene 발현 확인
# ============================================================


#########
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

# 1. 불러오기
seu <- readRDS("rds/seurat_with_singler.rds")
singler.res <- readRDS("rds/singler.res.rds")   

# 1. UMAP 좌표와 annotation label 추출
umap <- Embeddings(seu, "umap")
df <- data.frame(UMAP_1 = umap[,1], UMAP_2 = umap[,2], label = seu@meta.data$SingleR)


df <- df %>% filter(!is.na(label) & label != "NA")


# 2. label별 중심좌표 계산 (median 기준, 평균 써도 무방)
centers <- df %>%
  group_by(label) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))


# 3. UMAP plot with label (글씨 겹치지 않게!)
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = label)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_text_repel(
    data = centers,
    aes(label = label),
    color = 'black',
    size = 4,
    max.overlaps = 30,
    box.padding = 0.6
  ) +
  ggtitle("SingleR annotation") +
  theme_classic() +
  theme(plot.title = element_text(size=22, face='bold', hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size=5)))

pdf("umap_singler_seuratstyle.pdf", width=8, height=6)
print(p)
dev.off()


#########
# 원하는 cell만 색상, 나머지는 회색
# 주요 7개 cell만 색상, 나머지는 회색
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)

main_cells <- c("T_cells", "NK_cell", "Monocyte", "B_cell", "DC", "Neutrophils", "Pre-B_cell_CD34-", 'Pro-B_cell_CD34+')
df$label2 <- ifelse(df$label %in% main_cells, df$label, "Other")  # 기타는 Other로 통합

# 주요 cell color + 기타 회색
cell_types2 <- c(main_cells, "Other")
main_colors <- viridis_pal(option="D")(length(main_cells))
mycolors2 <- c(main_colors, "#D3D3D3")  # #D3D3D3 = light gray
names(mycolors2) <- cell_types2

# 중심좌표 (주요 cell만)
centers_main <- df %>%
  filter(label %in% main_cells) %>%
  group_by(label) %>%
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# Plot
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = label2)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_text_repel(
    data = centers_main,
    aes(x = UMAP_1, y = UMAP_2, label = label),
    color = "black", size = 4, max.overlaps = 30, box.padding = 0.6
  ) +
  scale_color_manual(values = mycolors2) +
  ggtitle("Main Cell Types Highlighted") + theme_classic() +
  theme(plot.title=element_text(size=22,face='bold', hjust=0.5)) +
  guides(color = guide_legend(override.aes = list(size=5)))

pdf("umap_main_gray_others.pdf", width=8, height=6)
print(p)
dev.off()



##########
# 모든 annotation labeled celltype에 대해 plot
library(Seurat)
library(ggplot2)
library(dplyr)

seu <- readRDS("rds/seurat_with_singler.rds")
umap <- Embeddings(seu, "umap")
df <- data.frame(
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2],
  label = seu@meta.data$SingleR
)

df <- df %>% filter(!is.na(label))
all_labels <- unique(df$label)

# PDF 저장 (혹시 이전 파일이 열려 있으면, 세션을 껐다가 다시 켜도 좋아요)
pdf("umap_single_highlight_each_label.pdf", width=8, height=6)

for (this_label in all_labels) {
  tryCatch({
    cat("Processing: ", this_label, "\n")  # 진단 메시지
    df$plot_label <- ifelse(df$label == this_label, this_label, "Other")
    mycolors <- c("#FF5252", "#D3D3D3")
    names(mycolors) <- c(this_label, "Other")
    centers <- df %>%
      filter(plot_label == this_label) %>%
      group_by(plot_label) %>%
      summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=plot_label)) +
      geom_point(size=0.5, alpha=0.7) +
      geom_text(
        data=centers,
        aes(x=UMAP_1, y=UMAP_2, label=plot_label),
        color="black", size=5
      ) +
      scale_color_manual(values=mycolors) +
      ggtitle(paste0(this_label, " only highlighted")) +
      theme_classic() +
      theme(
        plot.title = element_text(size=18,face='bold', hjust=0.5),
        legend.position = "none"
      )
    print(p)
  }, error=function(e){
    cat("Error with label: ", this_label, " -- ", e$message, "\n")
  })
}
dev.off()
cat("PDF writing complete.\n")




#########
# 특정 cluster에서 특정 gene 발현 확인하기
# 1. cluster 14에 속하는 cell의 이름(Barcode) 추출
cluster14_cells <- WhichCells(oo, idents = 14)

# 2. 해당 cell에서 CD38 유전자 발현값 추출
cd19_exp_in_cluster14 <- FetchData(oo, vars = "CD19", cells = cluster14_cells)

# 3. 결과 예시 확인 (head)
head(cd38_exp_in_cluster3)
