library(Seurat)

# v5에서 count matrix만 뽑아내기
obj <- readRDS("integrated_SCT_with_markers.rds")
print(str(obj@assays))

rna_counts <- obj@assays$RNA@layers[["counts"]]
rownames(rna_counts) <- rownames(features_obj@.Data)
colnames(rna_counts) <- rownames(cells_obj@.Data)

meta_data <- obj@meta.data
saveRDS(list(counts=rna_counts, meta=meta_data), "v4_aziref.rds")

# Seurat v4 object 만들기 
library(Seurat)

newobj <- readRDS('v4_aziref.rds')

v4_aziref_seurat <- CreateSeuratObject(newobj$counts, meta.data=newobj$meta)

saveRDS(v4_aziref_seurat, 'v4_seurat_aziref.rds')

# web에서 Azimuth 돌린 뒤
anno <- read.table("/Annotation/azimuth/azimuth_predI2.tsv", header=TRUE, sep="\t")

obj_barcodes <- rownames(obj@meta.data)
head(obj_barcodes)
all(obj_barcodes == anno$cell_id)

# Idents()로 cell type 지정
obj$azimuth_type <- anno$predicted.celltype.l2
Idents(obj) <- obj@meta.data$azimuth_type

library(ggplot2)

pdf('Azimuth_annotation_l2.pdf', width=8, height=8)
DimPlot(obj, group.by="azimuth_type", label=TRUE, repel=TRUE)+
  theme(legend.text = element_text(size=7),     # 텍스트 크기
        legend.key.size = unit(0.5, "lines"))+
        coord_fixed(ratio=1)
dev.off()



####개별 cluster plot
# Seurat object 불러오기
obj <- readRDS("integrated_SCT_with_markers.rds")
anno <- read.table("Annotation/azimuth/azimuth_predI2.tsv", header=TRUE, sep="\t")  # cell_id, predicted.celltype.l2 등

obj_barcodes <- rownames(obj@meta.data)
all(obj_barcodes == anno$cell_id)  # TRUE여야 바로 붙일 수 있음

obj@meta.data$azimuth_type <- anno$predicted.celltype.l2

df <- data.frame(
  UMAP_1 = Embeddings(obj, "umap")[,1],
  UMAP_2 = Embeddings(obj, "umap")[,2],
  label = obj@meta.data$azimuth_type
)

df <- df %>% filter(!is.na(label))
all_labels <- unique(df$label)

pdf("umap_azimuth_highlight_each_label.pdf", width=8, height=6)

for (this_label in all_labels) {
  df$plot_label <- ifelse(df$label == this_label, this_label, "Other")
  mycolors <- c("#1976D2", "#D3D3D3")  # 강조/회색
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
}
dev.off()

