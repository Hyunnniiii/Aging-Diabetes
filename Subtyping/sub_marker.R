suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

lineages    <- c("T","NK","Mono")
sub_root    <- "Subtyping"                         # subset_{LIN}_SCTintegrated.rds 위치
rds_name    <- "subset_%s_SCTintegrated.rds"
out_root    <- "Markers"; dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
CLUSTER_COL <- "integrated_snn_res.0.4"
RES_TAG <- "res0.4"                


for (lin in lineages) {
  message(sprintf("\n[%s] ---", lin))
  obj <- readRDS(file.path(sub_root, lin, sprintf(rds_name, lin)))

  # 1) 클러스터 지정 (integrated_snn_res.0.4)
  stopifnot(CLUSTER_COL %in% colnames(obj@meta.data))
  Idents(obj) <- obj[[CLUSTER_COL]][,1]

  ## 2) DE는 'RNA'로 (전체 유전자). v5 Assay5라면 data 레이어 생성
    if (any(grepl("^counts\\.", Layers(obj[["RNA"]])))) {
    obj <- JoinLayers(obj, assay = "RNA")}
    obj <- NormalizeData(obj, assay = "RNA", layer = "counts", verbose = FALSE)

   # ★ 3) MT 유전자 제외한 feature 목록
  all_genes <- rownames(GetAssayData(obj, assay = "RNA", slot = "data"))
  genes_use <- all_genes[!grepl("^(MT-|MTRNR|mt-|mt-Rnr)", all_genes, ignore.case = TRUE)]

  # 3) 클러스터별 마커 계산 (양의 마커만)
  markers <- FindAllMarkers(
    obj, assay = "RNA", slot = "data",
    features = genes_use,
    only.pos        = TRUE,
    test.use        = "wilcox",
    min.pct         = 0.10,
    logfc.threshold = 0.25
  ) %>% arrange(cluster, desc(avg_log2FC))

  # 4) 저장 (전체 + 클러스터별 top10)
  write.table(markers,
              file = file.path(out_root, sprintf("markers_%s_%s.tsv", lin, RES_TAG)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  top10 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
  write.table(top10,
              file = file.path(out_root, sprintf("markers_%s_%s_top10.tsv", lin, RES_TAG)),
              sep = "\t", quote = FALSE, row.names = FALSE)

 # DotPlot으로 확인
  top5 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
    arrange(cluster, desc(avg_log2FC))


    pdf(file.path(out_root, sprintf("DotPlot_%s_%s_top5.pdf", lin, RES_TAG)), width = 10, height = 6)
    print(DotPlot(obj, features = unique(top5$gene), assay = "RNA") + RotatedAxis())
    dev.off()  
  
  message(sprintf("[%s] clusters(%s): %s", lin, RES_TAG, paste(levels(Idents(obj)), collapse = ", ")))
}
message("\nDone.")


#### [cluster: GENE] 형식의 파일 생성
suppressPackageStartupMessages(library(dplyr))

lineages <- c("T","NK","Mono")
RES_TAG  <- "res0.4"
in_root  <- "Markers"

for (lin in lineages) {
  m <- read.delim(file.path(in_root, sprintf("markers_%s_%s.tsv", lin, RES_TAG)))
  lines <- m %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    summarise(genes = paste(gene, collapse = ", "), .groups = "drop") %>%
    arrange(as.numeric(cluster))
  writeLines(paste0(lines$cluster, ": ", lines$genes),
             file.path(in_root, sprintf("cluster_genes_%s_%s_top10.txt", lin, RES_TAG)))
}
