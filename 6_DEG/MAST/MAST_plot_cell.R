# ============================================================
## Date | 2025-08-29 | Yoon JiHyun
## Harmony로 추가 보정한 결과 파일로 이후 분석 진행

## Description | celltype별 MAST 결과를 모아 각 대비(HY vs HO vs DO)에 대한 log2FC 분포를 plot으로 확인
# ============================================================


suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(readr)
})

## === 설정 ===
dir <- "/home/tanya0721/PROJECT/DEG/MAST/celltypeDEG"   # RDS가 있는 폴더
fc_thr  <- log2(1.5)    # 0.585
p_thr   <- 0.05
pct_thr <- 0.10
label_top_per_ct <- 20

## === 로더: 파일명에서 celltype/contrast 추출 + log2FC 생성 ===
read_one <- function(f){
  df <- readRDS(f)
  if (!"gene" %in% names(df)) { df <- cbind(gene = rownames(df), df); rownames(df) <- NULL }
  df <- df[!grepl("^(MT-|MTRNR|mt-|mt-Rnr)", df$gene), , drop = FALSE]
  lfc_col <- if ("avg_log2FC" %in% names(df)) "avg_log2FC" else "avg_logFC"
  df$lfc <- df[[lfc_col]]

  base  <- sub("_MAST\\.rds$", "", basename(f))   # "<celltype>__G1_vs_G2"
  parts <- strsplit(base, "__", fixed = TRUE)[[1]]
  df$celltype <- parts[1]
  df$contrast <- parts[2]                          # "Healthy_Young_vs_Diabetes_Old"
  df
}

# celltype별 RDS만 모음(__가 들어간 것)
rds_files <- list.files(dir, pattern = "__.*_MAST\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) > 0)

all_df <- bind_rows(lapply(rds_files, read_one))

make_ct_plot <- function(con){
  d <- all_df %>% dplyr::filter(contrast == con) %>%
    dplyr::mutate(
      keep = (pct.1 >= pct_thr) | (pct.2 >= pct_thr),
      sig  = keep & p_val_adj < p_thr & abs(lfc) >= fc_thr,
      dir  = ifelse(lfc > 0, "Up", "Down")
    )
    
  d$celltype <- factor(d$celltype, levels = unique(d$celltype))
  
    # d <- ... mutate(...) 다음에
    deg <- d %>% dplyr::filter(sig) %>%
        dplyr::mutate(celltype = factor(celltype, levels = unique(d$celltype)))

    lab <- dplyr::bind_rows(
        deg %>% dplyr::filter(lfc > 0) %>%                      # Up: log2FC 큰 순
            dplyr::group_by(celltype) %>%
            dplyr::arrange(dplyr::desc(lfc), .by_group = TRUE) %>%
            dplyr::slice_head(n = 20),
        deg %>% dplyr::filter(lfc < 0) %>%                      # Down: 더 음수(절대값 큼)
            dplyr::group_by(celltype) %>%
            dplyr::arrange(lfc, .by_group = TRUE) %>%
            dplyr::slice_head(n = 20)
            ) %>% dplyr::ungroup()
  


  p <- ggplot(d, aes(x = celltype, y = lfc)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_hline(yintercept = c(-fc_thr, fc_thr), linetype = "dashed", linewidth = 0.3) +
    geom_point(color = "grey70", alpha = 0.5, size = 0.6,
               position = position_jitter(width = 0.22, height = 0)) +
    geom_point(data = deg, aes(color = dir),
           size = 1.2,
           position = position_jitter(width = 0.15, height = 0)) +

    # 모든 유의 DEG 라벨
        # 라벨: 셀타입당 최대 20개만
    { if (nrow(lab)) ggrepel::geom_text_repel(
        data = lab, aes(label = gene, color = dir),
        size = 2.2, max.overlaps = Inf, box.padding = 0.2,
        point.padding = 0.08, segment.size = 0.2
      ) } +

    scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4")) +
    labs(title = gsub("_", " ", sub("^(.+)_vs_(.+)$", "\\1 vs \\2", con)),
         x = "Cell type", y = "log2 Fold Change") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  width <- max(10, length(levels(d$celltype)) * 0.7)
  out <- file.path(dir, paste0("DEG_across_celltypes_", con, "_log2FC.pdf"))
  ggsave(out, p, width = width, height = 6, dpi = 200)
  message("Saved: ", out)
}


# === 3개 대비 각각 1개 PDF 생성 ===
invisible(lapply(unique(all_df$contrast), make_ct_plot))