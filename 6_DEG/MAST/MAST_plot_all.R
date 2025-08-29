# ============================================================
## Date | 2025-08-29 | Yoon JiHyun
## Harmony로 추가 보정한 결과 파일로 이후 분석 진행

## Description | HY vs HO vs DO의 MAST RDS 결과를 불러와 유의 유전자를 Volcano plot으로 확인 및 저장
# ============================================================


##### 전체 cell에서의 DEG volcano plot
suppressPackageStartupMessages({library(readr)
                                library(ggplot2)
                                library(ggrepel)})

dir <- "/home/tanya0721/PROJECT/DEG/MAST"          # 파일이 있는 폴더
fc_thr <- log2(1.5); p_thr <- 0.05; pct_thr <- 0.10
fdr_log10_thr <- -log10(p_thr)

# 로더 (RDS만)
load_one <- function(pfx, dir="/home/tanya0721/PROJECT/DEG/MAST", rm_mt=TRUE) {
  rds <- file.path(dir, paste0(pfx, "_MAST.rds"))
  if (!file.exists(rds)) stop("파일 없음: ", rds)
  df <- readRDS(rds)
  if (!"gene" %in% names(df)) { df <- cbind(gene = rownames(df), df); rownames(df) <- NULL }
  if (rm_mt) df <- df[!grepl("^(MT-|MTRNR|mt-|mt-Rnr)", df$gene), , drop = FALSE]  # 인간/마우스 둘 다 커버
  df
}

  # logFC 컬럼명 찾기
logfc_name <- function(df) if ("avg_log2FC" %in% names(df)) "avg_log2FC" else "avg_logFC"

# 볼케이노(PDF) 함수
make_volcano_pdf <- function(df, title, outfile_pdf) {
  lfc_col <- logfc_name(df)
  dd <- df
  dd$lfc <- dd[[lfc_col]]
  dd$neglog10padj <- -log10(dd$p_val_adj + 1e-300)

  # min.pct: 두 집단 중 하나라도 10% 이상 검출
  keep <- (dd$pct.1 >= pct_thr) | (dd$pct.2 >= pct_thr)
  
  dd$status <- ifelse(keep & dd$neglog10padj >= fdr_log10_thr & abs(dd$lfc) >= fc_thr,
                    ifelse(dd$lfc > 0, "Up", "Down"), "NS")
  dd$label  <- ifelse(dd$status != "NS", dd$gene, "") 

  p <- ggplot(dd, aes(x = lfc, y = neglog10padj, color = status)) +
    geom_point(size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed") +
    geom_hline(yintercept = fdr_log10_thr, linetype = "dashed") +    
    scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey70")) +
    labs(
      title = sprintf("%s (FC≥1.5, FDR<0.05, min.pct≥10%%)", title),
      x = "log2 fold change", y = "-log10(adj p-value)", color = ""
    ) +
    theme_classic(base_size = 12)

  if (any(nzchar(dd$label))) {
    p <- p + ggrepel::geom_text_repel(
      data = dd[nzchar(dd$label), ],
      aes(label = label),
      size = 2.8, max.overlaps = Inf, box.padding = 0.25, point.padding = 0.1
    )
  }
  ggsave(outfile_pdf, p, width = 7, height = 5, dpi = 200)  # 확장자 .pdf면 PDF로 저장
  message("Saved: ", outfile_pdf)

  invisible(dd)
}

# ----- 개별 파일 로드 -----
res_HY_DO <- load_one("HY_vs_DO", dir)
res_HY_HO <- load_one("HY_vs_HO", dir)
res_HO_DO <- load_one("HO_vs_DO", dir)

# ----- 각 비교별 PDF 저장 (새 폴더 생성 안 함) -----
d1 <- make_volcano_pdf(res_HY_DO, "HY vs DO (MAST)", file.path(dir, "HY_vs_DO_volcano.pdf"))
d2 <- make_volcano_pdf(res_HY_HO, "HY vs HO (MAST)", file.path(dir, "HY_vs_HO_volcano.pdf"))
d3 <- make_volcano_pdf(res_HO_DO, "HO vs DO (MAST)", file.path(dir, "HO_vs_DO_volcano.pdf"))
