suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
})

# ---- 설정 ----
alpha <- 0.05
fc_thr <- 1.5
log2_thr <- log2(fc_thr)   # ≈ 0.5849625
ln_thr   <- log(fc_thr)    # ≈ 0.4054651

# p_val_adj가 없으면 BH로 보정해서 추가
ensure_fdr <- function(df) {
  if (!"p_val_adj" %in% names(df)) {
    if ("p_val" %in% names(df)) {
      df$p_val_adj <- p.adjust(df$p_val, method = "BH")
    } else {
      stop("p_val_adj도 p_val도 없습니다.")
    }
  }
  df
}

# logFC 컬럼(Seurat 버전별로 다름) 정규화: 항상 log2 스케일로 비교
extract_log2fc <- function(df) {
  if ("avg_log2FC" %in% names(df)) {
    df$log2fc_std <- df$avg_log2FC
  } else if ("avg_logFC" %in% names(df)) {
    # avg_logFC는 보통 natural log(ln) 스케일 → log2로 변환
    df$log2fc_std <- df$avg_logFC / log(2)
  } else {
    stop("avg_log2FC/avg_logFC 컬럼을 찾지 못했습니다.")
  }
  df
}

# MT 유전자 제거 (사람: 'MT-'; 마우스 데이터면 'mt-' 대응)
drop_mt <- function(df, gene_col = "gene") {
  df %>% filter(!str_detect(.data[[gene_col]], "^(MT-|mt-)"))
}

# up/down 카운트 및 리스트 뽑기
count_updown <- function(df, label) {
  df <- df %>% ensure_fdr() %>% extract_log2fc() %>% drop_mt("gene")

  keep <- df %>% filter(!is.na(p_val_adj) & p_val_adj < alpha)
  up   <- keep %>% filter(log2fc_std >= log2_thr)
  down <- keep %>% filter(log2fc_std <= -log2_thr)

  # 저장: 유전자 리스트
  dir.create("MAST_DEG", showWarnings = FALSE)
  write_tsv(up  %>% select(gene, p_val_adj, log2fc_std),
            file = file.path("MAST_DEG", paste0(label, "_UP_genes.tsv")))
  write_tsv(down %>% select(gene, p_val_adj, log2fc_std),
            file = file.path("MAST_DEG", paste0(label, "_DOWN_genes.tsv")))

  tibble(
    contrast = label,
    n_up   = nrow(up),
    n_down = nrow(down),
    n_sig  = nrow(keep)
  )
}

# ---- 입력: RDS로 계산 (TSV로도 가능) ----
res_HY_DO <- readRDS("/data/project/diabetes_LYH/tanya/rds/MAST/HY_vs_DO_MAST.rds")
res_HY_HO <- readRDS("/data/project/diabetes_LYH/tanya/rds/MAST/HY_vs_HO_MAST.rds")
res_HO_DO <- readRDS("/data/project/diabetes_LYH/tanya/rds/MAST/HO_vs_DO_MAST.rds")

# ---- 집계 ----
sum_HY_DO <- count_updown(res_HY_DO, "HY_vs_DO")
sum_HY_HO <- count_updown(res_HY_HO, "HY_vs_HO")
sum_HO_DO <- count_updown(res_HO_DO, "HO_vs_DO")

summary_tbl <- bind_rows(sum_HY_DO, sum_HY_HO, sum_HO_DO) %>%
  arrange(contrast)

print(summary_tbl)
write_tsv(summary_tbl, "MAST_DEG/DEG_count_summary.tsv")
