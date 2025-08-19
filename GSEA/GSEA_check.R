# ========================================================
# GSEA 결과로 Figure 생성 (Reactome/GO:BP 분리)
# 색상 규칙:
#   - Diverging barplot:  Up=RED (#d73027), Down=BLUE (#4575b4)
#   - NES Heatmap:        NEGATIVE=BLUE → 0=GREY → POSITIVE=RED
#   - Leading-edge bars:  비교별 1장 (내부에 Reactome/GOBP × Up/Down 패널)
# 입력: out_dir 폴더의 *__GSEA_result.tsv
# 출력: fig_*.pdf (fig_dir 폴더)
# ========================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr); library(tibble)
  library(ggplot2); library(pheatmap); library(patchwork)
})

# ----- 설정 -----
out_dir  <- "GSEA_results_Reactome_GO/result"  # __GSEA_result.tsv들이 있는 폴더
fig_dir  <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

padj_cut_heatmap <- 0.25   # 히트맵 선정 기준(탐색용 넉넉하게)
padj_cut_bar     <- 0.05   # 바차트/리딩엣지 기준(보고용 보수)
top_terms_each   <- 20     # 히트맵 대표 용어 수(컬렉션별)
topN_bar_each    <- 20     # diverging barplot 총 개수(Up/Down 합→각 10씩)
top_genes_each   <- 20     # leading-edge 막대: 방향(Up/Down)별 상위 gene 수
wrap_labels_width <- 34    # 긴 용어 줄바꿈 폭(겹침 완화)

# ----- 유틸 -----
parse_AB <- function(bn){
  m <- stringr::str_match(bn, "(HY|HO|DO)_vs_(HY|HO|DO)")
  tibble(A = m[,2], B = m[,3], contrast = paste0(m[,2], "_vs_", m[,3]))
}
collection_from <- function(bn){
  tag <- stringr::str_match(bn, "__(C2REACTOME|C5BP)__")[,2]
  dplyr::recode(tag, "C2REACTOME"="Reactome", "C5BP"="GO:BP")
}
pretty_term <- function(x){
  x <- gsub("^(GOBP_|GO:BP_|GO_BP_|REACTOME_|R-HSA-\\d+_|R-HSA-)", "", x)
  stringr::str_to_sentence(gsub("_", " ", x))
}
read_cp_gsea <- function(fp){
  df <- readr::read_tsv(fp, show_col_types = FALSE)
  if (!"Description" %in% names(df))      df$Description      <- df$ID
  if (!"core_enrichment" %in% names(df))  df$core_enrichment  <- NA_character_
  if (!"setSize" %in% names(df))          df$setSize          <- NA_integer_
  df
}

# 라벨 겹침 최소화된 diverging barplot
make_diverging <- function(tbl, title, topn = 20, padj_cut = 0.05, wrap = 34){
  stopifnot(all(c("NES","padj") %in% names(tbl)))
  base <- dplyr::filter(tbl, !is.na(padj))
  sig  <- dplyr::filter(base, padj < padj_cut)
  pool <- if (nrow(sig) > 0) sig else base

  n_each <- max(1, floor(topn/2))
  up   <- pool %>% dplyr::filter(NES >= 0) %>% dplyr::arrange(padj, dplyr::desc(NES)) %>% dplyr::slice_head(n = n_each)
  down <- pool %>% dplyr::filter(NES < 0)  %>% dplyr::arrange(padj, NES)             %>% dplyr::slice_head(n = n_each)
  sel  <- dplyr::bind_rows(down, up)
  if (nrow(sel) == 0) return(NULL)

  lab <- sel$term %>% gsub("^(GOBP_|GO:BP_|GO_BP_|REACTOME_|R-HSA-\\d+_|R-HSA-)", "", .)
  lab <- stringr::str_to_sentence(gsub("_", " ", lab))
  lab <- stringr::str_wrap(lab, width = wrap)

  sel$label <- lab
  sel$dir   <- ifelse(sel$NES >= 0, "Up", "Down")
  lim <- max(abs(sel$NES), na.rm = TRUE)

  sub_txt <- if (nrow(sig) > 0) {
    paste0("Top ", nrow(sel), " (padj<", padj_cut, ")")
  } else {
    paste0("Top ", nrow(sel), " (no padj<", padj_cut, "; best available)")
  }

  p <- ggplot2::ggplot(sel, ggplot2::aes(x = reorder(label, NES), y = NES, fill = dir)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = c(0.05, 0.08))) +
    ggplot2::scale_fill_manual(values = c("Down" = "#4575b4", "Up" = "#d73027")) +
    ggplot2::labs(x = NULL, y = "NES", title = title, subtitle = sub_txt) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.95),
      plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5)
    )

  list(plot = p, n = nrow(sel))
}

# ----- 데이터 읽기 & 공통 테이블 -----
files <- list.files(out_dir, pattern="__GSEA_result\\.tsv$", full.names=TRUE)
stopifnot(length(files) > 0)

all_res <- dplyr::bind_rows(lapply(files, function(fp){
  bn  <- basename(fp)
  ab  <- parse_AB(bn)
  col <- collection_from(bn)
  df  <- read_cp_gsea(fp)
  df %>%
    dplyr::transmute(
      contrast   = ab$contrast,
      collection = col,
      term_id    = ID,
      term_raw   = Description,
      term       = pretty_term(Description),
      NES        = NES,
      pvalue     = pvalue,
      padj       = p.adjust,
      setSize    = setSize,
      core_enrichment = core_enrichment,
      side_up    = ifelse(NES >= 0, ab$A, ab$B),
      dir        = ifelse(NES >= 0, "Up", "Down")
    )
}))

# ========================================================
# 1) NES Heatmap (Reactome/GO:BP 각각)
# ========================================================
make_heatmap <- function(tag = "Reactome", top_each = 20, padj_cut = 0.25){
  dfc <- dplyr::filter(all_res, collection == tag)

  picks <- dfc %>%
    dplyr::filter(!is.na(padj) & padj < padj_cut) %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(best = min(padj, na.rm = TRUE), .groups="drop") %>%
    dplyr::arrange(best) %>%
    dplyr::slice_head(n = top_each)

  if (nrow(picks) == 0) {
    message("[heatmap] ", tag, " : padj<", padj_cut, " 항목이 없어 스킵")
    return(invisible(NULL))
  }

  mat <- dfc %>%
    dplyr::semi_join(picks, by = "term") %>%
    dplyr::select(term, contrast, NES) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = contrast, values_from = NES)

  mm <- as.data.frame(mat)
  rownames(mm) <- mm$term
  mm$term <- NULL
  m <- as.matrix(mm)

  # ★ NEGATIVE→BLUE, 0→GREY, POSITIVE→RED
  pal <- colorRampPalette(c("#4575b4","#f7f7f7","#d73027"))(101)

  ofn <- file.path(fig_dir, paste0("fig_overview_heatmap_", ifelse(tag=="GO:BP","GOBP","Reactome"), ".pdf"))
  pheatmap::pheatmap(m, color = pal, cluster_rows = TRUE, cluster_cols = TRUE,
                     border_color = NA, fontsize = 9,
                     filename = ofn, width = 7.5, height = max(6, nrow(m)*0.24))
  message("✔ saved: ", ofn)
}

make_heatmap("Reactome", top_each = top_terms_each, padj_cut = padj_cut_heatmap)
make_heatmap("GO:BP",   top_each = top_terms_each, padj_cut = padj_cut_heatmap)

# ========================================================
# 2) Per-contrast Diverging Barplot (비교별, Reactome & GO:BP 나란히)
# ========================================================
contrasts <- unique(all_res$contrast)

for (ct in contrasts) {
  df_ct <- dplyr::filter(all_res, contrast == ct)

  # Reactome
  dfr <- dplyr::filter(df_ct, collection == "Reactome")
  r   <- if (nrow(dfr)>0) make_diverging(dfr, paste0(ct, " · Reactome"),
                                         topn = topN_bar_each, padj_cut = padj_cut_bar, wrap = wrap_labels_width) else NULL
  # GO:BP
  dfg <- dplyr::filter(df_ct, collection == "GO:BP")
  g   <- if (nrow(dfg)>0) make_diverging(dfg, paste0(ct, " · GO:BP"),
                                         topn = topN_bar_each, padj_cut = padj_cut_bar, wrap = wrap_labels_width) else NULL

  if (is.null(r) && is.null(g)) {
    message("[barplot] ", ct, " : no data, skip")
    next
  }

  plt <- if (!is.null(r) && !is.null(g)) r$plot + g$plot + patchwork::plot_layout(ncol=2) else (r$plot %||% g$plot)

  # 높이: 항목 수 기반(겹침 완화)
  nmax <- max(c(if (!is.null(r)) r$n else 0, if (!is.null(g)) g$n else 0))
  ofn <- file.path(fig_dir, paste0("fig_summary_diverging_", ct, ".pdf"))
  ggsave(ofn, plt, width = 14, height = max(7, 0.5 * nmax))
  message("✔ saved: ", ofn)
}

# ========================================================
# 3) Leading-edge 핵심 유전자 막대 (비교별 1장: 내부에 Reactome/GOBP × Up/Down)
# ========================================================
make_le_bar <- function(dfc, dir_sel = c("Up","Down"), title = "", top_genes = 20){
  dir_sel <- match.arg(dir_sel)
  x <- dfc %>% dplyr::filter(dir == dir_sel, !is.na(core_enrichment), padj < padj_cut_bar)
  if (nrow(x) == 0) return(list(plot=NULL, n=0))
  lead <- x %>%
    dplyr::transmute(gene = strsplit(core_enrichment, "/", fixed = TRUE)) %>%
    tidyr::unnest(gene) %>%
    dplyr::filter(gene != "") %>%
    dplyr::count(gene, sort = TRUE) %>%
    dplyr::slice_head(n = top_genes) %>%
    dplyr::mutate(gene = factor(gene, levels = rev(gene)))
  col_fill <- if (dir_sel == "Up") "#d73027" else "#4575b4"
  p <- ggplot2::ggplot(lead, ggplot2::aes(x = gene, y = n)) +
    ggplot2::geom_col(fill = col_fill) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::labs(x = NULL, y = "Count in leading-edge",
                  title = paste0(title, " · ", dir_sel)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.95),
      plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5)
    )
  list(plot = p, n = nrow(lead))
}

for (ct in contrasts) {
  df_ct <- dplyr::filter(all_res, contrast == ct)

  # Reactome (Up/Down)
  dfR <- dplyr::filter(df_ct, collection == "Reactome")
  RU  <- make_le_bar(dfR, "Up",   title = "Reactome", top_genes = top_genes_each)
  RD  <- make_le_bar(dfR, "Down", title = "Reactome", top_genes = top_genes_each)

  # GOBP (Up/Down)
  dfG <- dplyr::filter(df_ct, collection == "GO:BP")
  GU  <- make_le_bar(dfG, "Up",   title = "GO:BP", top_genes = top_genes_each)
  GD  <- make_le_bar(dfG, "Down", title = "GO:BP", top_genes = top_genes_each)

  # 패널 모으기 (존재하는 것만)
  panels <- list(RU$plot, RD$plot, GU$plot, GD$plot)
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) == 0) { message("[leading-edge] ", ct, " : no padj<", padj_cut_bar, " terms"); next }

  # 2x2(또는 1x2/1x1 자동) 배치
  plt <- patchwork::wrap_plots(panels, ncol = 2)

  # 높이: 각 패널 항목 수 고려(겹침 완화)
  nmax <- max(c(RU$n, RD$n, GU$n, GD$n, 0))
  nrow_panels <- ceiling(length(panels) / 2)
  ofn <- file.path(fig_dir, paste0("fig_leading_edge_genes_", ct, ".pdf"))
  ggsave(ofn, plt, width = 14, height = max(6, 0.34 * nmax * nrow_panels))
  message("✔ saved: ", ofn)
}

message("\n=== Figures saved in: ", normalizePath(fig_dir), " ===")
