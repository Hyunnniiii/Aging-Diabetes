## =========================================
## Date | 2025-08-29 | Yoon JiHyun

## Description | T/NK/Mono subset에서 res=0.4 클러스터의 그룹별 조성 비율(HY/HO/DO) 계산 후 stack bar plot 확인
## =========================================


suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(readr); library(purrr); library(scales)
})

# ===== 설정 =====
lineages    <- c("T","NK","Mono")
GROUP_COL   <- "group"   # sub RDS에 있는 그룹 컬럼명
grp_lvls    <- c("Healthy_Young","Healthy_Old","Diabetes_Old")  # 실제 값
x_labels    <- c(Healthy_Young="HY", Healthy_Old="HO", Diabetes_Old="DO")  # 축 라벨(옵션)
CLUSTER_COL <- "integrated_snn_res.0.4"  # ★ res0.4 고정

sub_rds_path <- function(lin) file.path("/data/project/diabetes_LYH/tanya/rds", sprintf("subset_%s_SCTintegrated.rds", lin))
out_dir <- file.path("Subtyping", "Proportions_res0.4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

walk(lineages, function(lin){
  sobj <- readRDS(sub_rds_path(lin))
  md   <- sobj@meta.data

  # 필수 컬럼 체크
  stopifnot(GROUP_COL %in% colnames(md))
  stopifnot(CLUSTER_COL %in% colnames(md))

  # 준비
  md[[GROUP_COL]] <- factor(as.character(md[[GROUP_COL]]), levels = grp_lvls)
  md$cluster      <- as.factor(md[[CLUSTER_COL]])

  # 유효값만 사용
  d <- md %>% filter(!is.na(.data[[GROUP_COL]]), !is.na(cluster)) %>%
       transmute(group = .data[[GROUP_COL]], cluster)

  # 비율(각 group 내에서 합=1)
  by_grp <- d %>%
    count(group, cluster, name = "n") %>%
    complete(group = factor(grp_lvls, levels = grp_lvls), cluster, fill = list(n = 0)) %>%
    group_by(group) %>% mutate(prop = n / sum(n)) %>% ungroup()

  # 스택 순서: 해당 라인리지 전체 빈도 내림차순
  cluster_order <- by_grp %>% group_by(cluster) %>% summarise(total = sum(n), .groups="drop") %>%
                    arrange(desc(total)) %>% pull(cluster)
  by_grp$cluster <- factor(by_grp$cluster, levels = cluster_order)


  # 스택 바(비율)
  p <- ggplot(by_grp, aes(x = group, y = prop, fill = cluster)) +
    geom_col() +
    scale_x_discrete(labels = x_labels) +  # hy/ho/do로 축 표시(원하면 이 줄 제거)
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
    labs(title = sprintf("[%s] Subset (res=0.4) composition within each group", lin),
         x = NULL, y = "Proportion", fill = "subset(res0.4)") +
    theme_classic(base_size = 12)

  ggsave(file.path(out_dir, sprintf("stacked_%s_subset_res0.4_by_group.pdf", lin)), p, width = 7, height = 5)

  message(sprintf("[%s] done → %s", lin, normalizePath(out_dir)))
})
