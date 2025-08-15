library(Seurat)
library(edgeR)
library(limma)
library(Matrix)
library(dplyr)

## 메타데이터 준비 및 배치 병합
seu <- readRDS('/data/project/diabetes_LYH/tanya/rds/SCT&Harmony_annot.rds')
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
meta   <- seu@meta.data

# sample_id 칼럼 정의
meta$sample_id <- meta$patient

# Run2_1과 Run2_2를 Run2로 병합하여 batch 변수 생성
meta$batch <- meta$orig.ident
meta$batch[meta$batch %in% c("Run2_1", "Run2_2")] <- "Run2"
meta$batch <- factor(meta$batch, levels = c("Run1", "Run2"))

meta$sample_id <- interaction(meta$patient, meta$batch, drop = TRUE)  # patient__batch

# 최종 메타데이터 확인 (일부)
head(meta)

# 새로운 배치 정보 확인
table(meta$sample_id, meta$batch)

## Pseudobulk 카운트 집계
# 3. 'sample_id(patient__batch)'별로 세포를 분류하여 pseudobulk 카운트 집계
stopifnot(all(colnames(counts) == rownames(meta)))

idx <- split(seq_len(ncol(counts)), meta$sample_id, drop = TRUE)
pseudobulk_counts <- do.call(
  "cbind",
  lapply(idx, function(i) Matrix::rowSums(counts[, i, drop = FALSE]))
)

#-------------------------------------------------------------------------------
## DEG 분석 시작
# 1. Pseudobulk 샘플에 맞는 메타데이터 생성/정렬
samples_info <- meta |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    group = dplyr::first(group),
    batch = dplyr::first(batch),
    .groups = "drop"
  ) |>
  as.data.frame()

# 열 순서 맞추기
pseudobulk_counts <- pseudobulk_counts[, samples_info$sample_id]
rownames(samples_info) <- samples_info$sample_id

dir.create("deg_out", showWarnings = FALSE)
write.table(table(samples_info$group, samples_info$batch),
            file = "deg_out/group_by_batch.txt", quote = FALSE, col.names = NA)


# edgeR 파이프라인: 필터링/정규화
y <- DGEList(counts = pseudobulk_counts, samples = samples_info)

# 미토 유전자 제거 (대/소문자 무시)
#non_mito <- !grepl("^MT-", rownames(y), ignore.case = TRUE)
#y <- y[non_mito, , keep.lib.sizes = FALSE]

# 저발현 필터 (필요시 min.count 조정)
# keep <- filterByExpr(y, group = y$samples$group, min.count = 10, min.prop = 0.5)
# y <- y[keep, , keep.lib.sizes = FALSE]

# TMM 정규화
y <- calcNormFactors(y)



# 5. 모델 설계: group 효과
design <- model.matrix(~ group, data = y$samples)

# 설계 진단
writeLines(paste("fullrank:", limma::is.fullrank(design)),
           "deg_out/design_fullrank.txt")
ne <- limma::nonEstimable(design)
if (!is.null(ne)) writeLines(ne, "deg_out/design_nonEstimable.txt")
write.table(colnames(design), file = "deg_out/design_cols.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# 6. 분산추정 + QL fit
y <- estimateDisp(y, design, robust = TRUE)

pdf("deg_out/BCV.pdf", width = 7, height = 6)
plotBCV(y)
dev.off()

fit <- glmFit(y, design, robust = TRUE)


# 8. 대비 정의 및 테스트
con_list <- list(
  HO_vs_DO = makeContrasts(groupHealthy_Old, levels = design),                   # HO - DO
  HY_vs_DO = makeContrasts(groupHealthy_Young, levels = design),                 # HY - DO
  HO_vs_HY = makeContrasts(groupHealthy_Old - groupHealthy_Young, levels = design) # HO - HY
)

deg_tables <- list()

pdf("deg_out/MA_plots.pdf", width = 7, height = 6)
for (nm in names(con_list)) {
  lrt <- glmLRT(fit, contrast = con_list[[nm]])
  tt  <- as.data.frame(topTags(lrt, n = Inf))
  tt$Gene <- rownames(tt)

  # MA plot
  plotMD(lrt, main = paste0("MA: ", nm))
  abline(h = c(-1, 1), col = "blue", lty = "dashed")

  deg_tables[[nm]] <- tt
}
dev.off()

## ★ 추가: Run2 내부에서 HO 비교들만 별도 계산 (HO가 Run2에만 존재)
keep_r2 <- y$samples$batch == "Run2" &
           y$samples$group %in% c("Healthy_Old","Healthy_Young","Diabetes_Old")
if (sum(keep_r2) >= 2) {
  y_r2 <- y[, keep_r2, keep.lib.sizes = FALSE]
  design_r2 <- model.matrix(~ 0 + group, data = y_r2$samples)

  y_r2  <- estimateDisp(y_r2, design_r2)
  fit_r2 <- glmFit(y_r2, design_r2)

  # HO vs HY (Run2)
  lrt_HO_HY_r2 <- glmLRT(fit_r2,
    contrast = makeContrasts(groupHealthy_Old - groupHealthy_Young, levels = design_r2))
  tt_HO_HY_r2 <- as.data.frame(topTags(lrt_HO_HY_r2, n = Inf)); tt_HO_HY_r2$Gene <- rownames(tt_HO_HY_r2)
  deg_tables[["HO_vs_HY_Run2_only"]] <- tt_HO_HY_r2
  write.table(tt_HO_HY_r2, "deg_out/DEG_HO_vs_HY_Run2_only_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  # HO vs DO (Run2)
  lrt_HO_DO_r2 <- glmLRT(fit_r2,
    contrast = makeContrasts(groupHealthy_Old - groupDiabetes_Old, levels = design_r2))
  tt_HO_DO_r2 <- as.data.frame(topTags(lrt_HO_DO_r2, n = Inf)); tt_HO_DO_r2$Gene <- rownames(tt_HO_DO_r2)
  deg_tables[["HO_vs_DO_Run2_only"]] <- tt_HO_DO_r2
  write.table(tt_HO_DO_r2, "deg_out/DEG_HO_vs_DO_Run2_only_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

## Volcano 플롯(유의미한 점만 라벨)
pdf("deg_out/Volcano_plots.pdf", width = 7, height = 6)
for (nm in names(deg_tables)) {
  tt <- deg_tables[[nm]]

  tt$Expression <- "Non-DEG"
  tt$Expression[tt$logFC >  0.58 & tt$FDR < 0.05] <- "Up"
  tt$Expression[tt$logFC < -0.58 & tt$FDR < 0.05] <- "Down"

  yvals <- -log10(pmax(tt$FDR, .Machine$double.xmin))

  plot(x = tt$logFC, y = yvals,
       col = ifelse(tt$Expression == "Up", "red",
             ifelse(tt$Expression == "Down", "blue", "gray")),
       pch = 16, cex = 0.7,
       xlab = expression(log[2]~"FC"),
       ylab = expression(-log[10]~"FDR"),
       main = paste0("Volcano: ", nm))
  abline(v = c(-0.58, 0.58), col = "black", lty = "dashed")
  abline(h = -log10(0.05),   col = "black", lty = "dashed")

  if (!"Gene" %in% names(tt)) tt$Gene <- rownames(tt)
  sig <- which(tt$Expression != "Non-DEG" & is.finite(tt$FDR) & tt$FDR > 0)
  if (length(sig)) with(tt[sig, ], text(logFC, -log10(FDR),
                                        labels = Gene, cex = 0.6, pos = 3, offset = 0.3))
}
dev.off()






