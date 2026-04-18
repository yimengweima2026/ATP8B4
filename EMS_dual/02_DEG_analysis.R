# =============================================================================
# 02 差异表达基因分析 (DEG Analysis)
# - Limma 筛选 EMS 和 RPL 的 DEGs
# - EMS 采用 1:1 子集法 (5次随机抽样取交集)
# - Venn 图展示 5 个子集共有 DEGs
# - EMS 与 RPL DEGs 取交集
# - pheatmap 热图可视化
# =============================================================================

library(GEOquery)
library(limma)
library(sva)

library(pheatmap)
library(VennDiagram)

# --- 1. 设置路径与参数 ---
base_dir   <- "d:/learnspace/MyRSci"
ems_dir    <- file.path(base_dir, "EMS")
rpl_dir    <- file.path(base_dir, "RPL")
output_dir <- file.path(base_dir, "EMS_dual")

# DEG 筛选阈值
logFC_cutoff   <- 0.585
adjP_cutoff    <- 0.05
n_subsets      <- 5
set.seed(42)

# --- 2. 读取 ComBat 校正后数据 ---
ems_expr <- readRDS(file.path(output_dir, "EMS_combat_expr.rds"))
rpl_expr <- readRDS(file.path(output_dir, "RPL_combat_expr.rds"))
ems_meta <- read.csv(file.path(output_dir, "EMS_sample_meta.csv"),
                     stringsAsFactors = FALSE)
rpl_meta <- read.csv(file.path(output_dir, "RPL_sample_meta.csv"),
                     stringsAsFactors = FALSE)

# --- 3. 构建分组标签: 从原始GEO文件提取表型信息 ---

#' 为 EMS 样本标注 Disease/Control
assign_ems_group <- function(ems_dir, ems_meta) {
  ems_meta$group <- NA_character_

  ems_files <- list.files(ems_dir, pattern = "series_matrix\\.txt\\.gz$",
                          full.names = TRUE)
  ems_files <- ems_files[!grepl("GSE120103", ems_files)]

  for (f in ems_files) {
    gse_name <- gsub("_series_matrix\\.txt\\.gz$", "", basename(f))
    eset <- getGEO(filename = f, getGPL = FALSE)
    pd <- pData(eset)
    gsm_ids <- rownames(pd)

    if (gse_name == "GSE23339") {
      grp <- ifelse(grepl("control", pd$characteristics_ch1, ignore.case = TRUE),
                    "Control", "Disease")
    } else if (gse_name == "GSE51981") {
      col <- pd$`endometriosis/no endometriosis:ch1`
      if (is.null(col)) {
        char_col <- pd$characteristics_ch1.1
        grp <- ifelse(grepl("Non-Endometriosis", char_col), "Control", "Disease")
      } else {
        grp <- ifelse(grepl("Non", col, ignore.case = TRUE), "Control", "Disease")
      }
    } else if (gse_name == "GSE6364") {
      grp <- ifelse(grepl("Normal", pd$characteristics_ch1, ignore.case = TRUE),
                    "Control", "Disease")
    } else if (gse_name == "GSE7305") {
      grp <- ifelse(grepl("Normal", pd$title, ignore.case = TRUE),
                    "Control", "Disease")
    } else if (gse_name == "GSE7846") {
      grp <- ifelse(grepl("without", pd$title, ignore.case = TRUE),
                    "Control", "Disease")
    } else {
      next
    }

    idx <- match(gsm_ids, ems_meta$sample_id)
    idx_valid <- !is.na(idx)
    ems_meta$group[idx[idx_valid]] <- grp[idx_valid]
  }
  return(ems_meta)
}

#' 为 RPL 样本标注 RPL/Control (排除 UIF)
assign_rpl_group <- function(rpl_dir, rpl_meta) {
  rpl_meta$group <- NA_character_

  rpl_files <- list.files(rpl_dir, pattern = "series_matrix\\.txt\\.gz$",
                          full.names = TRUE)

  for (f in rpl_files) {
    gse_name <- gsub("_series_matrix\\.txt\\.gz$", "", basename(f))
    eset <- getGEO(filename = f, getGPL = FALSE)
    pd <- pData(eset)
    gsm_ids <- rownames(pd)

    if (gse_name == "GSE165004") {
      char <- pd$characteristics_ch1
      grp <- ifelse(grepl("Control", char), "Control",
             ifelse(grepl("RPL", char), "Disease", "UIF"))
    } else if (gse_name == "GSE26787") {
      src <- pd$source_name_ch1
      grp <- ifelse(grepl("fertile", src, ignore.case = TRUE), "Control",
             ifelse(grepl("recurrent", src, ignore.case = TRUE), "Disease", "UIF"))
    } else {
      next
    }

    idx <- match(gsm_ids, rpl_meta$sample_id)
    idx_valid <- !is.na(idx)
    rpl_meta$group[idx[idx_valid]] <- grp[idx_valid]
  }
  return(rpl_meta)
}

cat("正在构建分组标签...\n")
ems_meta <- assign_ems_group(ems_dir, ems_meta)
rpl_meta <- assign_rpl_group(rpl_dir, rpl_meta)

cat("\nEMS 分组:\n")
print(table(ems_meta$group, useNA = "ifany"))
cat("\nRPL 分组 (含UIF):\n")
print(table(rpl_meta$group, useNA = "ifany"))

# --- 4. Limma 差异分析函数 ---

run_limma_deg <- function(expr_mat, group_vector, logFC_cut, adjP_cut) {
  group <- factor(group_vector, levels = c("Control", "Disease"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  contrast_mat <- makeContrasts(Disease - Control, levels = design)

  fit <- lmFit(expr_mat, design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)

  results <- topTable(fit2, number = Inf, sort.by = "none")
  results$gene <- rownames(results)

  # 筛选 DEGs
  degs <- results[abs(results$logFC) > logFC_cut &
                    results$adj.P.Val < adjP_cut, ]

  return(list(full_results = results, degs = degs))
}

# --- 5. EMS: 1:1 子集法 ---
cat("\n", strrep("=", 60), "\n")
cat("EMS 差异表达分析 (1:1 子集法, 共", n_subsets, "次)\n")
cat(strrep("=", 60), "\n")

ems_disease_idx <- which(ems_meta$group == "Disease")
ems_control_idx <- which(ems_meta$group == "Control")

n_disease <- length(ems_disease_idx)
n_control <- length(ems_control_idx)
n_min     <- min(n_disease, n_control)

cat("Disease样本:", n_disease, "| Control样本:", n_control,
    "| 每次子集取:", n_min, ": ", n_min, "\n\n")

# 确定多数组和少数组
if (n_disease > n_control) {
  majority_idx <- ems_disease_idx
  minority_idx <- ems_control_idx
  majority_label <- "Disease"
} else {
  majority_idx <- ems_control_idx
  minority_idx <- ems_disease_idx
  majority_label <- "Control"
}

ems_subset_degs <- list()

for (i in seq_len(n_subsets)) {
  cat(">>> 子集", i, "\n")

  # 从多数组中随机抽取 n_min 个样本
  sampled_majority <- sample(majority_idx, n_min, replace = FALSE)

  # 子集 = 全部少数组 + 抽取的多数组
  subset_idx <- c(minority_idx, sampled_majority)
  subset_expr <- ems_expr[, subset_idx]
  subset_group <- ems_meta$group[subset_idx]

  cat("  样本数:", ncol(subset_expr),
      " (Disease:", sum(subset_group == "Disease"),
      ", Control:", sum(subset_group == "Control"), ")\n")

  res <- run_limma_deg(subset_expr, subset_group, logFC_cutoff, adjP_cutoff)
  deg_genes <- res$degs$gene
  ems_subset_degs[[paste0("Subset", i)]] <- deg_genes
  cat("  DEGs数量:", length(deg_genes), "\n\n")

  # 保存第一个子集的完整结果供后续使用
  if (i == 1) ems_full_results <- res$full_results
}

# 五个子集取交集
ems_common_degs <- Reduce(intersect, ems_subset_degs)
cat("五个子集共有 DEGs (交集):", length(ems_common_degs), "\n")

# --- 6. Venn 图: EMS 5 个子集 ---
cat("\n绘制 Venn 图...\n")

venn_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

venn.diagram(
  x = ems_subset_degs,
  category.names = paste0("Subset ", 1:n_subsets),
  filename = file.path(output_dir, "EMS_venn_5subsets.tiff"),
  output = TRUE,
  imagetype = "tiff",
  height = 3000,
  width = 3600,
  resolution = 300,
  lwd = 2,
  fill = venn_colors,
  alpha = 0.3,
  cat.cex = 1.0,
  cex = 1.0,
  main = "EMS DEGs across 5 subsets",
  main.cex = 1.5
)
cat("Venn图已保存: EMS_venn_5subsets.tiff\n")

# --- 7. RPL: 标准 Limma 分析 ---
cat("\n", strrep("=", 60), "\n")
cat("RPL 差异表达分析\n")
cat(strrep("=", 60), "\n")

# 排除 UIF 样本，仅保留 RPL (Disease) 和 Control
rpl_keep <- rpl_meta$group %in% c("Disease", "Control")
rpl_expr_sub <- rpl_expr[, rpl_keep]
rpl_meta_sub <- rpl_meta[rpl_keep, ]

cat("RPL 分析样本: Disease =", sum(rpl_meta_sub$group == "Disease"),
    ", Control =", sum(rpl_meta_sub$group == "Control"), "\n")

rpl_res <- run_limma_deg(rpl_expr_sub, rpl_meta_sub$group,
                         logFC_cutoff, adjP_cutoff)
rpl_degs <- rpl_res$degs$gene
cat("RPL DEGs数量:", length(rpl_degs), "\n")

# --- 8. EMS 与 RPL DEGs 取交集 ---
cat("\n", strrep("=", 60), "\n")
cat("EMS 与 RPL DEGs 交集\n")
cat(strrep("=", 60), "\n")

common_degs <- intersect(ems_common_degs, rpl_degs)
cat("EMS 共有DEGs:", length(ems_common_degs), "\n")
cat("RPL DEGs:", length(rpl_degs), "\n")
cat("交集 DEGs:", length(common_degs), "\n")

if (length(common_degs) > 0) {
  cat("交集基因:\n")
  cat(paste(sort(common_degs), collapse = ", "), "\n")
}

# Venn 图: EMS vs RPL
venn.diagram(
  x = list(EMS = ems_common_degs, RPL = rpl_degs),
  filename = file.path(output_dir, "EMS_RPL_venn.tiff"),
  output = TRUE,
  imagetype = "tiff",
  height = 2400,
  width = 2800,
  resolution = 300,
  lwd = 2,
  fill = c("#E41A1C", "#377EB8"),
  alpha = 0.3,
  cat.cex = 1.4,
  cex = 1.5,
  main = "DEGs: EMS vs RPL",
  main.cex = 1.5
)
cat("Venn图已保存: EMS_RPL_venn.tiff\n")

# --- 9. pheatmap 热图: 交集 DEGs ---
cat("\n绘制热图...\n")

if (length(common_degs) > 0) {

  # 合并 EMS 和 RPL 表达矩阵 (取共同基因)
  common_all_genes <- intersect(rownames(ems_expr), rownames(rpl_expr))

  # drop=FALSE 确保单基因时仍为矩阵
  # 单行时 scale="row" 会产生NaN，改用 scale="none" 并手动标准化
  use_scale <- ifelse(length(common_degs) > 1, "row", "none")

  # (a) EMS 热图
  ems_heatmap_expr <- ems_expr[common_degs, , drop = FALSE]
  ems_anno <- data.frame(
    Group = factor(ems_meta$group, levels = c("Control", "Disease")),
    Dataset = factor(ems_meta$dataset),
    row.names = ems_meta$sample_id
  )

  pheatmap(
    ems_heatmap_expr,
    scale = use_scale,
    cluster_rows = length(common_degs) > 1,
    clustering_method = "ward.D2",
    show_colnames = FALSE,
    annotation_col = ems_anno,
    annotation_colors = list(
      Group = c(Control = "#377EB8", Disease = "#E41A1C"),
      Dataset = c(GSE23339 = "#66C2A5", GSE51981 = "#FC8D62",
                  GSE6364 = "#8DA0CB", GSE7305 = "#E78AC3",
                  GSE7846 = "#A6D854")
    ),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = paste0("EMS: Common DEGs Expression (n=", length(common_degs), ")"),
    fontsize_row = max(4, min(10, 200 / length(common_degs))),
    filename = file.path(output_dir, "heatmap_EMS_common_degs.pdf"),
    width = 12, height = max(4, length(common_degs) * 0.15 + 2)
  )
  cat("EMS 热图已保存: heatmap_EMS_common_degs.pdf\n")

  # (b) RPL 热图 (仅保留 Disease/Control 样本)
  rpl_heatmap_expr <- rpl_expr_sub[common_degs, , drop = FALSE]
  rpl_anno <- data.frame(
    Group = factor(rpl_meta_sub$group, levels = c("Control", "Disease")),
    Dataset = factor(rpl_meta_sub$dataset),
    row.names = rpl_meta_sub$sample_id
  )

  pheatmap(
    rpl_heatmap_expr,
    scale = use_scale,
    cluster_rows = length(common_degs) > 1,
    clustering_method = "ward.D2",
    show_colnames = FALSE,
    annotation_col = rpl_anno,
    annotation_colors = list(
      Group = c(Control = "#377EB8", Disease = "#E41A1C"),
      Dataset = c(GSE165004 = "#66C2A5", GSE26787 = "#FC8D62")
    ),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = paste0("RPL: Common DEGs Expression (n=", length(common_degs), ")"),
    fontsize_row = max(4, min(10, 200 / length(common_degs))),
    filename = file.path(output_dir, "heatmap_RPL_common_degs.pdf"),
    width = 10, height = max(4, length(common_degs) * 0.15 + 2)
  )
  cat("RPL 热图已保存: heatmap_RPL_common_degs.pdf\n")

} else {
  cat("警告: EMS与RPL没有共同DEGs，无法绘制热图\n")
}

# --- 10. 保存结果 ---
cat("\n保存分析结果...\n")

# 保存 EMS DEGs 完整列表
write.csv(ems_full_results,
          file = file.path(output_dir, "EMS_limma_full_results.csv"),
          row.names = FALSE)

# 保存 RPL DEGs 完整结果
write.csv(rpl_res$full_results,
          file = file.path(output_dir, "RPL_limma_full_results.csv"),
          row.names = FALSE)

# 保存各子集DEGs列表
ems_subset_df <- data.frame(
  gene = ems_common_degs,
  in_EMS = TRUE,
  in_RPL = ems_common_degs %in% rpl_degs,
  stringsAsFactors = FALSE
)
write.csv(ems_subset_df,
          file = file.path(output_dir, "EMS_common_DEGs_5subsets.csv"),
          row.names = FALSE)

# 保存交集DEGs
if (length(common_degs) > 0) {
  # 提取logFC和adjP信息
  ems_info <- ems_full_results[ems_full_results$gene %in% common_degs,
                               c("gene", "logFC", "adj.P.Val")]
  colnames(ems_info) <- c("gene", "EMS_logFC", "EMS_adjP")

  rpl_info <- rpl_res$full_results[rpl_res$full_results$gene %in% common_degs,
                                    c("gene", "logFC", "adj.P.Val")]
  colnames(rpl_info) <- c("gene", "RPL_logFC", "RPL_adjP")

  common_df <- merge(ems_info, rpl_info, by = "gene")
  common_df <- common_df[order(common_df$EMS_adjP), ]
  write.csv(common_df,
            file = file.path(output_dir, "EMS_RPL_common_DEGs.csv"),
            row.names = FALSE)
}

# 保存带分组标签的元数据
write.csv(ems_meta, file = file.path(output_dir, "EMS_sample_meta_grouped.csv"),
          row.names = FALSE)
write.csv(rpl_meta, file = file.path(output_dir, "RPL_sample_meta_grouped.csv"),
          row.names = FALSE)

# --- 11. 汇总 ---
cat("\n", strrep("=", 60), "\n")
cat("分析完成汇总\n")
cat(strrep("=", 60), "\n")
cat("\n筛选阈值: |logFC| >", logFC_cutoff, ", adj.P.Val <", adjP_cutoff, "\n")
cat("\nEMS (1:1 子集法):\n")
for (i in seq_len(n_subsets)) {
  cat("  子集", i, ":", length(ems_subset_degs[[i]]), "DEGs\n")
}
cat("  五子集交集:", length(ems_common_degs), "DEGs\n")
cat("\nRPL:", length(rpl_degs), "DEGs\n")
cat("\nEMS ∩ RPL:", length(common_degs), "共同DEGs\n")

cat("\n输出文件:\n")
cat("  EMS_venn_5subsets.tiff       -> EMS 5子集 Venn图\n")
cat("  EMS_RPL_venn.tiff            -> EMS vs RPL Venn图\n")
cat("  heatmap_EMS_common_degs.pdf  -> EMS 交集DEGs热图\n")
cat("  heatmap_RPL_common_degs.pdf  -> RPL 交集DEGs热图\n")
cat("  EMS_limma_full_results.csv   -> EMS limma完整结果\n")
cat("  RPL_limma_full_results.csv   -> RPL limma完整结果\n")
cat("  EMS_common_DEGs_5subsets.csv -> EMS 5子集共有DEGs\n")
cat("  EMS_RPL_common_DEGs.csv      -> 共同DEGs详情\n")
