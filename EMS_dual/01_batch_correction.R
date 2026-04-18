# =============================================================================
# 批次效应校正脚本 (v2 - 改进数据质量控制)
# 使用 sva::ComBat 消除不同GEO数据集间的批次效应
# 分别针对 EMS (子宫内膜异位症) 和 RPL (复发性流产) 数据
# =============================================================================

# --- 1. 安装和加载所需R包 ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs_bioc <- c("GEOquery", "sva", "limma")
for (pkg in pkgs_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(GEOquery)
library(sva)
library(limma)

# --- 2. 设置路径 ---
base_dir   <- "d:/learnspace/MyRSci"
ems_dir    <- file.path(base_dir, "EMS")
rpl_dir    <- file.path(base_dir, "RPL")
output_dir <- file.path(base_dir, "EMS_dual")

# GSE120103 (GPL6480) 数据为 log2 ratio (中位数=0, 50%负值),
# 与其他单通道绝对表达量数据不兼容，排除
ems_exclude <- c("GSE120103")

# --- 3. 定义辅助函数 ---

#' 从ExpressionSet的fData中提取基因符号
extract_gene_symbols <- function(eset) {
  fd <- fData(eset)
  candidate_cols <- c(
    "Gene symbol", "Gene Symbol", "gene_symbol", "GENE_SYMBOL",
    "Symbol", "SYMBOL", "GeneName", "gene_name", "GENE"
  )
  for (col in candidate_cols) {
    if (col %in% colnames(fd)) {
      cat("    使用注释列:", col, "\n")
      return(as.character(fd[[col]]))
    }
  }
  idx <- grep("symbol|gene.?name", colnames(fd), ignore.case = TRUE)
  if (length(idx) > 0) {
    col_name <- colnames(fd)[idx[1]]
    cat("    使用注释列(模糊匹配):", col_name, "\n")
    return(as.character(fd[[col_name]]))
  }
  stop("无法找到基因符号列。可用列名:\n  ",
       paste(colnames(fd), collapse = "\n  "))
}

#' 判断表达矩阵是否需要log2转换
#' 综合判断 max、median、是否含负值
needs_log2 <- function(expr_mat) {
  val_max    <- max(expr_mat, na.rm = TRUE)
  val_median <- median(expr_mat, na.rm = TRUE)
  pct_neg    <- mean(expr_mat < 0, na.rm = TRUE)

  # 已log2的数据特征: max<30, 可能含少量负值
  # 原始强度数据特征: max>100, 无负值或极少负值
  if (val_max > 100) {
    return(TRUE)   # 最大值过大，必定为原始强度
  }
  if (val_median > 20) {
    return(TRUE)   # 中位数大于20，也是原始强度
  }
  return(FALSE)    # 其余情况视为已在log尺度
}

#' 探针 -> 基因: 多探针取平均表达最高者
probe_to_gene <- function(expr_mat, gene_symbols) {
  gene_symbols <- as.character(gene_symbols)

  valid <- !is.na(gene_symbols) &
    nchar(trimws(gene_symbols)) > 0 &
    gene_symbols != "---"
  expr_mat     <- expr_mat[valid, , drop = FALSE]
  gene_symbols <- gene_symbols[valid]

  # 多基因映射取第一个
  gene_symbols <- sapply(strsplit(gene_symbols, "\\s*///\\s*"), `[`, 1)

  avg_expr <- rowMeans(expr_mat, na.rm = TRUE)
  probe_df <- data.frame(
    probe = rownames(expr_mat),
    gene  = gene_symbols,
    avg   = avg_expr,
    stringsAsFactors = FALSE
  )
  probe_df <- probe_df[order(probe_df$gene, -probe_df$avg), ]
  probe_df <- probe_df[!duplicated(probe_df$gene), ]

  result <- expr_mat[probe_df$probe, , drop = FALSE]
  rownames(result) <- probe_df$gene
  return(result)
}

#' 处理一组GEO数据集: 读取 -> QC -> log2 -> 标准化 -> 注释 ->
#'                     合并 -> 过滤低表达 -> ComBat
process_disease_datasets <- function(data_dir, label, exclude = character(0)) {
  cat(strrep("=", 60), "\n")
  cat("处理", label, "数据集\n")
  cat(strrep("=", 60), "\n")

  files <- list.files(data_dir, pattern = "series_matrix\\.txt\\.gz$",
                      full.names = TRUE)

  # 排除不兼容数据集
  if (length(exclude) > 0) {
    exclude_pattern <- paste0("(", paste(exclude, collapse = "|"), ")")
    excluded_idx <- grep(exclude_pattern, basename(files))
    if (length(excluded_idx) > 0) {
      cat("排除不兼容数据集:", exclude, "\n")
      files <- files[-excluded_idx]
    }
  }
  cat("纳入分析数据集:", length(files), "个\n\n")

  gene_expr_list   <- list()
  batch_list       <- list()
  sample_meta_list <- list()

  for (f in files) {
    gse_name <- gsub("_series_matrix\\.txt\\.gz$", "", basename(f))
    cat(">>> 正在处理:", gse_name, "\n")

    eset <- getGEO(filename = f, getGPL = TRUE, AnnotGPL = TRUE)
    expr_mat <- exprs(eset)
    platform <- annotation(eset)
    n_probes <- nrow(expr_mat)
    n_samples <- ncol(expr_mat)
    cat("  平台:", platform, " | 探针数:", n_probes,
        " | 样本数:", n_samples, "\n")

    # (a) 移除全NA探针
    na_rows <- apply(expr_mat, 1, function(x) all(is.na(x)))
    if (any(na_rows)) {
      cat("  移除", sum(na_rows), "个全NA探针\n")
      expr_mat <- expr_mat[!na_rows, , drop = FALSE]
    }

    # (b) log2转换判断 (改进: 综合max和median)
    if (needs_log2(expr_mat)) {
      cat("  检测到原始强度数据 (max=", round(max(expr_mat, na.rm=TRUE), 1),
          ", median=", round(median(expr_mat, na.rm=TRUE), 1),
          "), 进行log2转换\n")
      expr_mat[expr_mat <= 0] <- min(expr_mat[expr_mat > 0], na.rm = TRUE) / 2
      expr_mat <- log2(expr_mat)
    } else {
      cat("  数据已在log2尺度 (max=", round(max(expr_mat, na.rm=TRUE), 1),
          ", median=", round(median(expr_mat, na.rm=TRUE), 1), ")\n")
    }

    # (c) 数据集内分位数标准化 (消除样本间技术偏差)
    expr_mat <- normalizeBetweenArrays(expr_mat, method = "quantile")
    cat("  完成分位数标准化\n")

    # (d) 过滤低表达探针: 至少在20%样本中表达值 > 全局25分位数
    threshold <- quantile(expr_mat, 0.25, na.rm = TRUE)
    keep <- rowSums(expr_mat > threshold, na.rm = TRUE) >= (n_samples * 0.2)
    n_removed <- sum(!keep)
    expr_mat <- expr_mat[keep, , drop = FALSE]
    cat("  过滤低表达探针:", n_removed, "个移除,",
        nrow(expr_mat), "个保留\n")

    # (e) 探针 -> 基因符号
    gene_symbols <- extract_gene_symbols(eset)
    gene_symbols_filtered <- gene_symbols[names(gene_symbols) %in% rownames(expr_mat)]
    if (length(gene_symbols_filtered) == 0) {
      gene_symbols_filtered <- gene_symbols[match(rownames(expr_mat),
                                                   names(gene_symbols))]
      if (all(is.na(gene_symbols_filtered))) {
        gene_symbols_filtered <- gene_symbols[seq_len(nrow(expr_mat))]
      }
    }
    gene_expr <- probe_to_gene(expr_mat, gene_symbols_filtered)
    cat("  映射后基因数:", nrow(gene_expr), "\n")

    gene_expr_list[[gse_name]] <- gene_expr
    batch_list[[gse_name]] <- rep(gse_name, n_samples)

    pd <- pData(eset)
    sample_meta_list[[gse_name]] <- data.frame(
      sample_id = colnames(expr_mat),
      dataset   = gse_name,
      platform  = platform,
      title     = pd$title,
      stringsAsFactors = FALSE
    )
    cat("  完成!\n\n")
  }

  # 取共同基因
  common_genes <- Reduce(intersect, lapply(gene_expr_list, rownames))
  cat("跨数据集共同基因数:", length(common_genes), "\n")

  # 合并表达矩阵
  merged_expr <- do.call(cbind, lapply(gene_expr_list, function(x) {
    x[common_genes, , drop = FALSE]
  }))
  batch_vector <- unlist(batch_list)

  # 过滤低方差基因 (方差接近零的基因无法有效校正)
  gene_var <- apply(merged_expr, 1, var, na.rm = TRUE)
  low_var <- gene_var < quantile(gene_var, 0.05, na.rm = TRUE)
  if (any(low_var)) {
    cat("过滤低方差基因:", sum(low_var), "个\n")
    merged_expr <- merged_expr[!low_var, , drop = FALSE]
  }

  cat("合并后矩阵维度:", nrow(merged_expr), "基因 x",
      ncol(merged_expr), "样本\n")
  cat("各批次样本数:\n")
  print(table(batch_vector))

  # ComBat
  cat("\n正在运行 ComBat 批次效应校正...\n")
  combat_expr <- ComBat(
    dat       = as.matrix(merged_expr),
    batch     = batch_vector,
    mod       = NULL,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  cat("ComBat 校正完成!\n")

  sample_meta <- do.call(rbind, sample_meta_list)
  rownames(sample_meta) <- NULL

  return(list(
    combat_expr  = combat_expr,
    sample_meta  = sample_meta,
    common_genes = rownames(combat_expr)
  ))
}

# --- 4. 分别处理EMS和RPL ---

ems_result <- process_disease_datasets(
  ems_dir, "EMS (子宫内膜异位症)", exclude = ems_exclude
)

rpl_result <- process_disease_datasets(
  rpl_dir, "RPL (复发性流产)"
)

# --- 5. 保存结果 ---
cat("\n", strrep("=", 60), "\n")
cat("保存结果\n")
cat(strrep("=", 60), "\n")

saveRDS(ems_result$combat_expr,
        file = file.path(output_dir, "EMS_combat_expr.rds"))
write.csv(ems_result$combat_expr,
          file = file.path(output_dir, "EMS_combat_expr.csv"))
write.csv(ems_result$sample_meta,
          file = file.path(output_dir, "EMS_sample_meta.csv"),
          row.names = FALSE)

saveRDS(rpl_result$combat_expr,
        file = file.path(output_dir, "RPL_combat_expr.rds"))
write.csv(rpl_result$combat_expr,
          file = file.path(output_dir, "RPL_combat_expr.csv"))
write.csv(rpl_result$sample_meta,
          file = file.path(output_dir, "RPL_sample_meta.csv"),
          row.names = FALSE)

cat("文件已保存至:", output_dir, "\n")

# --- 6. 汇总 ---
cat("\n", strrep("=", 60), "\n")
cat("处理完成汇总\n")
cat(strrep("=", 60), "\n")

cat("\n[EMS 子宫内膜异位症]\n")
cat("  排除数据集:", paste(ems_exclude, collapse=", "),
    "(log2 ratio数据, 与单通道数据不兼容)\n")
cat("  纳入数据集:", paste(unique(ems_result$sample_meta$dataset),
                           collapse = ", "), "\n")
cat("  校正后维度:", nrow(ems_result$combat_expr), "基因 x",
    ncol(ems_result$combat_expr), "样本\n")
cat("  表达值范围:", paste(round(range(ems_result$combat_expr), 2),
                           collapse = " ~ "), "\n")

cat("\n[RPL 复发性流产]\n")
cat("  纳入数据集:", paste(unique(rpl_result$sample_meta$dataset),
                           collapse = ", "), "\n")
cat("  校正后维度:", nrow(rpl_result$combat_expr), "基因 x",
    ncol(rpl_result$combat_expr), "样本\n")
cat("  表达值范围:", paste(round(range(rpl_result$combat_expr), 2),
                           collapse = " ~ "), "\n")

cat("\n输出文件:\n")
cat("  EMS_combat_expr.rds / .csv  -> EMS ComBat校正后表达矩阵\n")
cat("  RPL_combat_expr.rds / .csv  -> RPL ComBat校正后表达矩阵\n")
cat("  EMS_sample_meta.csv         -> EMS样本元数据\n")
cat("  RPL_sample_meta.csv         -> RPL样本元数据\n")
