# =============================================================================
# 03 ATP8B4 诊断模型 & ROC 曲线
# - 逻辑回归 (Logistic Regression) 构建诊断模型
# - pROC 绘制 ROC 曲线并计算 AUC
# - 分别对 EMS 和 RPL 进行建模
# - 箱线图展示 Disease vs Control 表达差异
# =============================================================================

library(pROC)
library(ggplot2)

# --- 1. 读取数据 ---
base_dir   <- "d:/learnspace/MyRSci"
output_dir <- file.path(base_dir, "EMS_dual")

ems_expr <- readRDS(file.path(output_dir, "EMS_combat_expr.rds"))
rpl_expr <- readRDS(file.path(output_dir, "RPL_combat_expr.rds"))
ems_meta <- read.csv(file.path(output_dir, "EMS_sample_meta_grouped.csv"),
                     stringsAsFactors = FALSE)
rpl_meta <- read.csv(file.path(output_dir, "RPL_sample_meta_grouped.csv"),
                     stringsAsFactors = FALSE)

gene <- "ATP8B4"

# --- 2. 构建分析数据框 ---

# EMS
ems_df <- data.frame(
  sample  = ems_meta$sample_id,
  group   = factor(ems_meta$group, levels = c("Control", "Disease")),
  dataset = ems_meta$dataset,
  expr    = as.numeric(ems_expr[gene, ]),
  stringsAsFactors = FALSE
)

# RPL: 排除 UIF，仅保留 Disease 和 Control
rpl_keep <- rpl_meta$group %in% c("Disease", "Control")
rpl_df <- data.frame(
  sample  = rpl_meta$sample_id[rpl_keep],
  group   = factor(rpl_meta$group[rpl_keep], levels = c("Control", "Disease")),
  dataset = rpl_meta$dataset[rpl_keep],
  expr    = as.numeric(rpl_expr[gene, rpl_keep]),
  stringsAsFactors = FALSE
)

cat("EMS 样本:", nrow(ems_df), " (Disease:", sum(ems_df$group == "Disease"),
    ", Control:", sum(ems_df$group == "Control"), ")\n")
cat("RPL 样本:", nrow(rpl_df), " (Disease:", sum(rpl_df$group == "Disease"),
    ", Control:", sum(rpl_df$group == "Control"), ")\n")

# --- 3. 逻辑回归建模 ---
cat("\n===== 逻辑回归建模 =====\n")

# EMS
ems_glm <- glm(group ~ expr, data = ems_df, family = binomial)
cat("\n[EMS] 逻辑回归摘要:\n")
print(summary(ems_glm))

# RPL
rpl_glm <- glm(group ~ expr, data = rpl_df, family = binomial)
cat("\n[RPL] 逻辑回归摘要:\n")
print(summary(rpl_glm))

# --- 4. ROC 曲线 ---
cat("\n===== ROC 曲线 =====\n")

ems_pred <- predict(ems_glm, type = "response")
rpl_pred <- predict(rpl_glm, type = "response")

ems_roc <- roc(ems_df$group, ems_pred, levels = c("Control", "Disease"),
               direction = "<", quiet = TRUE)
rpl_roc <- roc(rpl_df$group, rpl_pred, levels = c("Control", "Disease"),
               direction = "<", quiet = TRUE)

cat("EMS AUC:", round(auc(ems_roc), 4), "\n")
cat("RPL AUC:", round(auc(rpl_roc), 4), "\n")

# 95% 置信区间
ems_ci <- ci.auc(ems_roc, method = "delong")
rpl_ci <- ci.auc(rpl_roc, method = "delong")
cat("EMS AUC 95% CI:", round(ems_ci[1], 4), "-", round(ems_ci[3], 4), "\n")
cat("RPL AUC 95% CI:", round(rpl_ci[1], 4), "-", round(rpl_ci[3], 4), "\n")

# 最佳截断值 (Youden Index)
ems_best <- coords(ems_roc, "best", ret = c("threshold", "sensitivity",
                                              "specificity"), best.method = "youden")
rpl_best <- coords(rpl_roc, "best", ret = c("threshold", "sensitivity",
                                              "specificity"), best.method = "youden")
cat("\nEMS 最佳截断点: threshold =", round(ems_best$threshold, 4),
    ", Sensitivity =", round(ems_best$sensitivity, 4),
    ", Specificity =", round(ems_best$specificity, 4), "\n")
cat("RPL 最佳截断点: threshold =", round(rpl_best$threshold, 4),
    ", Sensitivity =", round(rpl_best$sensitivity, 4),
    ", Specificity =", round(rpl_best$specificity, 4), "\n")

# --- 5. 绘制 ROC 曲线 (合并图) ---
cat("\n绘制 ROC 曲线...\n")

pdf(file.path(output_dir, "ROC_ATP8B4_combined.pdf"), width = 7, height = 6)

plot(ems_roc, col = "#E41A1C", lwd = 2.5,
     main = paste0("ROC Curve: ", gene, " Diagnostic Model"),
     legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     print.auc = FALSE)

plot(rpl_roc, col = "#377EB8", lwd = 2.5, add = TRUE)

# 对角线
abline(a = 0, b = 1, lty = 2, col = "grey50")

# 图例
legend("bottomright",
       legend = c(
         paste0("EMS (AUC = ", sprintf("%.3f", auc(ems_roc)),
                ", 95%CI: ", sprintf("%.3f", ems_ci[1]), "-",
                sprintf("%.3f", ems_ci[3]), ")"),
         paste0("RPL (AUC = ", sprintf("%.3f", auc(rpl_roc)),
                ", 95%CI: ", sprintf("%.3f", rpl_ci[1]), "-",
                sprintf("%.3f", rpl_ci[3]), ")")
       ),
       col = c("#E41A1C", "#377EB8"),
       lwd = 2.5, cex = 0.85, bty = "n")

dev.off()
cat("ROC 合并图已保存: ROC_ATP8B4_combined.pdf\n")

# --- 6. 分病种独立 ROC 图 ---

# EMS 独立 ROC
pdf(file.path(output_dir, "ROC_ATP8B4_EMS.pdf"), width = 6, height = 6)
plot(ems_roc, col = "#E41A1C", lwd = 2.5,
     main = paste0("EMS: ", gene, " Diagnostic ROC"),
     legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)")
abline(a = 0, b = 1, lty = 2, col = "grey50")
text(0.4, 0.3,
     labels = paste0("AUC = ", sprintf("%.3f", auc(ems_roc)),
                     "\n95%CI: ", sprintf("%.3f", ems_ci[1]), "-",
                     sprintf("%.3f", ems_ci[3]),
                     "\nSens = ", sprintf("%.1f%%", ems_best$sensitivity * 100),
                     "\nSpec = ", sprintf("%.1f%%", ems_best$specificity * 100)),
     cex = 0.9, adj = 0)
dev.off()
cat("EMS 独立 ROC 已保存: ROC_ATP8B4_EMS.pdf\n")

# RPL 独立 ROC
pdf(file.path(output_dir, "ROC_ATP8B4_RPL.pdf"), width = 6, height = 6)
plot(rpl_roc, col = "#377EB8", lwd = 2.5,
     main = paste0("RPL: ", gene, " Diagnostic ROC"),
     legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)")
abline(a = 0, b = 1, lty = 2, col = "grey50")
text(0.4, 0.3,
     labels = paste0("AUC = ", sprintf("%.3f", auc(rpl_roc)),
                     "\n95%CI: ", sprintf("%.3f", rpl_ci[1]), "-",
                     sprintf("%.3f", rpl_ci[3]),
                     "\nSens = ", sprintf("%.1f%%", rpl_best$sensitivity * 100),
                     "\nSpec = ", sprintf("%.1f%%", rpl_best$specificity * 100)),
     cex = 0.9, adj = 0)
dev.off()
cat("RPL 独立 ROC 已保存: ROC_ATP8B4_RPL.pdf\n")

# --- 7. 箱线图: Disease vs Control 表达量对比 ---
cat("\n绘制箱线图...\n")

# EMS 箱线图
p_ems <- ggplot(ems_df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_fill_manual(values = c(Control = "#377EB8", Disease = "#E41A1C")) +
  labs(title = paste0("EMS: ", gene, " Expression"),
       x = NULL, y = "Expression (log2, ComBat corrected)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  # Wilcoxon 检验 p 值
  annotate("text", x = 1.5, y = max(ems_df$expr) + 0.3,
           label = paste0("Wilcoxon p = ",
                          format.pval(wilcox.test(expr ~ group, data = ems_df)$p.value,
                                      digits = 3)),
           size = 4)

ggsave(file.path(output_dir, "boxplot_ATP8B4_EMS.pdf"), p_ems,
       width = 4.5, height = 5)
cat("EMS 箱线图已保存: boxplot_ATP8B4_EMS.pdf\n")

# RPL 箱线图
p_rpl <- ggplot(rpl_df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_fill_manual(values = c(Control = "#377EB8", Disease = "#E41A1C")) +
  labs(title = paste0("RPL: ", gene, " Expression"),
       x = NULL, y = "Expression (log2, ComBat corrected)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  annotate("text", x = 1.5, y = max(rpl_df$expr) + 0.3,
           label = paste0("Wilcoxon p = ",
                          format.pval(wilcox.test(expr ~ group, data = rpl_df)$p.value,
                                      digits = 3)),
           size = 4)

ggsave(file.path(output_dir, "boxplot_ATP8B4_RPL.pdf"), p_rpl,
       width = 4.5, height = 5)
cat("RPL 箱线图已保存: boxplot_ATP8B4_RPL.pdf\n")

# --- 8. Kaplan-Meier 曲线说明 ---
cat("\n", strrep("=", 60), "\n")
cat("关于 Kaplan-Meier 生存曲线\n")
cat(strrep("=", 60), "\n")
cat("当前 GEO 数据集中没有生存/预后随访数据 (OS, PFS 等)。\n")
cat("EMS 和 RPL 属于非肿瘤性妇科疾病，GEO 转录组数据通常\n")
cat("仅包含病例/对照标签，不包含生存时间和事件信息。\n")
cat("因此无法构建 Kaplan-Meier 生存曲线。\n")
cat("\n如需生存分析，可考虑:\n")
cat("  1. 查找含随访数据的临床队列数据集\n")
cat("  2. 在 TCGA 等肿瘤数据库中验证 ATP8B4 的预后价值\n")
cat("     (如子宫内膜癌 UCEC 等相关癌种)\n")

# --- 9. 汇总 ---
cat("\n", strrep("=", 60), "\n")
cat("分析完成汇总\n")
cat(strrep("=", 60), "\n")
cat("\n模型基因:", gene, "\n")
cat("\nEMS 诊断模型:\n")
cat("  AUC =", sprintf("%.3f", auc(ems_roc)),
    " (95%CI:", sprintf("%.3f-%.3f", ems_ci[1], ems_ci[3]), ")\n")
cat("  最佳截断 - Sensitivity:", sprintf("%.1f%%", ems_best$sensitivity * 100),
    ", Specificity:", sprintf("%.1f%%", ems_best$specificity * 100), "\n")
cat("\nRPL 诊断模型:\n")
cat("  AUC =", sprintf("%.3f", auc(rpl_roc)),
    " (95%CI:", sprintf("%.3f-%.3f", rpl_ci[1], rpl_ci[3]), ")\n")
cat("  最佳截断 - Sensitivity:", sprintf("%.1f%%", rpl_best$sensitivity * 100),
    ", Specificity:", sprintf("%.1f%%", rpl_best$specificity * 100), "\n")
cat("\nKaplan-Meier: 数据中无生存/预后信息，无法绘制\n")
cat("\n输出文件:\n")
cat("  ROC_ATP8B4_combined.pdf  -> EMS+RPL 合并 ROC 曲线\n")
cat("  ROC_ATP8B4_EMS.pdf       -> EMS 独立 ROC 曲线\n")
cat("  ROC_ATP8B4_RPL.pdf       -> RPL 独立 ROC 曲线\n")
cat("  boxplot_ATP8B4_EMS.pdf   -> EMS 表达量箱线图\n")
cat("  boxplot_ATP8B4_RPL.pdf   -> RPL 表达量箱线图\n")
