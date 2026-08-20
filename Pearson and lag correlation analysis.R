# 安装并加载必要包
if (!require(corrplot)) install.packages("corrplot")
if (!require(ggplot2)) install.packages("ggplot2")
library(corrplot)
library(ggplot2)

# 数据准备
data_hc <- data.frame(
  Time = c(0.5, 1, 2, 4, 6),
  Reducing_Sugar = c(3085.1, 1351.11, 1822.19, 1180.72, 1220.81), # mg/L
  Amino_Acid = c(1220.55, 931.8, 1285.1, 2209.09, 1149.21), # mg/L
  X5.HMF = c(38.7, 22.74, 52.56, 25.64, 22.92), # mg/L
  Furfuryl_Alcohol = c(51.29, 44.12, 58.79, 72.5, 11.8), # mg/L
  Aldehyde = c(23.42, 31.26, 36.37, 28.79, 1.98), # %
  Ketone = c(2.81, 19.29, 4.08, 5.4, 8.12), # %
  Phenol = c(24.76, 9.25, 25.11, 31.52, 39.44), # %
  Amine = c(0, 0, 2.14, 3.12, 13.83), # %
  Thiophene = c(8.81, 5.34, 5.42, 6.89, 7.99), # %
  Furan = c(18.28, 2.89, 5.57, 2.18, 1.5), # %
  NCHC = c(18.37, 17.39, 16.63, 18.08, 24.81), # %
  pH = c(4.95, 4.99, 4.57, 4.64, 4.51),
  Humic_Acid = c(3.41, 5.66, 6.35, 6.76, 14.20) # %
)

data_phc <- data.frame(
  Time = c(0.5, 1, 2, 4, 6),
  Reducing_Sugar = c(408.94, 849.95, 829.91, 809.86, 609.4), # mg/L
  Amino_Acid = c(3571.3, 3890.62, 4576.82, 4559.83, 2063.01), # mg/L
  X5.HMF = c(20.28, 30.16, 32.81, 48.71, 12.42), # mg/L
  Furfuryl_Alcohol = c(1.27, 24.76, 117.48, 226.15, 74.35), # mg/L
  Aldehyde = c(21.23, 9.89, 3.42, 2.34, 3.05), # %
  Ketone = c(7.23, 6.17, 5.59, 7.61, 5.47), # %
  Phenol = c(18.59, 20.57, 27.35, 26.35, 23.76), # %
  Amine = c(24.5, 21.89, 17.59, 14.85, 28.64), # %
  Thiophene = c(6.07, 8.98, 17.79, 15.63, 11.7), # %
  Furan = c(0.32, 0.31, 0.57, 0.63, 0.48), # %
  NCHC = c(10.61, 14.61, 21.99, 27.29, 23.27), # %
  pH = c(8.87, 7.95, 8.13, 8.68, 7.61), 
  Humic_Acid = c(7.97, 10.34, 11.90, 20.19, 14.67) # %
)

# 检查列名
cat("HC Data Columns:\n")
print(colnames(data_hc))
cat("\nPHC Data Columns:\n")
print(colnames(data_phc))

# 选择分析变量
vars_to_analyze <- c("Reducing_Sugar", "Amino_Acid", "X5.HMF", "Furfuryl_Alcohol", 
                     "Aldehyde", "Ketone", "Phenol", "Amine", "Thiophene", 
                     "Furan", "NCHC", "pH", "Humic_Acid")

# 验证列名是否存在
if (!all(vars_to_analyze %in% colnames(data_hc))) {
  stop("Some variables in vars_to_analyze are not found in data_hc: ", 
       paste(vars_to_analyze[!vars_to_analyze %in% colnames(data_hc)], collapse = ", "))
}
if (!all(vars_to_analyze %in% colnames(data_phc))) {
  stop("Some variables in varstras_to_analyze are not found in data_phc: ", 
       paste(vars_to_analyze[!vars_to_analyze %in% colnames(data_phc)], collapse = ", "))
}

# 提取分析数据
hc_data <- data_hc[, vars_to_analyze]
phc_data <- data_phc[, vars_to_analyze]

# 皮尔逊相关性
pearson_hc <- cor(hc_data, method = "pearson")
pearson_phc <- cor(phc_data, method = "pearson")



# 以 PDF 格式输出相关性热图
pdf("D:/硕博期间/论文/论文1/相关性分析/correlation111.pdf", width = 12, height = 48, family = "Times")
# 设置图形参数
par(mfrow = c(4, 1), mar = c(2, 2, 2, 2) + 0.1)  # 调整边距以适应标签

corrplot(pearson_hc, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, 
         tl.cex = 2.0, cl.cex = 2.2, number.cex = 1.6, diag = FALSE, title = "HC-pearson", mar = c(0, 0, 2, 0),
         col.lim = c(-1, 1), col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
         tl.pos = "td", cl.pos = "r", cl.ratio = 0.2)
corrplot(pearson_phc, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 50, 
         tl.cex = 2.0, cl.cex = 2.2, number.cex = 1.6, diag = FALSE, title = "PHC-pearson", mar = c(0, 0, 2, 0),
         col.lim = c(-1, 1), col = colorRampPalette(c("navy", "white", "firebrick3"))(200),
         tl.pos = "td", cl.pos = "r", cl.ratio = 0.2)
dev.off()

# 时间滞后相关性分析
lag_correlation <- function(df, target_var, other_vars, max_lag = 2) {
  results <- list()
  for (var in other_vars) {
    results[[var]] <- numeric(max_lag + 1)
    for (lag in 0:max_lag) {
      lagged_var <- c(rep(NA, lag), df[[var]][1:(nrow(df) - lag)])
      corr <- cor(df[[target_var]], lagged_var, method = "pearson", use = "complete.obs")
      results[[var]][lag + 1] <- corr
    }
  }
  return(results)
}

# 针对PHC组，分析Humic_Acid与关键变量的滞后相关性
lag_vars <- c("Reducing_Sugar", "X5.HMF", "Furfuryl_Alcohol", "Amino_Acid", "Amine", "NCHC","Aldehyde","Ketone","Phenol","Thiophene","Furan")
lag_results_phc <- lag_correlation(phc_data, "Humic_Acid", lag_vars, max_lag = 2)
cat("\nPHC Group Lag Correlation (Humic_Acid):\n")
for (var in names(lag_results_phc)) {
  cat(sprintf("%s:\n", var))
  for (lag in 0:2) {
    cat(sprintf("  Lag %d: %.3f\n", lag, lag_results_phc[[var]][lag + 1]))
  }
}

# 针对HC组，分析Humic_Acid与关键变量的滞后相关性
lag_results_hc <- lag_correlation(hc_data, "Humic_Acid", lag_vars, max_lag = 2)
cat("\nHC Group Lag Correlation (Humic_Acid):\n")
for (var in names(lag_results_hc)) {
  cat(sprintf("%s:\n", var))
  for (lag in 0:2) {
    cat(sprintf("  Lag %d: %.3f\n", lag, lag_results_hc[[var]][lag + 1]))
  }
}
 
 