# Smooth extinction percentages and predict on a 0.1-Ma grid.
args <- gsub("~+~", " ", commandArgs(FALSE), fixed = TRUE)
file_path <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
dir.create(file.path(file_path, "generated"), showWarnings = FALSE)
# 1. 读取原始数据（保留原始列名）
data <- read.csv(file.path(file_path, "inputs", "N_extinction%_Lowery 2020.csv"), check.names = FALSE)

# 2. 按 Age(Ma) 升序排序
data <- data[order(data$`Age(Ma)`), ]

# 3. 平滑拟合（smooth.spline）
fit <- smooth.spline(x = data$`Age(Ma)`, y = data$`Extinction(%)`, spar = 0.1)
# spar 控制平滑程度，可调为 0.3 ~ 0.7 看哪种更自然

# 4. 构建每 0.1 Ma 的预测点
Ma_new <- seq(min(data$`Age(Ma)`), max(data$`Age(Ma)`), by = 0.1)

# 5. 从拟合曲线中获取预测值
predicted <- predict(fit, x = Ma_new)

# 6. 整理成结果数据框
result <- data.frame(age_Ma = predicted$x, extinction_percent = predicted$y)

# 7. 可视化：原始数据 + 拟合曲线
plot(data$`Age(Ma)`, data$`Extinction(%)`,
     pch = 16, col = "gray40", cex = 1.2,
     xlab = "Age (Ma)", ylab = "Extinction (%)",
     main = "Smoothed Extinction Curve (Calcareous Nannofossils)")
lines(result$age_Ma, result$extinction_percent, col = "blue", lwd = 2)
legend("topright", legend = c("Original Data", "Smoothed Curve"),
       col = c("gray40", "blue"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

# 8. 导出平滑结果
write.csv(result,
          file.path(file_path, "generated", "N_e% (100kyr predict).csv"),
          row.names = FALSE)
