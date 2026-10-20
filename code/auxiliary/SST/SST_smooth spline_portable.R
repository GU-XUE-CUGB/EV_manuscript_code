# Smooth the calibrated SST input table.
# 加载所需的包
library(ggplot2)

args <- gsub("~+~", " ", commandArgs(FALSE), fixed = TRUE)
file_path <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
dir.create(file.path(file_path, "generated"), showWarnings = FALSE)

# 读取数据
data <- read.csv(file.path(file_path, "inputs", "SST.csv"))

# 重命名列，方便使用
colnames(data) <- c("Age", "SST")

# 按年龄升序排序
data <- data[order(data$Age), ]

# 设置平滑参数
spar_values <- c(0.1, 0.4, 0.8)

# 生成平滑结果列表
smooth_list <- lapply(spar_values, function(s) {
  smooth.spline(data$Age, data$SST, spar = s)
})

# 将每条平滑结果转成数据框并添加标签
smooth_df <- do.call(rbind, lapply(1:3, function(i) {
  data.frame(Age = smooth_list[[i]]$x,
             SST = smooth_list[[i]]$y,
             Spar = paste0("spar = ", spar_values[i]))
}))

# 创建图对象
p <- ggplot() +
  geom_point(data = data, aes(x = Age, y = SST), color = "black", size = 1.5) +
  geom_line(data = smooth_df, aes(x = Age, y = SST, color = Spar, linetype = Spar), size = 1) +
  scale_x_reverse() +  # 年龄从左到右递减
  labs(title = "SST (45°) with Spline Smoothing",
       x = "Age (Ma)", y = "SST (°C)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("blue", "green", "red"))

# 保存图像为 PDF
ggsave(
  filename = file.path(file_path, "generated", "SST_spline.pdf"),
  plot = p,
  width = 8,
  height = 5,
  dpi = 300,
  device = "pdf"
)
