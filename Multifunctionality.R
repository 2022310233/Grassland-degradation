#################################################scale
#################################################
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
df <- read.csv("Multifunctionality.csv",row.names = 1)

# 0-1标准化函数
min_max_scale <- function(x) {
  if (length(unique(x)) == 1) {  # 处理常数列
    return(rep(0, length(x)))
  }
  (x - min(x)) / (max(x) - min(x))
}

# 对除第一列外的所有列应用标准化
scaled_data <- as.data.frame(lapply(df, min_max_scale))
# 合并样本名称和标准化后的数据
result_df <- cbind(rownames(df), scaled_data)
# 查看结果
print(result_df)

write.csv(result_df, "Multifunctionality_scaled.csv")

#################################################trade-off (evenness)
#################################################
library(vegan)
result_df <- read.csv("Multifunctionality_scaled.csv")

func_data <- result_df[, -1]  # 移除第一列（样本名称）
rownames(func_data) <- result_df[, 1]  # 设置行名为样本名称

# 2. 计算香农多样性指数
shannon <- diversity(func_data, index = "shannon")

# 3. 计算物种丰富度（每个样本中非零功能参数的数量）
richness <- rowSums(func_data > 0)  # 计算每行大于0的个数

# 4. 计算Pielou均匀度指数 (J = H'/ln(S))
evenness <- shannon / log(richness)

# 6. 创建结果数据框
result_evenness <- data.frame(
  样本 = rownames(func_data),
  香农指数 = shannon,
  丰富度 = richness,
  evenness = evenness
)
# 7. 查看结果
print(result_evenness)

#################################################multi-threshold
#################################################
df <- read.csv("Multifunctionality_scaled.csv")

# 1️⃣ 提取功能数据部分
function_data <- df[, 2:30]

# 2️⃣ 设置多个阈值
thresholds <- seq(0.1, 0.9, by = 0.05)

# 3️⃣ 初始化结果列表
multifunc_results <- list()

# 4️⃣ 遍历每一个阈值，计算每个样本达到该阈值的功能个数
for (thr in thresholds) {
  # 每个样本有多少功能 >= 当前阈值
  passed_count <- apply(function_data, 1, function(x) sum(x >= thr))
  
  # 存储结果：行名为样本，列为该阈值下的多功能性得分
  multifunc_results[[paste0("thr_", thr)]] <- passed_count
}
# 5️⃣ 合并结果为数据框
multifunc_df <- do.call(cbind, multifunc_results)
multifunc_df <- data.frame(Sample = df$X, multifunc_df)
# 查看结果
print(multifunc_df)

#################################################LMM_effect size
#################################################
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("Multifunctionality_results.csv")
dat <- data
library(lme4)
fm <- lmer(Multi_thres ~ SampleType + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result

#################################################LMM_correlation
#################################################
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("correlation.csv")
dat <- data
library(lme4)
library(performance)
fm <- lmer(Multi_mean ~ evenness + (1|TimePoint)+ (1|Position),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
####################plot
library(ggplot2)
dat$pred_Multi_mean <- predict(fm, re.form = NA)

ggplot(dat, aes(evenness, Multi_mean)) +
  geom_point(size=3, alpha=0.7) +
  geom_line(aes(y=pred_Multi_mean), color="blue", size=1.2) +
  theme_classic() +
  xlab("Evenness") +
  ylab("Multifunctionality")

#################################################beta diversity
#################################################
##################################jaccard
setwd("/Users/zhenchengye/Desktop/Phd thesis/Duolun/Microbial data/02_scaling pattern/STAR")
# 读取数据
otu_table <- read.csv('otutab_res.csv', check.names = FALSE)
sample_info <- read.csv('treat_res.csv')
rownames(sample_info) <- sample_info$SampleID

# 转置OTU表，使样本编号成为行
otu_table_transposed <- t(otu_table[-1]) # 移除OTU ID列并转置
colnames(otu_table_transposed) <- otu_table$OTUID # 设置OTU ID为列名
rownames(otu_table_transposed) <- colnames(otu_table)[-1] # 使用原始OTU表的列名（样本编号）作为行名

# 将转置后的OTU表转换为数据框，并添加样本编号列
otu_df <- as.data.frame(otu_table_transposed, stringsAsFactors = FALSE)
otu_df$SampleID <- rownames(otu_df)
# 合并OTU数据和样本信息
merged_data <- merge(sample_info,otu_df, by = "SampleID")
merged_data[,1:15]
rownames(merged_data) <- merged_data$SampleID
# 定义时间段
timepoints <- lapply(0:3, function(x) paste0("T", x))

# 初始化一个空的数据框来存储结果
results_df <- data.frame(
  SampleID = character(),
  Timepoint_Group = character(),
  Distance = numeric(),
  stringsAsFactors = FALSE
)

# 计算每个时间段和位置组合的beta diversity
for (time_group in timepoints) {
  
  time_data <- merged_data[merged_data$TimePoint %in% time_group, ]
  
  # 只保留OTU列
  otu_data <- time_data[, 8:ncol(time_data)]
  
  # 计算Jaccard距离
  dis1 <- ecodist::distance(otu_data, method = "jaccard")
  
  # 转换为矩阵
  dis_mat <- as.matrix(dis1)
  
  # 按行求平均
  dis_mean <- rowMeans(dis_mat, na.rm = TRUE)
  
  # 构建数据框
  temp_df <- data.frame(
    SampleID = rownames(dis_mat),
    Timepoint_Group = paste(time_group, collapse = ", "),
    Distance = dis_mean,
    stringsAsFactors = FALSE
  )
  
  # 合并结果
  results_df <- rbind(results_df, temp_df)
}

results_df
write.csv(results_df,"/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new/beta_jaccard_res.csv")

#################################################plot for manuscript
#################################################
#################################################Fig.1
###Fig. 1_effect size for each functional indicator
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("Multifunctionality_results.csv")
dat <- data

library(lme4)
library(car)

# 从第10列开始为所有指标
response_vars <- names(dat)[10:ncol(dat)]

# 建立空列表存储结果
result_list <- list()

for (resp in response_vars) {
  
  # 构建公式
  form <- as.formula(paste(resp, "~ SampleType + (1|TimePoint) + (1|Position)+(1|PlotID)"))
  
  # 拟合模型
  fm <- lmer(form, data = dat)
  
  # Anova结果
  presult <- car::Anova(fm, type = 2)
  
  # 提取系数
  coefs <- coef(summary(fm))[ , "Estimate"]
  names(coefs) <- paste0(names(coefs), ".mean")
  
  # 提取标准误
  SEvalues <- coef(summary(fm))[ , "Std. Error"]
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  
  # 提取P值
  chisqP <- presult[,3]
  names(chisqP) <- paste0(row.names(presult), ".P")
  
  # 合并结果
  result <- c(coefs, SEvalues, chisqP)
  
  # 转为data.frame
  result_list[[resp]] <- data.frame(
    Indicator = resp,
    t(result),
    check.names = FALSE
  )
}
# 合并所有结果
final_result <- do.call(rbind, result_list)
# 查看结果
final_result
# 可选：导出结果
write.csv(final_result, "results_LMM_functional indicators.csv", row.names = FALSE)

#################################################Fig. 1_effect size plot_按照不同的分组设置颜色
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(ggplot2)
# 读入数据
cohen_d_df1 <- read.csv("results_LMM_functional indicators.csv")
cohen_d_df <- cohen_d_df1[1:29,]

# 保持Indicator顺序
cohen_d_df$Indicator <- factor(cohen_d_df$Indicator, 
                               levels = cohen_d_df$Indicator)

# 将Type设为因子（保证颜色固定）
cohen_d_df$Type <- factor(cohen_d_df$Type)

jama_colors <- c(
  "#003366",
  "#800080",
  "#8B0000",
  "#8B4513",
  "#2F4F4F",
  "#4B0082",
  "#006400",
  "#FF8C00",
  "#00CED1",
  "#DC143C",
  "#1E90FF",
  "#228B22"
)

n_type <- length(unique(cohen_d_df$Type))
jama_colors <- jama_colors[1:n_type]

# 绘图
p <- ggplot(cohen_d_df, 
            aes(x = SampleType.mean, 
                y = Indicator, 
                fill = Type)) +
  
  geom_bar(stat = "identity") +
  
  geom_errorbar(
    aes(xmin = SampleType.mean - SampleType.se,
        xmax = SampleType.mean + SampleType.se),
    width = 0.2
  ) +
  
  labs(
    x = "Effect size",
    y = NULL,
    fill = "Type"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    
    # 统一四边框
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    plot.title = element_text(hjust = 0.5),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  
  scale_fill_manual(values = jama_colors)
p
ggsave("Fig1_effect size_functional indicators.tiff", p, dpi = 300, width = 7, height = 7)

#################################################Fig. 1_barplot_evenness, multifunctionality
###Fig. 1_evenness
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
# 提取处理和第19列的数据
data_19 <- data[, c(5, 38)]  # 假设第一列是处理，19列是要检验的变量

# 执行T检验（假设有两个处理组：A和B）
t_test_result <- t.test(data_19$evenness ~ data_19$SampleType1)
print(t_test_result)  # 打印T检验结果

# 计算每个处理组的均值和标准差
summary_stats <- data_19 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(evenness, na.rm = TRUE),
            se_value = sd(evenness, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制箱线图，添加散点
p <- ggplot(data_19, aes(x = SampleType1, y = evenness, fill = SampleType1)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 5) +
  
  scale_fill_manual(values = c(
    "Degradation" = "#003366", 
    "Undegradation" = "#8B0000"
  )) +
  
  labs(x = NULL, y = "Evenness of functions") +
  
  theme_minimal() +
  
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig1_evenness_box.png", p, dpi = 300, width = 2, height = 2.5)

#################################################Fig. 1_multifunctionality_mean
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
# 提取第18列的数据和处理信息
data_18 <- data[, c(5, 39)]  # 第一列为处理，18列为需要检验的变量
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(Multi_mean, na.rm = TRUE),
            sd_value = sd(Multi_mean, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType1, y = mean_value, fill = SampleType1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  
  scale_fill_manual(values = c(
    "Degradation" = "#003366",
    "Undegradation" = "#8B0000"
  )) +
  
  labs(x = NULL, y = "Multifunctionality") +
  
  theme_minimal() +
  
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )
p
ggsave("Fig1_multi_averaging.png", p, dpi = 300, width = 2, height = 2)

#################################################Fig. 1_linear between evenness and multifunctionality_mean
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(ggplot2)
library(lme4)
# 读取CSV文件
data <- read.csv("Multifunctionality_results.csv",row.names = 1)
# 进行线性回归
fm <- lmer(Multi_mean ~ evenness + (1|TimePoint)+ (1|Position),data=data)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(evenness, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = dat,
    aes(evenness, Multi_mean),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "Evenness of functions",
    y = "Multifunctionality"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig1_correlation1.png", p, dpi = 300, width = 2.5, height = 2.5)

#################################################Fig. 1_multifunctionality_multi-threshold
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
# 提取第18列的数据和处理信息
data_18 <- data[, c(5, 40)]  # 第一列为处理，18列为需要检验的变量
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(Multi_thres, na.rm = TRUE),
            sd_value = sd(Multi_thres, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType1, y = mean_value, fill = SampleType1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  
  scale_fill_manual(values = c(
    "Degradation" = "#003366",
    "Undegradation" = "#8B0000"
  )) +
  
  labs(x = NULL, y = "Multifunctionality") +
  
  theme_minimal() +
  
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )
p
ggsave("Fig1_multi_thres.tiff", p, dpi = 300, width = 2, height = 2)

#################################################Fig. 1_linear between evenness and multifunctionality_multi-threshold
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(dplyr)
# 读取数据
data <- read.csv("Multifunctionality_results.csv", row.names = 1)
dat <- data
# 响应变量
response_vars <- names(dat)[40:57]
# 存储预测结果
pred_list <- list()

for (resp in response_vars) {
  
  form <- as.formula(
    paste(resp, "~ evenness + (1|TimePoint) + (1|Position)")
  )
  
  fm <- lmer(form, data = dat)
  
  # 预测值
  pred <- data.frame(
    evenness = dat$evenness,
    pred = predict(fm, re.form = NA),
    Indicator = resp
  )
  
  pred_list[[resp]] <- pred
}
# 合并所有预测结果
pred_all <- bind_rows(pred_list)

# 绿色渐变
green_colors <- colorRampPalette(c(
  "#b7e4c7",
  "#52b788",
  "#2d6a4f",
  "#004b23"
))(length(response_vars)-1)

# 颜色向量
color_vec <- c(
  "Multi_thres" = "#8B0000",   # 深红色
  setNames(green_colors,
           response_vars[response_vars != "Multi_thres"])
)

p <- ggplot(pred_all, aes(evenness, pred, color = Indicator)) +
  
  geom_line(
    linewidth = 1
  ) +
  
  scale_color_manual(values = color_vec) +
  
  labs(
    x = "Evenness of functions",
    y = "Multifunctionality"
  ) +
  
  theme_minimal() +
  
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )
p
ggsave("Fig1_correlation2.png", p, dpi = 300, width = 2.5, height = 2.5)

data <- read.csv("Multifunctionality_results.csv", row.names = 1)
fm <- lmer(Multi_thres ~ evenness + (1|TimePoint)+ (1|Position),data=data)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2

#################################################Fig.3
#################################################Fig.3_NMDS
setwd("/Users/zhenchengye/Desktop/Phd thesis/Duolun/Microbial data/03_composition_diversity/beta")
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggsci)

sp <- read.csv("otutab_NMDS.csv",row.names = 1)
sp.nmds <- metaMDS(sp[,-1],distance = "jaccard",k=2)
sp.nmds$points
nmds_sp_site <- data.frame(sp$SampleType,sp.nmds$points)
colnames(nmds_sp_site) <- c("SampleType","NMDS1","NMDS2")
adonis <- adonis2(sp[,-1] ~ SampleType,data = sp) #Pr
anosim<-anosim(sp[,-1],sp$SampleType,permutations = 999,distance = "bray") #Significance
mrpp<-mrpp(sp[,-1],sp$SampleType,permutations = 999,distance = "bray") #Significance of delta

p <- ggplot(data = nmds_sp_site,aes(NMDS1,NMDS2))+
  geom_point(aes(color=SampleType),size=5,alpha = 0.7)+
  #stat_ellipse(aes(fill=Site),geom = "polygon",level = 0.95,alpha=0.3,show.legend = F)+
  theme_bw()+
  geom_vline(xintercept = 0,linetype=3,size=1)+
  geom_hline(yintercept = 0,linetype=3,size=1)+
  #annotate("text",x=-0.4,y=0.22,hjust=0,vjust=0,label=paste("Stress:",round(sp.nmds$stress,4)))+
  #annotate("text",x=-0.4,y=0.2,hjust=0,vjust=0,label=paste("Anosim_p:",round(anosim$signif,4)))+
  scale_color_manual(values = c("Bare" = "#003366", "Restoration" = "#8B0000"))+  # 自定义 SampleType 颜色
  theme(axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
        axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色为黑色
        panel.border = element_rect(color = "black", size = 1, fill = "transparent"))  # 增大Y轴标签字体大小
p # 一般要求NMDS的stress<0.2
#pairwise.anosim(sp[,-ncol(sp)], sp$Site, sim.method="bray", p.adjust.m= "fdr") #两两比较
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig5_NMDS.png", p, dpi = 300,width = 4,height = 3)

#################################################Fig.3_alpha diversity
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("correlation.csv",row.names = 1)  # 请替换为你的文件路径
dat <- data
library(lme4)
fm <- lmer(Richness ~ SampleType + (1|TimePoint)+ (1|Position),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
# 提取第18列的数据和处理信息
data_18 <- data[, c(5, 41)]  # 第一列为处理，18列为需要检验的变量
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(Richness, na.rm = TRUE),
            sd_value = sd(Richness, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制箱线图，添加散点
p <- ggplot(data_18, aes(x = SampleType1, y =Richness, fill = SampleType1)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 3) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("Degradation" = "#003366", "Undegradation" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Richness") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("Fig2_Richness.tiff", p, dpi = 300, width = 2.5, height = 2)


#################################################Fig.2_beta diversity
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("correlation.csv",row.names = 1)  # 请替换为你的文件路径
dat <- data
library(lme4)
fm <- lmer(Distance ~ SampleType + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
# 提取第18列的数据和处理信息
data_18 <- data[, c(5, 53)]  # 第一列为处理，18列为需要检验的变量
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(Distance, na.rm = TRUE),
            sd_value = sd(Distance, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制箱线图，添加散点
p <- ggplot(data_18, aes(x = SampleType1, y =Distance, fill = SampleType1)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 3) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("Degradation" = "#003366", "Undegradation" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Distance_jaccard") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("Fig2_distance.tiff", p, dpi = 300, width = 2.5, height = 2)


#################################################Fig.4
#################################################Fig.4_network
setwd("/Users/zhenchengye/Desktop/Phd thesis/Duolun/Microbial data/03_network/Bare")
dat <- read.csv("bareT0 Pearson Correlation.csv", header = TRUE, row.names = 1)
dat <- as.matrix(dat)
dat[1:10,1:10]
# 获取矩阵的行数和列数
n <- nrow(dat)
for (i in 2:n) {
  for (j in 1:(i-1)) {
    dat[i, j] <- dat[j, i]  # 用上三角填充下三角
  }
}

library(GO.db)
## 定义一些颜色
col_g <- "#C1C1C1"
cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")
occor.r <- dat
g <-  graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
# 删除自相关
g <- igraph::simplify(g)
# 删除孤立节点
g <- igraph::delete_vertices(g, which(igraph::degree(g)==0) )

#pdf(paste0("Example_1.pdf"), encoding="MacRoman", width=15, height=3)
#par(mfrow=c(1,2),mar=c(0,0,1,0),font.main=4)
g1 <- g
E(g1)$correlation <- E(g1)$weight
E(g1)$weight <- abs(E(g1)$weight)
set.seed(007)
V(g1)$modularity <- membership(cluster_fast_greedy(g1))

V(g1)$label <- V(g1)$name
V(g1)$label <- NA
modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T)

top_num <- 18
modu_name <- names(modu_sort[1:18])
modu_cols <- cols[1:length(modu_name)]
names(modu_cols) <- modu_name
V(g1)$color <- V(g1)$modularity
V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)],modu_name)]
V(g1)$frame.color <- V(g1)$color

E(g1)$color <- col_g
for ( i in modu_name){
  col_edge <- cols[which(modu_name==i)]
  otu_same_modu <-V(g1)$name[which(V(g1)$modularity==i)]
  E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
}

sub_net_layout <- layout_with_fr(g1, niter=999,grid = 'nogrid')
sub_net_layout <- layout_with_kk(g1)
plot(g1,layout=sub_net_layout, edge.color = E(g1)$color,vertex.size=3)
title(main = paste0('Nodes=',length(V(g1)$name),', ','Edges=',nrow(data.frame(as_edgelist(g1)))))

#################################################Fig.4_effect size_network property
#################################################LMM
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
data[,8:ncol(data)] <- scale(data[,8:ncol(data)])
dat <- data
library(lme4)
library(car)
# 从第8列开始作为响应变量
response_vars <- names(dat)[8:ncol(dat)]
# 建立空列表
result_list <- list()

for (resp in response_vars) {
  
  # 构建公式
  form <- as.formula(paste(resp, "~ SampleType + (1|TimePoint) + (1|Position)+(1|PlotID)"))
  
  # 拟合模型
  fm <- lmer(form, data = dat)
  
  # ANOVA
  presult <- car::Anova(fm, type = 2)
  
  # 提取系数
  coefs <- coef(summary(fm))[ , "Estimate"]
  names(coefs) <- paste0(names(coefs), ".mean")
  
  # 提取SE
  SEvalues <- coef(summary(fm))[ , "Std. Error"]
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  
  # 提取P值
  chisqP <- presult[,3]
  names(chisqP) <- paste0(row.names(presult), ".P")
  
  # 合并结果
  result <- c(coefs, SEvalues, chisqP)
  
  # 转为dataframe
  result_list[[resp]] <- data.frame(
    Indicator = resp,
    t(result),
    check.names = FALSE
  )
}

# 合并所有结果
final_result <- do.call(rbind, result_list)

# 查看结果
final_result

# 输出结果
write.csv(final_result, "results_Network_LMM.csv")

#################################################plot
library(ggplot2)
#plot
# 将Cohen's d值及其置信区间转换为数据框
cohen_d_df1 <- read.csv("results_Network_LMM.csv")
cohen_d_df <- cohen_d_df1[1:7,]
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = cohen_d_df$Column)

# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")
#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")
# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))

# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = SampleType.mean, y = Column, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  geom_errorbar(aes(xmin = SampleType.mean-SampleType.se, xmax = SampleType.mean+SampleType.se), width = 0.2) +  # 添加误差棒
  labs(x = "Effect size", y = NULL) +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p
ggsave("Fig3_effect size_network.tiff", p, dpi = 300, width = 6, height = 3)

#################################################Fig. 3_motifs plot
# 载入必要的包
library(ggplot2)
library(dplyr)
library(lme4)
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
dat <- data
fm <- lmer(trancom ~ SampleType + (1|TimePoint)+(1|Position),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result

# 提取第18列的数据和处理信息
data_18 <- data[, c(4, 26)]  # 第一列为处理，18列为需要检验的变量 26

# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType1) %>%
  summarise(mean_value = mean(n.p2, na.rm = TRUE),
            sd_value = sd(n.p2, na.rm = TRUE)/sqrt(24),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType1, y = mean_value, fill = SampleType1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Degradation" = "#003366", "Undegradation" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Proportion of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    # 添加横纵坐标刻度线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    # 刻度线长度
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("Fig3_trancom_Proportion.tiff", p, dpi = 300, width = 1.5, height = 1.9)

#################################################Fig. 3_correlation between motifs and soil nutrients
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(n.p ~ NO3 + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2

# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(TN, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = dat,
    aes(TN, n.p),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "TN",
    y = "n.p"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig3_correlation with TN.tiff", p, dpi = 300, width = 2.5, height = 2.5)


#################################################Fig. 3_assambly
# 载入必要的包
library(ggplot2)
library(dplyr)
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/iCAMP")
data <- read.csv("icamp.csv") 
fm <- lmer(HoS ~ Treat + (1|Time),data=data)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result

# 提取第18列的数据和处理信息
data_18 <- data[, c(1, 4)]  # 第一列为处理，18列为需要检验的变量
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(HoS, na.rm = TRUE),
            sd_value = sd(HoS, na.rm = TRUE),
            n = n())
##boxplot
data <- read.csv("icamp.csv") 
# 提取第18列的数据和处理信息
data_18 <- data[1:8, c(1, 4)] 
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(HoS, na.rm = TRUE),
            sd_value = sd(HoS, na.rm = TRUE),
            n = n())

# 绘制箱线图，添加散点
p <- ggplot(data_18, aes(x = Treat, y = HoS, fill = Treat)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 5) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("degradation" = "#003366", "undegradation" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Importance of HoS") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.2, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )

# 显示图形
print(p)
ggsave("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new/Fig3_HOS_box.png", p, dpi = 300, width = 1.6, height = 2)


#堆积柱状图
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取CSV数据
data <- read.csv("icamp.csv") 

# 提取处理列和第3到第7列的数据
data_sub <- data[9:10, c(1, 3:7)]

# 将数据从宽格式转换为长格式
data_long <- data_sub %>%
  gather(key = "Process", value = "Value", -Treat)  # -Treat表示不转换处理列

# 绘制堆积柱状图
p <- ggplot(data_long, aes(x = Treat, y = Value, fill = Process)) +
  geom_bar(stat = "identity") +  # 堆积柱状图
  labs(x = NULL, y = "Relative importance") +
  scale_fill_manual(values = c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")) +  # 设置自定义配色
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.2, fill = "transparent"),  # 黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
p
ggsave("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new/Fig3_assembly.tiff", p, dpi = 300, width = 2.4, height = 2)

#################################################Fig.4
#################################################Fig.4_pearson correlation_evenness
# 加载所需的包
# 提取数据
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(linkET)
library(dplyr)
library(ggplot2)
data <- read.csv("correlation.csv",row.names = 1)
evenness_data <- data[, 38]  # 第15列为Evenness
network_data <- data[, c(41,53, 42:45,47,51)]  # 前7列为网络性质

# 计算Evenness与网络性质之间的Pearson相关性
cor_results <- sapply(1:ncol(network_data), function(i) cor(evenness_data, network_data[, i], method = "pearson"))

# 计算相关性p值
cor_pvals <- sapply(1:ncol(network_data), function(i) cor.test(evenness_data, network_data[, i])$p.value)

# FDR校正p值
fdr_pvals <- p.adjust(cor_pvals, method = "fdr")

# 将相关性和FDR p值整理成数据框
cor_df <- data.frame(
  Variable = colnames(network_data),
  Pearson_r = cor_results,
  P_value = cor_pvals,
  FDR_P_value = fdr_pvals
)

# 为了方便图形展示，添加相关性分类和p值分类
cor_df <- cor_df %>%
  mutate(
    rd = cut(Pearson_r, breaks = c(-Inf, 0.3, 0.5, Inf), labels = c("< 0.3", "0.3 - 0.5", ">= 0.5")),
    pd = cut(FDR_P_value, breaks = c(-Inf, 0.05, 0.1, Inf), labels = c("< 0.05", "0.05 - 0.1", ">= 0.1"))
  )
cor_df$spec <- rep("Evenness",8)
cor_df <- cor_df[, c("spec", setdiff(names(cor_df), "spec"))]
cor_df <- cor_df %>%
  mutate(colour = ifelse(pd == ">= 0.1", "gray", "#8B0000"))  # 将pd为">= 0.05"的行设置为灰色，其他设置为暗红色
cor_df$colour[2] <- "darkblue"
cor_df$rd[2] <- "< 0.3"

# 使用qcorrplot绘制相关性热图
p <- qcorrplot(correlate(network_data), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = colour, size = rd), data = cor_df, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#8B0000","darkblue","#CCCCCC99")) +
  guides(
    size = guide_legend(title = "Pearson's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "FDR p-value", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )
p
ggsave("Fig4_correlation.tiff", p, dpi = 300, width = 5, height = 5)

#################################################Fig.4_pearson correlation_multi_mean
# 加载所需的包
# 提取数据
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(linkET)
library(dplyr)
library(ggplot2)
data <- read.csv("correlation.csv",row.names = 1)
evenness_data <- data[, 39]  # 第15列为Evenness
network_data <- data[, c(41,53, 42:45,47,51)]  # 前7列为网络性质

# 计算Evenness与网络性质之间的Pearson相关性
cor_results <- sapply(1:ncol(network_data), function(i) cor(evenness_data, network_data[, i], method = "pearson"))

# 计算相关性p值
cor_pvals <- sapply(1:ncol(network_data), function(i) cor.test(evenness_data, network_data[, i])$p.value)

# FDR校正p值
fdr_pvals <- p.adjust(cor_pvals, method = "fdr")

# 将相关性和FDR p值整理成数据框
cor_df <- data.frame(
  Variable = colnames(network_data),
  Pearson_r = cor_results,
  P_value = cor_pvals,
  FDR_P_value = fdr_pvals
)

# 为了方便图形展示，添加相关性分类和p值分类
cor_df <- cor_df %>%
  mutate(
    rd = cut(Pearson_r, breaks = c(-Inf, 0.3, 0.5, Inf), labels = c("< 0.3", "0.3 - 0.5", ">= 0.5")),
    pd = cut(FDR_P_value, breaks = c(-Inf, 0.05, 0.1, Inf), labels = c("< 0.05", "0.05 - 0.1", ">= 0.1"))
  )
cor_df$spec <- rep("Multi",8)
cor_df <- cor_df[, c("spec", setdiff(names(cor_df), "spec"))]
cor_df <- cor_df %>%
  mutate(colour = ifelse(pd == ">= 0.1", "gray", "#8B0000"))  # 将pd为">= 0.05"的行设置为灰色，其他设置为暗红色
cor_df$colour[2] <- "darkblue"
cor_df$rd[2] <- "0.3 - 0.5"

# 使用qcorrplot绘制相关性热图
p <- qcorrplot(correlate(network_data), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = colour, size = rd), data = cor_df, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#8B0000","darkblue","#CCCCCC99")) +
  guides(
    size = guide_legend(title = "Pearson's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "FDR p-value", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )
p
ggsave("Fig4_correlation_multi.tiff", p, dpi = 300, width = 5, height = 5)



#################################################Fig.4_multiple regression
#################################################evenness
data <- read.csv("correlation.csv",row.names = 1)
data[,38:53] <- scale(data[,38:53] )

# 模型
lm <- lm(evenness~Richness+avgK+Total.nodes+Total.links+Density+Distance,data=data)
summary(lm)
lms <- step(lm,direction = "both")
summary(lms)
aov <- anova(lms)
aovss <- aov$`Sum Sq`
result <- cbind(aov,exp=aovss/sum(aovss)*100)
result

#################################################plot
result$Column <- rownames(result)

cohen_d_df <- result[1:2,]
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = cohen_d_df$Column)

# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")
#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")
# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))

# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = Column, y = exp, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  labs(x = NULL, y = "Explained variation") +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p
ggsave("Fig4_step regression_evenness.tiff", p, dpi = 300, width = 2, height =3)

#################################################Fig4_correlation between evenness and n.p
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(evenness ~ n.p + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
#########plot
# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(n.p, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = dat,
    aes(n.p, evenness),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "n.p",
    y = "Evenness"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig4_evenness and np.tiff", p, dpi = 300, width = 2.5, height = 2.5)

#################################################Fig4_correlation between evenness and trancom之间关联
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(evenness ~ trancom + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
#################################################plot
# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(trancom, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = dat,
    aes(trancom, evenness),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "trancom",
    y = "Evenness"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig4_evenness and trancom.tiff", p, dpi = 300, width = 2.5, height = 2.5)


#################################################Fig.4_多元回归
#################################################multi_mean
data <- read.csv("correlation.csv",row.names = 1)
data[,38:53] <- scale(data[,38:53] )

# 模型
lm <- lm(Multi_mean~Richness+avgK+Total.nodes+Total.links+Density+Distance,data=data)
summary(lm)
lms <- step(lm,direction = "both")
summary(lms)
aov <- anova(lms)
aovss <- aov$`Sum Sq`
result <- cbind(aov,exp=aovss/sum(aovss)*100)
result

#################################################plot
result$Column <- rownames(result)

cohen_d_df <- result[1:3,]
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = cohen_d_df$Column)

# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#006400")
#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")
# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))

# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = Column, y = exp, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  labs(x = NULL, y = "Explained variation") +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p
ggsave("Fig4_step regression_multi.tiff", p, dpi = 300, width = 2, height =3)

#################################################Fig4_correlation between multi and n.p
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(Multi_mean ~ n.p + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
#########plot
# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(n.p, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = dat,
    aes(n.p, Multi_mean),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "n.p",
    y = "Multi_mean"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig4_Multi_mean and np.tiff", p, dpi = 300, width = 2.5, height = 2.5)

#################################################Fig4_correlation between multi and trancom
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(Multi_mean ~ trancom + (1|TimePoint)+ (1|Position)+(1|PlotID),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2
#########plot
# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(trancom, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = dat,
    aes(trancom, Multi_mean),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "trancom",
    y = "Multi_mean"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("Fig4_Multi_mean and trancom.tiff", p, dpi = 300, width = 2.5, height = 2.5)


#################################################Fig.4_SEM
#################################################SEM_evenness
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(lavaan)
library(haven)
library(Hmisc)
mydata11 <- read.csv("correlation.csv") #Multifunctionality_results.csv
#选择用于构建结构方程模型的变量
mydata11 <- mydata11[c("SampleType","Richness","avgK","trancom","evenness","SOC")]
#对数据进行标准化
mydata11[,2:6]<-scale(mydata11[,2:6])
mydata <- mydata11
head(mydata)
#模型构建
model3<-'Richness~ avgK + SOC + SampleType
SOC ~ SampleType
avgK ~ SampleType + SOC
trancom ~ avgK
evenness ~ Richness + avgK + trancom
evenness  ~~ SOC
evenness  ~~ SampleType
trancom ~~ SOC
trancom ~~ SampleType'
#拟合模型
fit1 <- sem(model3,data = mydata)
modificationIndices(fit1, standardized=F)
summary(fit1,standardized =TRUE,fit.measures =TRUE,rsquare = T)

#################################################SEM_Multi_mean
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(lavaan)
library(haven)
library(Hmisc)
mydata11 <- read.csv("correlation.csv")
#选择用于构建结构方程模型的变量
mydata11 <- mydata11[c("SampleType","Richness","avgK","trancom","Multi_mean","SOC")]
#对数据进行标准化
mydata11[,2:6]<-scale(mydata11[,2:6])
mydata <- mydata11
head(mydata)
#模型构建
model3<-'Richness~ avgK + SOC + SampleType
SOC ~ SampleType
avgK ~ SampleType + SOC
trancom ~ avgK
Multi_mean ~ Richness + trancom
Multi_mean  ~~ SOC
Multi_mean  ~~ SampleType
trancom ~~ SOC
trancom ~~ SampleType
Multi_mean ~~ avgK'
#拟合模型
fit1 <- sem(model3,data = mydata)
modificationIndices(fit1, standardized=F)
summary(fit1,standardized =TRUE,fit.measures =TRUE,rsquare = T)


#################################################Uncertainty analysis
#################################################multi_mean_effect size
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(vegan)
library(lme4)
library(lmerTest)
library(car)

## 1. 读入数据 --------------------------------------------------------------
dat <- read.csv("correlation.csv", header = TRUE)

## 2. 设置功能列 ------------------------------------------------------------
func_cols <- 10:38

## 3. 创建结果表 ------------------------------------------------------------
res_list <- list()
idx <- 1

for (X in 3:29) {
  
  start_pos <- 1:(length(func_cols) - X + 1)
  
  for (i in start_pos) {
    
    this_cols <- func_cols[i:(i + X - 1)]
    this_fun_names <- names(dat)[this_cols]
    
    ## ---- 计算 multifunctionality ----------------------------------------
    
    func_data <- dat[, this_cols, drop = FALSE]
    
    multi_avg <- rowMeans(func_data)
    multi_avg[!is.finite(multi_avg)] <- NA
    
    dat$multi_avg <- multi_avg
    
    ## ---- 去除 NA ---------------------------------------------------------
    
    cc <- complete.cases(dat[, c("multi_avg",
                                 "SampleType",
                                 "TimePoint",
                                 "Position")])
    
    dat_sub <- dat[cc, ]
    n_used <- nrow(dat_sub)
    
    ## ---- 样本太少跳过 ---------------------------------------------------
    
    if (n_used < 5) {
      
      res_list[[idx]] <- data.frame(
        X            = X,
        func_indices = paste(this_cols, collapse = ","),
        func_names   = paste(this_fun_names, collapse = "|"),
        n_used       = n_used,
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1
      next
    }
    
    ## ---- 线性混合模型 ---------------------------------------------------
    
    fm <- lmer(
      multi_avg ~ SampleType + 
        (1|TimePoint) + 
        (1|Position),
      data = dat_sub
    )
    
    ## ---- 提取 LMM 结果 ---------------------------------------------------
    
    presult <- car::Anova(fm, type = 2)
    
    coefs <- coef(summary(fm))[ , "Estimate"]
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm))[ , "Std. Error"]
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    chisqP <- c(presult[,3])
    names(chisqP) <- paste0(row.names(presult), ".P")
    
    result <- c(coefs, SEvalues, chisqP)
    
    ## ---- 保存结果 -------------------------------------------------------
    
    res_list[[idx]] <- data.frame(
      X            = X,
      func_indices = paste(this_cols, collapse = ","),
      func_names   = paste(this_fun_names, collapse = "|"),
      n_used       = n_used,
      t(result),
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
}

## 6. 合并结果 --------------------------------------------------------------

res_lmm <- do.call(rbind, res_list)
res_lmm
write.csv(
  res_lmm,
  "results_Uncertainty_LMM_multimean.csv",
  row.names = FALSE
)

#################################################evenness_effect size
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(vegan)
library(lme4)
library(lmerTest)
library(car)

## 1. 读入数据 --------------------------------------------------------------
dat <- read.csv("correlation.csv", header = TRUE)

## 2. 设置功能列 ------------------------------------------------------------
func_cols <- 10:38

## 3. 创建结果表 ------------------------------------------------------------
res_list <- list()
idx <- 1

for (X in 3:29) {
  
  start_pos <- 1:(length(func_cols) - X + 1)
  
  for (i in start_pos) {
    
    this_cols <- func_cols[i:(i + X - 1)]
    this_fun_names <- names(dat)[this_cols]
    
    ## ---- 计算 multifunctionality ----------------------------------------
    
    func_data <- dat[, this_cols, drop = FALSE]
    
    # 2. 计算香农多样性指数
    shannon <- diversity(func_data, index = "shannon")
    
    # 3. 计算物种丰富度（每个样本中非零功能参数的数量）
    richness <- rowSums(func_data > 0)  # 计算每行大于0的个数
    
    # 4. 计算Pielou均匀度指数 (J = H'/ln(S))
    multi_avg <- shannon / log(richness)
    multi_avg[!is.finite(multi_avg)] <- NA
    
    dat$multi_avg <- multi_avg
    
    ## ---- 去除 NA ---------------------------------------------------------
    
    cc <- complete.cases(dat[, c("multi_avg",
                                 "SampleType",
                                 "TimePoint",
                                 "Position")])
    
    dat_sub <- dat[cc, ]
    n_used <- nrow(dat_sub)
    
    ## ---- 样本太少跳过 ---------------------------------------------------
    
    if (n_used < 5) {
      
      res_list[[idx]] <- data.frame(
        X            = X,
        func_indices = paste(this_cols, collapse = ","),
        func_names   = paste(this_fun_names, collapse = "|"),
        n_used       = n_used,
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1
      next
    }
    
    ## ---- 线性混合模型 ---------------------------------------------------
    
    fm <- lmer(
      multi_avg ~ SampleType + 
        (1|TimePoint) + 
        (1|Position),
      data = dat_sub
    )
    
    ## ---- 提取 LMM 结果 ---------------------------------------------------
    
    presult <- car::Anova(fm, type = 2)
    
    coefs <- coef(summary(fm))[ , "Estimate"]
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm))[ , "Std. Error"]
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    chisqP <- c(presult[,3])
    names(chisqP) <- paste0(row.names(presult), ".P")
    
    result <- c(coefs, SEvalues, chisqP)
    
    ## ---- 保存结果 -------------------------------------------------------
    
    res_list[[idx]] <- data.frame(
      X            = X,
      func_indices = paste(this_cols, collapse = ","),
      func_names   = paste(this_fun_names, collapse = "|"),
      n_used       = n_used,
      t(result),
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
}

## 6. 合并结果 --------------------------------------------------------------

res_lmm <- do.call(rbind, res_list)
res_lmm
write.csv(
  res_lmm,
  "results_Uncertainty_LMM_evenness.csv",
  row.names = FALSE
)

#################################################Multi_mean_multiple regression
library(vegan)
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
## 1. 读入数据 --------------------------------------------------------------
dat <- read.csv("correlation.csv", header = TRUE)
## 2. 设置功能列 ------------------------------------------------------------
# 生态系统功能所在列
func_cols <- 10:38
## 3. 逐步回归中要保留的自变量（列名要和 dat 完全一致） -----------------
rich_col  <- "Richness"
avgK_col  <- "avgK"
nodes_col <- "Total.nodes"
links_col <- "Total.links"
dens_col  <- "Density"

predictors <- c("Richness", "avgK", "Total.nodes", "Total.links", "Density")

## 4. 创建结果表列表 --------------------------------------------------------
res_list <- list()
idx <- 1

for (X in 3:29) {
  
  # 相邻列滑动窗口
  start_pos <- 1:(length(func_cols) - X + 1)
  
  for (i in start_pos) {
    
    # 取相邻列
    this_cols <- func_cols[i:(i + X - 1)]
    
    this_fun_names <- names(dat)[this_cols]
    
    ## ---- 计算功能均匀度（Pielou Evenness） -----------------------------
    func_data <- dat[, this_cols, drop = FALSE]
    
    multi_avg <- rowMeans(func_data)
    multi_avg[!is.finite(multi_avg)] <- NA
    
    ## ---- 组装回归所需的数据框 -----------------------------------------
    dat$multi_avg <- multi_avg
    
    # 完整观测（响应变量 + 所有候选自变量都非 NA）
    cc <- complete.cases(dat[, c("multi_avg",
                                 rich_col, avgK_col, nodes_col,
                                 links_col, dens_col)])
    dat_sub <- dat[cc, ]
    
    n_used <- nrow(dat_sub)
    
    # 如果有效样本过少，直接给一行 NA 结果避免报错
    if (n_used < 5) {
      res_list[[idx]] <- data.frame(
        X            = X,
        func_indices = paste(this_cols, collapse = ","),
        func_names   = paste(this_fun_names, collapse = "|"),
        n_used       = n_used,
        R2           = NA,
        adjR2        = NA,
        
        Richness.F   = NA,
        Richness.P   = NA,
        Richness.exp = NA,
        
        avgK.F       = NA,
        avgK.P       = NA,
        avgK.exp     = NA,
        
        Total.nodes.F   = NA,
        Total.nodes.P   = NA,
        Total.nodes.exp = NA,
        
        Total.links.F   = NA,
        Total.links.P   = NA,
        Total.links.exp = NA,
        
        Density.F    = NA,
        Density.P    = NA,
        Density.exp  = NA,
        
        Residuals.exp = NA,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
      next
    }
    
    ## ---- 多元线性回归 + 逐步回归 ----------------------------------------
    lm_full <- lm(multi_avg ~ Richness + avgK + Total.nodes + Total.links + Density,
                  data = dat_sub)
    
    # 逐步回归，双向 step
    lms <- step(lm_full, direction = "both", trace = 0)
    
    sm <- summary(lms)
    aov_tab <- anova(lms)
    
    # 模型 R2、调整 R2
    R2    <- sm$r.squared
    adjR2 <- sm$adj.r.squared
    
    # 各项 Sum Sq 及其百分比
    aovss <- aov_tab$`Sum Sq`
    exp   <- aovss / sum(aovss) * 100
    names(exp) <- rownames(aov_tab)  # 如 "Richness","avgK",...,"Residuals"
    
    # 对应的 F 和 P
    F_val <- aov_tab$`F value`
    P_val <- aov_tab$`Pr(>F)`
    names(F_val) <- rownames(aov_tab)
    names(P_val) <- rownames(aov_tab)
    
    ## ---- 按固定变量名提取，缺失的填 NA ----------------------------------
    get_exp <- function(term) {
      if (term %in% names(exp)) exp[term] else NA
    }
    get_F <- function(term) {
      if (term %in% names(F_val)) F_val[term] else NA
    }
    get_P <- function(term) {
      if (term %in% names(P_val)) P_val[term] else NA
    }
    
    exp_Richness   <- get_exp("Richness")
    exp_avgK       <- get_exp("avgK")
    exp_nodes      <- get_exp("Total.nodes")
    exp_links      <- get_exp("Total.links")
    exp_density    <- get_exp("Density")
    exp_residuals  <- get_exp("Residuals")
    
    F_Richness   <- get_F("Richness")
    F_avgK       <- get_F("avgK")
    F_nodes      <- get_F("Total.nodes")
    F_links      <- get_F("Total.links")
    F_density    <- get_F("Density")
    
    P_Richness   <- get_P("Richness")
    P_avgK       <- get_P("avgK")
    P_nodes      <- get_P("Total.nodes")
    P_links      <- get_P("Total.links")
    P_density    <- get_P("Density")
    
    ## ---- 每个组合一行，缺失指标用 NA -----------------------------------
    res_list[[idx]] <- data.frame(
      X            = X,
      func_indices = paste(this_cols, collapse = ","),
      func_names   = paste(this_fun_names, collapse = "|"),
      n_used       = n_used,
      R2           = R2,
      adjR2        = adjR2,
      
      Richness.F   = F_Richness,
      Richness.P   = P_Richness,
      Richness.exp = exp_Richness,
      
      avgK.F       = F_avgK,
      avgK.P       = P_avgK,
      avgK.exp     = exp_avgK,
      
      Total.nodes.F   = F_nodes,
      Total.nodes.P   = P_nodes,
      Total.nodes.exp = exp_nodes,
      
      Total.links.F   = F_links,
      Total.links.P   = P_links,
      Total.links.exp = exp_links,
      
      Density.F    = F_density,
      Density.P    = P_density,
      Density.exp  = exp_density,
      
      Residuals.exp = exp_residuals,
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
}

## 6. 汇总结果并导出 -------------------------------------------------------
res_step <- do.call(rbind, res_list)

head(res_step)

write.csv(res_step,
          "results_Uncertainty_multimean.csv",
          row.names = FALSE)

#################################################evenness_multiple regression
library(vegan)
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
## 1. 读入数据 --------------------------------------------------------------
dat <- read.csv("correlation.csv", header = TRUE)
## 2. 设置功能列 ------------------------------------------------------------
# 生态系统功能所在列
func_cols <- 10:38
## 3. 逐步回归中要保留的自变量（列名要和 dat 完全一致） -----------------
rich_col  <- "Richness"
avgK_col  <- "avgK"
nodes_col <- "Total.nodes"
links_col <- "Total.links"
dens_col  <- "Density"

predictors <- c("Richness", "avgK", "Total.nodes", "Total.links", "Density")

## 4. 创建结果表列表 --------------------------------------------------------
res_list <- list()
idx <- 1

for (X in 3:29) {
  
  # 相邻列滑动窗口
  start_pos <- 1:(length(func_cols) - X + 1)
  
  for (i in start_pos) {
    
    # 取相邻列
    this_cols <- func_cols[i:(i + X - 1)]
    
    this_fun_names <- names(dat)[this_cols]
    
    ## ---- 计算功能均匀度（Pielou Evenness） -----------------------------
    func_data <- dat[, this_cols, drop = FALSE]
    
    # 2. 计算香农多样性指数
    shannon <- diversity(func_data, index = "shannon")
    
    # 3. 计算物种丰富度（每个样本中非零功能参数的数量）
    richness <- rowSums(func_data > 0)  # 计算每行大于0的个数
    
    # 4. 计算Pielou均匀度指数 (J = H'/ln(S))
    multi_avg <- shannon / log(richness)
    
    multi_avg[!is.finite(multi_avg)] <- NA
    
    ## ---- 组装回归所需的数据框 -----------------------------------------
    dat$multi_avg <- multi_avg
    
    # 完整观测（响应变量 + 所有候选自变量都非 NA）
    cc <- complete.cases(dat[, c("multi_avg",
                                 rich_col, avgK_col, nodes_col,
                                 links_col, dens_col)])
    dat_sub <- dat[cc, ]
    
    n_used <- nrow(dat_sub)
    
    # 如果有效样本过少，直接给一行 NA 结果避免报错
    if (n_used < 5) {
      res_list[[idx]] <- data.frame(
        X            = X,
        func_indices = paste(this_cols, collapse = ","),
        func_names   = paste(this_fun_names, collapse = "|"),
        n_used       = n_used,
        R2           = NA,
        adjR2        = NA,
        
        Richness.F   = NA,
        Richness.P   = NA,
        Richness.exp = NA,
        
        avgK.F       = NA,
        avgK.P       = NA,
        avgK.exp     = NA,
        
        Total.nodes.F   = NA,
        Total.nodes.P   = NA,
        Total.nodes.exp = NA,
        
        Total.links.F   = NA,
        Total.links.P   = NA,
        Total.links.exp = NA,
        
        Density.F    = NA,
        Density.P    = NA,
        Density.exp  = NA,
        
        Residuals.exp = NA,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
      next
    }
    
    ## ---- 多元线性回归 + 逐步回归 ----------------------------------------
    lm_full <- lm(multi_avg ~ Richness + avgK + Total.nodes + Total.links + Density,
                  data = dat_sub)
    
    # 逐步回归，双向 step
    lms <- step(lm_full, direction = "both", trace = 0)
    
    sm <- summary(lms)
    aov_tab <- anova(lms)
    
    # 模型 R2、调整 R2
    R2    <- sm$r.squared
    adjR2 <- sm$adj.r.squared
    
    # 各项 Sum Sq 及其百分比
    aovss <- aov_tab$`Sum Sq`
    exp   <- aovss / sum(aovss) * 100
    names(exp) <- rownames(aov_tab)  # 如 "Richness","avgK",...,"Residuals"
    
    # 对应的 F 和 P
    F_val <- aov_tab$`F value`
    P_val <- aov_tab$`Pr(>F)`
    names(F_val) <- rownames(aov_tab)
    names(P_val) <- rownames(aov_tab)
    
    ## ---- 按固定变量名提取，缺失的填 NA ----------------------------------
    get_exp <- function(term) {
      if (term %in% names(exp)) exp[term] else NA
    }
    get_F <- function(term) {
      if (term %in% names(F_val)) F_val[term] else NA
    }
    get_P <- function(term) {
      if (term %in% names(P_val)) P_val[term] else NA
    }
    
    exp_Richness   <- get_exp("Richness")
    exp_avgK       <- get_exp("avgK")
    exp_nodes      <- get_exp("Total.nodes")
    exp_links      <- get_exp("Total.links")
    exp_density    <- get_exp("Density")
    exp_residuals  <- get_exp("Residuals")
    
    F_Richness   <- get_F("Richness")
    F_avgK       <- get_F("avgK")
    F_nodes      <- get_F("Total.nodes")
    F_links      <- get_F("Total.links")
    F_density    <- get_F("Density")
    
    P_Richness   <- get_P("Richness")
    P_avgK       <- get_P("avgK")
    P_nodes      <- get_P("Total.nodes")
    P_links      <- get_P("Total.links")
    P_density    <- get_P("Density")
    
    ## ---- 每个组合一行，缺失指标用 NA -----------------------------------
    res_list[[idx]] <- data.frame(
      X            = X,
      func_indices = paste(this_cols, collapse = ","),
      func_names   = paste(this_fun_names, collapse = "|"),
      n_used       = n_used,
      R2           = R2,
      adjR2        = adjR2,
      
      Richness.F   = F_Richness,
      Richness.P   = P_Richness,
      Richness.exp = exp_Richness,
      
      avgK.F       = F_avgK,
      avgK.P       = P_avgK,
      avgK.exp     = exp_avgK,
      
      Total.nodes.F   = F_nodes,
      Total.nodes.P   = P_nodes,
      Total.nodes.exp = exp_nodes,
      
      Total.links.F   = F_links,
      Total.links.P   = P_links,
      Total.links.exp = exp_links,
      
      Density.F    = F_density,
      Density.P    = P_density,
      Density.exp  = exp_density,
      
      Residuals.exp = exp_residuals,
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
}

## 6. 汇总结果并导出 -------------------------------------------------------
res_step <- do.call(rbind, res_list)

head(res_step)

write.csv(res_step,
          "results_Uncertainty_evenness.csv",
          row.names = FALSE)

#################################################Fig.S1
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(Hmisc)
library(ggcorrplot)
# 读取数据
data <- read.csv("correlation.csv", header = TRUE)
# 提取第10-38列
dat <- data[,10:38]
# 计算相关系数
cor_matrix <- cor(dat, method = "pearson", use = "pairwise.complete.obs")
# 计算p值
p_matrix <- rcorr(as.matrix(dat), type = "pearson")$P
# FDR矫正 (Benjamini-Hochberg)
p_matrix_fdr <- matrix(
  p.adjust(p_matrix, method = "fdr"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix)
)

# 保持行列名一致
rownames(p_matrix_fdr) <- rownames(p_matrix)
colnames(p_matrix_fdr) <- colnames(p_matrix)

# 绘图
p <- ggcorrplot(
  cor_matrix,
  method = "square",
  type = "upper",
  p.mat = p_matrix_fdr,   # 使用FDR矫正后的p值
  sig.level = 0.05,
  insig = "blank",
  colors = c("#003366", "white", "#B2182B"),
  lab = FALSE,
  outline.color = "white",
  ggtheme = theme_bw()
)

p
ggsave("FigS1_correlation.tiff", p, dpi = 300,width = 9,height = 9)


############################Fig.S_correlation between motifs and soil nutrients
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
# 读取CSV文件
data <- read.csv("correlation.csv",row.names = 1)
dat <- data
library(lme4)
library(performance)
fm <- lmer(n.p ~ NO3 + (1|TimePoint)+ (1|Position),data=dat)
presult<-car::Anova(fm,type=2)
coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
names(coefs)<-paste0(names(coefs),".mean")
SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
names(SEvalues)<-paste0(names(SEvalues),".se")
# tvalues<-coef(summary(fm))[ , "t value"]#t values
# names(tvalues)<-paste0(names(tvalues),".t")
# chisqP<-c(presult[,1],presult[,3])
# names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
chisqP<-c(presult[,3])
names(chisqP)<-c(paste0(row.names(presult),".P"))
# result<-c(coefs,tvalues,SEvalues,chisqP)
result<-c(coefs,SEvalues,chisqP)
result
r2 <- performance::r2(fm)
r2

# 预测值
dat$pred_Multi_mean <- predict(fm, re.form = NA)

# 计算95%置信区间
pred_se <- predict(fm, re.form = NA, se.fit = TRUE)

dat$conf.low  <- pred_se$fit - 1.96 * pred_se$se.fit
dat$conf.high <- pred_se$fit + 1.96 * pred_se$se.fit

# 绘图
p <- ggplot(dat, aes(NO3, pred_Multi_mean)) +
  
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "grey",
    alpha = 0.3
  ) +
  
  geom_line(
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = dat,
    aes(NO3, n.p),
    color = "darkgreen",
    size = 5,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "NO3",
    y = "n.p"
  ) +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",
      linewidth = 1.2,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    
    legend.position = "none"
  )

p
ggsave("FigS_correlation with NO3.tiff", p, dpi = 300, width = 2.5, height = 2.5)

#################################################FigS1
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(tidyverse)
# 读取CSV文件
data <- read.csv("correlation.csv")

# 指标列（第10-38列）
vars <- colnames(data)[10:38]

# 转成长格式
data_long <- data %>%
  select(SampleType1, all_of(vars)) %>%
  pivot_longer(
    cols = -SampleType1,
    names_to = "Indicator",
    values_to = "Value"
  )

# 计算均值和标准误
data_sum <- data_long %>%
  group_by(SampleType1, Indicator) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 绘制柱状图
p <- ggplot(data_sum, aes(x = Indicator, y = mean, fill = SampleType1)) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  scale_fill_manual(
    values = c("#D55E00", "#0072B2")
  ) +
  
  labs(
    x = "",
    y = "Mean value",
    fill = "Sample type"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    
    legend.position = "top",
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    )
  )
p

#################################################FigS_CNP cycling
###Fig. 1_effect size for CNP cycling
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
data <- read.csv("CNP cycling.csv")
dat <- data

library(lme4)
library(car)

# 从第10列开始为所有指标
response_vars <- names(dat)[10:ncol(dat)]

# 建立空列表存储结果
result_list <- list()

for (resp in response_vars) {
  
  # 构建公式
  form <- as.formula(paste(resp, "~ SampleType + (1|TimePoint) + (1|Position)"))
  
  # 拟合模型
  fm <- lmer(form, data = dat)
  
  # Anova结果
  presult <- car::Anova(fm, type = 2)
  
  # 提取系数
  coefs <- coef(summary(fm))[ , "Estimate"]
  names(coefs) <- paste0(names(coefs), ".mean")
  
  # 提取标准误
  SEvalues <- coef(summary(fm))[ , "Std. Error"]
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  
  # 提取P值
  chisqP <- presult[,3]
  names(chisqP) <- paste0(row.names(presult), ".P")
  
  # 合并结果
  result <- c(coefs, SEvalues, chisqP)
  
  # 转为data.frame
  result_list[[resp]] <- data.frame(
    Indicator = resp,
    t(result),
    check.names = FALSE
  )
}
# 合并所有结果
final_result <- do.call(rbind, result_list)
# 查看结果
final_result
# 可选：导出结果
write.csv(final_result, "results_LMM_CNP cycling.csv", row.names = FALSE)

#################################Fig. 1_effect size plot_按照不同的分组设置颜色
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/all treat_new")
library(ggplot2)
# 读入数据
cohen_d_df1 <- read.csv("results_LMM_CNP cycling.csv")
cohen_d_df <- cohen_d_df1[1:3,]

# 保持Indicator顺序
cohen_d_df$Indicator <- factor(cohen_d_df$Indicator, 
                               levels = cohen_d_df$Indicator)

# 将Type设为因子（保证颜色固定）
cohen_d_df$Type <- factor(cohen_d_df$Type)

jama_colors <- c(
  "#003366",
  "#800080",
  "#8B0000",
  "#8B4513",
  "#2F4F4F",
  "#4B0082",
  "#006400",
  "#FF8C00",
  "#00CED1",
  "#DC143C",
  "#1E90FF",
  "#228B22"
)

n_type <- length(unique(cohen_d_df$Type))
jama_colors <- jama_colors[1:n_type]

# 绘图
p <- ggplot(cohen_d_df, 
            aes(x = SampleType.mean, 
                y = Indicator, 
                fill = Type)) +
  
  geom_bar(stat = "identity") +
  
  geom_errorbar(
    aes(xmin = SampleType.mean - SampleType.se,
        xmax = SampleType.mean + SampleType.se),
    width = 0.2
  ) +
  
  labs(
    x = "Effect size",
    y = NULL,
    fill = "Type"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    
    # 统一四边框
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.1, "cm"),
    
    plot.title = element_text(hjust = 0.5),
    
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  
  scale_fill_manual(values = jama_colors)
p
ggsave("Fig1_effect size_CNP cycling.tiff", p, dpi = 300, width = 7, height = 3.5)




