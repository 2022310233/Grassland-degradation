##############################scale the functions
##############################
##############################
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

##############################trade-off (the evenness calculation)
##############################
##############################
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
# 5. 创建结果数据框
result_evenness <- data.frame(
  样本 = rownames(func_data),
  香农指数 = shannon,
  丰富度 = richness,
  evenness = evenness
)
# 6. 查看结果
print(result_evenness)

##############################Multifunctionality based on Multi-threshold method
##############################
##############################
df <- read.csv("Multifunctionality_scaled.csv")
# 1. 提取功能数据部分
function_data <- df[, -1]
# 2. 设置多个阈值
thresholds <- seq(0.1, 0.9, by = 0.05)
# 3. 初始化结果列表
multifunc_results <- list()
# 4. 遍历每一个阈值，计算每个样本达到该阈值的功能个数
for (thr in thresholds) {
  # 每个样本有多少功能 >= 当前阈值
  passed_count <- apply(function_data, 1, function(x) sum(x >= thr))
  # 存储结果：行名为样本，列为该阈值下的多功能性得分
  multifunc_results[[paste0("thr_", thr)]] <- passed_count
}
# 5. 合并结果为数据框
multifunc_df <- do.call(cbind, multifunc_results)
multifunc_df <- data.frame(Sample = df$X, multifunc_df)
# 6. 查看结果
print(multifunc_df)

##############################linear mixed model to test the effect size of grassland degradation on bacterial network
##############################the network properties were generated from MENA pipeline
##############################
library(lme4)
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
data[,c(1:14,21:30)] <- scale(data[,c(1:14,21:30)])
dat <- data
fm <- lmer(p.n2 ~ SampleType + (1|TimePoint)+ (1|Position),data=dat)
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

##############################rrn copy number
##############################
##############################
otu_table <- read.csv("otutab_rarefied.csv", row.names = 1)
otu_table1 <- read.csv("otutab_rarefied.csv")
copy_number_table <- read.csv("copy number.csv") ###the copy number.csv was generated from ribosomal RNA operons (rrn) DataBase
# 计算每个OTU在每个样本中的相对丰度
otu_relative_abundance <- apply(otu_table, 2, function(x) x / sum(x))
otu_relative_abundance[1:10, 1:20]
write.csv(otu_relative_abundance,"otutab_rel.csv")
otu_phylum_copy <- merge(otu_table1, copy_number_table, by = "OTU.ID")
head(otu_phylum_copy)
nrow(otu_phylum_copy)
community_copy_numbers <- list()
community_copy_numbers1 <- list()
otu_data <- otu_phylum_copy

sample_results <- list()  
# 遍历每个样本
for (sample in colnames(otu_relative_abundance)) {
    
    # 提取该样本的相对丰度和copy number
    Si <- as.numeric(otu_relative_abundance[otu_data$OTU.ID, sample])
    ni <- as.numeric(otu_data$copy.number)
    
    # 计算分子和分母
    numerator <- sum(Si)
    denominator <- sum(Si / ni)
    
    # 计算community-level copy number
    community_level_copy_number <- numerator / denominator
    
    # 存储该样本的结果
    sample_results[[sample]] <- community_level_copy_number
}
copy_number <- unlist(sample_results)
copy_number <- as.data.frame(copy_number)

#####test the effect of grassland degradation on rrn copy number
library(lme4)
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
data[,c(1:14,21:30)] <- scale(data[,c(1:14,21:30)])
dat <- data
library(lme4)
fm <- lmer(copy_number ~ SampleType + (1|TimePoint)+ (1|Position),data=dat)
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


#########################################################################Figures in manuscript
############################################Fig.1
############################################Fig. 1_effect size for each function
library(effsize)  # 用于计算 Cohen's d
library(ggplot2)  # 用于绘制条形图
# 读取数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为实际文件路径
# 提取第11-17列
data_sub <- data[, 12:18]

# 初始化Cohen's d值的向量
cohen_d_values <- c()
lower_ci_values <- c()
upper_ci_values <- c()

# 对每一列进行处理
for (i in 1:ncol(data_sub)) {
  # 提取当前列的数据
  column_data <- data_sub[, i]
  
  # 提取前6行和后6行的数据
  first_6 <- column_data[1:6]
  last_6 <- column_data[(length(column_data)-5):length(column_data)]
  
  # 计算Cohen's d及其95%置信区间
  cohen_d_result <- cohen.d(last_6, first_6)
  cohen_d <- cohen_d_result$estimate
  lower_ci <- cohen_d_result$conf.int[1]
  upper_ci <- cohen_d_result$conf.int[2]
  
  # 将结果添加到列表中
  cohen_d_values <- c(cohen_d_values, cohen_d)
  lower_ci_values <- c(lower_ci_values, lower_ci)
  upper_ci_values <- c(upper_ci_values, upper_ci)
  
  # 进行T检验并打印p值
  t_test_result <- t.test(first_6, last_6)
  print(paste("T-test p-value for column", colnames(data_sub)[i], ": ", t_test_result$p.value))
}

# 将Cohen's d值及其置信区间转换为数据框
cohen_d_df <- data.frame(
  Column = colnames(data_sub),
  Cohen_d = cohen_d_values,
  Lower_CI = lower_ci_values,
  Upper_CI = upper_ci_values
)
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = unique(colnames(data_sub)))
# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")

#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")

# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))

# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = Cohen_d, y = Column, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2) +  # 添加误差棒
  labs(x = "Effect size", y = NULL) +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black")  # 设置坐标轴标题颜色为黑色
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p

############################################Fig. 1_effect of grassland degradation on functional evenness
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
# 提取处理和第19列的数据
data_19 <- data[, c(1, 19)]  # 假设第一列是处理，19列是要检验的变量
# 执行T检验（假设有两个处理组：A和B）
t_test_result <- t.test(data_19$Evenness ~ data_19$Treat)
print(t_test_result)  # 打印T检验结果
# 计算每个处理组的均值和标准差
summary_stats <- data_19 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Evenness, na.rm = TRUE),
            sd_value = sd(Evenness, na.rm = TRUE),
            n = n())
# 绘制箱线图，添加散点
p <- ggplot(data_19, aes(x = Treat, y = Evenness, fill = Treat)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 5) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("Bare land" = "#003366", "Grassland" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Evenness of functions") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
# 显示图形
print(p)

############################################Fig. 1_effect of grassland degradation on multifunctionality (averaging)
library(ggplot2)
library(dplyr)
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
data_18 <- data[, c(1, 20)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
t_test_result <- t.test(data_18$Multifunctionality_mean ~ data_18$Treat)
print(t_test_result)
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Multifunctionality_mean, na.rm = TRUE),
            sd_value = sd(Multifunctionality_mean, na.rm = TRUE),
            n = n())

p <- ggplot(summary_stats, aes(x = Treat, y = mean_value, fill = Treat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Bare land" = "#003366", "Grassland" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Multifunctionality") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig. 1_linear regression between evenness and multifunctionality_averaging
data <- read.csv("Multifunctionality_results.csv",row.names = 1)
x <- data[, 19]
y <- data[, 20]
# 进行线性回归
lm_model <- lm(y ~ x)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "Evenness of functions",
       y = "Multifunctionality") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig. 1_effect of grassland degradation on multifunctionality (Multi-threshold)
library(ggplot2)
library(dplyr)
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径
data_18 <- data[, c(1, 21)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$Multifunctionality_Multi.threshold ~ data_18$Treat)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Multifunctionality_Multi.threshold, na.rm = TRUE),
            sd_value = sd(Multifunctionality_Multi.threshold, na.rm = TRUE),
            n = n())
# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = Treat, y = mean_value, fill = Treat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Bare land" = "#003366", "Grassland" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Multifunctionality") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig. 1_linear regression between evenness and multifunctionality_Multi-threshold
data <- read.csv("Multifunctionality_results.csv",row.names = 1)
x_data <- data[, 19]
y_data <- data[, 23:39]
plot <- ggplot() 
# 遍历第23-39列，为每列进行线性回归并绘制回归直线
for (i in 1:ncol(y_data)) {
  y_col <- y_data[, i]
  
  # 拟合线性回归模型
  model <- lm(y_col ~ x_data)
  
  # 添加回归直线到图中
  plot <- plot + geom_smooth(data = data.frame(x = x_data, y = y_col), 
                             aes(x = x, y = y), 
                             method = "lm", 
                             se = FALSE, 
                             color = "darkgreen")  # 设定线条为黑色
}
# 设置图形标签和主题
p <- plot + labs(x = "Evenness of functions", y = "Multifunctionality") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig.2
############################################Fig. 2_effect of grassland degradation on bacterial alpha diversity
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("alpha_diversity.csv",row.names = 1)  # 请替换为你的文件路径
dat <- data
library(lme4)
fm <- lmer(Chao1 ~ SampleType + (1|TimePoint),data=dat)
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
###################################plot
data_18 <- data[, c(3, 7)]
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$Chao1 ~ data_18$SampleType)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType) %>%
  summarise(mean_value = mean(Chao1, na.rm = TRUE),
            sd_value = sd(Chao1, na.rm = TRUE),
            n = n())
# 绘制箱线图，添加散点
p <- ggplot(data_18, aes(x = SampleType, y = Chao1, fill = SampleType)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 3) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("Bare" = "#003366", "Restoration" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Chao1") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
                                
############################################Fig.3
############################################Fig.3_network
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

############################################Fig.3_effect of grassland degradation on network properties
cohen_d_df <- read.csv("Network_LMM_result.csv") ###linear mixed model results
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = cohen_d_df$Column)
# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")
#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")
# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))
# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = Cohen_d, y = Column, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2) +  # 添加误差棒
  labs(x = "Effect size", y = NULL) +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black")  # 设置坐标轴标题颜色为黑色
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p

############################################Fig.3_effect of grassland degradation on motif numbers
library(ggplot2)
library(dplyr)
library(lme4)
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality/network")
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
dat <- data
fm <- lmer(n.p2 ~ SampleType + (1|TimePoint),data=dat)
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
#################################plot                                
# 提取第18列的数据和处理信息
data_18 <- data[, c(17, 29)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$n.p2 ~ data_18$SampleType)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType) %>%
  summarise(mean_value = mean(n.p2, na.rm = TRUE),
            sd_value = sd(n.p2, na.rm = TRUE),
            n = n())
# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType, y = mean_value, fill = SampleType)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Degradation" = "#003366", "Undegradation" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Proportion of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig.3_linear between motifs and soil nutrients
data <- read.csv("Multifunctionality_results.csv",row.names = 1)
# 提取第19列和第20列
x <- data[, 6]
y <- data[, 41]
# 进行线性回归
lm_model <- lm(y ~ x)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "TC (g/kg)",
       y = "Proportion of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig.3_effect of grassland degradation on assembly of networked bacterial community
library(ggplot2)
library(dplyr)
library(lme4)                                
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
###############################boxplot                                
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
  scale_fill_manual(values = c("bare" = "#003366", "res" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Importance of HoS") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
print(p)

############################################Fig.3_assembly of networked bacterial community(堆积柱状图)
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
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black")  # 设置坐标轴标题颜色为黑色
  )
p

############################################Fig.4
############################################Fig.4_pearson correlation
library(linkET)
data <- read.csv("correlation.csv",row.names = 1)
evenness_data <- data[, 15]  # 第15列为Evenness
network_data <- data[, c(1:5,19:20)]  # 前7列为网络性质
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
    rd = cut(Pearson_r, breaks = c(-Inf, 0.5, 0.7, Inf), labels = c("< 0.5", "0.5 - 0.7", ">= 0.7")),
    pd = cut(FDR_P_value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )
cor_df$spec <- rep("Evenness",7)
cor_df <- cor_df[, c("spec", setdiff(names(cor_df), "spec"))]
cor_df <- cor_df %>%
  mutate(colour = ifelse(pd == ">= 0.05", "gray", "#8B0000"))  # 将pd为">= 0.05"的行设置为灰色，其他设置为暗红色
cor_df$colour[6] <- "darkblue"
cor_df$rd[6] <- "0.5 - 0.7"
# 使用qcorrplot绘制相关性热图
p <- qcorrplot(correlate(network_data), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = colour, size = rd), data = cor_df, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")) +
  scale_size_manual(values = c(0.5, 1, 1.5)) +
  scale_colour_manual(values = c("#8B0000","darkblue","#CCCCCC99")) +
  guides(
    size = guide_legend(title = "Pearson's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "FDR p-value", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )
p

############################################Fig.4_random forest
library(randomForest)
library(rfPermute)
data <- read.csv("correlation.csv",row.names = 1)
dat <- data[,c(1:5,15)]
####################################100次random forest
set.seed(123)
# 初始化空列表以保存结果
mse_list <- list()
pval_list <- list()
# 循环运行100次模型
for (i in 1:100) {
  rf <- rfPermute(Evenness ~ ., data = dat, importance = TRUE, proximity = TRUE)
  imp <- importance(rf)
  mse_list[[i]] <- imp[, "%IncMSE"]
  pval_list[[i]] <- imp[, "%IncMSE.pval"]
}
# 将列表转换为数据框矩阵并计算平均值
mse_matrix <- do.call(cbind, mse_list)
pval_matrix <- do.call(cbind, pval_list)
mse_mean <- rowMeans(mse_matrix, na.rm = TRUE)
pval_mean <- rowMeans(pval_matrix, na.rm = TRUE)
# 组合结果为数据框
result <- data.frame(
  Mean_IncMSE = mse_mean,
  Mean_IncMSE_pval = pval_mean
)
# 输出结果
print(result)
result$Column <- rownames(result)
cohen_d_df <- result
cohen_d_df$Column <- factor(cohen_d_df$Column, levels = cohen_d_df$Column)
# 定义JAMA期刊的配色方案，可以选择多种蓝色
jama_colors <- c("#003366", "#800080", "#8B0000", "#8B4513", "#2F4F4F", "#4B0082", "#006400")
#c("#4B0082", "#8B0000", "#8E3D59", "#F1A22D", "#7C6A3D", "#346751", "#2866A1")
# 如果有更多条形图，可以通过调整颜色数量来进行循环分配
cohen_d_df$Color <- rep(jama_colors, length.out = nrow(cohen_d_df))
# 绘制Cohen's d条形图，带有置信区间的误差棒，并使用不同颜色
p <- ggplot(cohen_d_df, aes(x = Column, y = Mean_IncMSE, fill = Color)) +
  geom_bar(stat = "identity") +  # 使用JAMA期刊的配色
  labs(x = NULL, y = "Relative importance") +
  theme_minimal(base_size = 14) +  # 使用最小化主题并增加基础字体大小
  theme(
    panel.background = element_rect(fill = "white", color = "black",size = 1),  # 设置白色背景和黑色边框
    panel.grid = element_blank(),  # 去除网格线
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色为黑色
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_identity()  # 确保填充颜色按照Color列进行填充
p
##############计算解释量
# 初始化一个向量用于存储每次模型的 % Var explained
var_explained <- numeric(100)
# 运行 100 次 randomForest 模型
for (i in 1:100) {
  rf <- randomForest(Evenness ~ ., data = dat,
                     importance = TRUE, proximity = TRUE)
  var_explained[i] <- tail(rf$rsq, 1) * 100  # 转换为百分比
}
# 输出平均值
mean_var_explained <- mean(var_explained)
print(mean_var_explained)

############################################Fig.4_regression between evenness and positive (or negative) associations
data <- read.csv("correlation.csv",row.names = 1)
# 进行线性回归
lm_model <- lm(Evenness~Proportion.of.negative,data=data)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = Proportion.of.negative, y = Evenness)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "Negative association",
       y = "Evenness of functions") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p

############################################Fig.5
############################################Fig.5_relative abundance of the taxa involved in potential competition
dat <- read.csv("otutab_bare_T3_rel.csv", row.names = 1)
otu <- read.csv("keystone species.csv", header = T)  #***note:col 2-9:keystone species;col 10-26: taxa involved in potential competition
otu_names <- otu[,17,drop = F]
# 筛选出OTU名称
selected_otus <- otu_names$resT3.1
# 根据选定的OTU筛选相关性矩阵
filtered_otu <- dat[rownames(dat) %in% selected_otus, ]
# 读取筛选后的OTU表
filtered_otu_table <- filtered_otu
filtered_otu_table$OTUID <- rownames(filtered_otu_table)
# 读取OTU编号与菌属名称的对应表
otu_taxonomy <- otu[1:87,17:18,drop = F]  # 假设第一列是OTU编号，第二列是菌属名称
rownames(otu_taxonomy) <- otu_taxonomy$resT3.1
otu_taxonomy <- otu_taxonomy[,-1, drop = F]
otu_taxonomy$OTUID <- rownames(otu_taxonomy)
# 合并OTU表与菌属名称表，根据OTU编号进行匹配
merged_data <- merge(filtered_otu_table, otu_taxonomy, by = "OTUID")
# 按菌属名称进行分组，并对每个分组进行求和
grouped_data <- merged_data %>%
  group_by(Genus3) %>%
  dplyr::summarize(across(starts_with("X"), sum, na.rm = TRUE))  # 假设OTU表中的样本列以X开头
grouped_data <- as.data.frame(grouped_data)
grouped_data

############################################Fig.5_keystone species
setwd("/Users/zhenchengye/Desktop/博士期间项目/多伦数据/Multifunctionality")
library(ggplot2)
dat <- read.csv("keystone species.csv", header = TRUE)
competition_OTUs <- dat[,15,drop = F]
data <- read.csv("resT2_keystone.csv",row.names = 1)
data$Competition_Status <- ifelse(data$ID %in% competition_OTUs$resT2.1, "Competing", "Non-Competing")
# 使用ifelse()函数根据给定的条件为每个节点分配分类
data$node.Classification <- ifelse(data$node.zi > 2.5 & data$node.Pi < 0.62, "Module_hubs",
                                   ifelse(data$node.zi < 2.5 & data$node.Pi > 0.62, "Connectors",
                                          ifelse(data$node.zi > 2.5 & data$node.Pi > 0.62, "Network_hubs", 
                                                 ifelse(data$node.zi < 2.5 & data$node.Pi < 0.62, "Peripherals", NA))))
# 创建颜色映射：暗红色用于潜在竞争OTU，其他颜色用于其他OTU
data$Color <- ifelse(data$Competition_Status == "Competing" & 
                       data$node.Classification %in% c("Module_hubs", "Connectors", "Network_hubs"), 
                     "#8B0000",  # 暗红色
                     ifelse(data$node.Classification == "Peripherals", "#808080",  # 灰色
                            "#443B84FF"))  # 默认其他颜色
p1 <- ggplot(data, aes(node.Pi, node.zi)) +
  geom_point(aes(color = Color), alpha = 0.6, size = 4, shape = 17) +
  scale_color_identity() +  # 使用自定义颜色
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'node.Pi', y = 'node.Zi', color = '') +
  geom_vline(xintercept = 0.61, linetype = 2, linewidth = 1) +  # 垂直线
  geom_hline(yintercept = 2.5, linetype = 2, linewidth = 1) +  # 水平线  
  theme_bw() +
  theme(axis.ticks.length = unit(-0.25, "cm"), 
        axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")), 
        axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        panel.border = element_rect(color = "black", size = 1, fill = "transparent"),
        axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
        axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色为黑色
        axis.title.x = element_text(size = 14),  # 增大X轴标签字体大小
        axis.title.y = element_text(size = 14))  # 增大Y轴标签字体大小
print(p1)

############################################Fig.5_identify the taxa involved in potential competition
library(ieggr)
occor.r <- read.csv("resT3 Pearson Correlation.csv",row.names = 1, header = T)
occor.r[is.na(occor.r)] <- 0
diag(occor.r) <- 0
occor.r <- as.matrix(occor.r)
sum(occor.r != 0)
nrow(occor.r)

trancom<- function(mat) {
  # 计算每行和每列的 1 的数量
  row_sums <- rowSums(abs(mat))
  col_sums <- colSums(abs(mat))
  # 按照 1 的数量对行和列进行排序（降序）
  sorted_rows <- order(row_sums, decreasing = TRUE)
  sorted_cols <- order(col_sums, decreasing = TRUE)
  # 重新排序矩阵
  sorted_mat <- mat[sorted_rows, sorted_cols]
  otu_names <- rownames(sorted_mat)  # OTU 名称
  diag(sorted_mat) <- 0
  sorted_mat[!lower.tri(sorted_mat, diag = TRUE)] <- 0
  nsp <- nrow(sorted_mat)
  sorted_mat[is.na(sorted_mat)] <- 0
  matp <- sorted_mat; matp[matp < 0] <- 0; matp[matp > 0] <- 1
  matn <- sorted_mat; matn[matn > 0] <- 0; matn[matn < 0] <- 1
  matt <- sorted_mat; matt[matt != 0] <- 1
  ntrip <- 0
  triplet_list <- list()
  for (i in 1:nsp) {
    nei <- sum(matt[, i])
    idnei <- which(matt[, i] == 1)
    
    if (nei >= 2) {
      for (k in 1:(nei - 1)) {
        for (z in (k + 1):nei) {
          if ((matn[idnei[z], idnei[k]] + matn[idnei[z], i] + matn[idnei[k], i] == 2) &&
              (matp[idnei[z], idnei[k]] + matp[idnei[z], i] + matp[idnei[k], i] == 0)) {
        
            ntrip <- ntrip + 1
            triplet <- sort(c(otu_names[idnei[z]], otu_names[idnei[k]], otu_names[i]))
            triplet_list[[length(triplet_list) + 1]] <- triplet
          }
        }
      }
    }
  }
  return(list(
    ntrip = ntrip,
    triplets = unique(triplet_list)
  ))
}
                    
result <- trancom(occor.r)
# 提取所有 triplet 组成的 OTU
otu_triplets <- result$triplets
# 将 list 中所有三元组打平成一个向量
otu_vector <- unlist(otu_triplets)
# 去重
otu_unique <- unique(otu_vector)
# 转换为一列的数据框
otu_df <- data.frame(OTU = otu_unique)
# 显示结果
print(otu_df)

#######################################icamp for each plot each timepoint
#######################################
#######################################
# 设置根目录路径
base_dir <- "/vd04/home/YeZhencheng/C_input/icamp"
# 获取所有文件夹路径
folders <- paste0(base_dir, "/bare_T", 0:3)
# 设置通用文件路径
tree.file <- paste0(base_dir, "/rooted_tree.nwk")
clas.file <- paste0(base_dir, "/tax.csv")
treat.file <- paste0(base_dir, "/treat.csv")
# 设置参数
rand.time <- 1000
nworker <- 80
memory.G <- 350
prefix <- "Bacteria"
# 加载所需的R包
library(iCAMP)
library(ape)
library(ieggr)
                    
# 遍历每个文件夹
for (folder in folders) {
  setwd(folder)
  # 提取文件夹中的时间点信息，例如 "T0"
  timepoint <- sub(".*_(T\\d+)$", "\\1", folder)
  
  # 设置当前文件的路径
  com.file <- paste0("otutab_bare_", timepoint, ".csv")
  
  # 运行代码
  save.wd <- folder  # 将保存路径设为当前文件夹
  
  # 读取OTU、分类和处理数据
  comm <- t(read.csv(com.file, header = TRUE, row.names = 1, as.is = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE))
  tre <- read.tree(file = tree.file)
  clas <- read.csv(clas.file, header = TRUE, row.names = 1, as.is = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  treat <- read.csv(treat.file, header = TRUE, row.names = 1)
  treat$ID <- rownames(treat)
  
  # 样本ID匹配
  sampid.check <- match.name(rn.list = list(comm = comm, treat = treat))
  treat <- sampid.check$treat
  comm <- sampid.check$comm
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  
  # OTU ID匹配
  spid.check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tre = tre))
  comm <- spid.check$comm
  clas <- spid.check$clas
  tree <- spid.check$tre
  
  tree$root.edge <- 0
  tree <- root(tre, 1)
  
  # 计算成对的系统发育距离矩阵
  if(!file.exists("pd.desc")) 
  {
    pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
    # output files:
    # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
    # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
    # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
    # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
  }else{
    # if you already calculated the phylogenetic distance matrix in a previous run
    pd.big=list()
    pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
    pd.big$pd.wd=save.wd
    pd.big$pd.file="pd.desc"
    pd.big$pd.name.file="pd.taxon.name.csv"
  }
  
  # iCAMP分析
  bin.size.limit = 24 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
  sig.index="SES.RC" # see other options in help document of icamp.big.
  icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                         prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                         phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                         phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                         nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                         qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                         correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                         ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
  # there are quite a few parameters in this function, please check the help document of "icamp.big".
  # output files:
  # Test.iCAMP.detail.rda: the object "icres" saved in R data format. it is a list object. The first element bNRIiRCa is the result of relative importance of each assembly process in each pairwise comparison (each turnover). The second element "detail" including binning information (named taxabin), phylogenetic and taxonomic metrics results in each bin (named like bNRIi, RCa, etc.), relative abundance of each bin (bin.weight), relative importance of each process in each turnover between communities (processes), input settings (setting), and input community data matrix (comm). See help document of the function icamp.big for more details.
  
  write.csv(icres$.united,file = paste0(prefix,".icres.detail.csv"),row.names = FALSE)
  
  # 10 # iCAMP bin level statistics
  
  treat = treat[,4,drop = F]
  
  icbin=icamp.bins(icamp.detail = icres$detail,treat = treat,
                   clas=clas,silent=FALSE, boot = TRUE,
                   rand.time = rand.time,between.group = F)
  save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
  write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
  write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
  write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
  write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
  write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
            file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
  write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)
  # 输出进度提示
  cat("Processed:", timepoint, "\n")
}
# 所有文件夹处理完成的提示
cat("All folders processed and saved.\n")                    
                    
