setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
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

#####################trade-off
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


#######################significance
df <- read.csv("Multifunctionality_results.csv")
df1 <- df[1:6,]
df2 <- df[7:12,]

t.test(df1$Multifunctionality_Multi.threshold,df2$Multifunctionality_Multi.threshold)
t.test(df1$Evenness,df2$Evenness)

lm <- lm(Multifunctionality_Multi.threshold~Evenness,data = df)
summary(lm)

t.test(df1$productivity,df2$productivity) #+
t.test(df1$Stability.1,df2$Stability.1) # 
t.test(df1$Pathogen.control,df2$Pathogen.control) #+
t.test(df1$carbon.stocks,df2$carbon.stocks) #+
t.test(df1$maintenance.of.soil.nutrients,df2$maintenance.of.soil.nutrients) #+
t.test(df1$biological.activity,df2$biological.activity) #+
t.test(df1$SM,df2$SM) #+


################################多阈值法
df <- read.csv("Multifunctionality_scaled.csv")

# 1️⃣ 提取功能数据部分
function_data <- df[, -1]

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

#######################################pearson correlation
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 加载必要的包
library(ggplot2)
library(reshape2)
library(pheatmap)

df <- read.csv("correlation.csv",row.names = 1)

# 第1-7列是指标组A，第8-16列是指标组B
groupA <- df[, 8:16]
groupB <- df[, 1:7]

# 初始化相关性和p值矩阵
cor_matrix <- matrix(NA, nrow = ncol(groupB), ncol = ncol(groupB))
p_matrix <- matrix(NA, nrow = ncol(groupB), ncol = ncol(groupB))
rownames(cor_matrix) <- colnames(groupB)
colnames(cor_matrix) <- colnames(groupB)
rownames(p_matrix) <- colnames(groupB)
colnames(p_matrix) <- colnames(groupB)

# 计算 Pearson 相关性和 p 值
for (i in 1:ncol(groupB)) {
  for (j in 1:ncol(groupB)) {
    test <- cor.test(groupB[[i]], groupB[[j]], method = "pearson")
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# FDR 校正
p_matrix_adj <- matrix(p.adjust(p_matrix, method = "fdr"),
                       nrow = nrow(p_matrix),
                       dimnames = dimnames(p_matrix))

# 添加星号表示显著性
stars_matrix <- ifelse(p_matrix_adj < 0.001, "***",
                       ifelse(p_matrix_adj < 0.01, "**",
                              ifelse(p_matrix_adj < 0.05, "*", "")))

# 使用 pheatmap 画图
pheatmap(cor_matrix,
         display_numbers = stars_matrix,  # 显著性标记
         number_color = "black",
         cluster_rows = FALSE,         # ❗ 保持行顺序
         cluster_cols = FALSE,         # ❗ 保持列顺序
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         fontsize = 12,
         main = "Pearson Correlation with Significance Stars")

#############################random forest
library(randomForest)
library(rfPermute)
data <- read.csv("correlation.csv",row.names = 1)
dat <- data[,c(1:7,16)]

rf <- randomForest(Multifunctionality_mean~., data=dat, importance=TRUE, proximity=TRUE)
print(rf)
importance(rf)

rf <- rfPermute(Multifunctionality_mean~., data=dat, importance=TRUE, proximity=TRUE)
print(rf)
importance(rf)

lm <- lm(dat$Multifunctionality_mean~dat$Positive.links)
summary(lm)

####################################100次random forest
set.seed(123)
# 初始化空列表以保存结果
mse_list <- list()
pval_list <- list()

# 循环运行100次模型
for (i in 1:100) {
  rf <- rfPermute(Multifunctionality_mean ~ ., data = dat, importance = TRUE, proximity = TRUE)
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


# 初始化一个向量用于存储每次模型的 % Var explained
var_explained <- numeric(100)
# 运行 100 次 randomForest 模型
for (i in 1:100) {
  rf <- randomForest(Multifunctionality_mean ~ ., data = dat,
                     importance = TRUE, proximity = TRUE)
  var_explained[i] <- tail(rf$rsq, 1) * 100  # 转换为百分比
}

# 输出平均值
mean_var_explained <- mean(var_explained)
print(mean_var_explained)


####################################################################evenness
data <- read.csv("correlation.csv",row.names = 1)
dat <- data[,c(1:5,15)]

rf <- randomForest(Evenness~., data=dat, importance=TRUE, proximity=TRUE)
print(rf)
importance(rf)

rf <- rfPermute(Evenness~., data=dat, importance=TRUE, proximity=TRUE)
print(rf)
importance(rf)

lm <- lm(dat$Evenness~dat$Links)
summary(lm)

lm <- lm(dat$Evenness~.)
summary(lm)
aov <- anova(lm)
aovss <- aov$`Sum Sq`
result <- cbind(aov,exp=aovss/sum(aovss)*100)
result

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

#############################positive or negative
set.seed(123)
# 初始化空列表以保存结果
mse_list <- list()
pval_list <- list()

# 循环运行100次模型
for (i in 1:100) {
  rf <- rfPermute(Evenness ~ Positive.links + Negative.links, data = data, importance = TRUE, proximity = TRUE)
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



##############################################dispersion
###position
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\01_raw_data")
df1 <- read.csv("otutab_rarefied.csv",row.names = 1, check.names = T, header = T)
dis1 <- vegdist(t(df1),method = "bray")
mean(dis1)
groups <- factor(rep(1:12,each = 4),
                 labels = c("11","12","17","18","22","23","28","30","34","36","4","5"))
mod <- betadisper(dis1, groups)
anova(mod)
per <- permutest(mod, pairwise = TRUE, permutations = 99)
boxplot(mod)
mean(mod$distances)


######################################################################tNST
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_network\\Bare")
dat <- read.table("bareT2 0.98 fast_greedy modularity.txt", header = T, sep = "\t",row.names = 1)
dat1 <- dat[,1,drop = F]

# 加载必要的包
library(tidyverse)
# 步骤1：读取数据
# 假设文件位于工作目录中，请替换为实际路径
otu_table <- read.csv("bareT2.csv", row.names = 1, check.names = FALSE)
selected_otus <- dat1

# 步骤2：数据预处理
# 确保OTU编号格式一致（转换为字符型）
selected_otus$ID <- as.character(selected_otus$ID)
rownames(otu_table) <- as.character(rownames(otu_table))

# 步骤3：筛选OTU表
filtered_otu <- otu_table %>%
  rownames_to_column(var = "ID") %>%       # 将行名转换为列
  filter(ID %in% selected_otus$ID) %>% # 筛选目标OTU
  column_to_rownames(var = "ID")           # 将OTU编号恢复为行名
dim(filtered_otu)

# 步骤4：保存结果
write.csv(filtered_otu, "D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST\\otutab_bare_T2.csv")

##################################NST
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST")
# 加载必要的包
library(iCAMP)
library(ape)
library(ieggr)
library(NST)
# 设置参数
otu_dir <- "D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST"  # OTU表的文件夹路径
output_dir <- "D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST"  # 输出文件路径
treat.file <- "D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST\\treat.csv"

rand.time <- 1000
nworker <- 10
memory.G <- 10
prefix.m <- "Bacteria"

# 获取所有时间点的OTU表文件
otu_files <- list.files(otu_dir, pattern = "^otutab_bare_T\\d+\\.csv$", full.names = TRUE)

# 遍历每个时间点的OTU表文件
for (otu_file in otu_files) {
  # 获取时间点信息，例如 "T0"
  timepoint <- sub(".*_(T\\d+)\\.csv$", "\\1", basename(otu_file))
  
  # 读取OTU表数据
  comm <- t(read.csv(otu_file, header = TRUE, row.names = 1, check.names = FALSE))
  treat <- read.csv(treat.file, header = TRUE, row.names = 1)
  treat$ID <- rownames(treat)
  
  # 匹配分类和处理信息
  sampid.check <- match.name(rn.list = list(comm = comm, treat = treat))
  treat <- sampid.check$treat
  comm <- sampid.check$comm
  
  # 删除丰度全为零的列
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  
  # 使用处理信息
  treat.use <- treat[,4,drop = F]
  
  # 计算tNST
  tnstout <- tNST(comm = comm, group = treat.use, dist.method = "bray", 
                  abundance.weighted = FALSE, rand = rand.time,  
                  nworker = nworker, null.model = "PF", output.rand = TRUE,
                  dirichlet = TRUE, SES = FALSE, RC = FALSE)
  
  # 保存tNST输出结果，以时间点命名文件
  output_file <- file.path(output_dir, paste0(prefix.m, ".PF.tNST.pairwise.", timepoint, ".csv"))
  write.csv(tnstout$index.pair.grp, file = output_file)
}


############icamp_循环
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


##############################识别motif中的物种
library(ieggr)
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_network\\Res")
occor.r <- read.csv("resT3 Pearson Correlation.csv",row.names = 1, header = T)
occor.r[is.na(occor.r)] <- 0
diag(occor.r) <- 0
occor.r <- as.matrix(occor.r)
sum(occor.r != 0)
nrow(occor.r)

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


#########################icamp
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\iCAMP")
data <- read.csv("icamp.csv")

library(lme4)
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

#####################################################################plot for manuscript
############################################Fig.1
###Fig. 1_effect size for each function
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 加载必要的库
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
ggsave("Fig1_effect size_degradation.png", p, dpi = 300, width = 5.5, height = 3)

###Fig. 1_evenness
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径

# 提取第18列的数据和处理信息
data_18 <- data[, c(1, 19)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$Evenness ~ data_18$Treat)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Evenness, na.rm = TRUE),
            sd_value = sd(Evenness, na.rm = TRUE),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = Treat, y = mean_value, fill = Treat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Bare land" = "#003366", "Grassland" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Evenness of functions") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
ggsave("Fig1_evenness.png", p, dpi = 300, width = 2, height = 2.5)

###boxplot
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
ggsave("Fig1_evenness_box.png", p, dpi = 300, width = 2, height = 2.5)

###Fig. 1_multifunctionality_mean
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径

# 提取第18列的数据和处理信息
data_18 <- data[, c(1, 20)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$Multifunctionality_mean ~ data_18$Treat)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Multifunctionality_mean, na.rm = TRUE),
            sd_value = sd(Multifunctionality_mean, na.rm = TRUE),
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
ggsave("Fig1_multi_averaging.png", p, dpi = 300, width = 2, height = 2)

###Fig. 1_linear between evenness and multifunctionality_mean
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 读取CSV文件
data <- read.csv("Multifunctionality_results.csv",row.names = 1)

# 提取第19列和第20列
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

ggsave("Fig1_correlation1.png", p, dpi = 300, width = 2.2, height = 2.2)

###Fig. 1_multifunctionality_multiple
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("Multifunctionality_results.csv",row.names = 1)  # 请替换为你的文件路径

# 提取第18列的数据和处理信息
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
ggsave("Fig1_multi_multiple.png", p, dpi = 300, width = 2, height = 2)


###Fig. 1_linear between evenness and multifunctionality_multiple
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 读取CSV文件
data <- read.csv("Multifunctionality_results.csv",row.names = 1)

# 提取第19列与第23-39列的数据
x_data <- data[, 19]
y_data <- data[, 23:39]

# 创建一个空的ggplot对象
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
ggsave("Fig1_correlation2.png", p, dpi = 300, width = 2.2, height = 2.2)

############################################Fig.2
###Fig.2_network
######################################network plot
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_network\\Bare")
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

###Fig.2_network property
#LMM
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\network")
data <- read.csv("Network_property_LMM.csv",row.names = 1) 
#data <- read.csv("Network_property_timepoint.csv",row.names = 1) 
#Average.degree..avgK.,Average.clustering.coefficient..avgCC.,Connectedness..Con.,Module,Modularity,Robustness,Vulnerability
#Average.path.distance..GD., Density..D.  Krackhardt.Connectedness..Con.
data[,c(1:14,21:27)] <- scale(data[,c(1:14,21:27)])
dat <- data
library(lme4)
fm <- lmer(Total.nodes ~ SampleType + (1|TimePoint),data=dat)
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

#plot
# 将Cohen's d值及其置信区间转换为数据框
cohen_d_df <- read.csv("Network_LMM_result.csv")
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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig2_effect size_degradation.png", p, dpi = 300, width = 6, height = 3)

###Fig. 2_motifs
# 载入必要的包
library(ggplot2)
library(dplyr)
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\network")
data <- read.csv("Network_property_LMM.csv",row.names = 1) 

# 提取第18列的数据和处理信息
data_18 <- data[, c(17, 26)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$trancom ~ data_18$SampleType)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(SampleType) %>%
  summarise(mean_value = mean(trancom, na.rm = TRUE),
            sd_value = sd(trancom, na.rm = TRUE),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType, y = mean_value, fill = SampleType)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Bare" = "#003366", "Res" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Number of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig2_trancom.png", p, dpi = 300, width = 1.5, height = 2)

###Fig. 2_linear between motifs and soil nutrients
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 读取CSV文件
data <- read.csv("Multifunctionality_results.csv",row.names = 1)

# 提取第19列和第20列
x <- data[, 9]
y <- data[, 22]

# 进行线性回归
lm_model <- lm(y ~ x)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "TN (g/kg)",
       y = "Number of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("Fig2_correlation2.png", p, dpi = 300, width = 2.2, height = 2.2)

###Fig. 2_assambly
# 载入必要的包
library(ggplot2)
library(dplyr)
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\iCAMP")
data <- read.csv("icamp.csv") 

# 提取第18列的数据和处理信息
data_18 <- data[, c(1, 4)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$HoS ~ data_18$Treat)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(HoS, na.rm = TRUE),
            sd_value = sd(HoS, na.rm = TRUE),
            n = n())

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = Treat, y = mean_value, fill = Treat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("bare" = "#003366", "res" = "#8B0000")) +  # 设置自定义配色
  labs(x = NULL, y = "Importance of HoS") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1,fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig2_HoS.png", p, dpi = 300, width = 1.6, height = 2)

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

# 显示图形
print(p)
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig2_HOS_box.png", p, dpi = 300, width = 1.6, height = 2)


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
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black")  # 设置坐标轴标题颜色为黑色
  )

ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig2_assembly.png", p, dpi = 300, width = 2.4, height = 2)

#########################################################Fig.3
###Fig.3_pearson correlation
# 加载所需的包
# 提取数据
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
data <- read.csv("correlation.csv",row.names = 1)
evenness_data <- data[, 15]  # 第15列为Evenness
network_data <- data[, 1:7]  # 前7列为网络性质

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

# 使用qcorrplot绘制相关性热图
p <- qcorrplot(correlate(network_data), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = cor_df, curvature = nice_curvature()) +
  scale_fill_gradientn(colours = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#8B0000","#CCCCCC99")) +
  guides(
    size = guide_legend(title = "Pearson's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "FDR p-value", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  )
p
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig3_correlation.png", p, dpi = 300, width = 5, height = 5)

###Fig.3_random forest
library(randomForest)
library(rfPermute)
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig3_randomforest.png", p, dpi = 300, width = 3, height =3)


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

#############################positive or negative
set.seed(123)
# 初始化空列表以保存结果
mse_list <- list()
pval_list <- list()

# 循环运行100次模型
for (i in 1:100) {
  rf <- rfPermute(Evenness ~ Positive + Negative, data = data, importance = TRUE, proximity = TRUE)
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

cohen_d_df <- result
cohen_d_df$Column <- rownames(cohen_d_df)
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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig3_randomforest1.png", p, dpi = 300, width = 2, height = 3)


##############计算解释量
# 初始化一个向量用于存储每次模型的 % Var explained
var_explained <- numeric(100)
# 运行 100 次 randomForest 模型
for (i in 1:100) {
  rf <- randomForest(Evenness ~ Positive + Negative, data = data,
                     importance = TRUE, proximity = TRUE)
  var_explained[i] <- tail(rf$rsq, 1) * 100  # 转换为百分比
}

# 输出平均值
mean_var_explained <- mean(var_explained)
print(mean_var_explained)


#######################################################Fig.4
###Fig.4_network
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST")
# 读取OTU表
otu_table <- read.csv("otutab_bare_T3.csv", row.names = 1,header = T)  # 假设OTU表的行名是OTU编号，列名是样本编号
# 计算每个样本的总丰度
sample_sums <- colSums(otu_table)
# 计算相对丰度
relative_abundance <- otu_table / sample_sums
# 查看相对丰度结果
head(relative_abundance)
# 如果你想将结果保存为CSV文件
write.csv(relative_abundance, "otutab_bare_T3_rel.csv")

dat <- read.csv("otutab_bare_T3_rel.csv", row.names = 1)
otu <- read.csv("keystone species.csv", header = T)  # OTU名称表
otu_names <- otu[,17,drop = F]
# 筛选出OTU名称
selected_otus <- otu_names$resT3.1  # 这里假设第二个表中的OTU名称存储在第一列
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
write.csv(grouped_data, "resT2_networked_genus.csv")

#plot
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\NST")
# 加载必要的库
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# 读取CSV文件
data <- read.csv("resT0_networked_genus.csv",row.names = 1)

# 计算每个genus的平均相对丰度和标准误
data_long <- data[,1:6] %>%
  gather(key = "Sample", value = "Relative_Abundance", -Genus)  # 将数据从宽格式转换为长格式
anova_result <- aov(Relative_Abundance ~ Genus, data = data_long)
# 输出方差分析的结果
summary(anova_result)

# 计算平均值和标准误
summary_stats <- data_long %>%
  group_by(Genus) %>%
  summarise(
    Mean_Abundance = mean(Relative_Abundance, na.rm = TRUE),
    SE_Abundance = sd(Relative_Abundance, na.rm = TRUE) / sqrt(n())
  )

# 按照Mean_Abundance降序排序并选取前10的genus
top_10_genus <- summary_stats %>%
  arrange(desc(Mean_Abundance)) %>%
  head(10)

# 绘制柱状图
p <- ggplot(top_10_genus, aes(x = reorder(Genus, Mean_Abundance), y = Mean_Abundance)) +
  geom_bar(stat = "identity", fill = "#003366") +  # 使用柱状图
  geom_errorbar(aes(ymin = Mean_Abundance - SE_Abundance, ymax = Mean_Abundance + SE_Abundance), width = 0.2) +  # 添加误差棒
  labs(x = NULL, y = "Mean relative abundance") +  # 添加标签
  theme_minimal() +  # 使用最小化主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签以避免重叠
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),  # 设置坐标轴标题颜色为黑色
    panel.grid = element_blank(),  # 去除背景网格
    panel.border = element_rect(color = "black", size = 1, fill = "transparent")  # 添加黑色边框
  )
p
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\FigS7_2021.png", p, dpi = 300, width = 3, height =3)

###Fig.4_keystone
################################################keystone species
##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
library(ggplot2)
# 读取潜在竞争OTU的数据
dat <- read.csv("keystone species.csv", header = TRUE)
competition_OTUs <- dat[,15,drop = F]
data <- read.csv("resT2_keystone.csv",row.names = 1)

# 假设您的数据框（data）中的OTU编号列为"OTU", 根据"OTU"列检查是否属于潜在竞争OTUs
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

# 作图
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

# 打印图形
print(p1)
ggsave("FigS7_keystone_resT2.png", p1, dpi = 300,width = 3.5,height = 3.5)


####################################################Fig.5
###Fig.5_NMDS
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_composition_diversity\\beta")
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggsci)

sp <- read.csv("otutab_NMDS.csv",row.names = 1)
sp.nmds <- metaMDS(sp[,-1],distance = "bray",k=4)
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

###Fig.5_diversity
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_composition_diversity\\alpha")
# 载入必要的包
library(ggplot2)
library(dplyr)
# 读取CSV数据
data <- read.csv("alpha_diversity.csv",row.names = 1)  # 请替换为你的文件路径

# 提取第18列的数据和处理信息
data_18 <- data[, c(3, 7)]  # 第一列为处理，18列为需要检验的变量
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

# 绘制柱状图，显示每个处理的均值及其标准差
p <- ggplot(summary_stats, aes(x = SampleType, y = mean_value, fill = SampleType)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # 柱状图
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
                width = 0.2, position = position_dodge(0.7)) +  # 添加误差棒
  scale_fill_manual(values = c("Bare" = "#003366", "Restoration" = "#8B0000")) +  # 设置自定义配色
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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\Fig5_Chao1.png", p, dpi = 300, width = 2.5, height = 2)

#################################Fig.S4
###Fig.S4_linear between motifs and other soil nutrients
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 读取CSV文件
data <- read.csv("Multifunctionality_results.csv",row.names = 1)

# 提取第19列和第20列
x <- data[, 11]
y <- data[, 22]

# 进行线性回归
lm_model <- lm(y ~ x)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "NO3- (mg/kg)",
       y = "Number of trancom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = "transparent"),  # 添加黑色边框
    panel.grid = element_blank(),  # 去除网格线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色为黑色
    axis.title = element_text(color = "black"),
    legend.position = "none"  # 去除图例
  )
p
ggsave("FigS4_NO3.png", p, dpi = 300, width = 2.2, height = 2.2)

#################################Fig.S5_deterministic
# 载入必要的包
library(ggplot2)
library(dplyr)
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\iCAMP")
data <- read.csv("icamp.csv") 

# 提取第18列的数据和处理信息
data_18 <- data[, c(1, 9)]  # 第一列为处理，18列为需要检验的变量
# 进行T检验，假设你有两个处理组：处理A和处理B
# 例如，如果处理列中的值为 "A" 和 "B"，我们进行两组比较
t_test_result <- t.test(data_18$Sto ~ data_18$Treat)
# 打印T检验结果
print(t_test_result)
# 计算每个处理组的均值和标准差
summary_stats <- data_18 %>%
  group_by(Treat) %>%
  summarise(mean_value = mean(Sto, na.rm = TRUE),
            sd_value = sd(Sto, na.rm = TRUE),
            n = n())

# 绘制箱线图，添加散点
p <- ggplot(data_18, aes(x = Treat, y = Sto, fill = Treat)) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.9) +  # 设置箱线图的透明度为0.8
  geom_jitter(width = 0.1, height = 0, color = "black", alpha = 0.5, size = 5) +  # 添加散点，透明度和大小可调整
  scale_fill_manual(values = c("bare" = "#003366", "res" = "#8B0000")) +  # 设置处理组的配色
  labs(x = NULL, y = "Importance of Stochastic") +
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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\FigS5_Sto_box.png", p, dpi = 300, width = 2, height = 2)


#堆积柱状图
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取CSV数据
data <- read.csv("icamp.csv") 

# 提取处理列和第3到第7列的数据
data_sub <- data[9:10, c(1, 8:9)]

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
ggsave("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality\\FigS5_assembly.png", p, dpi = 300, width = 2.4, height = 2)


#################################Fig.S6_multiple regression
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
data <- read.csv("correlation.csv",row.names = 1)
lm <- lm(Evenness~Positive + Negative,data=data)
summary(lm)
lms <- step(lm,direction = "both")
aov <- anova(lms)
aovss <- aov$`Sum Sq`
result <- cbind(aov,exp=aovss/sum(aovss)*100)
result

# 进行线性回归
lm_model <- lm(Evenness~Negative,data=data)
summary(lm_model)
# 绘制散点图并添加线性回归线和95%置信区间
p <- ggplot(data, aes(x = Negative, y = Evenness)) +
  geom_point(color = "darkgreen", size = 5, alpha = 0.5) +  # 增大散点大小，设置为深绿色，透明度为0.5
  geom_smooth(method = "lm", color = "black", size = 1, fill = "grey") +  # 设置回归线颜色为黑色，粗细为1.5，置信区间为灰色
  labs(x = "Number of trancom",
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
ggsave("FigS6_multiregression.png", p, dpi = 300, width = 2.5, height = 2.5)


##############################TableS1_LMM
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_composition_diversity\\pathogenic")
data <- read.csv("pathogen_abundance_LMM.csv",row.names = 1) 
data[,1] <- scale(data[,1])
dat <- data
library(lme4)
fm <- lmer(Pathogen ~ SampleType + (1|TimePoint),data=dat)
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

setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\03_soil property")
data <- read.csv("plant biomass.csv",row.names = 1) 
data[,3] <- scale(data[,3])
dat <- data
library(lme4)
fm <- lmer(Plant.biomass.2017 ~ SampleType + (1|TimePoint),data=dat)
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


setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
data <- read.csv("correlation.csv",row.names = 1) 
dat <- data[,9,drop = F]
t.test(dat[1:6,],dat[7:12,])


#####################################功能之间的相关性
setwd("D:\\桌面\\Phd thesis\\Duolun\\Microbial data\\Multifunctionality")
# 加载包
library(Hmisc)      # rcorr()
library(corrplot)   # 绘图
library(reshape2)   # 数据整理
# 读取 CSV 文件（替换为你的路径）
df <- read.csv("Multifunctionality_scaled.csv", header = TRUE, row.names = 1)
df_numeric <- df
# 计算 Pearson 相关性和 P 值
corr_result <- rcorr(as.matrix(df_numeric), type = "pearson")
cor_matrix <- corr_result$r     # 相关系数矩阵
p_matrix   <- corr_result$P     # P值矩阵
# 展开矩阵并 FDR 校正
p_melt <- melt(p_matrix, na.rm = TRUE)
p_melt$FDR <- p.adjust(p_melt$value, method = "fdr")
# 重新构建 FDR 校正后的矩阵
p_adj_matrix <- acast(p_melt, Var1 ~ Var2, value.var = "FDR")
# 对称性处理（因为相关矩阵是对称的）
p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
# 构造符号矩阵
sig_matrix <- ifelse(p_adj_matrix < 0.001, "***",
                     ifelse(p_adj_matrix < 0.01, "**",
                            ifelse(p_adj_matrix < 0.05, "*", 
                                   ifelse(p_adj_matrix < 0.1, "#", ""))))
# 设置颜色
col_scheme <- colorRampPalette(c("blue", "white", "red"))(200)
# 绘图
corrplot(cor_matrix,
         method = "color",              # 用颜色块展示相关性
         col = col_scheme,              # 红-白-蓝配色
         type = "upper",                # 仅显示上三角
         order = "hclust",              # 聚类排序
         p.mat = p_adj_matrix,          # FDR 校正后的显著性矩阵
         sig.level = 0.1,              # 显著性水平
         insig = "blank",               # 不显著位置空白
         addCoef.col = "black",         # 显示相关系数数字
         tl.col = "black",              # 标签颜色
         tl.srt = 45,                   # 标签旋转角度
         addgrid.col = "black")         # 添加黑色网格线作为边框







