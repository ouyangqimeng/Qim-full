# Qimeng R

## Markdown语法 **[Markdown官网](https://markdown.com.cn/basic-syntax/ "最好的Markdown语法教程")**.
Markdown语法讲解了markdown语句的基本用法包括:<br>
**[标题语法](https://markdown.com.cn/basic-syntax/headings.html)** <br>
**[段落语法](https://markdown.com.cn/basic-syntax/paragraphs.html)** <br>
**[换行语法](https://markdown.com.cn/basic-syntax/line-breaks.html)** <br>
**[强调语法](https://markdown.com.cn/basic-syntax/emphasis.html)** <br>
**[引用语法](https://markdown.com.cn/basic-syntax/blockquotes.html)** <br>
**[列表语法](https://markdown.com.cn/basic-syntax/lists.html)** <br>
**[代码语法](https://markdown.com.cn/basic-syntax/code.html)** <br>
**[分割线语法](https://markdown.com.cn/basic-syntax/horizontal-rules.html)** <br>
**[链接语法](https://markdown.com.cn/basic-syntax/links.html)** <br>
**[图片语法](https://markdown.com.cn/basic-syntax/images.html)** <br>
**[转义字符语法](https://markdown.com.cn/basic-syntax/escaping-characters.html)** <br>

## R语言PCA & PLS-DA 
```   
#------ PCA and PLS-DA ------
BiocManager::install('mixOmics')
library(mixOmics)

data(srbct)
X = srbct$gene  #the gene expression data
dim(X)
X[1:6,1:10]
summary(srbct$class)
#------ PCA分析 ------
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)

# 输出每个主成分可解释方差(变异)
plot(pca.srbct) 
#在下面的样本图中，样本用前两个主成分表示，并根据肿瘤类型着色。
#这里我们观察到变异的主要来源可能不能用肿瘤类型来解释。
#注意，由于PCA是无监督的，为了可视化目的，我们只考虑PCA之后的样本类型信息。

plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

#------ PLS-DA分析 ------

#对于判别分析，我们设置因子Y来表示每个样本的类别隶属度。
#在PLS-DA过程中，将Y因子转化为一个虚拟矩阵。
Y = srbct$class 
summary(Y) 

#PLS-DA模型采用10个成分来评估最终模型所需的性能和成分数量(见下文)。
##样品图:将样本投影到前两个成分所形成的子空间中
srbct.plsda <- plsda(X, Y, ncomp = 10) 
plotIndiv(srbct.plsda , comp = 1:2,
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')
#从样本图中可以看出，与无监督的PCA样本图相比，四种肿瘤类型明显分离。
#绘制每个类的置信椭圆以突出区分的强度(置信水平默认设置为95%，参数ellipse.level)。

#在覆盖样本图之前，可以通过计算背景面来可视化预测区域。
#with background
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
#optional: xlim = c(-40,40), ylim = c(-30,30))

plotIndiv(srbct.plsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

#------ 性能 ------

#PLS-DA模型的分类性能通过重复10次的5折交叉验证来评估。
#重复次数对于确保分类错误率的良好估计是必要的(因为cv -fold是以随机方式确定的)。
#从性能结果中我们可以决定选择最终的PLS模型的成分数量。
set.seed(2543) # 设置种子，便于重复
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = FALSE, auc = TRUE, nrepeat = 10)

#perf.plsda.srbct$error.rate  错误率
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
#从性能图中，可以看出总体错误率和平衡错误率(BER)相似，从1个成分急剧下降到3个成分。
#6个成分后错误率趋于稳定。ncomp = 6时，BER和最大距离足以实现良好的性能(0.06错误率)。

auc.plsda = auroc(srbct.plsda, roc.comp = 6)

#------ 以数据iris为例复现PCA ------
#导入iris数据集
data<-iris
head(data)

#对原数据进行z-score归一化；
dt<-as.matrix(scale(data[,1:4])) #不含Species列
head(dt)

#计算协方差
rm1<-cor(dt)
rm1

#特征值与特征向量均为矩阵分解的结果。
#特征值表示标量部分，一般为某个主成分的方差，其相对比例可理解为方差解释度或贡献度 ；
#特征值从第一主成分会逐渐减小。
#特征向量为对应主成分的线性转换向量（线性回归系数），特征向量与原始矩阵的矩阵积为主成分得分。
#特征向量是单位向量，其平方和为1。
#特征分解
rs1<-eigen(rm1)
rs1

val <- rs1$values #转换成标准差Standard deviation
Standard_deviation <- sqrt(val)
Standard_deviation
#计算方差贡献率和累积贡献率；
Proportion_of_Variance <- val/sum(val)
Proportion_of_Variance
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
Cumulative_Proportion

#碎石图绘制
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",     
                        cex=2,     
                        cex.lab=2,     
                        cex.axis=2,     
                        lty=2,     
                        lwd=2,     
                        xlab = "PC",     
                        ylab="Proportion_of_Variance")

#计算主成分得分
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
U<-as.matrix(rs1$vectors)
U
#进行矩阵乘法，获得PC score；
PC <-dt %*% U
colnames(PC) <- c("PC1","PC2","PC3","PC4")
head(PC)

#可视化
#合并Species列
df<-data.frame(PC,iris$Species)
head(df) 

library(ggplot2)
#提取主成分的方差贡献率，生成坐标轴标题
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
#绘制散点图并添加置信椭圆
p1<-ggplot(data = df,aes(x=PC1,y=PC2,color=iris.Species))+
  stat_ellipse(aes(fill=iris.Species),type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point(size = 2)+
  labs(x=xlab,y=ylab,color="")+
  guides(fill=F) +
  theme_classic(base_line_size = 1) +
  theme(axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  ) 
p1
```




 
