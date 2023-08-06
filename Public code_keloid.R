##合并样本，矫正批次效应
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")
#引用包
library(limma)
library(sva)
outFile="merge.txt"       #输出文件

#获取目录下所有".txt"结尾的文件
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
	if(file==outFile){next}
    rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
    geneNames=as.vector(rt[,1])      #提取基因名称
    uniqGene=unique(geneNames)       #基因取unique
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #读取输入文件，并对输入文件进行整理
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)

    #对数值大的数据取log2
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    
    #数据合并
    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)




##数据标准化
library(limma)               #引用包
expFile="merge.txt"     #表达数据文件
conFile="control.txt"             #对照组样品信息文件
treatFile="keloid.txt"           #实验组样品信息文件
#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#读取对照组样品信息文件,提取对照组的表达数据
s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]
#读取实验组样品信息文件,提取实验组的表达数据
s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]
#数据合并
rt=cbind(conData, treatData)
#如果数据没有取log2,会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
#输出矫正后的表达数据, 同时在样品名字后面加上样品的分组信息
conNum=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)



#铜死亡基因的表达数据
library(limma)              #引用包
expFile="normalize.txt"     #表达数据文件
geneFile="CuproptosisGene.txt"         #基因列表文件
#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#读取基因列表文件，获取铜死亡相关基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
#输出铜死亡基因的表达数据
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="CupGeneExp.txt", sep="\t", quote=F, col.names=F)




#铜死亡相关基因组间差异
#引用包
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)
expFile="CupGeneExp.txt"      #表达数据文件
#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
exp=data
#提取样品的分组信息
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
#基因差异分析
sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ Type)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	if(pvalue<0.05){
	sigVec=c(sigVec, paste0(i, Sig))
	sigGeneVec=c(sigGeneVec, i)}
}
#输出差异基因的表达量
data=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec
#对差异基因进行可视化，绘制热图
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=7, height=4.5)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
#把表达数据转换成ggplot2输入文件
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")
#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     add="point",
	     width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")
#输出箱线图
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()



#绘制铜死亡基因染色体定位圈图
#install.packages("RCircos")
library("RCircos")       #引用包
#初始化圈图
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
#设置圈图的参数
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)
#输出文件
pdf(file="RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
#读取基因注释文件，标注基因的名称
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()




#铜死亡基因相关性分析
#引用包
library(corrplot)
library(circlize)
inputFile="CupGeneExp.txt"    #输入文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)
#计算基因间相关系数
cor1=cor(rt)
#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))
#绘制圈图
pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
#绘制图例
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()

#绘制相关性图形
pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col=colorRampPalette(c("green", "white", "red"))(50)
         )
dev.off()








#求致病基因和分型基因交集
#install.packages("VennDiagram")
library(VennDiagram)       #引用包
diseaseFile="hubGenes_MMblack&brown.txt"                 #疾病共表达分析的结果文件
clusterFile="cluster.hubGenes_MMbrown&turquoise.txt"         #分型共表达分析的结果文件
geneList=list()
#读取疾病共表达分析的结果文件
rt=read.table(diseaseFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["Disease WGCNA"]]=uniqGene     #将疾病WGCNA的核心基因保存到geneList里面
#读取分型共表达分析的结果文件
rt=read.table(clusterFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["Cluster WGCNA"]]=uniqGene     #把分型WGCNA的核心基因保存到geneList里面
#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex=1)
pdf(file="venn.pdf", width=3, height=3)
grid.draw(venn.plot)
dev.off()
#输出交集核心基因的列表
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)



#机器学习筛选特征基因
#引用包
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)
set.seed(123)      #设置种子
inputFile="normalize.txt"      #表达数据文件
geneFile="interGenes.txt"      #基因列表文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
write.table(data, file="interExp.txt", sep="\t", quote=F, col.names=T,row.names=T)

row.names(data)=gsub("-", "_", row.names(data))
#获取样品分组信息
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group
#对数据进行分组
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]
#RF随机森林树模型
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)
#SVM机器学习模型
mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)
#XGB模型
mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)
#GLM模型
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)
#定义预测函数
p_fun=function(object, newdata){
	predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)
#RF随机森林树模型预测结果
explainer_rf=explain(mod_rf, label = "RF",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#SVM机器学习模型预测结果
explainer_svm=explain(mod_svm, label = "SVM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm=model_performance(explainer_svm)
#XGB模型预测结果
explainer_xgb=explain(mod_xgb, label = "XGB",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)
#GLM模型预测结果
explainer_glm=explain(mod_glm, label = "GLM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_glm=model_performance(explainer_glm)
#绘制四种方法的残差反向累计分布图
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()
#绘制四种方法的残差箱线图
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()
#绘制ROC曲线
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
	   c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
	     paste0('SVM: ',sprintf("%.03f",roc2$auc)),
	     paste0('XGB: ',sprintf("%.03f",roc3$auc)),
	     paste0('GLM: ',sprintf("%.03f",roc4$auc))),
	   col=c("red","blue","green","yellow"), lwd=2, bty = 'n')
dev.off()

#对四种方法进行基因的重要性分析,得到四种方法基因重要性评分
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
#绘制基因重要性图形
pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()
#输出重要性评分最高的基因
geneNum=113     #设置基因的数目
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)





#验证组
inputFile="GSE7890.normalize.txt"      #表达数据文件
geneFile="AUC》0.65.txt"      #基因列表文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
write.table(data, file="GSE7890Exp.txt", sep="\t", quote=F, col.names=T,row.names=T)
#ROC曲线
library(pROC)              
inputFile="input4.txt"      
outFile="ROC4.pdf"         
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)  
y=colnames(rt)[1]
#定义颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
	bioCol=rainbow(ncol(rt))}
#绘制
pdf(file=outFile,width=5,height=5)
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(rt)){
	roc1=roc(rt[,y], as.vector(rt[,i]))
	lines(roc1, col=bioCol[i-1])
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()
aucText

inputFile="GSE92566.normalize.txt"      #表达数据文件
geneFile="AUC》0.65.txt"      #基因列表文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
write.table(data, file="GSE92566Exp.txt", sep="\t", quote=F, col.names=T,row.names=T)
#ROC曲线
library(pROC)              
inputFile="input5.txt"      
outFile="ROC5.pdf"         
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)  
y=colnames(rt)[1]
#定义颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
	bioCol=rainbow(ncol(rt))}
#绘制
pdf(file=outFile,width=5,height=5)
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(rt)){
	roc1=roc(rt[,y], as.vector(rt[,i]))
	lines(roc1, col=bioCol[i-1])
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()
aucText

inputFile="GSE121618.normalize.txt"      #表达数据文件
geneFile="AUC》0.65.txt"      #基因列表文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#读取基因列表文件,提取交集核心基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
write.table(data, file="GSE121618Exp.txt", sep="\t", quote=F, col.names=T,row.names=T)
#ROC曲线
library(pROC)              
inputFile="input6.txt"      
outFile="ROC6.pdf"         
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)  
y=colnames(rt)[1]
#定义颜色
bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
	bioCol=rainbow(ncol(rt))}
#绘制
pdf(file=outFile,width=5,height=5)
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(rt)){
	roc1=roc(rt[,y], as.vector(rt[,i]))
	lines(roc1, col=bioCol[i-1])
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()
aucText


#PCA分析

##PCA分析
#引用包
library(limma)
library(ggplot2)
clusterFile="input2-1.txt"      #手动转置并在尾列添加分型信息
#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
#设置图形的颜色
bioCol=c("#00AFBB","#E7B800")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
#绘制PCA图形
pdf(file="PCA2.pdf", width=4, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


 
##PCA分析
#引用包
library(limma)
library(ggplot2)
clusterFile="input3-1.txt"      #手动转置并在尾列添加分型信息
#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
#设置图形的颜色
bioCol=c("#00AFBB","#E7B800")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
#绘制PCA图形
pdf(file="PCA3.pdf", width=4, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


 
 
##PCA分析
#引用包
library(limma)
library(ggplot2)
clusterFile="input4-1.txt"      #手动转置并在尾列添加分型信息
#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
#设置图形的颜色
bioCol=c("#00AFBB","#E7B800")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
#绘制PCA图形
pdf(file="PCA4.pdf", width=4, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


 
##PCA分析
#引用包
library(limma)
library(ggplot2)
clusterFile="input5-1.txt"      #手动转置并在尾列添加分型信息
#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
#设置图形的颜色
bioCol=c("#00AFBB","#E7B800")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
#绘制PCA图形
pdf(file="PCA5.pdf", width=4, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


 
##PCA分析
#引用包
library(limma)
library(ggplot2)
clusterFile="input6-1.txt"      #手动转置并在尾列添加分型信息
#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
Cluster=as.vector(rt[,ncol(rt)])
#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)
PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)
#设置图形的颜色
bioCol=c("#00AFBB","#E7B800")
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]
#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))), Cluster=g))
}
#绘制PCA图形
pdf(file="PCA6.pdf", width=4, height=3)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
	scale_colour_manual(name="Cluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()




##GO, KEGG富集分析
#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件
#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)     #读取输入文件
#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
#GO富集分析
go=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(go)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}
#柱状图
pdf(file="GO_barplot.pdf", width=9, height=7)
bar=barplot(go, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()	
#气泡图
pdf(file="GO_bubble.pdf", width=9, height=7)
bub=dotplot(go, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件
#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)     #读取输入文件
#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
#KEGG=as.data.frame(kk)
#KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
#KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(kk, file="kk.txt", sep="\t", quote=F, row.names = F)
#定义显示通路的数目
showNum=9
if(nrow(kk)<showNum){
	showNum=nrow(kk)
}
#柱状图
pdf(file="KEGG_barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()
#气泡图
pdf(file="KEGG_bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, color=colorSel)
dev.off()



#GSEA_KEGG分析
#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
inputFile="all.txt"         #输入文件
gmtFile="c2.cp.kegg.v7.4.symbols.gmt"      #基因集文件
#读取文件,并对输入文件进行整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)
#读入基因集文件
gmt=read.gmt(gmtFile)
#富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)
#输出实验组富集的图形
termNum=5      #展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")
	pdf(file="GSEA.treat.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}
#输出正常组富集的图形
termNum=3      #展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA.con.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

#GSEA_GO分析
#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
inputFile="all.txt"         #输入文件
gmtFile="c5.go.v7.4.symbols.gmt"      #基因集文件
#读取文件,并对输入文件进行整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)
#读入基因集文件
gmt=read.gmt(gmtFile)
#富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA_GO.result.txt",sep="\t",quote=F,row.names = F)
#输出实验组富集的图形
termNum=5      #展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")
	pdf(file="GSEA_GO.treat.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}
#输出正常组富集的图形
termNum=5      #展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA_GO.con.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}




#CIBERSORT免疫浸润分析
inputFile="merge.txt"      #表达数据文件
source("CIBERSORT.R")       #引用包
#免疫细胞浸润分析
outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)



#免疫细胞的组间差异
#引用包
library(reshape2)
library(ggpubr)
inputFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
#读取免疫细胞浸润文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#对样品进行分组
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))
#绘制柱状图
pdf(file="barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()

##################绘制箱线图##################
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Type))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()





#铜死亡基因的表达数据
library(limma)              #引用包
expFile="merge_group.txt"     #表达数据文件
geneFile="gene.txt"         #基因列表文件
#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#读取基因列表文件，获取铜死亡相关基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
#输出铜死亡基因的表达数据
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="CupGeneExp.txt", sep="\t", quote=F, col.names=F)



##免疫细胞相关性分析
library(corrplot)           
inputFile="CIBERSORT-Results.txt"       
rt=read.table(inputFile,sep="\t",header=T,row.names=1)      #读取文件
#rt=t(rt)      #数据转置
M=cor(rt)     #相关型矩阵
#绘制相关性图形

pdf(file="corpot2.pdf",width=10,height=10)
corrplot(M,
         order="original",
         method = "color",
		 number.cex = 0.7, #相关系数
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("#00AFBB", "white", "#E7B800"))(50))
dev.off()



##免疫细胞和铜死亡基因的相关性
#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
expFile="CupGeneExp.txt"           #表达数据
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)
#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
#相关性分析
outTab=data.frame()
for(cell in colnames(immune)){
	if(sd(immune[,cell])==0){next}
	for(gene in colnames(data)){
		x=as.numeric(immune[,cell])
		y=as.numeric(data[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
	}
}
#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #去掉背景
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 8, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()




#特征基因的表达数据
library(limma)              #引用包
expFile="merge_group.txt"     #表达数据文件
geneFile="AUC》0.65.txt"         #基因列表文件
#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#读取基因列表文件，获取特征基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]
#输出铜死亡基因的表达数据
out=rbind(ID=colnames(geneExp), geneExp)
write.table(out, file="AUC》0.65Exp.txt", sep="\t", quote=F, col.names=F)



##免疫细胞和特征基因的相关性
#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
expFile="AUC》0.65Exp.txt"           #表达数据
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)
#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
#相关性分析
outTab=data.frame()
for(cell in colnames(immune)){
	if(sd(immune[,cell])==0){next}
	for(gene in colnames(data)){
		x=as.numeric(immune[,cell])
		y=as.numeric(data[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
	}
}
#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="AUC》0.65Exp.cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #去掉背景
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 8, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()





