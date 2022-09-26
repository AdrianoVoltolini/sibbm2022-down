
raw <- read.csv("NormalizedCount.txt",sep="\t")
newdf <- data.frame(raw,row.names = 1)

nome_controlli <- c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10",
                    "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
                    "C21","C22","C23")
nome_down <- c("D01","D02","D03","D04","D05","D06","D07","D08","D09","D10",
               "D11","D12","D13","D14","D15","D16","D17")

colnames(newdf) <- c(nome_controlli, nome_down)

ex2 <- log2(newdf+0.001)

boxplot(ex2,main = "Dataset Boxplot", las=2, ylab = "Log(Signal)")

newdf <- newdf[,-24]
ex2 <- log2(newdf+0.001)

boxplot(ex2,main = "Dataset Boxplot", las=2, ylab = "Log(Signal)")


#PCA

library(scatterplot3d)

pca <- prcomp(t(ex2))
summary(pca)
grpcol <- c(rep("blue",23), rep("red",16))

source('addgrids3d.r')

tmp <- scatterplot3d(pca$x[,1], pca$x[,3], pca$x[,2], xlab="PCA1", ylab="PCA3",
       zlab="PCA2",pch = "",color=grpcol, box = FALSE, grid = FALSE)
addgrids3d(pca$x[,1], pca$x[,3], pca$x[,2],grid = c("xy", "xz", "yz"))

tmp$points3d(pca$x[,1], pca$x[,3], pca$x[,2], pch = 16,col=grpcol)

legend("topright",tmp$xyz.convert(7.5, 3, 4.5), legend = c("Control","Down"),
       col =  c("blue", "red"), pch = 16)


#K-means

library(factoextra)

km.out <- kmeans(t(ex2), 2, nstart=100)
km.out$cluster

fviz_cluster(km.out, data=t(newdf), 
             palette="jco", # this loads one of the many predefined palettes
             xlab="",
             ylab="", main = "K-Means Clustering with K = 2")

#Hierarchical Clustering

dist.matrix <- dist(t(newdf))
hc.complete <- hclust(dist.matrix, method="complete")

fviz_dend(hc.complete, main="Complete Linkage Hierarchical Clustering with K = 2", cex=0.7,
          k=2,type = "rectangle",
          palette="jco")


#filtraggio geni

library(genefilter)
set.seed(80085)
f <- factor(c(rep("CTRL",23), rep("DOWN",16)))
tt.newdf <- rowttests(as.matrix(newdf),f) #confronta la media dei due gruppi per ogni gene, con t test
adjusted_p.value <- p.adjust(tt.newdf$p.value,"fdr")
tt.newdf[,"adjusted p.value"] <- adjusted_p.value
keepers.newdf <- which(tt.newdf$p.value<0.1)
ex.newdf <- newdf[keepers.newdf,]# prendi geni con diff significativa tra medie dei gruppi (pvalue < 0.1)
dat.newdf <- cbind(as.data.frame(t(ex.newdf)),f) 
colnames(dat.newdf)[ncol(dat.newdf)] <- "DOWN"


#SUPERVISED LEARNING

library(caret)
set.seed(80085)

ctrl <- trainControl(method = "repeatedcv", number = 13, repeats = 3, savePredictions = TRUE,
                     classProbs = TRUE) #13-fold cross validation, repeated 3 times

#random forest
set.seed(80085)
fit.rf.newdf  <- train(DOWN ~ ., data = dat.newdf, method = "rf",ntree = 1000,
                       trControl = ctrl, metric = "Accuracy")

#lda
set.seed(80085)
fit.lda.newdf <- train(DOWN ~ ., data = dat.newdf, method = "lda",
                       trControl = ctrl, metric = "Accuracy")   

library(pROC)
plot.roc(as.numeric(fit.lda.newdf$pred$pred),as.numeric(fit.lda.newdf$pred$obs), 
         main = "LDA AUC-ROC CURVE")

#lasso
set.seed(80085)
fit.lasso.newdf <- train(DOWN ~ ., data = dat.newdf, method = "glmnet",              
                         family = "binomial",trControl = ctrl,
                         tuneGrid = expand.grid(alpha = 1,lambda = seq(0,1,by=0.05)),
                         metric = "Accuracy")                   


plot(fit.lasso.newdf, xlab = "Lambda")
plot(fit.lasso.newdf$finalModel,"lambda")

#scudo
library(caret)
library(rScudo)
library(igraph)
set.seed(80085)

model.scudo <- scudoModel(nTop = (9:13)*5, nBottom = (13:17)*5, N = 0.5)

fit.scudo.newdf  <- train(DOWN ~ ., data = dat.newdf, method = model.scudo,
                       trControl = ctrl, metric = "Accuracy")

#sample validation with optimal parameters:
set.seed(80085)
inTrain <- createDataPartition(f, list = FALSE)
trainData <- dat.newdf[inTrain,-4787]
testData <- dat.newdf[-inTrain,-4787]


trainRes <- scudoTrain(t(trainData), groups = f[inTrain],
                       nTop = 50, nBottom = 80, alpha = 0.05)
trainNet <- scudoNetwork(trainRes, N = 0.5)
scudoPlot(trainNet, vertex.label = NA)
testRes <- scudoTest(trainRes, t(testData), f[-inTrain],
                     nTop = 50, nBottom = 80)

testNet <- scudoNetwork(testRes, N = 0.5)
scudoPlot(testNet, vertex.label = NA)
#variabili importanti

#random forest
ImpMeasure.rf<-data.frame(varImp(fit.rf.newdf,scale = FALSE)$importance)
ImpMeasure.rf$Vars<-row.names(ImpMeasure.rf)
top.rf <- ImpMeasure.rf[order(-ImpMeasure.rf$Overall),]
top200.rf <- top.rf[1:200,2]
write(top200.rf,"top200.rf.txt")

plot(top.rf[1:200,1],type = "h",ylab = "MeanDecreaseGini",
     xlab ="Top 200 most important RNAs for Random Forest Classification",
     main = "Random Forest Variable Importance")
axis(1,labels=rownames(top.rf[1:200,]),at=1:200,las=2, cex.axis=0.5)

#lda
ImpMeasure.lda<-data.frame(varImp(fit.lda.newdf, scale = FALSE)$importance)
ImpMeasure.lda$Vars<-row.names(ImpMeasure.lda)
top.lda <- ImpMeasure.lda[order(-ImpMeasure.lda$DOWN),]
top200.lda <- top.lda[1:200,3]
write(top200.lda,"top200.lda.txt")

plot(top.lda[1:200,1],type = "h",ylab = "Importance",xlab ="",
     main = "Top 200 LDA Variable Importance")
axis(1,at=1:200, labels=rownames(top.lda[1:200,]),las=2, cex.axis=0.5)

#lasso
ImpMeasure.lasso<-data.frame(varImp(fit.lasso.newdf)$importance)
ImpMeasure.lasso$Vars<-row.names(ImpMeasure.lasso)
top.lasso <- ImpMeasure.lasso[order(-ImpMeasure.lasso$Overall),]
top6.lasso <- top.lasso[top.lasso$Overall != 0.0,2]
write(top6.lasso,"top6.lasso.txt")

length(intersect(top200.lda,top6.lasso))#  6 in comune
length(intersect(top200.lda,top200.rf))  # 48 in comune
length(intersect(top200.rf,top6.lasso)) #  3 in comune



#david

david <- read.table("david.txt",sep = "\t", header = TRUE)
david.pulito <- david[1:3,c("Category","Term","Count","PValue","FDR")]

#pathfindr
library(pathfindR)
input_pathfindr <- data.frame(row.names(tt.newdf),tt.newdf[,"statistic"]*-1,
                              tt.newdf[,"adjusted p.value"])

output_df <- run_pathfindR(input_pathfindr)
output_cluster <- cluster_enriched_terms(output_df) #usa kegg. 

term_gene_graph(output_df)

#enrichnet

geni.id.50 <- read.table("top50.lda.csv",sep = ",", header = T)
input.enrichnet.50 <- geni.id.50[which(geni.id.50[,"converted_alias"] != "None"),]
write(input.enrichnet.50[,"converted_alias"],"ensembl.50.lda.txt")

XD.raw <- read.table("tissue_200.txt",sep = "\t")
XD <- XD.raw[,c("atrioventricular.node","heart")]
XD.filtered <- XD[which(XD[,"atrioventricular.node"] > 2.5 | XD[,"heart"] > 2.5),]

#heatmap

library(pheatmap)
library(grid)
library(RColorBrewer)

heat.df <- newdf[top200.lda,]
group <- c(rep("CTRL",23),rep("DOWN",16))
colnames(heat.df) <- group
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)

setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(as.matrix(heat.df),scale = "row", color=hmcol,treeheight_row = 0,
         show_rownames = F, fontsize = 6)
setHook("grid.newpage", NULL, "replace")
grid.text("Top 200 RNAs for LDA classification", x=-0.03, rot=90, gp=gpar(fontsize=14))

