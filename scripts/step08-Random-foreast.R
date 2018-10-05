## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-08-10 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-10  First version
### Update Log: 2018-10-10  second version
###
### ---------------


### https://github.com/jmzeng1314/ML

rm(list=ls())
options(stringsAsFactors = F)
Sys.setenv(R_MAX_NUM_DLLS=999)

Rdata_dir='../Rdata/'
Figure_dir='../figures/'
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
table(group_list)

exprSet=na.omit(expr)

load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
head(phe)
exprSet[1:4,1:4]
 

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)

x=t(log2(exprSet+1))
y=phe$event
table(y)

tmp = as.vector(table(y))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)
tmp_rf=file.path(Rdata_dir,'TCGA_KIRC_miRNA_rf_output.Rdata')

if(!file.exists(tmp_rf)){
  rf_output=randomForest(x=x, y=y,importance = TRUE, ntree = 10001, proximity=TRUE )
  save(rf_output,file = tmp_rf)
}
load(file = tmp_rf)
rf_output
str(rf_output)


rf_importances=importance(rf_output, scale=FALSE)
head(rf_importances)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, 
           main="Variable Importance (Gini) for top 30 predictors",cex =.7)
target_labels=as.vector(y)
MDSplot(rf_output, y, k=2, xlab="", ylab="", 
        pch=target_labels, palette=c("red", "blue"), main="MDS plot")

print(rf_output) 
choose_gene=rownames(tail(rf_importances[order(rf_importances[,2]),],50))

library(pheatmap) 
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]
choose_matrix=t(scale(t(log2(choose_matrix+1)))) 
## http://www.bio-info-trainee.com/1980.html
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)
pheatmap(choose_matrix,show_colnames = F,annotation_col = annotation_col,
         filename = 'rf_genes.heatmap.png')


library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('rf_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()

 







