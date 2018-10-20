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

rm(list=ls())
options(stringsAsFactors = F)

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
dim(exprSet)
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
dim(exprSet) ## remove the nomral
head(phe)
exprSet[1:4,1:4]
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)


library(lars) 
library(glmnet) 
x=t(log2(exprSet+1))
y=phe$event

model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)
print(model_lasso)
# 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。
# 它在0和1之间，越接近1说明模型的表现越好，
# 如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。
head(coef(model_lasso, s=c(model_lasso$lambda[29],0.009)))
plot.glmnet(model_lasso, xvar = "norm", label = TRUE)
plot(model_lasso, xvar="lambda", label=TRUE)
cv_fit <- cv.glmnet(x=x, y=y, alpha = 1, nlambda = 1000)
plot.cv.glmnet(cv_fit)
# 两条虚线分别指示了两个特殊的λ值:
c(cv_fit$lambda.min,cv_fit$lambda.1se) 
 
model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)
dat=as.data.frame(re[,1:2])
colnames(dat)=c('event','prob')
dat$event=as.factor(dat$event)
library(ggpubr) 
p <- ggboxplot(dat, x = "event", y = "prob",
               color = "event", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()

library(ROCR)
library(glmnet)
library(caret)
# calculate probabilities for TPR/FPR for predictions
pred <- prediction(re[,2], re[,1])
perf <- performance(pred,"tpr","fpr")
performance(pred,"auc") # shows calculated AUC for model
plot(perf,colorize=FALSE, col="black") # plot ROC curve
lines(c(0,1),c(0,1),col = "gray", lty = 4 )


fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
head(fit$beta)
#一倍SE内的更简洁的模型,是22个miRNA
#fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
#head(fit$beta)# 这里是40个miRNA
choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
length(choose_gene)
myexpr=x[,choose_gene]
mysurv=phe[,c("days","event")]
mysurv$days[mysurv$days< 1] = 1 
# 详细代码参见这个网站https://github.com/jeffwong/glmnet/blob/master/R/coxnet.R#
fit <- glmnet( myexpr, Surv(mysurv$days,mysurv$event), 
              family = "cox") 
#用包自带的函数画图
plot(fit, xvar="lambda", label = TRUE)
plot(fit, label = TRUE)
## 如果需要打印基因名，需要修改函数，这里不展开。

library(pheatmap) 
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]
choose_matrix=t(scale(t(log2(choose_matrix+1)))) 
## http://www.bio-info-trainee.com/1980.html
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)
pheatmap(choose_matrix,show_colnames = F,annotation_col = annotation_col,
         filename = 'lasso_genes.heatmap.png')


library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('lasso_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()
