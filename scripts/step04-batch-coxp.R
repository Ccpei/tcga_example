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

### https://github.com/jmzeng1314/GEO/blob/master/GSE11121/step5-surivival.R

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

load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
head(phe)
exprSet[1:4,1:4]



library(survival)
library(survminer)

## 批量生存分析 使用 coxph 回归方法
# http://www.sthda.com/english/wiki/cox-proportional-hazards-model
colnames(phe)
mySurv=with(phe,Surv(time, event))
cox_results <-apply(exprSet , 1 , function(gene){
  # gene= exprSet[1,]
  group=ifelse(gene>median(gene),'high','low') 
  survival_dat <- data.frame(group=group,stage=phe$stage,age=phe$age,
                             gender=phe$gender,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ gender + age + stage+ group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results[cox_results[,4]<0.05,]

## 批量生存分析 使用  logrank test 方法
mySurv=with(phe,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  # gene=exprSet[1,]
  phe$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

require("VennDiagram")
VENN.LIST=list(cox=rownames(cox_results[cox_results[,4]<0.05,]),
             log=names(log_rank_p[log_rank_p<0.05]))
venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          fill=c("darkmagenta", "darkblue"), 
                          alpha=c(0.5,0.5), cex = 2, 
                          cat.fontface=4,  
                          main="overlap of coxph and log-rank test")
grid.draw(venn.plot) 

save(log_rank_p,cox_results ,
     file = 
       file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_results.Rdata')
)
 
library(pheatmap)
choose_gene=rownames(cox_results[cox_results[,4]<0.05,])
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]
choose_matrix=t(scale(t(log2(choose_matrix+1)))) 
## http://www.bio-info-trainee.com/1980.html
annotation_col = data.frame( group_list=group_list  )
rownames(annotation_col)=colnames(expr)
pheatmap(choose_matrix,show_colnames = F,annotation_col = annotation_col,
         filename = 'cox_genes.heatmap.png')


library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('cox_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()














