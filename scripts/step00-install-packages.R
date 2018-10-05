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
###
### ---------------

# https://github.com/jmzeng1314/biotrainee 

## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 强调，不是所有的R包都需要安装成功的。
## 失败就失败，大不了从头再来，卸载R语言，从新开始。

## 强调，中国大陆的粉丝务必注意下载镜像。
## 强调，管是什么电脑，都请务必安装好R及Rstudio哦
# 所有的软件都安装在c盘哦，然后系统用户名最好是不要用中文，写代码最怕中文字符串哦！
# 生信0基础第一步，下载R和Rstudio并且安装在自己的电脑上面。官网链接是 
# - R: https://mirrors.tuna.tsinghua.edu.cn/CRAN/
# - RStudio：https://www.rstudio.com/products/rstudio/download/#download 
# 如果你的网络不好，可以从我整理的网盘下载，链接：https://share.weiyun.com/5hW6VAA  密码：3fuhrm

Sys.setenv(R_MAX_NUM_DLLS=999)
# 首先需要更改一些镜像配置
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
source("http://bioconductor.org/biocLite.R") 
## 如果你的网络实在是太差，试试看：
# install.packages("BiocInstaller",repos="http://bioconductor.org/packages/3.7/bioc")  
## 很可惜你在中国大陆，不得不承受这个痛苦。
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
library('BiocInstaller')

## 检查镜像是否修改成功。
options()$BioC_mirror
options()$repos
if(F){
  install.packages('shiny')
  install.packages(c("devtools","ggplot2","pheatmap"))
  install.packages(c("ggpubr","ggstatsplot"))
  
  #source("http://bioconductor.org/biocLite.R") 
  library('BiocInstaller') 
  biocLite(c('airway','DESeq2','edgeR','limma')) 
  biocLite('clusterProfiler')
  ## 到这里，R里面已经有了 234 个包
}
install.packages(c("devtools","reshape2","pheatmap",
                   "ggplot2","ggfortify","stringr",
                   "survival","survminer","lars",
                   "glmnet","timeROC","ggpubr",
                   "randomForest","ROCR","Hmisc",
                   "caret","ggstatsplot","tableone", 
                   "devtools","reshape2","randomForest"))

library(devtools) 

if(! require('edgeR')){
  biocLite(c('airway','DESeq2','edgeR','limma'))
}
if(! require("maftools")) biocLite("maftools")
if(! require("genefilter")) biocLite("genefilter")

if(! require("CLL")) biocLite("CLL")
if(! require("org.Hs.eg.db")) biocLite('org.Hs.eg.db')
library(BiocInstaller)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
if(! require("maftools")) biocLite("maftools")
if(! require("RTCGA")) biocLite("RTCGA")
if(! require("RTCGA.clinical")) biocLite("RTCGA.clinical")
# https://bioconductor.org/packages/3.6/data/experiment/src/contrib/RTCGA.clinical_20151101.8.0.tar.gzn
if(! require("RTCGA.miRNASeq")) biocLite("RTCGA.miRNASeq")


# Then from : https://github.com/ShixiangWang 
# You don't need run the codes below, I will explain to you face to face.

if(F){
  source("http://bioconductor.org/biocLite.R") 
  packs = c("devtools", "reshape2", "ggplot2", "pheatmap", "ggfortify", "stringr", "survival",
            "survminer", "lars", "glmnet", "timeROC", "ggpubr", "randomForest", "ROCR", "genefilter",
            "Hmisc", "caret", "airway","DESeq2","edgeR","limma", "CLL", "org.Hs.eg.db", "maftools")
  if(! require(pacman)) install.packages("pacman", dependencies = TRUE)
  pacman::p_load(packs, dependencies=TRUE, character.only = TRUE)
  # check
  pacman::p_loaded(packs, character.only = TRUE)
  all(pacman::p_loaded(packs, character.only = TRUE))
}






