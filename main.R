### ---------------
###
### Create: Jianming Zeng
### Date: 2018-10-05 15:06:55
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log:  2018-10-05 15:06:55  First version
###
### ---------------
setwd('scripts/')

### ---------------
###
### install the packages 
###
### ---------------
#source('step00-install-packages.R')
# we don't need to run it by source.

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R') 
dim(expr)
dim(meta)

### ---------------
###
### Do DEG for miRNA expression matrix by DESeq2,edgeR,limma
###
### ---------------
source('step02-DEG-3-packages.R')

### ---------------
###
###  
###
### ---------------
source('step03-batch-logRank.R')
dim(expr)
dim(meta)


### ---------------
###
###  
###
### ---------------
source('step04-batch-coxp.R')
 

### ---------------
###
### 
###
### ---------------
source('step05-lasso.R') 


### ---------------
###
###  
###
### ---------------
source('step06-coxph-forest.R') 

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)

### ---------------
###
### Get miRNA expression matrix from RTCGA package 
###
### ---------------
source('step01-getData-from-RTCGA.R')
dim(expr)
dim(meta)




