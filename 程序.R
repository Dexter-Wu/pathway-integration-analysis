######
##分析
##j-GRP

compUp<-function(d1,d2,ad=0.9){#d1-实验组矩阵，d2-对照组；上调控
  up=apply(d1,2,function(x){rowSums(x>=d2)/ncol(d2)>=ad});
  up=cbind(up,apply(d2,2,function(x){rowSums(x<d1)/ncol(d1)>=ad}));
  return(up)
}
compDown<-function(d1,d2,ad=0.9){##d1-实验组矩阵，d2-对照组；下调控
  down=apply(d1,2,function(x){rowSums(x<=d2)/ncol(d2)>=ad});
  down=cbind(down,apply(d2,2,function(x){rowSums(x>d1)/ncol(d1)>=ad}));
  return(down)
}

regulation<-function(expr,samples,ad=0.9,miss=0){##miss=1随机置换实验
  n=nrow(expr[[1]])  ## number of genes (features)
  s=nrow(samples)  ## number of studies
  
  up=NULL;
  down=NULL;  
  upg=array(NA,c(n,1))#上调控比例
  downg=array(NA,c(n,1))
  if(miss==0){
    for(i in 1:s){#研究
      s1=samples[i,1]
      s2=samples[i,1]+samples[i,2]
      upt=compUp(expr[[i]][,1:s1],expr[[i]][,(s1+1):s2],ad)
      downt=compDown(expr[[i]][,1:s1],expr[[i]][,(s1+1):s2],ad)
      up=cbind(up,upt); 
      down=cbind(down,downt);
    }
    upg=rowSums(up)/ncol(up);
    downg=rowSums(down)/ncol(down);
    pd=upg-downg
  }else{
    for(i in 1:s){#研究
      s1=samples[i,1]
      s2=samples[i,1]+samples[i,2]
      label=c(1:s2)
      label=sample(label)###打乱标签
      upt=compUp(expr[[i]][,label[c(1:s1)]],expr[[i]][,label[c((s1+1):s2)]],ad)
      downt=compDown(expr[[i]][,label[c(1:s1)]],expr[[i]][,label[c((s1+1):s2)]],ad)
      up=cbind(up,upt); 
      down=cbind(down,downt);
    }
    
    upg=rowSums(up)/ncol(up);
    downg=rowSums(down)/ncol(down);
    pd=upg-downg
  }
  return(pd)
}



#list结构
jgrp<-function(expr,lab,ad=0.9,nper=1000){##lab: a factor variable of {"control","treatment"}
  studys=NULL
  for(i in 1:length(expr)){
    studys[[i]]=cbind(expr[[i]][,c(which(lab[[i]]=="control"))],expr[[i]][,c(which(lab[[i]]=="treatment"))])
  }
  samples=array(0,c(length(expr),2))
  for(i in 1:length(expr)){
    samples[i,1]=length(which(lab[[i]]=="control"))
    samples[i,2]=length(which(lab[[i]]=="treatment"))
  }
  expr=studys
  pd=regulation(expr,samples,ad)#正常序列做的真实pd
  pm=NULL;
  pn=list();
  for(j in 1:nper){
    a=regulation(expr,samples,ad,1)
    pm=c(pm,a)#打乱序列做的pd
    pn[[j]]=a
    row.names(pn[[j]])=rownames(pm)
  }
  pm=abs(pm);
  pval=sapply(abs(pd),function(x){length(pm[pm>=x])/length(pm)})
  pval_BH=p.adjust(pval,method = "BH",n=length(pval))
  return(list(pd=pd,pval=pval,pval_BH=pval_BH,pn=pn))
}


Sort_H <- function(RWH){
  #输入参数为RWH,即经过非负矩阵分解获得的W和H矩阵组成的list;
  #此函数为提取非负矩阵分解后的H矩阵;
  ###
  len=length(RWH);
  H=list();
  for(i in 1:len){
    H[[i]]=as.data.frame(RWH[[i]][[2]])  #提取H矩阵
    rownames(H[[i]])=c('a','b','c','d','e')
  }
  return(H)
}


Test_addcltfunc<-function(input,obsu,cutoff=0.05){
  ###obsu:真实结果
  ###input:P-value
  ###cutoff:为不同的临界值(threshold)
  ##
  ##此函数计算混淆矩阵
  ##计算Accracy,Sensitivity,Specificity
  ##计算TPR和FNR
  test_result=list()
  tp <- sum(input <= cutoff & obsu == 1)
  fp <- sum(input <= cutoff & obsu == 0)
  fn <- sum(input > cutoff & obsu == 1)
  tn <- sum(input > cutoff & obsu == 0)
  
  ###构建混淆矩阵
  mi_matrix=rbind(cbind(tp,fp),cbind(fn,tn))
  colnames(mi_matrix)=c("positive","negetive")
  rownames(mi_matrix)=c("Ture","False")
  ####
  n=length(input)
  ###准确率
  Accracy = (tp+tn)/n
  ####Sensitivity
  Sensitivity = tp/(tp+fn)
  ####Specificity
  Specificity = tn/(fp+tn)
  ##假阳性率
  FPR=fp/(fp+tn)
  ##真阳性率
  FNR=fn/(tp+fn)
  ###合并成列表
  test_result=list(mi_matrix=mi_matrix,Accracy=Accracy,Sensitivity=Sensitivity,Specificity=Specificity,FPR=FPR,FNR=FNR);
  return(test_result)
}


####=========================
#####=============
##########数据导入
load("C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/模拟数据集/NMF/NMF_a=0.5_n1=20_seed=8_lamb=0.3.Rdata")
####提取H矩阵
H=Sort_H(RWH)
pathname=c('a','b','c','d','e')
###合并Lab
Lab[Lab=="2"]="treatment"
Lab[Lab=="1"]="control"
Lab_GRP <- list(Lab,Lab,Lab,Lab,Lab,Lab,Lab,Lab,Lab,Lab)


###========================
###===========
#j-RGP通路分析
jGRP_0.7 <- jgrp (H,Lab_GRP,ad=0.7,nper=1000)
save(jGRP_0.7,file = "C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/j_GRP结果/jGRP_0.7_a=0.5_n1=20_seed=8_lamb=0.3.Rdata")
#load("C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/j_GRP结果/jGRP_0.7_a=0.5_n1=20_seed=10_lamb=0.1.Rdata")


#########构造真实label,即真实差异通路,差异表达通路记为1
obsu <- c(1,1,1,1,0)
input=jGRP_0.7[["pval"]]
###计算各个评估量与混淆矩阵
test_result=Test_addcltfunc(input,obsu,cutoff = 0.05)
save(test_result,file = "C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/test结果/test_a=0.5_n1=20_seed=8_lamb=0.3.Rdata")
#load("C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/test结果/test_a=0.5_n1=20_seed=10_lamb=0.1.Rdata")
test_result


#use package PORC to analyze the result of j-GRP
#install.packages("pROC")
###计算AUC
library(pROC)
############
###############
#draw picture 2
load("C:/Users/Administrator/Desktop/模拟数据结果分析/第二次分析/j_GRP结果/jGRP_0.7_a=0.5_n1=20_seed=1_lamb=0.1.Rdata")
input=jGRP_0.7[["pval"]]
roc1 <- roc(obsu, input, main="Statistical comparison", col="1")
#显示auc
#roc1[["auc"]]
