library(SIMLR)
library(igraph)
library(randomForest)
library(dplyr)
library(Seurat)
library(SCMarker)
library(M3Drop)

rnaseq_processing<-function(data,eratio=0.06,log2=T,median=T,expressionThreshold=2){
  #data: genes in rows and samples in column
  
  # rm cnst genes
  kvar = which(apply(data,1,sd)>0)
  data=data[kvar,]
  
  # control for expression ratio
  if (eratio>0){
    minsample = floor(ncol(data)*eratio)
    kgood = rowSums(data > expressionThreshold) > minsample
    data = data[kgood,]
  }
  
  # log2-transform
  if (log2){data=log2(data+1)}
  
  # median-center
  # if (median){data=sweep(data,1,apply(data,1,median))}
  
  return(data)
}



RFCell <- function(mydata){
  # 扰动
  x=mydata
  x=t(x)
  x=as.data.frame(x)
  
  
  z<-matrix(nrow=nrow(x),ncol=1)
  for(i in 1:ncol(x))
  {
    
    t<- sample(x[[i]],length(x[[i]]))
    t=as.data.frame(t)
    
    z=cbind(z,t)
  }
  z=z[,-1]
  
  # 分正负样本
  xlabel=rep(1,time=nrow(x))
  zlabel=rep(0,time=nrow(x))
  xlabel=as.data.frame(xlabel)
  zlabel=as.data.frame(zlabel)
  
  
  names(xlabel)=c('1')
  names(zlabel)=c('1')
  y=rbind(xlabel,zlabel)
  
  names(z)=c(1:ncol(x))
  names(x)=c(1:ncol(x))
  xnew=rbind(x,z)
  
  cnames=paste("x",1:ncol(x),sep="")
  colnames(xnew)=cnames
  
  colnames(y)=c('y')
  
  train=cbind(xnew,y)
  train$y=as.factor(train$y) 
  #train$y
  
  # 随机森林选基因
  
  #cnames=paste("x",1:13961,sep="")
  #colnames(train)=cnames
  rf=randomForest(y ~ ., data=train, importance=TRUE, proximity=TRUE)
  #rf['type']
  
  rf1=rf$importance
  rf1=as.data.frame(rf1)
  c=rf1[order(-rf1$MeanDecreaseAccuracy),]

  ca <- subset(c, MeanDecreaseAccuracy > 0)
  gene=rownames(ca)
  #gene=as.data.frame(gene)
  #gene=gene[[1]]
  
  cnames=paste("x",1:ncol(x),sep="")
  colnames(x)=cnames
  newdata=x[gene]
  newdata=t(newdata)
  newdata=as.data.frame(newdata)
  ex=newdata
  return(ex)
}

# SCMarker

SCMarker <- function(mydata){
  
  Res <- ModalFilter(data = mydata, geneK = 10, cellK = 10)
  Res <- GeneFilter(obj = Res)
  Res <- getMarker(obj = Res, k = 300, n = 30)
  scMarker_genes <- Res$marker
  x=mydata
  x=t(x)
  x=as.data.frame(x)
  ex=x[scMarker_genes]
  ex=t(ex)
  ex=as.data.frame(ex)
  return(ex)
}


# Seurat （HVG）

HVG <- function(mydata){
  s=nrow(mydata)*0.1
  s=floor(s)
  pbmc_small <- CreateSeuratObject(mydata)
  pbmc_small <- NormalizeData(object = pbmc_small)
  pbmc_small <- FindVariableFeatures(object = pbmc_small,nfeatures = s)
  pbmc_small <- ScaleData(
    object = pbmc_small
  )
  ex=pbmc_small@assays[["RNA"]]@scale.data
  ex=as.data.frame(ex)
  return(ex)
}


# M3Drop
M3Drop <- function(mydata){
  uso_list <- M3Drop::M3DropCleanData(
    mydata,
    labels = colnames(mydata),
    min_detected_genes = 100,
    is.counts = FALSE
  )
  expr_matrix <- uso_list$data
  ex=M3Drop::M3DropFeatureSelection(expr_matrix,mt_method = "fdr",mt_threshold = 0.01)
  ex=ex[[1]]
  x=mydata
  x=t(x)
  x=as.data.frame(x)
  ex=x[ex]
  ex=t(ex)
  ex=as.data.frame(ex)
  return(ex)
}


# Expr
Expr <- function(mydata){
  meanExp = apply(mydata,1,mean)
  meanExp = meanExp[order(meanExp,decreasing = TRUE)]
  s=nrow(mydata)*0.1
  s=floor(s)
  HAEgenes = names(meanExp)[1:s]
  x=mydata
  x=t(x)
  x=as.data.frame(x)
  HAEgenes=x[HAEgenes]
  HAEgenes=t(HAEgenes)
  ex=as.data.frame(HAEgenes)
  return(ex)
}


evalcluster<-function(truelabel,predlabel){
  if(length(truelabel)!=length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)
  #Mutual information
  MI = 0.0
  for(idx in x_ids){
    for(idy in y_ids){
      idxOccur = which(truelabel==idx)
      idyOccur = which(predlabel==idy)
      idxyOccur = intersect(idxOccur,idyOccur)
      if(length(idxyOccur)>0){
        MI = MI + (length(idxyOccur)/total)*log2((length(idxyOccur)*total)/(length(idxOccur)*length(idyOccur)));
      }
    }
  }
  
  #Normalized Mutual information
  Hx = 0; #Entropies
  for(idx in x_ids){
    idxOccurCount = length(which(truelabel==idx));
    Hx = Hx - (idxOccurCount/total) * log2(idxOccurCount/total);
  }
  Hy = 0;#Entropies
  for(idy in y_ids){
    idyOccurCount = length(which(predlabel==idy));
    Hy = Hy - (idyOccurCount/total) * log2(idyOccurCount/total);
  }
  nmi = 2 * MI / (Hx+Hy)
  
  #(adjusted) Rand Index
  tab = table(truelabel,predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
  ari=c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2)
  
  out = c(nmi,ri,ari)
  names(out)=c("NMI","RI","ARI")
  return(out)
}



setwd("E:\\2020\\cluster\\code\\a_rfcell\\a_rfcell\\data\\")
files <- list.files('.', pattern="*.Rdata")
res = data.frame()

for (file in files){
  load(file)
  
  # 1.数据预处理
  mydata=rnaseq_processing(data,eratio=0.06,log2=T,median=T,expressionThreshold=2)
  
  # 2.基因xuanze
  newdata4 = RFCell(mydata)
  newdata3 = SCMarker(mydata)
  # newdata3 = M3Drop(mydata)
  newdata2 = HVG(mydata)
  newdata1 = Expr(mydata)
  
  data_all = list(newdata1, newdata2, newdata3, newdata4)
  # 3.最终聚类
  set.seed(11111)
  ARI = c()
  NMI = c()
  dataset = c()
  methods = c("Expr", "HVG", "SCMarker", "RFCell")
  for (i in c(1:4)){
    newdata = as.data.frame(data_all[i])
    example = SIMLR(X = newdata, c = length(table(label)),normalize = TRUE)
    # "NMI","RI","ARI"
    cluster_res = evalcluster(as.matrix(label), example$y$cluster)
    ARI = c(ARI, cluster_res[3])
    NMI = c(NMI, cluster_res[1])
    dataset = c(dataset, file)
  }
  tmp = as.data.frame(ARI)
  tmp['NMI'] = NMI
  tmp['dataset'] = dataset
  tmp['methods'] = methods
  
  res = rbind(tmp,res)
}

write.csv(res,"res.csv")


