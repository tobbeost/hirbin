
runAnalysis<-function(groups,reference,filename){
  groups<-as.factor(groups)
  groups<-relevel(groups,reference)
  #calculate differential expression based on overdispersed quasi-poisson glm.
  x<-read.table(filename,sep='\t',header=TRUE,row.names=1)
  x<-as.matrix(x)
  ngroup<-length(groups)
  
  N <- apply(x,2,sum) #normalizing factor: total number of reads for each sample.
  #xnorm<-apply(x,1,function(xrow){xrow/N})
  #xnorm<-t(xnorm)
  outputlist <- oGlmAnovaTest4(x,groups,N) #run glm
  pvals<-outputlist[[1]]
  logFC<-outputlist[[2]]
  stderr<-outputlist[[3]]
  names(pvals)<-rownames(x)
  names(logFC)<-rownames(x)
  names(stderr)<-rownames(x)
  fdr <- p.adjust(pvals,method="BH") #adjust pvalues
  idlist<-strsplit(rownames(x),'_')
  idlist2<-matrix(nrow=length(idlist),ncol=3)
  for(i in 1:length(idlist)){
    idlist2[i,1]<-idlist[[i]][1]
    idlist2[i,2]<-idlist[[i]][2]
    idlist2[i,3]<-idlist[[i]][3]
  }
  output<-cbind(idlist2,logFC,stderr,pvals,fdr)
  
  colnames(output)<-c("name","identity","clusterID","logFC","std. error","p.value","fdr")
  output<-as.data.frame(output)
  output<-output[!is.na(output$fdr),]
  newfilename<-gsub("abundance_matrix","results",filename)
  write.table(output,newfilename,sep='\t',quote=FALSE,row.names=FALSE)
  return(output)
}

oGlmAnovaTest4 <- function(Y,group,N,nR=dim(Y)[1]) {
  #A function to fit a overdispersed poisson glm model
  #and perform a likelihood ratio test.Takes a matrix of data as
  #input and returns the pvalue for each row
  #
  #
  #Y m by n matrix
  #group 1 by m vetor detailing which samples belong to witch group
  #ex group <- c(0,0,0,1,1,1) for a 3 vs 3 comparison
  #N normalizing constant usualy total number of reads
  
  
  pVals <- rep(NA,nR)
  logFC <- rep(NA,nR)
  stderr <- rep(NA,nR)
  
  for(i in 1:nR){
    irow<-Y[i,!is.na(Y[i,])]
    igroup<-group[!is.na(Y[i,])]
    iN<-N[!is.na(Y[i,])]
    fitPoisson <- glm(irow~igroup+offset(log(iN)),data = as.data.frame(cbind(irow,igroup)),family=quasipoisson)
    
    fitPoisson2 <- glm(irow~1+offset(log(iN)),data = as.data.frame(cbind(irow,igroup)),family=quasipoisson)
    tempAnova <- anova(fitPoisson2,fitPoisson,test="F")
    
    pVals[i] <- tempAnova[2,6]
    logFC[i] <- summary(fitPoisson)$coefficients[2,1]
    stderr[i] <- summary(fitPoisson)$coefficients[2,2]
  }
  output=list(pVals,logFC,stderr)
  return(output)
}


groups<-commandArgs(TRUE)[1]
groups<-read.table(text=groups,sep=',')
groups<-as.factor(as.matrix(groups))
reference<-commandArgs(TRUE)[2]
groups<-relevel(groups,reference)
filename<-commandArgs(TRUE)[3]
output<-runAnalysis(groups,reference,filename)

