groups<-as.factor(c(rep("Diabetic",15),rep("Control",15)))
groups<-relevel(groups,"Control")
#calculate differential expression based on overdispersed quasi-poisson glm.

#load the function for the GLMi
source('TobbeGLM4.R')
#read data from files
calculate.stats<-function(identityCutoff){
  x<-read.table(paste('results_filtered2/abundance_matrix_all_clust',identityCutoff,'.txt',sep=''),sep='\t',header=TRUE,row.names=1)
  x<-as.matrix(x)
  ngroup<-length(groups)
  cutoff<-round(ngroup*0.75)
  
  #apply filter
  n1<-apply(x,1,function(xrow) {return(sum(xrow>0))}) #check if number of nonzero entries is above cutoff
  test1<-(n1>=cutoff)
  mean1<-apply(x,1,mean) #test if avg count is above 10.
  test2<-(mean1>=10)
  test <- (test1 & test2)
  x.selected<-x[test,]
  
  N <- apply(x.selected,2,sum) #normalizing factor: total number of reads for each sample.
  xnorm<-x.selected/N
  outputlist <- oGlmAnovaTest4(x.selected,groups,N) #run glm
  pvals<-outputlist[[1]]
  logFC<-outputlist[[2]]
  stderr<-outputlist[[3]]
  names(pvals)<-rownames(x.selected)
  names(logFC)<-rownames(x.selected)
  names(stderr)<-rownames(x.selected)
  fdr <- p.adjust(pvals,method="BH") #adjust pvalues
  descr<-read.table('TIGRFAM_domain_info.txt',sep='\t',header=TRUE,quote='\"')
  

  idlist<-strsplit(rownames(x.selected),'_')
  idlist2<-matrix(nrow=length(idlist),ncol=3)
  for(i in 1:length(idlist)){
    idlist2[i,1]<-idlist[[i]][1]
    idlist2[i,2]<-idlist[[i]][2]
    idlist2[i,3]<-idlist[[i]][3]
  }
  output<-merge(cbind(idlist2,logFC,stderr,pvals,fdr),descr[,1:2],cbind(idlist2,as.numeric(fc),as.numeric(fdr)),by.x=1,by.y=1,all.x=TRUE)
  colnames(output)<-c("name","identity","clusterID","logFC","std. error","p.value","fdr","description")
  #write.table(output,'Results_Tigrfam_DanishvsSpanish_140625.txt',sep='\t',quote=FALSE,row.names=FALSE)
  #plot(fc,fdr)
  output<-output[!is.na(output$fdr),]
  write.table(output,paste('results_t2d/Results_Tigrfam_t2d_',identityCutoff,'_150422.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
  return(output)
}
a<-commandArgs(TRUE)[1]
output<-calculate.stats(a)

