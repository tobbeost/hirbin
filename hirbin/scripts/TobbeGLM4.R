

oGlmAnovaTest4 <- function(Y,group,N,nR=dim(Y)[1]) {
  #A simple wrapper to fit a overdispersed poisson glm model
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