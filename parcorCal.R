#' Calculation of the partial correlation between coding genes and epigenetic regulators.
#'
#' @param exp1 A numeric matrix containing the expression of mRNA (without 'interested.ER') with
#'   rownames and colnames.
#' @param exp2 A numeric matrix containing the expression of 'interested.ER' with
#'   rownames and colnames.
#' @param tum.pur A numeric vector indicating the samples' tumor purity with sample
#' names.
#' @param is.adjusted Whether or not the 'exp.profile' is preprocessed. If the expression values
#' were log-transformed and the rows containing many 0 values were removed, is.adjusted = TRUE.
#' 
#' @return A list including the partial correlation coefficients and p values for epigenetic regulator-gene pairs.  


parcorCal <- function(exp1,exp2,tum.pur,is.adjusted){
 
  inter_samples <- intersect(intersect(colnames(exp1),colnames(exp2)),names(tum.pur))
  tumpur <- tum.pur[inter_samples]
  if(length(tumpur)<1){stop("no same samples")}
  if (is.adjusted){
    expre1 <- exp1[,inter_samples]
    expre2 <- exp2[,inter_samples]
  }else{
    exp1_inter <- exp1[,inter_samples]
    exp2_inter <- exp2[,inter_samples]
    #non-ER expression
    filterlocs1 <- which(apply(exp1_inter,1,function(x){return((sum(x==0)/length(x))>=0.3)}))
    if(length(filterlocs1)==0){
      exp1_filter <- exp1_inter
    }else{
      exp1_filter <- exp1_inter[-(filterlocs1),]
    }
    expre1 <- log2(exp1_filter+0.001)
	  #ER expression
    filterlocs2 <- which(apply(exp2_inter,1,function(x){return((sum(x==0)/length(x))>=0.3)}))
    if(length(filterlocs2)==0){
      exp2_filter <- exp2_inter
    }else{
      exp2_filter <- exp2_inter[-(filterlocs2),]
    }
    expre2 <- log2(exp2_filter+0.001)
  }
  n <- length(inter_samples)
  gn <- 1
  pccCalcu <- function(x,y,z){
    r12=cor(t(x),t(y))
    r13=cor(t(x),z)
    r23=cor(z,t(y))
    r123=r13%*%r23
    r1=r12-r123
    r2a=sqrt(1-r13*r13)
    r2b=sqrt(1-r23*r23)
    r2=r2a%*%r2b
    r=r1/r2
    return(r)
  }
  pcor <- pccCalcu(expre2,expre1,tumpur)
  statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  p.value <- 2*pnorm(-abs(statistic))
  rownames(pcor) <- rownames(expre2)  
  colnames(pcor) <- rownames(expre1) 
  rownames(p.value) <- rownames(expre2)
  colnames(p.value) <- rownames(expre1)
  ParCorRes <- list(pcor,p.value)
  names(ParCorRes) <- c("pcor.value","p.value")
  return(ParCorRes)
}
