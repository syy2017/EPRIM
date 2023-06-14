#' Calculation enrichment results for epigenetic regulators in immune-related pathways.
#'   
#' @param exp.profile A numeric matrix containing the expression of mRNA with
#' rownames and colnames,  row as genes, column as samples.
#' @param interested.ER A character vector of interested epigenetic regulator genes. 
#' @param is.adjusted Whether or not the 'exp.profile' is preprocessed. If the expression values
#' were log-transformed and the rows containing many 0 values were removed, is.adjusted = TRUE.
#' @param tum.put A numeric vector containing the samples' tumor purity with sample names.
#' @param pathways A list containing immune gene sets with pathway names. 
#' @param minSize The least immune signature genes matched in the expression profile. Defaults 
#' to 1. 
#' @param nperm Number of permutation times used for enrichment analysis.   
#'   
#' @return Detailed correlation results for epigenetic regulator-signature pairs.  


ImmuEpiReg <- function(exp.profile,interested.ER,is.adjusted,tum.put,pathways,minSize=1,nperm=1000){
  #ER expression 
  ERexpre <- exp.profile[intersect(rownames(exp.profile),interested.ER),] 
  #non-ER expression
  NERexpre <- exp.profile[setdiff(rownames(exp.profile),interested.ER),] 

  ERexpre <- as.matrix(ERexpre)
  NERexpre <- as.matrix(NERexpre)
  
  #PCC calculation               
  pccRes <- parcorCal(exp1 = NERexpre,exp2 = ERexpre,tum_pur = tum.put,is.adjusted = is.adjusted)
  pcc_pValue <- pccRes$p.value
  pcc_pcorValue <- pccRes$pcor.value
  #Rank list
  Rankscore <- -log10(pcc_pValue)*sign(pcc_pcorValue)
  #Imm-ER pairs
  ImmERCorRes_all <- c()
  for(m in 1:nrow(Rankscore)){
    #remove the rows contain Infs
    if(sum( is.infinite(Rankscore[m,])) != 0){
      next()
    }
    ranks <- Rankscore[m,]
    #GSEA
    ImmERCorRes <- fgsea::fgsea(pathways, ranks, minSize, maxSize=5000, nperm)
    #Immune-ER correlations
    ImmERcorscore <- c()
    for(n in 1:nrow(ImmERCorRes)){
      if(ImmERCorRes$ES[n]>0){
        sig_mn <- 1 - 2*ImmERCorRes$pval[n]
      }else{
        sig_mn <- 2*ImmERCorRes$pval[n] - 1
      }
      ImmERcorscore <- c(ImmERcorscore,sig_mn)
      
    }
 
    ER <- rownames(Rankscore)[m]
    ImmERCorRes_m <- cbind(ER,ImmERCorRes,ImmERcorscore)
    ImmERCorRes_all <- rbind(ImmERCorRes_all,ImmERCorRes_m)
  }

  return(ImmERCorRes_all)
}