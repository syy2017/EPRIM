#' Identify epigenetic regulators associated with immune gene sets
#'
#' @param exp.profile A numeric matrix containing the expression of mRNA with
#' rownames and colnames,  row as genes, column as samples.
#' @param exp.profile.file A character string representing the pathname of mRNA expression 
#' file, gene symbol as rownames.
#' @param interested.ER A character vector of interested epigenetic regulator genes. 
#' @param signature.list A list containing immune gene sets. 
#' @param cor.cutoff The cutoff of correlation coefficient used for immune signature
#' refinement, cor.cutoff = NULL by default. 
#' @param cor.method The method used for correlation calculation, either "spearman", "pearson"
#' or "kendall", "spearman" by default. 
#' @param platform Character string indicating platform type, either "affymetrix", "agilent"
#' or "illumina". Defaults to "illumina".  
#' @param is.adjusted Whether or not the 'exp.profile' is preprocessed. If the expression values
#'   were log-transformed and the rows containing many 0 values were removed, is.adjusted = TRUE. 
#' @param min.sz The least immune signature genes matched in the expression profile. Defaults to 
#' 1. 
#' @param perm.times Number of permutation times used for enrichment analysis. 
#'
#' @return A list including detailed correlation results for epigenetic regulator-signature pairs, 
#' and refined immune gene sets.

#' @examples
#' # test for "EPRIM" with example data
#' exp.profile <- load(file=paste0(workDir,"/data/PCG_tpm2.RData"))
#' exp.profile.file <- paste0(workDir,"/data/PCG_tpm2Name.txt")
#' interested.ER <- load(file=paste0(workDir,"/data/ERgenes.RData")) 
#' signature.list <- load(file=paste0(workDir,"/data/signature.list.RData")) 
#' test_res <- EPRIM(exp.profile,exp.profile.file,interested.ER, signature.list)
#' test_res$ImmERres[1:3,] # showing the GSEA results
#' 
#' 

 
EPRIM <-  function(exp.profile, exp.profile.file, interested.ER, signature.list, cor.cutoff = NULL, cor.method = "spearman",
                   platform ="illumina", is.adjusted = FALSE, min.sz = 1, perm.times = 100){
  

  # Refine immune pathways
  signature.df0 <- sapply(1:length(signature.list),function(k){
    temp <- cbind(signature.list[[k]],names(signature.list)[k])
    return(temp)
  })               
  signature.df <- do.call(rbind.data.frame,signature.df0)
  colnames(signature.df) <- c("Gene","setName")
  
  genesetMeancor.res <- genesetMeanFilter(exp.profile = exp.profile, cor.method = cor.method, signature.df = signature.df, cor.cutoff = cor.cutoff)



  # Calculate tumor purity
  filename.split <- unlist(strsplit(exp.profile.file,"/"))
  input.dir <- paste(filename.split[1:(length(filename.split)-1)],collapse = "/")
  file.name <- filename.split[length(filename.split)]
  
  tum.put <- tumpurCal(input.dir=input.dir,file.name=file.name,platform=platform)
  
                         

  #Identify epigenetic regulators associated with immune pathways
  pathways <- genesetMeancor.res$list
  ImmERres <- ImmuEpiReg(exp.profile=exp.profile,interested.ER = interested.ER, is.adjusted=is.adjusted,tum.put=tum.put,pathways=pathways,
                            minSize=min.sz, nperm=perm.times)


  # return results
  return (list(ImmERres=ImmERres,genesetMeancor.res=genesetMeancor.res))
}













