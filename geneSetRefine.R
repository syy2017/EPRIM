#' Calculation for the correlation between expression of each gene and the mean expression
#' level of the gene set
#'
#' @param  geneset.exp Expression profile with row as genes, column as samples.
#' @param  cor.method The method used for correlation calculation, either "spearman", "pearson"
#' or "kendall", "spearman" by default.
#' @param  cor.cutoff The cutoff of correlation coefficient used for immune signature
#' refinement, cor.cutoff = NULL by default.
#'
#' @return A data frame indicating the detailed correlation results.

genesetMeanCor <- function(geneset.exp, cor.method = "spearman", cor.cutoff = NULL) {
  # Calculate the mean expression level of the gene set
  geneset.mean <- colMeans(geneset.exp, na.rm = TRUE)
  
  # Calculate the correlation between expression of each gene and the mean
  # expression level of the gene set
  cor.withMean <- apply(as.matrix(geneset.exp), 1, function(x) {
    cor.res <- cor.test(x, geneset.mean, method = cor.method) 
    return(data.frame(cor = signif(cor.res$estimate, 3), p.val = signif(cor.res$p.value, 3)))
  })
  
  cor.withMean <- do.call(rbind.data.frame, cor.withMean)
  cor.withMean$gene <- rownames(geneset.exp)
  # significant level
  cor.withMean$signif <- "non-sig"
  cor.withMean$signif[which(cor.withMean$p.val < 0.05)] <- "sig"
  # correlation direction
  cor.withMean$direction <- "positive"
  cor.withMean$direction[which(cor.withMean$cor < 0)] <- "negative"
  
  # return filtered result
  if (!is.null(cor.cutoff)) cor.withMean <- cor.withMean[which(abs(cor.withMean$cor) >= cor.cutoff), ]
  
  return(cor.withMean)
}







#' Refinement of specific immune gene sets
#' @param  exp.profile A numeric matrix containing the expression of mRNA with
#' rownames and colnames,  row as genes, column as samples.
#' @param  signature.df A data frame containing immune signature, the first column (named 'Gene')
#' is the gene ID, the second column (named 'setName') is the immune pathway name.
#' @param  cor.cutoff The cutoff of correlation coefficient used for immune signature
#' refinement.
#'
#' @return A list containing the detailed information of refined immune signature.

genesetMeanFilter <- function(exp.profile, signature.df, cor.method = "spearman", cor.cutoff = 0.5) {
  
  ## Remove not matched genes  
 
  idx <- match(signature.df$Gene, rownames(exp.profile))
  rm.genes <- signature.df$Gene[which(is.na(idx))]
  if (length(rm.genes) != 0) signature.df <- signature.df[-which(is.na(idx)), ]
  if (dim(signature.df)[1] == 0) stop("All signature genes are not in the expression profile")
  
  
  ## Compute the correlation

  cor.res <- tapply(signature.df$Gene, signature.df$setName, function(geneset) {
    signature.exp <- exp.profile[geneset, ]
    # Only one signature gene in a gene set, make expression as a matrix
    if(is.null(dim(signature.exp))){
      signature.expMat <- matrix(as.numeric(signature.exp),nrow=1)
      colnames(signature.expMat) <- names(signature.exp)
      rownames(signature.expMat) <- geneset
      signature.exp <- signature.expMat
    }
    geneset.MeanCor <- genesetMeanCor(geneset.exp = signature.exp, cor.method = cor.method, cor.cutoff = cor.cutoff)
    return(geneset.MeanCor)
  })
  
  cor.res.df <- do.call(rbind.data.frame,cor.res)
  # Add gene set names
  cor.res0 <- unlist(sapply(cor.res,function(k){temp <- dim(k)[1];return(temp) }))
  setAnno <- rep(names(cor.res),cor.res0)
  cor.res.df$setAnno <- setAnno
  # Correlation results of all genes
  cor.res.dfAll <- cor.res.df
  # Genes with significant correlations
  cor.res.df <- cor.res.df[which(cor.res.df$signif == "sig"), ]
  
  ## Output
 
  # Retain significant correlated genes
  signatureGene.list <- split(cor.res.df$gene, cor.res.df$setAnno)
  
  # Add label of the correlation direction
  cor.res.df$label <- paste(cor.res.df$setAnno, cor.res.df$direction, sep = "_")
  
  signatureGene.list <- list(df = cor.res.df, list = signatureGene.list, dfAll = cor.res.dfAll)
  return(signatureGene.list)
}


