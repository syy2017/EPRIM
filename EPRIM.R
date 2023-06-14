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








ImmERPancan <-  function(workDir="/home/data/t150343/IJob",canType){
  
  
  outDir <- paste0(workDir,"/Result/",canType)
  if(!dir.exists(outDir)){
    dir.create(outDir)
  }
  
  
  outPath <- paste0(outDir,"/ImmLnc")
  if(!dir.exists(outPath)){
    dir.create(outPath)
  }
  
  ###########---提取肿瘤样本mRNA表达谱
  load(paste0(workDir,"/Data/TCGA/OSF_",canType,"_All_tpm.RData"))
  #提取表达谱中的肿瘤样本数据(表达谱谱中既有正常样本又有肿瘤样本)
  sampletype <- substring(colnames(PCG_tpm), first=14, last=15) 
  names(sampletype) <- substring(colnames(PCG_tpm), first=1, last=12)#向量名字为样本名 
  nn <- is.element(sampletype,paste0("0",1:9))#有37个正常样本
  PCG_tpm1 <- PCG_tpm[,which(nn)] 
  sampletype1 <- sampletype[which(nn)]#416肿瘤样本
  #将表达谱中TCGA样本的命名格式和浸润谱中保持一致
  colnames(PCG_tpm1) <- substring(colnames(PCG_tpm1), first=1, last=12) #57110,416
  #去除一个样本多个值的情况(肿瘤样本不同样本类型（比如原发，转移等）)
  multi_sample <- duplicated(colnames(PCG_tpm1)) | duplicated(colnames(PCG_tpm1), fromLast = TRUE) #有0个样本重复
  PCG_tpm1 <- PCG_tpm1[, !multi_sample]#416
  pcgInfo <- get(load(file=paste0(workDir,"/Data/gencode_pcg.RData")))#pcgInfo,20296
  #length(intersect(rownames(pcgInfo),rownames(PCG_tpm1)))#20074
  PCG_tpm2 <- PCG_tpm1[intersect(rownames(pcgInfo),rownames(PCG_tpm1)),]#20
  save(PCG_tpm2,file=paste0(outPath,"/PCG_tpm2.RData"))
  PCG_tpm2Name <- PCG_tpm2
  rownames(PCG_tpm2Name) <- pcgInfo[intersect(rownames(pcgInfo),rownames(PCG_tpm1)),"geneName"]
  write.table(PCG_tpm2Name,sep="\t",row.names=T,col.names=T,quote=F,file=paste0(outPath,"/PCG_tpm2Name.txt"))
  
  
  
  ###########---免疫通路基因Refine
  ###输入参数
  cor.cutoff = NULL
  cor.method = "spearman"
  #row as genes, column as samples
  TPM.profile <- get(load(file=paste0(outPath,"/PCG_tpm2.RData"))) #20074,416 
  #A data frame containing immunophenotype signature: the first column
  #' (named 'Gene') is the gene ID, the second column (named 'setAnno') is the immunophenotype name.
  pathways <- get(load(file=paste0(workDir,"/Data/ImmLnc17Pathways.RData")))
  signature.df0 <- sapply(1:length(pathways),function(k){
    temp <- cbind(pathways[[k]],names(pathways)[k])
    return(temp)
  })               
  signature.df <- do.call(rbind.data.frame,signature.df0)
  colnames(signature.df) <- c("Gene","setAnno")
  
  
  ###运行函数
  source(paste0(workDir,"/Code/Rscript/main/geneSetRefine.R"))
  genesetMeancor.res <- genesetMeanFilter(exp.profile = TPM.profile, cor.method = cor.method, signature.df = signature.df, cor.cutoff = cor.cutoff)
  save(genesetMeancor.res,file=paste0(outPath,"/genesetMeancor.res.RData")) #genesetMeancor.res$df、genesetMeancor.res$list、genesetMeancor.res$dfAll
  
  
  
  
  ###########---计算肿瘤纯度
  source(paste0(workDir,"/Code/Rscript/main/turpur.est.R"))
  dir_test <- outPath
  file_test <-"PCG_tpm2Name.txt"
  res_test <- turpur.est(dir_test,file_test)
  save(res_test,file=paste0(outPath,"/tumoPurity.RData"))
  
  
  ###########---计算ER-pathway关联
  load(file=paste0(outPath,"/PCG_tpm2.RData")) #PCG_tpm2
  load(file=paste0(workDir,"/Data/ERgenes.RData")) #ERgenes
  ###提取ER表达谱
  ERexpre <- PCG_tpm2[intersect(rownames(PCG_tpm2),ERgenes[,"ENSEMBL.Gene.ID"]),] 
  ###提取非ER表达谱
  NERexpre <- PCG_tpm2[setdiff(rownames(PCG_tpm2),ERgenes[,"ENSEMBL.Gene.ID"]),] 
  ###肿瘤纯度
  load(file=paste0(outPath,"/tumoPurity.RData")) #res_test
  ###Refine后通路
  genesetMeancor.res <- get(load(file=paste0(outPath,"/genesetMeancor.res.RData"))) #refine后的pathway
  pathways <- genesetMeancor.res$list #refine后的pathway   
  ###ER-pathway关联
  turpur_ori <- res_test 
  names(turpur_ori) <- gsub("[.]","-",names(res_test))##肿瘤样本名字需要更改(TCGA.CG.4462>>TCGA-CG-4462)
  lncRNA_exp1 <- as.matrix(ERexpre) 
  mRNA_exp1 <- as.matrix(NERexpre)
  k=0.995
  
  source(paste0(workDir,"/Code/Rscript/main/par.cor.R"))
  source(paste0(workDir,"/Code/Rscript/main/immu.LncRNA.R"))
  test_resAll <- immu.LncRNA(mRNA_exp1,lncRNA_exp1,adjusted=F,turpur_ori,pathways,k)
  save(test_resAll,file=paste0(outPath,"/test_resAllRefine.RData"))
  
  return (NULL)
}




