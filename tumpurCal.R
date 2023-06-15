#' Calculation of samples' tumor purity
#' 
#' @param input.dir Character string indicating the storage directory of mRNA expression file,
#' gene symbol as rownames. 
#' @param file.name Character string indicating the name of mRNA expression file.
#' @param platform Character string indicating platform type, either "affymetrix", "agilent"
#' or "illumina". Defaults to "illumina".  
#'   
#' @return A numeric vector containing samples' tumor purity.

tumpurCal <- function(input.dir,file.name,platform="illumina"){
  library(estimate)
  dirOri <- getwd()
  setwd(input.dir)
  estimate::filterCommonGenes(input.f=file.name, output.f="expreFilter.gct", id="GeneSymbol")
  estimate::estimateScore("expreFilter.gct", "estiScore.gct", platform=platform)
  estScore <- readLines("estiScore.gct")
  estScoresplit <- unlist(strsplit(estScore[grep("ESTIMATEScore",estScore)],"\t"))
  estScores <- as.numeric(estScoresplit[3:length(estScoresplit)])
  tumPurity = cos(0.6049872018+0.0001467884*estScores)
  sampSplit <- unlist(strsplit(estScore[grep("NAME",estScore)],"\t"))
  sampName <- sampSplit[3:length(sampSplit)]
  names(tumPurity) <- sampName
  setwd(dirOri)
  return(tumPurity)
  
}


