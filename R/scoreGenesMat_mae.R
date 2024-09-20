#' scoreGenesMat_mae
#'
#' \code{scoreGenesMat} scores scatterplots using a binary and a numeric schemes on a row-wise basis.
#'
#' @param mae MultiAssayExperiment object containing methylation and expression matrices.
#' @param aReqPercentsMat Matrix of minimum maximum percentage of counts to have in a given cell
#' @param aWeightMifL A matrix of weights to score the previous counts if the scatterplot has been classified as L.
#' @param aWeightMifNonL A matrix of weights to score the previous counts if the scatterplot has been classified as non-L
#' @param x1,x2 Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1 the default values are 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1,percY2 Values used to act as default for `y1`and `y2` when these are set to `NULL`
#' @export scoreGenesMat_mae
#' @return A data frame with two columns: "logicSc" (a logical score indicating whether the gene is considered "active") and "numericSc" (a numerical score).
#'
#' @examples 
#'
#' # Score genes based on example data
#' # Methylation data 
#' methylData = matrix(runif(50), nrow=10)
#' colnames(methylData) <- paste0("samp", 1:ncol(methylData))
#' rownames(methylData) <- paste0("gene", 1:nrow(methylData))
#' # Expression data
#' expresData = matrix(rnorm(50), nrow=10)
#' colnames(expresData) <- paste0("samp", 1:ncol(methylData))
#' rownames(expresData) <- paste0("gene", 1:nrow(methylData))
#' # ColData
#' colDat <- data.frame(sampleID=colnames(methylData), 
#'                      name=letters[1:ncol(methylData)])
#' rownames(colDat)<- colDat$sampleID
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#' experiments = list(methylation = methylData,
#'   expression = expresData),
#'   colData =colDat
#' )
#' sampleSize <- dim(mets)[2]
#' numGenes <-   dim(mets)[1]
#' reqPercentages <- matrix (c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow=3, byrow=TRUE)
#' (theWeightMifL=matrix (c(2,-2,-sampleSize/5,1,0,-2,1,1,2), nrow=3, byrow=TRUE))
#' (theWeightMifNonL=matrix (c(0,-2,-sampleSize/5,0,0,-2,0,0,0), nrow=3, byrow=TRUE))
#' scoreGenesMat_mae (mae,
#'              x1=1/3, x2=2/3,
#'              y1=NULL, y2=NULL, percY1=1/3, percY2=2/3,
#'              aReqPercentsMat = reqPercentages,
#'              aWeightMifL= theWeightMifL,
#'              aWeightMifNonL= theWeightMifNonL)
#'

scoreGenesMat_mae <- function(mae, 
                          x1=1/3, x2=2/3,
                          y1=NULL, y2=NULL, 
                          percY1=1/3, percY2=2/3,
                          aReqPercentsMat, 
                          aWeightMifL=0.5, 
                          aWeightMifNonL=0.25)
{
  if(sum(aReqPercentsMat)!=100) 
    stop("Error: Percentages must add up to 100")
  if (prod(names(mae@ExperimentList) == c("methylation" ,"expression" ))!=1)
    stop("Error: Names of layers must be 'methylation' and 'expression'")
  mets <- MultiAssayExperiment::assay(mae, "methylation")
  expres <- MultiAssayExperiment::assay(mae, "expression")
  N <- dim(mets)[2]
  Ngenes <-nrow(mets)
  scores <- data.frame(logicSc=rep(FALSE, Ngenes), numericSc=rep(0,Ngenes))
  rownames(scores)<- rownames(mets)
  minmaxCounts <- toReqMat(N, aReqPercentMat=aReqPercentsMat)
  indexes <- seq(from=1, to=Ngenes, by=1)
  for (gene in indexes){
    theGene <- rownames(expres)[gene]
    xVec<- mets[theGene,]
    yVec<- expres[theGene,]
    geneGrid <- calcFreqs(xMet=xVec, yExp=yVec, x1=x1, x2=x2,
                          y1=y1, y2=y2, percY1=percY1, percY2=percY2)
    binSc <- binScore(geneGrid, minmaxCounts)
    scores[gene, "logicSc"] <- binSc
    numSc <- numScore (geneGrid, LShaped=binSc,
                       aWeightMifL=aWeightMifL,
                       aWeightMifNonL=aWeightMifNonL)
    scores[gene, "numericSc"] <- numSc
  }
  return (scores)
}
