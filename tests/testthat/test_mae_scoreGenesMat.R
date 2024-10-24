# Score genes based on example data
library(MultiAssayExperiment)
# Methylation data 
methylData = matrix(runif(50), nrow=10)
colnames(methylData) <- paste0("samp", 1:ncol(methylData))
rownames(methylData) <- paste0("gene", 1:nrow(methylData))
# Expression data
expresData = matrix(rnorm(50), nrow=10)
colnames(expresData) <- paste0("samp", 1:ncol(methylData))
rownames(expresData) <- paste0("gene", 1:nrow(methylData))
# ColData
colDat <- data.frame(sampleID=colnames(methylData), 
                     name=letters[1:ncol(methylData)])
                     
rownames(colDat)<- colDat$sampleID
mae <- MultiAssayExperiment::MultiAssayExperiment(
experiments = list(methylation = methylData,
  expression = expresData),
  colData =colDat
)
sampleSize <- ncol(mae@ExperimentList[[1]])
reqPercentages <- matrix (c(3, 20, 5, 5, 40, 20, 4, 1, 2), nrow=3, byrow=TRUE)
(theWeightMifL=matrix (c(2,-2,-sampleSize/5,1,0,-2,1,1,2), nrow=3, byrow=TRUE))
(theWeightMifNonL=matrix (c(0,-2,-sampleSize/5,0,0,-2,0,0,0), nrow=3, byrow=TRUE))

mae_scoreGenesMat (mae, 
             x1=1/3, x2=2/3,
             y1=NULL, y2=NULL, percY1=1/3, percY2=2/3,
             aReqPercentsMat = reqPercentages,
             aWeightMifL= theWeightMifL,
             aWeightMifNonL= theWeightMifNonL)

             
