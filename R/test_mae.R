library(MultiAssayExperiment)

# Methylation data

methylData = matrix(runif(50), nrow=10)
colnames(methylData) <- paste0("samp", 1:ncol(methylData))
rownames(methylData) <- paste0("gene", 1:nrow(methylData))

colDat <- data.frame(sampleID=colnames(methylData), 
                     name=letters[1:ncol(methylData)])
rownames(colDat)<- colDat$sampleID

# Expression data

expresData = matrix(rnorm(50), nrow=10)
colnames(expresData) <- paste0("samp", 1:ncol(methylData))
rownames(expresData) <- paste0("gene", 1:nrow(methylData))
                    
mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(methylation = methylData,
                     expression = expresData),
  colData =colDat
  )

mae@colData
mae@ExperimentList
mae@sampleMap
prod(names(mae@ExperimentList) == c("methylation" ,"expression" ))==1

