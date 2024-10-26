#' plotGenesMat
#'
#' \code{plotGenesMat} wrapper function for plotting the scatterplots associated
#' with two matrices.
#'
#' @param mets A matrix containing methylation data.
#' @param expres A matrix containing expression data.
#' @param fileName The name of the file used to save the results as a PDF.
#' If NULL, the plot is displayed on the screen.
#' @param text4Title An optional title for the plot, incorporating the gene name
#' and L-shape score. Defaults to NULL.
#' @param x1 The x-coordinate of vertical points on the X-axis. 
#' Expected to contain methylation values ranging between 0 and 1,
#' with default values set to 1/3 and 2/3.
#' @param y1 The y-coordinate of vertical points on the Y-axis.
#' If NULL, these are set to the percentiles of `yVec` 
#' defined by `percY1` and `percY2`.
#' @param percY1 Values used as defaults for `y1` when it is set to NULL.
#' @param x2 The x-coordinate of vertical points on the X-axis.
#' Expected to contain methylation values ranging between 0 and 1,
#' with default values set to 1/3 and 2/3.
#' @param y2 The y-coordinate of vertical points on the Y-axis.
#' If NULL, these are set to the percentiles of `yVec` 
#' defined by `percY1` and `percY2`.
#' @param percY2 Values used as defaults for `y2` when it is set to NULL.
#' @param plotGrid A logical value; defaults to TRUE, indicating 
#' whether to plot gridlines on the graph.
#' @param logicSc A numeric score representing the L-shape score. 
#' Defaults to NULL.
#'
#' @return a pdf with scatterplots for all genes
#'
#' @keywords scatterplot gene plot matrix
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export plotGenesMat
#' @examples
#' mets <- matrix(runif(1000), nrow = 100)
#' expres <- matrix(rnorm(1000), nrow = 100)
#' rownames(mets) <- paste0("Gene", 1:nrow(mets))
#' rownames(expres) <- paste0("Gene", 1:nrow(expres))
#' # plotGenesMat (mets=mets, expres=expres, fileName = 'PlotAllGenes.pdf')
#'
plotGenesMat <- function(mets, expres, fileName = NULL, 
    text4Title = NULL, x1 = 1/3,
    x2 = 2/3, y1 = NULL, y2 = NULL, percY1 = 1/3, percY2 = 2/3, 
    plotGrid = TRUE,
    logicSc = NULL) {
    if (!is.null(fileName)) {
        grDevices::pdf(fileName)
    }
    if (!is.null(text4Title)) {
        text4Title <- paste(rownames(expres), text4Title, sep = ",")
    } else {
        if (is.null(logicSc)) {
            text4Title <- rownames(expres)
        } else {
            text4Title <- paste(rownames(expres), "\n L-shaped = ", 
                logicSc, sep = " ")
        }
    }
    # opt<-par(mfrow=c(2,2))
    for (gene in seq_len(nrow(expres))) {
        xVec <- as.numeric(mets[gene, ])
        yVec <- as.numeric(expres[gene, ])
        plotGeneSel(
            xMet = xVec, yExp = yVec, titleText = text4Title[gene],
            x1 = x1,
            x2 = x2, percY1 = percY1, percY2 = percY2, plotGrid = plotGrid
        )
    }
    # par(opt)
    if (!is.null(fileName)) {
        grDevices::dev.off()
    }
}
