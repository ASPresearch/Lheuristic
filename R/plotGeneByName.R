#' plotGeneByName
#'
#' \code{plotGeneByName} plots points on a scatterplot with a 3x3 grid
#' superimposed. The name of a the gene is provided jointly with the matrix
#' and used to select the row to be plotted.
#'
#' @param geneName The name of the gene to be plotted.
#' @param mets A matrix containing the methylation data for the specified gene.
#' @param expmat A matrix containing the expression data for the specified gene.
#' @param filename If provided, the name of the file to save the results
#' as a PDF;  defaults to NULL.
#' @param text4Title A string used as the main title for the plot. Defaults to 
#' `geneName` if set to NULL.
#' @param plotGrid A boolean parameter indicating whether to pass
#' the grid option to the `plotGeneSel` function.
#' @param figs A two-component vector defining the 2-dimensional structure of
#' the plots to be generated.
#'
#' @return a pdf with scatterplots of selected by gene name
#'
#' @keywords plot gene name
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export plotGeneByName
#'
#' @examples
#' # Plot gene by name based on example data
#' mets <- matrix(runif(1000), nrow = 100)
#' expres <- matrix(rnorm(1000), nrow = 100)
#' rownames(mets) <- paste0("Gene", 1:nrow(mets))
#' rownames(expres) <- paste0("Gene", 1:nrow(expres))
#' plotGeneByName(gene = "Gene40", mets = mets, expmat = expres)
#'
plotGeneByName <- function(geneName, mets, expmat, filename = NULL,
    text4Title = NULL,
    plotGrid = TRUE, figs = c(2, 2)) {
    if (!is.null(filename)) {
        grDevices::pdf(filename)
    }
    if (!is.null(text4Title)) {
        text4Title <- paste(geneName, text4Title, sep = ", ")
    } else {
        text4Title <- geneName
    }
    if (geneName %in% rownames(expmat)) {
        genePos <- which(rownames(expmat) == geneName)
    } else {
        genePos <- NULL
    }
    if (!(is.null(genePos))) {
        xVec <- as.numeric(mets[genePos, ])
        yVec <- as.numeric(expmat[genePos, ])
        plotGeneSel(
            xMet = xVec, yExp = yVec, titleText = text4Title, 
            x1 = 1/3, x2 = 2/3,
            plotGrid = plotGrid
        )
    }
    if (!is.null(filename)) {
        grDevices::dev.off()
    }
}
