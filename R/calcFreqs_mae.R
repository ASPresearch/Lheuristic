#' calcFreqs_mae
#'
#' \code{calcFreqs_mae} Given a MultiAssayExperiment object with methylation and expression matrices, 
#' this function overimposes a grid on the scatterplot defined by YMet~Xmet layers
#' and returns a 3X3 matrix of counts with the number of points in each cell of the grid
#' for given vertical and horizontal lines.
#'
#' @param mae MultiAssayExperiment object containing methylation and expression matrices.
#' @param x1,x2 Coordinates of vertical points in the X axis. Because it is expected to contain methylation values that vary between 0 and 1, the default values are 1/3 and 2/3.
#' @param y1,y2 Coordinates of vertical points in the Y axis. Leaving them as NULL assigns them the percentiles of yVec defined by `percY1` and `percY2`.
#' @param percY1,percY2 Values used to act as default for `y1` and `y2` when these are set to `NULL`
#'
#' @return a matrix with calculated frequencies
#'
#' @keywords calculation frequencies
#' @export calcFreqs_mae
#'
#' @examples
#' 
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'   experiments = list(methylation = matrix(runif(1000), nrow=100), 
#'   expression = matrix(rnorm(1000), nrow=100))
#' )
#' x1 <- 1/3
#' x2 <- 2/3
#' y1 <- NULL
#' y2 <- NULL
#' percY1 <- 1/3
#' percY2 <- 2/3
#' calcFreqs_mae(mae, x1, x2, y1, y2, percY1, percY2)
#'
calcFreqs_mae <- function (mae, x1, x2, y1=NULL, y2=NULL,
                       percY1=1/3, percY2=2/3)
{
  xMet <- as.numeric(assay(mae, "methylation"))
  yExp <- as.numeric(assay(mae, "expression"))
  freqsMat <- matrix(0, nrow=3, ncol=3)
  xVals <- c(x1, x2)
  minExp <- min(yExp)
  maxExp <- max(yExp)
  delta <- maxExp - minExp
  if (is.null(y1)) y1 <- minExp + percY1 * delta
  if (is.null(y2)) y2 <- minExp + percY2 * delta
  yVals <- c(y1, y2)
  condX <- c("(xMet<=x1)", "((xMet>x1) & (xMet<=x2))", "(xMet>x2)")
  condY <- c("(yExp>y2)", "((yExp<=y2) & (yExp>y1))", "(yExp<=y1)")
  for (i in seq_len(3)) {
    for (j in seq_len(3)) {
      condij <- paste(condX[j], condY[i], sep="&")
      freqsMat[i, j] <- sum(eval(parse(text=condij)))
    }
  }
  return(freqsMat)
}
