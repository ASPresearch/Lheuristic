% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcFreqs.R
\name{calcFreqs}
\alias{calcFreqs}
\title{calcFreqs}
\usage{
calcFreqs(
  mae,
  geneNum,
  x1,
  x2,
  y1 = NULL,
  y2 = NULL,
  percY1 = 1/3,
  percY2 = 2/3
)
}
\arguments{
\item{mae}{MultiAssayExperiment object containing methylation and expression
matrices.}

\item{geneNum}{row of expression/methylation matrix for which the frequencies
will be computed.}

\item{x1, x2}{Coordinates of vertical points in the X axis. Because it is
expected to contain.
methylation values that vary between 0 and 1, default values are 1/3 and 2/3.}

\item{y1, y2}{Coordinates of vertical points on the Y-axis. If set to NULL,
they will be automatically assigned as the percentiles of `y` defined 
by `percY1` and `percY2`.}

\item{percY1, percY2}{Default values for `y1` and `y2` when set to `NULL`}
}
\value{
a matrix with calculated frequencies
}
\description{
\code{calcFreqs} Given a MultiAssayExperiment with methylation and
expression data, the function overlays a grid on a YMet ~ Xmet scatterplot
and returns a 3x3 matrix with point counts per grid cell based on  the
vertical and horizontal lines.
}
\examples{

# Methylation data
methylData <- matrix(runif(50), nrow = 10)
colnames(methylData) <- paste0(
    "samp",
    1:ncol(methylData)
)
rownames(methylData) <- paste0(
    "gene",
    1:nrow(methylData)
)
# Expression data
expresData <- matrix(rnorm(50), nrow = 10)
colnames(expresData) <- paste0(
    "samp",
    1:ncol(methylData)
)
rownames(expresData) <- paste0(
    "gene",
    1:nrow(methylData)
)
# ColData
colDat <- data.frame(
    sampleID = colnames(methylData),
    name = letters[1:ncol(methylData)]
)

rownames(colDat) <- colDat$sampleID
mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(
        methylation = methylData,
        expression = expresData
    ),
    colData = colDat
)
geneRow <- 1
x1 <- 1 / 3
x2 <- 2 / 3
y1 <- NULL
y2 <- NULL
percY1 <- 1 / 3
percY2 <- 2 / 3

calcFreqs(
    mae = mae, geneNum = geneRow,
    x1, x2, y1, y2, percY1, percY2
)
}
\keyword{calculation}
\keyword{frequencies}
