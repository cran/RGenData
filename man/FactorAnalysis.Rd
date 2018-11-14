\name{FactorAnalysis}
\alias{FactorAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FactorAnalysis
}
\description{Analyzes comparison data with known factorial structures
}
\usage{
FactorAnalysis(data, corr.matrix = FALSE, max.iteration = 50,n.factors = 0,
corr.type = "pearson")
}
\arguments{
  \item{data}{Matrix to store the simulated data (matrix).}
  \item{corr.matrix}{Correlation matrix (default is FALSE)}
  \item{max.iteration}{Maximum number of iterations (scalar, default is 50).}
  \item{n.factors}{Number of factors (scalar, default is 0).}
  \item{corr.type}{Type of correlation (character, default is "pearson", user can also call "spearman").}
}
\value{

 \item{$loadings}{Factor loadings (vector, if one factor. matrix, if multiple factors)}
 \item{$factors}{Number of factors (scalar).}
}
\references{
Ruscio & Roche (2011)
}
\author{
John Ruscio
}
\examples{
# create data matrix x with n = 200 cases, k = 9 variables
# 3 variables load onto each of 3 orthogonal factors
# all marginal distributions are highly skewed
x <- matrix(nrow = 200, ncol = 9)
for (i in 1:3) {
  shared <- rchisq(200, 1)
  for (j in 1:3) {
    x[, (i - 1) * 3 + j] <- shared + rchisq(200, 1)
  }
}
# perform factor analysis of data matrix x
FactorAnalysis(x)
}
