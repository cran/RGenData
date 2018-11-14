\name{GenDataPopulation}
\alias{GenDataPopulation}
\title{{GenDataPopulation}
}
\description{Simulates multivariate nonnormal data using an iterative algorithm
}
\usage{
GenDataPopulation(supplied.data, n.factors, n.cases, max.trials = 5,
                  initial.multiplier = 1, corr.type = "pearson", seed = 0)
}
\arguments{
  \item{supplied.data}{Data supplied by user.}
  \item{n.factors}{Number of factors (scalar).}
  \item{n.cases}{Number of cases (scalar).}
  \item{max.trials}{Maximum number of trials (scalar, default is 5).}
  \item{initial.multiplier}{Value of initial multiplier (scalar, default is 1).}
  \item{corr.type}{Type of correlation (character, default is "pearson", user can also call "spearman").}
   \item{seed}{seed value (scalar, default is 0).}
}
\value{
{data}{Population of data}
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
# generate (finite) population of data reproducing distributions and correlations in x
GenDataPopulation(x, n.factors = 3, n.cases = 10000)
}
