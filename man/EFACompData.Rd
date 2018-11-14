\name{EFACompData}
\alias{EFACompData}
\title{EFACompData
}
\description{
Comparison data
}
\usage{
EFACompData(data, f.max, n.pop = 10000, n.samples = 500, alpha = .30, graph = FALSE,
corr.type = "pearson")
}
\arguments{
  \item{data}{Matrix to store the simulated data (matrix).}
  \item{f.max}{Largest number of factors to consider (scalar).}
  \item{n.pop}{Size of finite populations of comparison data (scalar, default is 10,000 cases).}
  \item{n.samples}{Number of samples drawn from each population (scalar, default is 500).}
  \item{alpha}{Alpha level when testing statistical significance of improvement with additional factor (scalar, default is .30)}
  \item{graph}{Whether to plot the fit of eigenvalues to those for comparison data (default is FALSE).}
  \item{corr.type}{Type of correlation (character, default is "pearson", user can also call "spearman").}
}
\value{
Nothing, displays number of factors on screen.
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
# empirically determine number of factors in data matrix x
EFACompData(x, f.max = 5)
}
