GenDataPopulation <- function(supplied.data, n.factors, n.cases, max.trials = 5,
                              initial.multiplier = 1, corr.type = "pearson",
                              seed = 0 ) {
# Simulates multivariate nonnormal data using an iterative algorithm
#
# Args:
#   supplied.data      : Data supplied by user
#   n.factors          : Number of factors (scalar).
#   n.cases            : Number of cases (scalar).
#   max.trials         : Maximum number of trials (scalar).
#   initial.multiplier : Value of initial multiplier (scalar).
#   corr.type          : Type of correlation (character, default is "pearson",
#                        user can also call "spearman").
#   seed               : seed value (scalar, default is 0).
#
# Returns:
#   data : Population of data
#
  set.seed(seed)
  n.variables <- dim(supplied.data)[2]
  data <- matrix(0, nrow = n.cases, ncol = n.variables)
  distributions <- matrix(0, nrow = n.cases, ncol = n.variables)
  iteration <- 0
  best.rmsr <- 1
  trials.without.improvement <- 0
  for (i in 1:n.variables) {
    if (min(supplied.data[, i]) == max(supplied.data[, i]))
      supplied.data[, i] <- supplied.data[, i] + rnorm(dim(supplied.data)[1], 0,
                           .0001)
    v <- FALSE
    while (!v) {
      distributions[, i] <- sort(sample(supplied.data[, i], size = n.cases,
                                replace = TRUE))
      v <- min(distributions[, i]) != max(distributions[, i])
    }
  }
  target.corr <- cor(supplied.data, method = corr.type)
  intermediate.corr <- target.corr
  shared.comp <- matrix(rnorm(n.cases * n.factors, 0, 1), nrow = n.cases,
                 ncol = n.factors)
  unique.comp <- matrix(rnorm(n.cases * n.variables, 0, 1), nrow = n.cases,
                 ncol = n.variables)
  shared.load <- matrix(0, nrow = n.variables, ncol = n.factors)
  unique.load <- matrix(0, nrow = n.variables, ncol = 1)
  while (trials.without.improvement < max.trials) {
    iteration <- iteration + 1
    factor.analysis <- FactorAnalysis(intermediate.corr, corr.matrix = TRUE,
                       max.iteration = 50, n.factors, corr.type)
    if (n.factors == 1) {
      shared.load[, 1] <- factor.analysis$loadings
    } else {
        for (i in 1:n.factors)
          shared.load[, i] <- factor.analysis$loadings[, i]
    }
      shared.load[shared.load > 1] <- 1
      shared.load[shared.load < -1] <- -1
      if (shared.load[1, 1] < 0)
        shared.load <- shared.load * -1
      for (i in 1:n.variables)
        if (sum(shared.load[i, ] * shared.load[i, ]) < 1) {
          unique.load[i, 1] <- (1 - sum(shared.load[i, ] * shared.load[i, ]))
        } else {
          unique.load[i, 1] <- 0
        }
      unique.load <- sqrt(unique.load)
      for (i in 1:n.variables)
        data[, i] <- (shared.comp %*% t(shared.load))[, i] + unique.comp[, i] *
                     unique.load[i, 1]
      for (i in 1:n.variables) {
        data <- data[sort.list(data[, i]), ]
        data[, i] <- distributions[, i]
      }
      reproduced.corr <- cor(data, method = corr.type)
      residual.corr <- target.corr - reproduced.corr
      rmsr <- sqrt(sum(residual.corr[lower.tri(residual.corr)] *
        residual.corr[lower.tri(residual.corr)]) /
          (.5 * (n.variables * n.variables - n.variables)))
      if (rmsr < best.rmsr) {
        best.rmsr <- rmsr
        best.corr <- intermediate.corr
        best.res <- residual.corr
        intermediate.corr <- intermediate.corr + initial.multiplier *
                             residual.corr
        trials.without.improvement <- 0
      } else {
        trials.without.improvement <- trials.without.improvement + 1
        current.multiplier <- initial.multiplier *
                              .5 ^ trials.without.improvement
        intermediate.corr <- best.corr + current.multiplier * best.res
      }
   }
   factor.analysis <- FactorAnalysis(best.corr, corr.matrix = TRUE,
                                     max.iteration = 50, n.factors,
                                     corr.type)
   if (n.factors == 1) {
     shared.load[, 1] <- factor.analysis$loadings
   } else {
     for (i in 1:n.factors)
        shared.load[, i] <- factor.analysis$loadings[, i]
   }
   shared.load[shared.load > 1] <- 1
   shared.load[shared.load < -1] <- -1
   if (shared.load[1, 1] < 0)
     shared.load <- shared.load * -1
   for (i in 1:n.variables)
     if (sum(shared.load[i, ] * shared.load[i, ]) < 1) {
       unique.load[i, 1] <- (1 - sum(shared.load[i, ] * shared.load[i, ]))
     } else {
       unique.load[i, 1] <- 0
     }
   unique.load <- sqrt(unique.load)
   for (i in 1:n.variables)
     data[, i] <- (shared.comp %*% t(shared.load))[, i] + unique.comp[, i] *
                  unique.load[i, 1]
   data <- apply(data, 2, scale) # standardizes each variable in the matrix
   for (i in 1:n.variables) {
     data <- data[sort.list(data[, i]), ]
     data[, i] <- distributions[, i]
   }
   return(data)
}
################################################################################
GenDataSample <- function(supplied.data, n.factors = 0, max.trials = 5,
                          initial.multiplier = 1, corr.type = "pearson",
                          seed = 0) {
# Bootstraps each variable's score distribution from a supplied data set.
#
# Args:
#   supplied.data      : Data supplied by user.
#   n.factors          : Number of factors (scalar, default is 0).
#   max.trials         : Maximum number of trials (scalar, default is 5).
#   initial.multipiler : Value of initial multiplier (scalar, default is 1).
#   corr.type          : Type of correlation (character, default is "pearson",
#                        user can also call "spearman").
#   seed               : seed value (scalar, default is 0).
#
# Returns:
#   data : Sample of data

  set.seed(seed)
  n.cases <- dim(supplied.data)[1]
  n.variables <- dim(supplied.data)[2]
  data <- matrix(0, nrow = n.cases, ncol = n.variables)
  distribution.matrix <- matrix(0, nrow = n.cases, ncol = n.variables)
  iteration.count <- 0
  best.rmsr.corr <- 1
  trial.counter <- 0
  if (seed != 0)
    set.seed(seed) # If user specified a nonzero seed, set it
  for (i in 1:n.variables)
    distribution.matrix[, i] <- sort(sample(supplied.data[, i], n.cases,
    replace = TRUE))
    # Modify this block of the program, as needed, to generate the desired
    # distribution(s).
  target.correlation <- cor(supplied.data)
  intermediate.correlation <- target.correlation
  if (n.factors == 0) {
    eigenvalues.observed <- eigen(intermediate.correlation)$values
    eigenvalues.random <- matrix(0, nrow = 100, ncol = n.variables)
    random.data <- matrix(0, nrow = n.cases, ncol = n.variables)
    for (i in 1:100) {
      for (j in 1:n.variables)
        random.data[, j] <- sample(distribution.matrix[, j], size = n.cases,
        replace = TRUE)
      eigenvalues.random[i, ] <- eigen(cor(random.data))$values
    }
    eigenvalues.random <- apply(eigenvalues.random, 2, mean)
    # calculate mean eigenvalue for each factor
    n.factors <- max(1, sum(eigenvalues.observed > eigenvalues.random))
  }
  shared.components <- matrix(rnorm(n.cases * n.factors, 0, 1), nrow = n.cases,
                                    ncol = n.factors)
  unique.components <- matrix(rnorm(n.cases * n.variables, 0, 1), nrow = n.cases,
                                    ncol = n.variables)
  shared.load <- matrix(0, nrow = n.variables, ncol = n.factors)
  unique.load <- matrix(0, nrow = n.variables, ncol = 1)
  while (trial.counter < max.trials) {
    iteration.count <- iteration.count + 1
    factor.analysis <- FactorAnalysis(intermediate.correlation,
                                      corr.matrix = TRUE, max.iteration = 50,
                                      n.factors, corr.type)
    if (n.factors == 1) {
      shared.load[, 1] <- factor.analysis$loadings
    } else {
      shared.load <- factor.analysis$loadings
    }
    shared.load[shared.load > 1] <- 1
    shared.load[shared.load < -1] <- -1
    if (shared.load[1, 1] < 0)
      shared.load <- shared.load * -1
    shared.load.sq <- shared.load * shared.load
    for (i in 1:n.variables)
      if (sum(shared.load.sq[i, ]) < 1) {
        unique.load[i, 1] <- (1 - sum(shared.load.sq[i, ]))
      } else {
        unique.load[i, 1] <- 0
      }
    unique.load <- sqrt(unique.load)
    for (i in 1:n.variables)
      data[, i] <- (shared.components %*% t(shared.load))[, i] +
                          unique.components[, i] * unique.load[i, 1]
    for (i in 1:n.variables) {
      data <- data[sort.list(data[, i]), ]
      data[, i] <- distribution.matrix[, i]
    }
    reproduced.corr <- cor(data)
    residual.corr <- target.correlation - reproduced.corr
    rmsr <- sqrt(sum(residual.corr[lower.tri(residual.corr)] *
                     residual.corr[lower.tri(residual.corr)]) /
                     (.5 * (n.variables * n.variables - n.variables)))
    if (rmsr < best.rmsr.corr) {
      best.rmsr.corr <- rmsr
      best.corr <- intermediate.correlation
      best.res <- residual.corr
      intermediate.correlation <- intermediate.correlation +
                                  initial.multiplier * residual.corr
      trial.counter <- 0
    } else {
      trial.counter <- trial.counter + 1
      current.multiplier <- initial.multiplier * .5 ^ trial.counter
      intermediate.correlation <- best.corr + current.multiplier * best.res
    }
  }
  factor.analysis <- FactorAnalysis(best.corr, corr.matrix = TRUE,
                                    max.iteration = 50, n.factors,
                                    corr.type)
  if (n.factors == 1) {
    shared.load[, 1] <- factor.analysis$loadings
  } else {
    shared.load <- factor.analysis$loadings
  }
  shared.load[shared.load > 1] <- 1
  shared.load[shared.load < -1] <- -1
  if (shared.load[1, 1] < 0)
    shared.load <- shared.load * -1
  shared.load.sq <- shared.load * shared.load
  for (i in 1:n.variables)
    if (sum(shared.load.sq[i, ]) < 1) {
      unique.load[i, 1] <- (1 - sum(shared.load.sq[i, ]))
    } else {
      unique.load[i, 1] <- 0
  }
  unique.load <- sqrt(unique.load)
  for (i in 1:n.variables)
    data[, i] <- (shared.components %*% t(shared.load))[, i] +
                        unique.components[, i] * unique.load[i, 1]
  data <- apply(data, 2, scale)
  for (i in 1:n.variables) {
    data <- data[sort.list(data[, i]), ]
    data[, i] <- distribution.matrix[, i]
  }
  data <- data[sample(1:n.cases, n.cases, replace = FALSE), ]
  return(data)
}
################################################################################
FactorAnalysis <- function(data, corr.matrix = FALSE, max.iteration = 50,
                            n.factors = 0, corr.type = "pearson") {
# Analyzes comparison data with known factorial structures
#
# Args:
#   data          : Matrix to store the simulated data.
#   corr.matrix   : Correlation matrix (default is FALSE)
#   max.iteration : Maximum number of iterations (scalar, default is 50).
#   n.factors     : Number of factors (scalar, default is 0).
#   corr.type     : Type of correlation (character, default is "pearson",
#                   user can also call "spearman").
#
# Returns:
#   $loadings : Factor loadings (vector, if one factor. matrix, if multiple
#               factors)
#   $factors  : Number of factors (scalar).
#
  data <- as.matrix(data)
  n.variables <- dim(data)[2]
  if (n.factors == 0) {
    n.factors <- n.variables
    determine <- TRUE
  } else {
    determine <- FALSE
  }
  if (!corr.matrix) {
    corr.matrix <- cor(data, method = corr.type)
  } else {
    corr.matrix <- data
  }
  criterion <- .001
  old.h2 <- rep(99, n.variables)
  h2 <- rep(0, n.variables)
  change <- 1
  iteration <- 0
  factor.loadings <- matrix(nrow = n.variables, ncol = n.factors)
  while ((change >= criterion) & (iteration < max.iteration)) {
    iteration <- iteration + 1
    eigenvalue <- eigen(corr.matrix)
    l <- sqrt(eigenvalue$values[1:n.factors])
    for (i in 1:n.factors)
      factor.loadings[, i] <- eigenvalue$vectors[, i] * l[i]
    for (i in 1:n.variables)
      h2[i] <- sum(factor.loadings[i, ] * factor.loadings[i, ])
    change <- max(abs(old.h2 - h2))
    old.h2 <- h2
    diag(corr.matrix) <- h2
  }
  if (determine) n.factors <- sum(eigenvalue$values > 1)
  return(list(loadings = factor.loadings[, 1:n.factors],
              factors = n.factors))
}

################################################################################

EFACompData <- function(data, f.max, n.pop = 10000, n.samples = 500,
                        alpha = .30, graph = FALSE, corr.type = "pearson") {
# Comparison data
#
# Args:
#   data      : Matrix to store the simulated data (matrix).
#   f.max     : Largest number of factors to consider (scalar).
#   n.pop     : Size of finite populations of comparison data (scalar, default
#               is 10,000 cases).
#   n.samples : Number of samples drawn from each population (scalar, default is
#               500)
#   alpha     : Alpha level when testing statistical significance of improvement
#               with additional factor (scalar, default is .30)
#   graph     : Whether to plot the fit of eigenvalues to those for comparison
#               data (default is FALSE)
#   corr.type : Type of correlation (character, default is "pearson", user can
#               also call "spearman").
#
# Returns:
#   Nothing, displays number of factors on screen.
#
  n <- dim(data)[1]
  n.variables <- dim(data)[2]
  corr.data <- cor(data, method = corr.type)
  eigs.data <- eigen(corr.data)$values
  rmsr.eigs <- matrix(0, nrow = n.samples, ncol = f.max)
  sig <- TRUE
  f.cd <- 1
  while ((f.cd <= f.max) & (sig)) {
    pop <- GenDataPopulation (data, n.factors = f.cd, n.cases = n.pop,
                              corr.type = corr.type)
    for (j in 1:n.samples) {
      v <- 0
      while (v == 0) {
        samp <- pop[sample(1:n.pop, size = n, replace = T), ]
        for (i in 1:n.variables)
          v <- v + min(samp[, i]) != max(samp[, i])
      }
      corr.samp <- cor(samp, method = corr.type)
      eigs.samp <- eigen(corr.samp)$values
      rmsr.eigs[j, f.cd] <- sqrt(sum((eigs.samp - eigs.data) *
                                (eigs.samp - eigs.data)) / n.variables)
    }
    if (f.cd > 1)
      sig <- (wilcox.test(rmsr.eigs[, f.cd], rmsr.eigs[, (f.cd - 1)],
              "less")$p.value < alpha)
    if (sig)
      f.cd <- f.cd + 1
  }
  cat("Number of factors to retain: ", f.cd - 1, "\n")
  if (graph) {
    if (sig) {
      x.max <- f.cd - 1
    } else {
      x.max <- f.cd
    }
    ys <- apply(rmsr.eigs[, 1:x.max], 2, mean)
    plot(x = 1:x.max, y = ys, ylim = c(0, max(ys)), xlab = "Factor",
         ylab = "RMSR Eigenvalue", type = "b",
         main = "Fit to Comparison Data")
    abline(v = f.cd - 1, lty = 3)
  }
}
