#' @useDynLib studentMetaMean
#' @importFrom Rcpp sourceCpp

dstudent <- function(x, nu, mean, sigma, log = FALSE) {
  if (log) {
    dt((x-mean)/sigma, nu, log=T)-log(sigma)
  } else {
    dt((x-mean)/sigma, nu, log=F)/sigma
  }
}

getLF <- function(row) {

  means <- row[1]
  sigmas <- row[2]
  nu <- row[3]

  list(
    # Exact density function
    exact = function(x) -dstudent(x, nu, means, sigmas, log=T),
    # Exact calculation of product of Student's t and normal distribution
    exact_sum = function(x, mean, sigma) -dstudent(x, nu, means, sigmas, log=T)-dnorm(x, mean, sigma, log=T),
    means = means,
    sigmas = sigmas,
    nu = nu
  )
}


INF <- .Machine$double.xmax/100

lfMetaMean <- function(P, likelihood_functions) {
  mean_total <- P[1]
  sigma_total <- P[2]

  if (sigma_total <= 0) return(INF)

  sum(sapply(likelihood_functions,
             function(l) {
               x <- optimise_overlap(mean_total, sigma_total, l$nu, l$means, l$sigmas)
               l$exact_sum( x, mean_total, sigma_total )
             }
  )) - log(sigma_total)
}

lfCombineMean <- function(mean_total, likelihood_functions) {
  sum(sapply(likelihood_functions, function(l) l$exact(mean_total)))
}


#' @export
studentMetaMean <- function(design, min_sigma = 0.001) {

  design_matrix <- as.matrix( design[ , c("mean", "sigma", "nu") ])
  likelihood_functions <- apply(design_matrix, 1, getLF)

  # first optimisation using both free mean and sigma
  o <- optim(c(mean(design_matrix[,1]), sd(design_matrix[,1])),
        lfMetaMean,
        likelihood_functions = likelihood_functions)

  l <- list()

  # if sigma is very small, compute weighted average
  if (o$par[2] < min_sigma) {
    l$mean <- optimise(lfCombineMean, range(design_matrix[,1]))$minimum
    l$sigma <- 0
  } else {
    l$mean <- o$par[1]
    l$sigma <- o$par[2]
  }
}
