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
    # Optmized version of the product, having dropped all factors that do not depend on x
    fast = function(x, mean, sigma) (x - mean)^2 + sigma^2 * (1 + nu) * log(((x - means)/sigmas)^2 + nu),
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
               means <- l$means
               if (means < mean_total) {
                 interval <- c(means, mean_total)
               } else {
                 interval <- c(mean_total, means)
               }

               x <- optimise(l$fast, interval, mean = mean_total, sigma = sigma_total)$minimum
               l$exact_sum( x, mean_total, sigma_total )
             }
  )) - log(sigma_total)
}

lfCombineMean <- function(mean_total, likelihood_functions) {
  sum(sapply(likelihood_functions, function(l) l$exact(mean_total)))
}


#' @export
studentMetaMean <- function(design) {

  design_matrix <- as.matrix( design[ , c("mean", "sigma", "nu") ])
  likelihood_functions <- apply(design_matrix, 1, getLF)

  optim(c(0, 1), lfMetaMean, likelihood_functions = likelihood_functions)
}
