#' @useDynLib studentMetaMean
#' @importFrom Rcpp sourceCpp

dstudent <- function(x, nu, mean, sigma, log = FALSE) {
  if (log) {
    dt((x-mean)/sigma, nu, log=T)-log(sigma)
  } else {
    dt((x-mean)/sigma, nu, log=F)/sigma
  }
}

qstudent <- function(p, nu, mean, sigma) (qt(p, nu)*sigma+mean)

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
    nu = nu,
    xmin = qstudent(0.005, nu, means, sigmas),
    xmax = qstudent(1-0.005, nu, means, sigmas)
  )
}


lfMetaMean <- function(P, likelihood_functions, INF, min_sigma) {
  mean_total <- P[[1]]
  sigma_total <- P[[2]]

  if (sigma_total < min_sigma) return(INF)

  colSums(do.call(rbind, (lapply(likelihood_functions,
             function(l) {
               x <- optimise_overlap(mean_total, sigma_total, l$nu, l$means, l$sigmas)
               l$exact_sum( x, mean_total, sigma_total )
             }
  )))) - log(sigma_total)
}

lfCombineMean <- function(mean_total, likelihood_functions) {
  colSums(do.call(rbind, (lapply(likelihood_functions, function(l) l$exact(mean_total)))))
}


#' @export
studentMetaMean <- function(design, min_sigma = 1) {

  if (max(design$nu) > 1000) message("Numerical instabilities may occur for large values of nu")

  design_matrix <- as.matrix( design[ , c("mean", "sigma", "nu") ])

  likelihood_functions <- apply(design_matrix, 1, getLF)

  INF <- .Machine$double.xmax/(1+length(likelihood_functions))

  # first optimisation using both free mean and sigma to find best sigma
  o1 <- optim(c(mean(design_matrix[,1]), sd(design_matrix[,1])),
        lfMetaMean,
        likelihood_functions = likelihood_functions,
        INF = INF,
        min_sigma = min_sigma
      )

  sigma_total <- o1$par[2]

  pMerged <- function(mean_total) {
    exp(-lfMetaMean(list(mean_total, sigma_total), likelihood_functions, INF, min_sigma)) +
      exp(-lfCombineMean(mean_total, likelihood_functions))
  }

  # second optimisation:
  o2 <- optimise(pMerged, range(design_matrix[,1]), maximum = T)

  l <- list(
    xmin = min(sapply(likelihood_functions, function(l) l$xmin)),
    xmax = max(sapply(likelihood_functions, function(l) l$xmax)),
    mean = o2$maximum,
    likelihood_functions = likelihood_functions,
    lfMetaMean = lfMetaMean,
    lfCombineMean = lfCombineMean,
    pMerged = pMerged
  )
  l
}

#' @export
findHDI <- function(l, p = 0.95, stepsize = 0.01) {

  xmin <- floor(l$xmin/stepsize)*stepsize
  xmax <- ceiling(l$xmax/stepsize)*stepsize
  xmid <- round(l$mean/stepsize)*stepsize

  X <- seq(xmin, xmax, stepsize)
  Y <- l$pMerged(X)

  target <- p*(sum(Y, na.rm=T))

  mid <- (xmid-xmin)/stepsize+1

  i <- mid
  j <- mid
  vi <- Y[i]
  vj <- Y[j]

  s <- 0

  while (s < target) {
    if (vi > vj) {
      i <- i-1
      s <- s + vi/2
      vi <- Y[i]
      s <- s + vi/2
    } else {
      j <- j+1
      s <- s + vj/2
      vj <- Y[j]
      s <- s + vj/2
    }
  }

  if (j-i < 10) message("Step size may be too small!")

  c( xmin+stepsize*(i-0.5), xmin + stepsize*(j-0.5))
}

#' @export
plotMetaMean <- function(l, p = 0.95, stepsize=0.1) {

  likelihood_functions <- l$likelihood_functions
  xmin <- floor(l$xmin/stepsize)*stepsize
  xmax <- ceiling(l$xmax/stepsize)*stepsize

  X <- seq(xmin, xmax, stepsize)

  d <- rbind(
    data.frame(x=X, y=l$pMerged(X), kind="metaMean", n = 0),
    do.call(rbind, lapply(1:length(likelihood_functions),
         function(i) data.frame(x=X, y=exp(-likelihood_functions[[i]]$exact(X)), kind="input", n = i))
    )
  )

  hdi <- findHDI(l, p, stepsize)

  d$n <- factor(d$n)
  dm <- d[ d$n == 0, ]
  dm <- dm[ dm$x >= hdi[1] & dm$x <= hdi[2], ]

  dm0 <- dm[1,]
  dm0$y <- 0
  dmN <- dm[nrow(dm),]
  dmN$y <- 0

  dm <- rbind(dm0, dm, dmN)

  p <- ggplot2::ggplot(d, aes(x, y, color=n, fill=n)) + geom_line() + facet_grid(kind~., scales = "free")
  p <- p + geom_polygon(data=dm)
  p + geom_vline(xintercept=l$mean)
}

