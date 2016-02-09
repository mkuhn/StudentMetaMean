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

safe_range <- function(v) {
  r <- range(v)
  if (r[1] == r[2]) {
    r[1] <- r[1] - 0.01
    r[2] <- r[2] + 0.01
  }
  r
}

getLF <- function(row) {

  means <- row[1]
  sigmas <- row[2]
  nu <- row[3]

  D0 <- -dstudent(means, nu, means, sigmas, log=T)

  list(
    # Exact density function
    exact = function(x) -dstudent(x, nu, means, sigmas, log=T),
    # Exact calculation of product of Student's t and normal distribution
    exact_sum = function(mean, sigma) D0-dnorm(means, mean, sigma, log=T),
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

  L <- lapply(likelihood_functions, function(l) l$exact_sum( mean_total, sigma_total ))

  Reduce("+", L) - log(sigma_total)
}

lfCombineMean <- function(mean_total, likelihood_functions) {
  Reduce("+", lapply(likelihood_functions, function(l) l$exact(mean_total)))
}


#' @export
studentMetaMean <- function(design = NULL, means = NULL, sigmas = NULL, nus = NULL, min_sigma = 1) {

  if (is.null(design)) {
    if (any(is.null(means), is.null(sigmas), is.null(nus)))
      stop("Missing either 'design' or 'means'/'simgas'/'nus' argument")

    design_matrix <- matrix( c(means, sigmas, nus), ncol=3)
  } else {
    design_matrix <- as.matrix( design[ , c("mean", "sigma", "nu") ])
  }

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
  o2 <- optimise(pMerged, safe_range(design_matrix[,1]), maximum = T)

  l <- list(
    xmin = min(sapply(likelihood_functions, function(l) l$xmin)),
    xmax = max(sapply(likelihood_functions, function(l) l$xmax)),
    mean = o2$maximum,
    sigma = sigma_total,
    likelihood_functions = likelihood_functions,
    pMerged = pMerged
  )
  l
}

#' @export
findHDI <- function(l, p = 0.95, stepsize = 0.01, density_resolution = 0.001) {

  xmid <- round(l$mean/stepsize)*stepsize

  # first, check if we need to expand the considered interval
  ytarget <- density_resolution*l$pMerged(xmid)
  while (l$pMerged(l$xmin) > ytarget) l$xmin <- l$mean - 1.5 * (l$mean - l$xmin)
  while (l$pMerged(l$xmax) > ytarget) l$xmax <- l$mean + 1.5 * (l$xmax - l$mean)

  xmin <- floor(l$xmin/stepsize)*stepsize
  xmax <- ceiling(l$xmax/stepsize)*stepsize

  # compute whole interval at once, thanks to the vectorization this is quite fast
  X <- seq(xmin, xmax, stepsize)
  Y <- l$pMerged(X)
  N <- length(Y)

  target <- p*(sum(Y, na.rm=T))

  mid <- as.integer((xmid-xmin)/stepsize+1)

  i <- mid
  j <- mid
  vi <- Y[i]
  vj <- Y[j]

  s <- 0

  while (s < target) {
    if (vi > vj) {
      i <- i-1

      # extend the range if needed
      if (i == 0) {
        ext <- as.integer(N/2)
        xminext <- xmin - stepsize * ext
        Xext <- seq(xminext, xmin-stepsize, stepsize)
        X <- c(Xext, X)
        Y <- c(l$pMerged(Xext), Y)
        xmin <- xminext
        i <- i + ext
        j <- j + ext
        N <- N + ext
        target <- p*(sum(Y, na.rm=T))
      }

      s <- s + vi/2
      vi <- Y[i]
      s <- s + vi/2
    } else {
      j <- j+1

      # extend the range if needed
      if (j > N) {
        ext <- as.integer(N/2)
        xmaxext <- xmax + stepsize * ext
        Xext <- seq(xmax+stepsize, xmaxext, stepsize)
        X <- c(X, Xext)
        Y <- c(Y, l$pMerged(Xext))
        xmax <- xmaxext
        N <- N + ext
        target <- p*(sum(Y, na.rm=T))
      }

      s <- s + vj/2
      vj <- Y[j]
      s <- s + vj/2
    }
  }

  if (j-i < 10) message("Step size may be too small!")

  c( xmin+stepsize*(i-0.5), xmin + stepsize*(j-0.5))
}

#' @export
plotMetaMean <- function(l, p = 0.95, stepsize=NULL) {

  if (is.null(stepsize)) {
    stepsize <- (l$xmax - l$xmin) / 1000
  }

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

  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, color=n, fill=n)) + ggplot2::geom_line() + ggplot2::facet_grid(kind~., scales = "free")
  p <- p + ggplot2::geom_polygon(data=dm)
  p + ggplot2::geom_vline(xintercept=l$mean)
}

