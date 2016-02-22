#' @useDynLib studentMetaMean
#' @importFrom Rcpp sourceCpp

log_dstudent <- function(x, nu, mean, sigma) {
  dt((x-mean)/sigma, nu, log=T)-log(sigma)
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

  m <- row[1]
  sigma <- row[2]
  nu <- row[3]

  if (is.na(nu)) {
    # normal distribution
    list(
      density = function(x) dnorm(x, m, sigma, log=T),
      xmin = qnorm(0.005, m, sigma),
      xmax = qnorm(1-0.005, m, sigma)
    )

  } else {
    # Student's t distribution
    list(
      density = function(x) log_dstudent(x, nu, m, sigma),
      xmin = qstudent(0.005, nu, m, sigma),
      xmax = qstudent(1-0.005, nu, m, sigma)
    )
  }
}

# Likelihood that the distribution means are part of a normal distribution
lfMetaMean <- function(sigma_total, mean_total, means, INF, min_sigma) {
  if (sigma_total < min_sigma) return(-INF)
  dnorms(mean_total, means, sigma_total)
  # sapply(mean_total, function(m) sum(dnorm(m, means, sigma_total, log=T))) + log(sigma_total)
}

# Likelihood that the distributions have one common mean
lfCombineMean <- function(mean_total, likelihood_functions) {
  Reduce("+", lapply(likelihood_functions, function(l) l$density(mean_total)))
}


#' @export
studentMetaMean <- function(design = NULL, means = NULL, sigmas = NULL, nus = NULL, min_sigma = 0.01) {

  if (is.null(design)) {
    if (any(is.null(means), is.null(sigmas), is.null(nus)))
      stop("Missing either 'design' or 'means'/'simgas'/'nus' argument")

    design_matrix <- matrix( c(means, sigmas, nus), ncol=3)
  } else {
    design_matrix <- as.matrix( design[ , c("mean", "sigma", "nu") ])
  }

  likelihood_functions <- apply(design_matrix, 1, getLF)
  means <- design_matrix[, 1]
  sigmas <- design_matrix[, 2]
  nus <- design_matrix[, 3]
  log_beta_precomp <- log(beta(nus/2, 0.5))

  nus[ is.na(nus) ] <- 0

  INF <- .Machine$double.xmax/(1+length(likelihood_functions))

  sigma_total <- max(min_sigma, sd(means))

  pMerged <- function(mean_total) {
    exp(lfMetaMean(sigma_total, mean_total, means, INF, min_sigma)) +
      exp(dstud(mean_total, means, nus, sigmas, log_beta_precomp))
      # exp(lfCombineMean(mean_total, likelihood_functions))
  }

  # find best mean
  o <- optimise(pMerged, safe_range(design_matrix[,1]), maximum = T)

  l <- list(
    xmin = min(sapply(likelihood_functions, function(l) l$xmin)),
    xmax = max(sapply(likelihood_functions, function(l) l$xmax)),
    mean = o$maximum,
    sigma = sigma_total,
    pMerged = pMerged,
    likelihood_functions = likelihood_functions
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
pMetaMean <- function(l, p) {
  tol <- l$pMerged(p)
  A <- integrate( l$pMerged, -Inf, p, abs.tol = tol  )$value
  B <- integrate( l$pMerged, p, Inf, abs.tol = tol  )$value
  A / (A+B)
}

#' @export
plotMetaMean <- function(l, p = 0.95, stepsize=NULL, xmin=NULL, xmax=NULL) {

  if (is.null(xmin)) xmin <- l$xmin
  if (is.null(xmax)) xmax <- l$xmax

  if (is.null(stepsize)) {
    stepsize <- (xmax - xmin) / 1000
  }

  likelihood_functions <- l$likelihood_functions
  xmin <- floor(xmin/stepsize)*stepsize
  xmax <- ceiling(xmax/stepsize)*stepsize

  X <- seq(xmin, xmax, stepsize)

  d <- rbind(
    data.frame(x=X, y=l$pMerged(X), kind="metaMean", n = 0),
    do.call(rbind, lapply(1:length(likelihood_functions),
         function(i) data.frame(x=X, y=exp(likelihood_functions[[i]]$density(X)), kind="input", n = i))
    )
  )

  hdi <- findHDI(l, p, stepsize)

  have_hdi <- hdi[1] < hdi[2]

  d$n <- factor(d$n)

  if (have_hdi) {
    dm <- d[ d$n == 0, ]
    dm <- dm[ dm$x >= hdi[1] & dm$x <= hdi[2], ]

    dm0 <- dm[1,]
    dm0$y <- 0
    dmN <- dm[nrow(dm),]
    dmN$y <- 0

    dm <- rbind(dm0, dm, dmN)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, color=n, fill=n)) + ggplot2::geom_line() + ggplot2::facet_grid(kind~., scales = "free")
  if (have_hdi) p <- p + ggplot2::geom_polygon(data=dm)
  p + ggplot2::geom_vline(xintercept=l$mean)
}
