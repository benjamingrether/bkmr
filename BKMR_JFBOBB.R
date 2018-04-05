# makeKpart <- function(r, Z) {
# Kpart <- as.matrix(dist(sqrt(matrix(r, byrow=TRUE, nrow(Z), ncol(Z)))*Z))^2
# Kpart
# }
makeKpart <- function(r, Z1, Z2 = NULL) {
  Z1r <- sweep(Z1, 2, sqrt(r), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  Kpart
}
makeVcomps <- function(r, lambda, Z, data.comps) {
  Kpart <- makeKpart(r, Z)
  V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
  if (data.comps$nlambda == 2) {
    V <- V + lambda[2]*data.comps$crossTT
  }
  cholV <- chol(V)
  Vinv <- chol2inv(cholV)
  logdetVinv <- -2*sum(log(diag(cholV)))
  Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  
  Vcomps
}

#' Fit Bayesian kernel machine regression
#'
#' Fits the Bayesian kernel machine regression (BKMR) model using Markov chain Monte Carlo (MCMC) methods.
#'
#' @export
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for \code{gaussian} and \code{binomial} families.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param Znew optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indictors for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#' @return an object of class "bkmrfit", which has the associated methods:
#' \itemize{
#'   \item \code{\link{print}} (i.e., \code{\link{print.bkmrfit}}) 
#'   \item \code{\link{summary}} (i.e., \code{\link{summary.bkmrfit}})
#' }
#' 
#' @seealso For guided examples, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' @references Bobb, JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, Godleski JJ, Coull BA (2015). Bayesian Kernel Machine Regression for Estimating the Health Effects of Multi-Pollutant Mixtures. Biostatistics 16, no. 3: 493-508.
#' @references Banerjee S, Gelfand AE, Finley AO, Sang H (2008). Gaussian predictive process models for large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(4), 825-848.
#' @import utils
kmbayes <- function(y, Z, X = NULL, iter = 1000, family = "gaussian", id = NULL, verbose = TRUE, Znew = NULL, starting.values = NULL, control.params = NULL, varsel = FALSE, groups = NULL, rmethod = "varying", est.h = FALSE) {
  
  missingX <- is.null(X)
  if (missingX) X <- matrix(0, length(y), 1)
  hier_varsel <- !is.null(groups)
  
  ##Argument check 1, required arguments without defaults
  ##check vector/matrix sizes
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  if (inherits(class(Z), "matrix") == FALSE)  Z <- as.matrix(Z)
  stopifnot (is.numeric(Z), nrow(Z) == length(y), anyNA(Z) == FALSE)
  if (inherits(class(X), "matrix") == FALSE)  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE) 
  
  ##Argument check 2: for those with defaults, write message and reset to default if invalid
  if (iter < 1) {
    message ("invalid input for iter, resetting to default value 1000")
    nsamp <- 1000
  } else {
    nsamp <- iter
  }
  if (!family %in% c("gaussian")) {
    stop("family", family, "not yet implemented; must specify either 'gaussian' or 'binomial'")
  }
  if (rmethod != "varying" & rmethod != "equal" & rmethod != "fixed") {
    message ("invalid value for rmethod, resetting to default varying")
    rmethod <- "varying"
  }
  
  ##Argument check 3: the rest id (below) znew, groups
  if (!is.null(id)) { 
    stopifnot(length(id) == length(y), anyNA(id) == FALSE)
  }
  if (!is.null(Znew)) { 
    if (class(Znew) != "matrix")  Znew <- as.matrix(Znew)
    stopifnot(is.numeric(Znew), ncol(Znew) == ncol(Z), anyNA(Znew) == FALSE)
  }
  if (!is.null(groups)) { 
    if (varsel == FALSE) {
      message ("groups should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(groups), length(groups) == ncol(Z), anyNA(groups) == FALSE)
    }
  }
  
  
  ## start JB code
  if (!is.null(id)) { ## for random intercept model
    randint <- TRUE
    id <- as.numeric(as.factor(id))
    nid <- length(unique(id))
    nlambda <- 2
    
    ## matrix that multiplies the random intercept vector
    TT <- matrix(0, length(id), nid)
    for (i in 1:nid) {
      TT[which(id == i), i] <- 1
    }
    crossTT <- tcrossprod(TT)
    print(TT)
    print(nid)
    rm(TT, nid)
  } else {
    randint <- FALSE
    nlambda <- 1
    crossTT <- 0
  }
  data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT)
  print(crossTT)
  rm(randint, nlambda, crossTT)
  
  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, nsamp, nrow(Z)),
                beta = matrix(0, nsamp, ncol(X)),
                lambda = matrix(NA, nsamp, data.comps$nlambda),
                sigsq.eps = rep(NA, nsamp),
                r = matrix(NA, nsamp, ncol(Z)),
                acc.r = matrix(0, nsamp, ncol(Z)),
                acc.lambda = matrix(0, nsamp, data.comps$nlambda),
                delta = matrix(1, nsamp, ncol(Z))
  )
  if (varsel) {
    chain$acc.rdelta <- rep(0, nsamp)
    chain$move.type <- rep(0, nsamp)
  }
  
  
  ## components to predict h(Znew)
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if (class(Znew) == "data.frame") Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
    ##Kpartall <- as.matrix(dist(rbind(Z,Znew)))^2
    chain$hnew <- matrix(0,nsamp,nrow(Znew))
    colnames(chain$hnew) <- rownames(Znew)
  }
  
  ## components if model selection is being done
  if (varsel) {
    ztest <- 1:ncol(Z)
    rdelta.update <- rdelta.comp.update
  } else {
    ztest <- NULL
  }
  
  ## control parameters
  control.params.default <- list(lambda.jump = rep(10, data.comps$nlambda), mu.lambda = rep(10, data.comps$nlambda), sigma.lambda = rep(10, data.comps$nlambda), a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3, b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1, r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
  } else {
    control.params <- control.params.default
  }
  
  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))
  
  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    if (!varsel) {
      stop("if doing grouped variable selection, must set varsel = TRUE")
    }
    rdelta.update <- rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }
  
  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)
  
  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
  }
  if (family == "gaussian") {
    if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
      lmfit0 <- lm(y ~ Z + X)
      if (is.null(starting.values$beta)) {
        coefX <- coef(lmfit0)[grep("X", names(coef(lmfit0)))]
        starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
      }
      if (is.null(starting.values$sigsq.eps)) {
        starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
      }
    } 
  } 
  
  print (starting.values)
  
  
  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1] <- starting.values$sigsq.eps
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }
  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h
  
  ## components
  Vcomps <- makeVcomps(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)
  
  ## start sampling ####
  chain$time1 <- Sys.time()
  for (s in 2:nsamp) {
    
    ## continuous version of outcome 
    if (family == "gaussian") {
      ycont <- y
    } 
    
    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)
    
    ## beta
    if (!missingX) {
      chain$beta[s,] <- beta.update(X = X, Vinv = Vcomps$Vinv, y = ycont, sigsq.eps = chain$sigsq.eps[s - 1])
    }
    
    ## \sigma_\epsilon^2
    if (family == "gaussian") {
      chain$sigsq.eps[s] <- sigsq.eps.update(y = ycont, X = X, beta = chain$beta[s,], Vinv = Vcomps$Vinv, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)
    }
    
    ## lambda
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- lambda.update(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = ycont, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim
    
    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- r.update(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update(r = rSim, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, ztest = ztest, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim
    
    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)
    
    if (est.h) {
      hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }
    
    ###################################################
    ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)
    
    if (!is.null(Znew)) {
      chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps)
    }
    
    ###################################################
    ## print details of the model fit so far
    opts <- set_verbose_opts(
      verbose_freq = control.params$verbose_freq, 
      verbose_digits = control.params$verbose_digits,
      verbose_show_ests = control.params$verbose_show_ests
    )
    print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = nsamp, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)
    
  }
  control.params$r.params <- NULL
  chain$time2 <- Sys.time()
  chain$iter <- nsamp
  chain$family <- family
  chain$starting.values <- starting.values
  chain$control.params <- control.params
  chain$X <- X
  chain$Z <- Z
  chain$y <- y
  chain$ztest <- ztest
  chain$data.comps <- data.comps
  if (!is.null(Znew)) chain$Znew <- Znew
  if (!is.null(groups)) chain$groups <- groups
  chain$varsel <- varsel
  class(chain) <- c("bkmrfit", class(chain))
  chain
}

#' Print basic summary of BKMR model fit
#'
#' \code{print} method for class "bkmrfit"
#'
#' @param x an object of class "bkmrfit"
#' @param digits the number of digits to show when printing
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
print.bkmrfit <- function(x, digits = 5, ...) {
  cat("Fitted object of class 'bkmrfit'\n")
  cat("Iterations:", x$iter, "\n")
  cat("Outcome family:", x$family, ifelse(x$family == "binomial", "(probit link)", ""), "\n")
  cat("Model fit on:", as.character(x$time2), "\n")
}

#' Summarizing BKMR model fits
#'
#' \code{summary} method for class "bkmrfit"
#'
#' @param object an object of class "bkmrfit"
#' @param q quantiles of posterior distribution to show
#' @param digits the number of digits to show when printing
#' @param show_ests logical; if \code{TRUE}, prints summary statistics of posterior distribution
#' @param show_MH logical; if \code{TRUE}, prints acceptance rates from the Metropolis-Hastings algorithm
#' @param ...	further arguments passed to or from other methods.
#'  
#' @export
summary.bkmrfit <- function(object, q = c(0.025, 0.975), digits = 5, show_ests = TRUE, show_MH = TRUE, ...) {
  x <- object
  elapsed_time <- difftime(x$time2, x$time1)
  
  print(x, digits = digits)  
  cat("Running time: ", round(elapsed_time, digits), attr(elapsed_time, "units"), "\n")
  cat("\n")
  
  if (show_MH) {
    cat("Acceptance rates for Metropolis-Hastings algorithm:\n")
    accep_rates <- data.frame()
    ## lambda
    nm <- "lambda"
    rate <- colMeans(x$acc.lambda[2:x$iter, ,drop = FALSE])
    if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
    accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    ## r_m
    if (!x$varsel) {
      nm <- "r"
      rate <- colMeans(x$acc.r[2:x$iter, , drop = FALSE])
      if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
    } else {
      nm <- "r/delta (overall)"
      rate <- mean(x$acc.rdelta[2:x$iter])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 1)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 1])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ##
      nm <- "r/delta  (move 2)"
      rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 2])
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      if (!is.null(x$groups)) {
        nm <- "r/delta  (move 3)"
        rate <- mean(x$acc.rdelta[2:x$iter][x$move.type[2:x$iter] == 3])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      }
    }
    print(accep_rates)
  }
  if (show_ests) {
    sel <- with(x, seq(floor(iter/2) + 1, iter))
    cat("\nParameter estimates (based on iterations ", min(sel), "-", max(sel), "):\n", sep = "")
    ests <- ExtractEsts(x, q = q, sel = sel)
    if (!is.null(ests$h)) {
      ests$h <- ests$h[c(1,2,nrow(ests$h)), ]
    }
    if (!is.null(ests$ystar)) {
      ests$ystar <- ests$ystar[c(1,2,nrow(ests$ystar)), ]
    }
    summ <- with(ests, rbind(beta, sigsq.eps, r, lambda))
    if (!is.null(ests$h)) {
      summ <- rbind(summ, ests$h)
    }
    if (!is.null(ests$ystar)) {
      summ <- rbind(summ, ests$ystar)
    }
    summ <- data.frame(param = rownames(summ), round(summ, digits))
    rownames(summ) <- NULL
    print(summ)
    if (x$varsel) {
      cat("\nPosterior inclusion probabilities:\n")
      pips <- ExtractPIPs(x)
      pips[, -1] <- round(pips[, -1], digits)
      print(pips)
    }
  }
}



make_r_params_comp <- function(r.params, whichcomp) {
  for(i in seq_along(r.params)) {
    if(length(r.params[[i]]) > 1) {
      r.params[[i]] <- r.params[[i]][whichcomp]
    }
  }
  r.params
}

set.r.params <- function(r.prior, comp, r.params) {
  if(r.prior == "gamma") {
    if(length(r.params$mu.r) > 1) r.params$mu.r <- r.params$mu.r[comp]
    if(length(r.params$sigma.r) > 1) r.params$sigma.r <- r.params$sigma.r[comp]
    if(length(r.params$r.jump1) > 1) r.params$r.jump1 <- r.params$r.jump1[comp]
    if(length(r.params$r.jump2) > 1) r.params$r.jump2 <- r.params$r.jump2[comp]
  }
  if(r.prior %in% c("unif", "invunif")) {
    if(length(r.params$r.a) > 1) r.params$r.a <- r.params$r.a[comp]
    if(length(r.params$r.b) > 1) r.params$r.b <- r.params$r.b[comp]
    if(length(r.params$r.jump2) > 1) r.params$r.jump2 <- r.params$r.jump2[comp]
  }
  r.params
}

set.r.MH.functions <- function(r.prior) {
  
  if(r.prior == "invunif") {
    # r.params <- list(r.a, r.b, r.jump2)
    rprior.logdens <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      ifelse(1/r.b <= x & x <= 1/r.a, -2*log(x) - log(r.b - r.a), log(0))
    }
    rprop.gen1 <- function(r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      1/runif(1, r.a, r.b)
    }
    rprop.logdens1 <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      ifelse(1/r.b <= x & x <= 1/r.a, -2*log(x) - log(r.b - r.a), log(0))
    }
    rprop.gen2 <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump)
    }
    rprop.logdens2 <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump))
    }
    rprop.gen <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump)
    }
    rprop.logdens <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a, mean = current, sd = r.jump))
    }
  }
  
  list(rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
}




beta.update <- function(X, Vinv, y, sigsq.eps) {
  XVinv <- crossprod(X, Vinv)
  Vbeta <- chol2inv(chol(XVinv %*% X))
  cholVbeta <- chol(Vbeta)
  betahat <- Vbeta %*% XVinv %*% y
  n01 <- rnorm(ncol(X))
  betahat + crossprod(sqrt(sigsq.eps)*cholVbeta, n01)
}

sigsq.eps.update <- function(y, X, beta, Vinv, a.eps=1e-3, b.eps=1e-3) {
  mu <- y - X%*%beta
  prec.y <- rgamma(1, shape=a.eps + nrow(X)/2, rate=b.eps + 1/2*crossprod(mu, Vinv)%*%mu)
  1/prec.y
}

ystar.update <- function(y, X, beta, h) {
  mu <-  drop(h + X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}
#' @importFrom tmvtnorm rtmvnorm
ystar.update.noh <- function(y, X, beta, Vinv, ystar) {
  mu <-  drop(X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Vinv, lower = lower, upper = upper, algorithm = "gibbs", start.value = ystar)
  #samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}

r.update <- function(r, whichcomp, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, data.comps, control.params, rprop.gen, rprop.logdens, rprior.logdens, ...) {
  # r.params <- set.r.params(r.prior = control.params$r.prior, comp = whichcomp, r.params = control.params$r.params)
  r.params <- make_r_params_comp(control.params$r.params, whichcomp)
  rcomp <- unique(r[whichcomp])
  if(length(rcomp) > 1) stop("rcomp should only be 1-dimensional")
  
  ## generate a proposal
  rcomp.star <- rprop.gen(current = rcomp, r.params = r.params)
  lambda.star <- lambda
  delta.star <- delta
  move.type <- NA
  
  ## part of M-H ratio that depends on the proposal distribution
  negdifflogproposal <- -rprop.logdens(rcomp.star, rcomp, r.params = r.params) + rprop.logdens(rcomp, rcomp.star, r.params = r.params)
  
  ## prior distribution
  diffpriors <- rprior.logdens(rcomp.star, r.params = r.params) - rprior.logdens(rcomp, r.params = r.params)
  
  r.star <- r
  r.star[whichcomp] <- rcomp.star
  
  ## M-H step
  return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

rdelta.comp.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen2, rprop.logdens1, rprior.logdens, rprior.logdens2, rprop.logdens2, rprop.gen1, ...) { ## individual variable selection
  r.params <- control.params$r.params
  a.p0 <- control.params$a.p0
  b.p0 <- control.params$b.p0
  delta.star <- delta
  r.star <- r
  
  move.type <- ifelse(all(delta[ztest] == 0), 1, sample(c(1,2),1))
  move.prob <- ifelse(all(delta[ztest] == 0), 1, 1/2)
  if(move.type == 1) {
    comp <- ifelse(length(ztest) == 1, ztest, sample(ztest, 1))
    r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)
    
    delta.star[comp] <- 1 - delta[comp]
    move.prob.star <- ifelse(all(delta.star[ztest] == 0), 1, 1/2)
    r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))
    
    diffpriors <- (lgamma(sum(delta.star[ztest]) + a.p0) + lgamma(length(ztest) - sum(delta.star[ztest]) + b.p0) - lgamma(sum(delta[ztest]) + a.p0) - lgamma(length(ztest) - sum(delta[ztest]) + b.p0)) + ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))
    
    negdifflogproposal <- -log(move.prob.star) + log(move.prob) - ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params))
    
  } else if(move.type == 2) {
    comp <- ifelse(length(which(delta == 1)) == 1, which(delta == 1), sample(which(delta == 1), 1))
    r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)
    
    r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)
    
    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)
    
    negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
  }
  
  lambda.star <- lambda
  
  ## M-H step
  return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

rdelta.group.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen1, rprior.logdens, rprop.logdens1, rprop.gen2, rprop.logdens2, ...) { ## grouped variable selection
  r.params <- control.params$r.params
  a.p0 <- control.params$a.p0
  b.p0 <- control.params$b.p0
  groups <- control.params$group.params$groups
  sel.groups <- control.params$group.params$sel.groups
  neach.group <- control.params$group.params$neach.group
  delta.star <- delta
  r.star <- r
  
  # if(length(mu.r) == 1) mu.r <- rep(mu.r, nz)
  # if(length(sigma.r) == 1) sigma.r <- rep(sigma.r, nz)
  
  delta.source <- sapply(sel.groups, function(x) ifelse(any(delta[which(groups == groups[x])] == 1), 1, 0))
  delta.source.star <- delta.source
  
  ## randomly select move type
  if(all(delta.source == 0)) {
    move.type <- 1
    move.prob <- 1
  } else if(length(which(neach.group > 1 & delta.source == 1)) == 0) {
    move.type <- sample(c(1, 3), 1)
    move.prob <- 1/2
  } else {
    move.type <- sample(1:3, 1)
    move.prob <- 1/3
  }
  # move.type <- ifelse(all(delta.source == 0), 1, ifelse(length(which(neach.group > 1 & delta.source == 1)) == 0, sample(c(1, 3), 1), sample(1:3, 1)))
  
  # print(move.type)
  
  if(move.type == 1) { ## randomly select a source and change its state (e.g., from being in the model to not being in the model)
    
    source <- sample(seq_along(delta.source), 1)
    source.comps <- which(groups == source)
    
    # r.params <- set.r.params(r.prior = control.params$r.prior, comp = source.comps, r.params = r.params)
    
    delta.source.star[source] <- 1 - delta.source[source]
    delta.star[source.comps] <- rmultinom(1, delta.source.star[source], rep(1/length(source.comps), length(source.comps)))
    move.prob.star <- ifelse(all(delta.source.star == 0), 1, ifelse(length(which(neach.group > 1 & delta.source.star == 1)) == 0, 1/2, 1/3))
    
    ## which component got switched
    comp <- ifelse(delta.source[source] == 1, source.comps[which(delta[source.comps] == 1)], source.comps[which(delta.star[source.comps] == 1)])
    r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)
    
    r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))
    
    # diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprior.logdens(x = r.sel, r.params = r.params))
    diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))
    
    # negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprop.logdens1(x = r.sel, r.params = r.params)))
    negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params)))
    
  } else if(move.type == 2) { ## randomly select a multi-component source that is in the model and change which component is included
    
    tmp <- which(neach.group > 1 & delta.source == 1)
    source <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
    source.comps <- which(groups == source)
    
    oldcomp <- source.comps[delta[source.comps] == 1]
    tmp <- source.comps[delta[source.comps] == 0]
    comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
    
    r.params.oldcomp <- set.r.params(r.prior = control.params$r.prior, comp = oldcomp, r.params = r.params)
    r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)
    
    delta.star[oldcomp] <- 0
    delta.star[comp] <- 1
    
    r.star[oldcomp] <- 0
    r.star[comp] <- rprop.gen1(r.params = r.params)
    
    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[oldcomp], r.params = r.params.oldcomp)
    
    negdifflogproposal <- -rprop.logdens1(r.star[comp], r.params = r.params) + rprop.logdens1(r[oldcomp], r.params = r.params.oldcomp)
    
  } else if(move.type == 3) { ## randomly select a component that is in the model and update it
    tmp <- which(delta == 1)
    comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
    
    r.params <- set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)
    
    r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)
    
    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)
    
    negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
  }
  
  lambda.star <- lambda
  
  ## M-H step
  return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

lambda.update <- function(r, delta, lambda, whichcomp=1, y, X, Z = Z, beta, sigsq.eps, Vcomps, data.comps, control.params) {
  lambda.jump <- control.params$lambda.jump[whichcomp]
  mu.lambda <- control.params$mu.lambda[whichcomp]
  sigma.lambda <- control.params$sigma.lambda[whichcomp]
  lambdacomp <- lambda[whichcomp]
  
  ## generate a proposal
  lambdacomp.star <- rgamma(1, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2)
  r.star <- r
  delta.star <- delta
  move.type <- NA
  
  ## part of M-H ratio that depends on the proposal distribution
  negdifflogproposal <- -dgamma(lambdacomp.star, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2, log=TRUE) + dgamma(lambdacomp, shape=lambdacomp.star^2/lambda.jump^2, rate=lambdacomp.star/lambda.jump^2, log=TRUE)
  
  ## prior distribution
  diffpriors <- dgamma(lambdacomp.star, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE) - dgamma(lambdacomp, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE)
  
  lambda.star <- lambda
  lambda.star[whichcomp] <- lambdacomp.star
  
  ## M-H step
  return(MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

MHstep <- function(r, lambda, lambda.star, r.star, delta, delta.star, y, X, Z, beta, sigsq.eps, diffpriors, negdifflogproposal, Vcomps, move.type, data.comps) {
  ## compute log M-H ratio
  Vcomps.star <- makeVcomps(r.star, lambda.star, Z, data.comps)
  mu <- y - X%*%beta
  diffliks <- 1/2*Vcomps.star$logdetVinv - 1/2*Vcomps$logdetVinv - 1/2/sigsq.eps*crossprod(mu, Vcomps.star$Vinv - Vcomps$Vinv)%*%mu
  logMHratio <- diffliks + diffpriors + negdifflogproposal
  logalpha <- min(0,logMHratio)
  
  ## return value
  acc <- FALSE
  if( log(runif(1)) <= logalpha ) {
    r <- r.star
    delta <- delta.star
    lambda <- lambda.star
    Vcomps <- Vcomps.star
    acc <- TRUE
  }
  return(list(r=r, lambda=lambda, delta=delta, acc=acc, Vcomps=Vcomps, move.type=move.type))
}

h.update <- function(lambda, Vcomps, sigsq.eps, y, X, beta, r, Z, data.comps) {
  if (is.null(Vcomps)) {
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
  }
  if(is.null(Vcomps$Q)) {
    Kpart <- makeKpart(r, Z)
    K <- exp(-Kpart)
    Vinv <- Vcomps$Vinv
    lambda <- lambda[1] ## in case with random intercept (randint==TRUE), where lambda is 2-dimensional
    lamKVinv <- lambda*K%*%Vinv
    h.postmean <- lamKVinv%*%(y-X%*%beta)
    ##h.postvar <- sigsq.eps*lamKVinv
    h.postvar <- sigsq.eps*lambda*(K - lamKVinv%*%K)
    h.postvar.sqrt <- try(chol(h.postvar), silent=TRUE)
    if(class(h.postvar.sqrt) == "try-error") {
      sigsvd <- svd(h.postvar)
      h.postvar.sqrt <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    hsamp <- h.postmean + crossprod(h.postvar.sqrt, rnorm(length(h.postmean)))
    hcomps <- list(hsamp = hsamp)
  } else {
    h.star.postvar.sqrt <- sqrt(sigsq.eps*lambda)*forwardsolve(t(Vcomps$cholR), Vcomps$Q)
    h.star.postmean <- lambda[1]*Vcomps$Q %*% Vcomps$Rinv %*% Vcomps$K10 %*% (y - X %*% beta)
    hsamp.star <- h.star.postmean + crossprod(h.star.postvar.sqrt, rnorm(length(h.star.postmean)))
    hsamp <- t(Vcomps$K10) %*% Vcomps$Qinv %*% hsamp.star
    hcomps <- list(hsamp = hsamp, hsamp.star = hsamp.star)
  }
  hcomps
}

newh.update <- function(Z, Znew, Vcomps, lambda, sigsq.eps, r, y, X, beta, data.comps) {
  
  n0 <- nrow(Z)
  n1 <- nrow(Znew)
  nall <- n0 + n1
  # Kpartall <- makeKpart(r, rbind(Z, Znew))
  # Kmat <- exp(-Kpartall)
  # Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
  # Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
  # Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
  Kmat1 <- exp(-makeKpart(r, Znew))
  Kmat10 <- exp(-makeKpart(r, Znew, Z))
  
  if(is.null(Vcomps)) {
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
  }
  Vinv <- Vcomps$Vinv
  
  lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
  Sigma.hnew <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
  mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
  root.Sigma.hnew <- try(chol(Sigma.hnew), silent=TRUE)
  if(class(root.Sigma.hnew) == "try-error") {
    sigsvd <- svd(Sigma.hnew)
    root.Sigma.hnew <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  hsamp <- mu.hnew + crossprod(root.Sigma.hnew, rnorm(n1))
  
  
  hsamp
}

## function to obtain posterior samples of h(znew) from fit of Bayesian kernel machine regression
predz.samps <- function(fit, Znew, verbose = TRUE) {
  if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
  if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
  Z <- fit$Z
  if(ncol(Z) != ncol(Znew)) {
    stop("Znew must have the same number of columns as Z")
  }
  
  hnew.samps <- sapply(1:fit$nsamp, function(s) {
    if(s%%(fit$nsamp/10)==0 & verbose) print(s)
    newh.update(Z = Z, Znew = Znew, Vcomps = NULL, lambda = fit$lambda[s], sigsq.eps = fit$sigsq.eps[s], r = fit$r[s,], y = fit$y, X = fit$X, beta = fit$beta[s,], data.comps = fit$data.comps)
  })
  rownames(hnew.samps) <- rownames(Znew)
  t(hnew.samps)
}

## function to approximate the posterior mean and variance as a function of the estimated tau, lambda, beta, and sigsq.eps
newh.postmean <- function(fit, Znew, sel) {
  if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
  if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
  
  Z <- fit$Z
  X <- fit$X
  y <- fit$y
  data.comps <- fit$data.comps
  lambda <- colMeans(fit$lambda[sel, ,drop = FALSE])
  sigsq.eps <- mean(fit$sigsq.eps[sel])
  r <- colMeans(fit$r[sel,])
  beta <- colMeans(fit$beta[sel, ,drop=FALSE])
  
  n0 <- nrow(Z)
  n1 <- nrow(Znew)
  nall <- n0 + n1
  Kpartall <- makeKpart(r, rbind(Z, Znew))
  Kmat <- exp(-Kpartall)
  Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
  Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
  Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
  
  Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
  Vinv <- Vcomps$Vinv
  
  lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
  Sigma.hnew <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
  mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
  
  ret <- list(postmean = drop(mu.hnew), postvar = Sigma.hnew)
  ret
}






#' Options for printing summary of model fit to the console
#'
#' Set options for what will be printed to the console when verbose = TRUE in the main kmbayes function
#'
#' @param verbose_freq After this percentage of iterations has been completed the summary of the model fit so far will be printed to the console 
#' @param verbose_show_ests TRUE or FALSE: flag indicating whether to print out summary statistics of all posterior samples obtained up until this point, for select parameters
#' @param verbose_digits Number of digits to be printed to the console 
#'
#' @export
#'
library(magrittr)
set_verbose_opts <- function(verbose_freq = NULL, verbose_show_ests = NULL, verbose_digits = NULL) {
  if (is.null(verbose_freq)) verbose_freq <- 10
  if (is.null(verbose_digits)) verbose_digits <- 5
  if (is.null(verbose_show_ests)) verbose_show_ests <- FALSE
  opts <- list(
    verbose_freq = verbose_freq,
    verbose_digits = verbose_digits,
    verbose_show_ests = verbose_show_ests
  )
  opts
}

print_diagnostics <- function(verbose, opts, curr_iter, tot_iter, chain, varsel, hier_varsel, ztest, Z, groups) {
  verbose_freq <- opts$verbose_freq
  verbose_digits <- opts$verbose_digits
  verbose_show_ests <- opts$verbose_show_ests
  s <- curr_iter
  nsamp <- tot_iter
  perc_iter_completed <- round(100*curr_iter/tot_iter, 1)
  
  all_iter <- 100*(1:nsamp)/nsamp
  sel_iter <- seq(verbose_freq, 100, by = verbose_freq)
  print_iter <- sapply(sel_iter, function(x) min(which(all_iter >= x)))
  
  elapsed_time <- difftime(Sys.time(), chain$time1)
  
  if (s %in% print_iter) {
    #if (verbose) message("------------------------------------------")
    if (verbose) cat("\n")
    message("Iteration: ", s, " (", perc_iter_completed, "% completed; ", round(elapsed_time, verbose_digits), " ", attr(elapsed_time, "units"), " elapsed)")
    
    if (verbose) {
      cat("Acceptance rates for Metropolis-Hastings algorithm:\n")
      accep_rates <- data.frame()
      ## lambda
      nm <- "lambda"
      rate <- colMeans(chain$acc.lambda[2:s, ,drop = FALSE])
      if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
      accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      ## r_m
      if (!varsel) {
        nm <- "r"
        rate <- colMeans(chain$acc.r[2:s, , drop = FALSE])
        if (length(rate) > 1) nm <- paste0(nm, seq_along(rate))
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
      } else {
        nm <- "r/delta (overall)"
        rate <- mean(chain$acc.rdelta[2:s])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        ##
        nm <- "r/delta  (move 1)"
        rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 1])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        ##
        nm <- "r/delta  (move 2)"
        rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 2])
        accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        if (hier_varsel) {
          nm <- "r/delta  (move 3)"
          rate <- mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 3])
          accep_rates %<>% rbind(data.frame(param = nm, rate = rate))
        }
      }
      print(accep_rates)
      
      ## extra information
      if (verbose_show_ests) {
        cat("\nCurrent parameter estimates:\n")
        chain$varsel <- varsel
        class(chain) <- c("bkmrfit", class(chain))
        chain$Z <- Z
        if (hier_varsel) chain$groups <- groups
        ests <- ExtractEsts(chain, q = c(0.025, 0.975), sel = 2:s)
        #ests$h <- ests$h[c(1,2,nrow(ests$h)), ]
        summ <- with(ests, rbind(beta, sigsq.eps, r, lambda))
        summ <- data.frame(param = rownames(summ), round(summ, verbose_digits))
        rownames(summ) <- NULL
        print(summ)
        if (varsel) {
          cat("\nCurrent posterior inclusion probabilities:\n")
          pips <- ExtractPIPs(chain, sel = 2:s)
          pips[, -1] <- round(pips[, -1], verbose_digits)
          print(pips)
        }
      }
    }
  }
  
}





# Validate control params list
##components of list
##lambda.jump             / default=10
##mu.lambda, sigma.lambda / default=10
##a.p0, b.p0              / default=1
##r.prior                 / default = "gamma", alt=invunif, unif
##a.sigsq, b.sigsq        / default=0.001
##mu.r, sigma.r           / default=5
##r.muprop                / default=1
##r.jump                  / default=0.2
##r.jump1, r.jump2        / default=2, 0.2
##r.a, r.b                / default=0, 100

# 
validateStartingParams <- function(varsel, Ylength, Xwidth, Zwidth, starting.params) {
  message ("Validating starting.params...")
  print(starting.params)
  stopifnot(starting.params$b.sigsq.eps > 0, starting.params$lambda > 0)
  ##messages only for scalar where vector required, expansion happens in main function
  ##beta length, any values
  if (length(starting.params$beta) != Xwidth) {
    message("beta should be a vector of length equal to the number of columns of X.  A vector will be created of repeating the input value.")
  }
  ##h.hat length and values
  if (length(starting.params$h.hat) != Ylength) {
    message("h.hat should be a vector of length equal to number of rows in Y.  A vector will be created of repeating the input value.")
  }
  for (i in 1:length(starting.params$h.hat)) {
    stopifnot(starting.params$h.hat > 0) 
  }
  ##delta length, 0 or 1 are valid values
  if (length(starting.params$delta) != Zwidth) {
    message("delta should be a vector of length equal to the number of columns of Z.  A vector will be created of repeating the input value.")
  }
  for (i in 1:length(starting.params$delta)) {
    stopifnot(starting.params$delta == 1 || starting.params$delta == 0) 
  }
  ## r length depends on varsel, truncate here but expand if necessary in main function
  if (varsel == TRUE) {
    if (length(starting.params$r) != Zwidth) {
      message("r should be a vector of length equal to the number of columns of Z.  A vector will be created of repeating the input value.")
    }
    else {
      if (length(starting.params$r) > 1) {
        message("r should a scalar.  Vector input will be truncated.")
        starting.params$r=starting.params$r[1]
      }
    }
  }
  for (i in 1:length(starting.params$r)) {
    stopifnot(starting.params$delta > 0) 
  }
}
