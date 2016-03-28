# Custom implementation for speedglm of extractAIC
library(speedglm)

extractAIC.speedglm <- function(fit, scale, k = 2, ...)
{
    edf <- fit$n - fit$df
    c(edf,  fit$deviance + k * edf)
}

nobs.speedglm <- function(fit, ...) {
  return (fit$ngoodobs)
}


mydropterm <-
  function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, sorted = FALSE, trace = FALSE, ...)
  {
    x <- model.matrix(object)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
      if(!is.character(scope))
        scope <- attr(terms(update.formula(object, scope)), "term.labels")
      if(!all(match(scope, tl, 0L)))
        stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ndrop <- match(scope, tl)
    rdf <- object$df.residual
    chisq <- object$deviance
    dfs <- numeric(ns)
    dev <- numeric(ns)
    y <- object$y
    if(is.null(y)) {
      y <- model.response(model.frame(object))
      if(!is.factor(y)) storage.mode(y) <- "double"
    }
    wt <- object$prior.weights
    if(is.null(wt)) wt <- rep.int(1, n)
    for(i in seq_len(ns)) {
      if(trace) {
        message(gettextf("trying - %s", scope[i]), domain = NA)
        utils::flush.console()
      }
      ii <- seq_along(asgn)[asgn == ndrop[i]]
      jj <- setdiff(seq(ncol(x)), ii)
      
      start <- NULL
      mustart <- NULL
      etastart <- NULL
      if (!is.null(object$call$start )) {
        #origstart <- eval(object$call$start)  # pass in start values.  These are original start values.
        # Now only use the start values that corresond to the actual new formula
        
        #start <- origstart[1:length(jj)]  
        start <- rep(1, length(jj))
        mustart <- rep(1, nrow(object$data))
        etastart <- rep(1, nrow(object$data))
        #print (c(jj, start))
      } 
      
      #print (start)
      
      new_formula <-  as.formula(paste0("~. -", names(coefficients(object))[ii]))
      print(new_formula)
      print(summary (mustart))
      print(summary (etastart))
      
      # hack for Gamma GLM on full sim data.  complains invalid start values.  Start at 1, far from zero.
      #z <- update(object, formula = new_formula, start=start, mustart=mustart, etastart=etastart)
#       z <-  glm.fit(x[, jj, drop = FALSE], y, wt, offset=object$offset,
#                     #start=start,  # pass in start values
#                     start=start,  # pass in start values
#                     etastart=etastart, #object$linear.predictors,
#                     mustart=mustart, #rep(mean(object$fitted.values, na.rm=TRUE), length(start)),
#                     family=object$family, control=object$control)  
#       
      z <- glm(new_formula, data=object$data, family=object$family,
                             start=start, mustart=mustart, etastart=etastart)
               
      dfs[i] <- z$rank
      dev[i] <- z$deviance
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    dev <- c(chisq, dev)
    dispersion <- if (is.null(scale) || scale == 0)
      summary(object, dispersion = NULL)$dispersion
    else scale
    fam <- object$family$family
    loglik <-
      if(fam == "gaussian") {
        if(scale > 0) dev/scale - n else n * log(dev/n)
      } else dev/dispersion
    aic <- loglik + k * dfs
    dfs <- dfs[1L] - dfs
    dfs[1L] <- NA
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic,
                      row.names = scope, check.names = FALSE)
    o <- if(sorted) order(aod$AIC) else seq_along(aod$AIC)
    if(all(is.na(aic))) aod <- aod[, -3]
    test <- match.arg(test)
    if(test == "Chisq") {
      dev <- pmax(0, loglik - loglik[1L])
      dev[1L] <- NA
      nas <- !is.na(dev)
      LRT <- if(dispersion == 1) "LRT" else "scaled dev."
      aod[, LRT] <- dev
      dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail=FALSE)
      aod[, "Pr(Chi)"] <- dev
    } else if(test == "F") {
      if(fam == "binomial" || fam == "poisson")
        warning(gettextf("F test assumes 'quasi%s' family", fam),
                domain = NA)
      dev <- aod$Deviance
      rms <- dev[1L]/rdf
      dev <- pmax(0, dev - dev[1L])
      dfs <- aod$Df
      rdf <- object$df.residual
      Fs <- (dev/dfs)/rms
      Fs[dfs < 1e-4] <- NA
      P <- Fs
      nas <- !is.na(Fs)
      P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail=FALSE)
      aod[, c("F value", "Pr(F)")] <- list(Fs, P)
    }
    aod <- aod[o, ]
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
    if(scale > 0)
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
  }


myaddterm.glm <-
  function(object, scope, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, sorted = FALSE, trace = FALSE, ...)
  {
    Fstat <- function(table, rdf) {
      dev <- table$Deviance
      df <- table$Df
      diff <- pmax(0, (dev[1L] - dev)/df)
      Fs <- diff/(dev/(rdf-df))
      Fs[df < .Machine$double.eps] <- NA
      P <- Fs
      nnas <- !is.na(Fs)
      P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas], lower.tail=FALSE)
      list(Fs=Fs, P=P)
    }
    if(missing(scope) || is.null(scope)) stop("no terms in scope")
    if(!is.character(scope))
      scope <- add.scope(object, update.formula(object, scope))
    if(!length(scope))
      stop("no terms in scope for adding to object")
    oTerms <- attr(terms(object), "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    dfs <- dev <- numeric(ns+1)
    names(dfs) <- names(dev) <- c("<none>", scope)
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(object, add.rhs)
    oc <- object$call
    Terms <- terms(new.form)
    oc$formula <- Terms
    ## model.frame.glm looks at the terms part for the environment
    fob <- list(call = oc, terms=Terms)
    class(fob) <- class(object)
    x <- model.matrix(Terms, model.frame(fob, xlev = object$xlevels),
                      contrasts = object$contrasts)
    n <- nrow(x)
    oldn <- length(object$residuals)
    y <- object$y
    newn <- length(y)
    if(newn < oldn)
      warning(sprintf(ngettext(newn,
                               "using the %d/%d row from a combined fit",
                               "using the %d/%d rows from a combined fit"),
                      newn, oldn), domain = NA)
    wt <- object$prior.weights
    if(is.null(wt)) wt <- rep(1, n)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
    if(int) ousex[1L] <- TRUE
    X <- x[, ousex, drop = FALSE]
    
    # Hack to get glm working for gamma()
    z <- update(object, formula=new.form, start=rep(1, length(ousex)))
#     z <-  glm.fit(X, y, wt, offset=object$offset,
#                   family=object$family, control=object$control)
    dfs[1L] <- z$rank
    dev[1L] <- z$deviance
    ## workaround for PR#7842. terms.formula may have flipped interactions
    sTerms <- sapply(strsplit(Terms, ":", fixed=TRUE),
                     function(x) paste(sort(x), collapse=":"))
    for(tt in scope) {
      if(trace) {
        message(gettextf("trying + %s", tt), domain = NA)
        utils::flush.console()
      }
      stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse=":")
      usex <- match(asgn, match(stt, sTerms), 0L) > 0L
      X <- x[, usex|ousex, drop = FALSE]
      # hack to get Gamma() glm working
      new.form <- as.formula(paste0("~. +", names(coefficients(object))[usex|ousex][2:length(ousex)]))
      z <- update(object, formula=new.form, start=rep(1, length(ousex)))
#       z <-  glm.fit(X, y, wt, offset=object$offset,
#                     family=object$family, control=object$control)
      dfs[tt] <- z$rank
      dev[tt] <- z$deviance
    }
    if (is.null(scale) || scale == 0)
      dispersion <- summary(object, dispersion = NULL)$dispersion
    else dispersion <- scale
    fam <- object$family$family
    if(fam == "gaussian") {
      if(scale > 0) loglik <- dev/scale - n
      else loglik <- n * log(dev/n)
    } else loglik <- dev/dispersion
    aic <- loglik + k * dfs
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L]) # same baseline for AIC
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic,
                      row.names = names(dfs), check.names = FALSE)
    o <- if(sorted) order(aod$AIC) else seq_along(aod$AIC)
    if(all(is.na(aic))) aod <- aod[, -3]
    test <- match.arg(test)
    if(test == "Chisq") {
      dev <- pmax(0, loglik[1L] - loglik)
      dev[1L] <- NA
      LRT <- if(dispersion == 1) "LRT" else "scaled dev."
      aod[, LRT] <- dev
      nas <- !is.na(dev)
      dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail=FALSE)
      aod[, "Pr(Chi)"] <- dev
    } else if(test == "F") {
      if(fam == "binomial" || fam == "poisson")
        warning(gettextf("F test assumes 'quasi%s' family", fam),
                domain = NA)
      rdf <- object$df.residual
      aod[, c("F value", "Pr(F)")] <- Fstat(aod, rdf)
    }
    aod <- aod[o, ]
    head <- c("Single term additions", "\nModel:", deparse(formula(object)))
    if(scale > 0)
      head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
  }


# stepAIC doesn't push in start values for cofficient estimation from GLM.
# But glm Gamma() complains when we don't give it start values.  Hack it here.
mystepAIC <- function (object, scope, scale = 0, direction = c("both", "backward", 
                                                  "forward"), trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
          k = 2, ...) 
{
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else extractAIC(x, k = 0)[2L]
  }
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste("\n", string[-1L], sep = "")
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:", 
                 deparse(formula(fit)), "\n")
    aod <- if (usingCp) 
      data.frame(Step = change, Df = ddf, Deviance = dd, 
                 `Resid. Df` = rdf, `Resid. Dev` = rd, Cp = AIC, 
                 check.names = FALSE)
    else data.frame(Step = change, Df = ddf, Deviance = dd, 
                    `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC, 
                    check.names = FALSE)
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme")) 
    object$call$fixed <- Terms
  else if (inherits(object, "gls")) 
    object$call$model <- Terms
  else object$call$formula <- Terms
 
    
    
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factors")
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep)) 
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  origstart <- eval(object$call$start)
  bAIC <- extractAIC(fit, scale, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
  if (bAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
  nm <- 1
  Terms <- terms(fit)
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))), 
        "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                         edf, change = "", AIC = bAIC)
  if (!is.null(keep)) 
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)) 
      ffac <- ffac[-st, ]
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- mydropterm(fit, scope$drop, scale = scale, trace = max(0, 
                                                                  trace - 1), k = k, ...)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
        if (any(is.finite(ch) & ch)) {
          warning("0 df terms are changing AIC")
          zdf <- zdf[!ch]
        }
        if (length(zdf) > 0L) 
          change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- myaddterm.glm(fit, scope$add, scale = scale, 
                        trace = max(0, trace - 1), k = k, ...)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], 
                                           sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) {
        print(aod[o, ])
        utils::flush.console()
      }
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    if (use.start) {
      #warning("'use.start' cannot be used with R's version of 'glm'")
      
      start <- origstart[1: length(coefficients(fit)) -1]  # NB:  this doesn't care which order the start values are. assumes they're all the same.
    } else {
      start <- NULL
    }
    
    usingCp <- match("Cp", names(aod), 0) > 0
    #fit <- update(fit, paste("~ .", change), start=start, evaluate = FALSE)
    fit <- update(fit, paste("~ .", change), start=start)
    #fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n) 
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      utils::flush.console()
    }
    if (bAIC >= AIC + 1e-07) 
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
                           edf, change = change, AIC = bAIC)
    if (!is.null(keep)) 
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep)) 
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}