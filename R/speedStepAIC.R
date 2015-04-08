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
