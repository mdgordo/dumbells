#' Function for taking the inverse hyperbolic sine of a number
#' 
#' @param x a double
#' @return a double
#' @export

ihs <- function(x){log(x + sqrt(x ^ 2 + 1))}

#' Calculates the weighted variance of a vector
#' 
#' @param x a vector
#' @param w a vector of weights same length as x
#' @return the variance
#' @export

weighted.var <- function(x, w){
  nas <- !is.na(x) & !is.na(w)
  x <- x[nas]
  w <- w[nas]
  xm <- weighted.mean(x, w)
  wv <- weighted.mean((x-xm)^2)
  return(wv)
}

#' Calculates robust standard errors from a regression model
#' 
#' @param x a model from lm or similar
#' @param c the clustering variable
#' @param setype standard error types
#' @return the variance
#' @export

robustses <- function(x, c, setype = "HC1"){
  r <- sqrt(diag(vcovPL(x, type = setype, cluster = ~ c)))
  return(r)
}

#' Creates a summary statistics from a data frame with weights
#' 
#' @param df a data frame
#' @param wts vector of weights same length as the data frame
#' @param labs vector of variable names - defaults to column names of df
#' @return the variance
#' @export

sumstats.wtd <- function(df, wts, labs = NULL, d = 2){
  if(is.null(labs)) labs <- colnames(df)
  ns <- apply(df, 2, function (x) sum(!is.na(x)))
  means <- apply(df, 2, function(x) weighted.mean(x, wts, na.rm = TRUE))
  sds <- apply(df, 2, function(x) sqrt(weighted.var(x, wts)))
  mins <- apply(df, 2, min, na.rm = TRUE)
  q25s <- apply(df, 2, function(x) quantile(x, .25, na.rm = TRUE))
  q75s <- apply(df, 2, function(x) quantile(x, .75, na.rm = TRUE))
  maxs <- apply(df, 2, max, na.rm = TRUE)
  tab <- data.frame("Statistic" = labs,
                    "N" = formatC(ns, format="d", big.mark=","),
                    "Mean" = formatC(means, format="f", big.mark=",", digits = 2),
                    "Std.Dev." = formatC(sds, format="f", big.mark=",", digits = 2),
                    "Min" = formatC(mins, format="d", big.mark=","),
                    "Pctle(25)" = formatC(q25s, format="d", big.mark=","),
                    "Pctle(75)" = formatC(q75s, format="d", big.mark=","),
                    "Max" = formatC(maxs, format="d", big.mark=","))
  rownames(tab) <- NULL
  xtab <- xtable(tab)
  return(xtab)
}