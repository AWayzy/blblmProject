#' @import purrr
#' @import furrr
#' @import stats
#' @import parallel
#' @import future
#' @import utils
#' @importFrom magrittr %>%
#' @aliases blblm-package
#' NULL
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))



#' @title Bag of Little Bootstraps Regression
#' @details
#' Uses Little Bag of Bootstraps to Compute Linear Regression with a specified number of subsamples and replications
#' @param formula Regression Formula
#' @param data Data Frame
#' @param m Integer (Number of Subsamples)
#' @param B Integer (Number of Replications)
#' @param parallel Boolean - If true, will use all available CPU cores
#' @param seed Optional seed, will allow for testing output.
#' @export
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE, seed = NaN) {
  if (parallel == FALSE) {
    if (!is.nan(seed)) {
      data_list <- split_data(data, m, seed)
    }
    else {
      data_list <- split_data(data, m)
    }
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, seed = seed)
    )
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
    return(res)
  }
  else if (parallel == TRUE) {
    # got this if / else chunk from https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions
    # R only allows 2 cores to be required for a package, and thus the test fails if I attempt to run with all 4.

    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      num_workers <- 2
    } else {
      num_workers <- parallel::detectCores()
    }
    data_list <- split_data(data, m)
    n <- nrow(data)
    cl <- parallel::makeCluster(num_workers)
    parallel::clusterExport(cl, c("data_list", "B", "n", "formula"), envir = environment())
    if (!is.nan(seed)) {
      parallel::clusterSetRNGStream(seed)
      estimates <- parallel::parLapply(cl, data_list, function(x) {
        lm_each_subsample(formula = formula, data = x, n = n, B = B, parallel = TRUE, seed = seed)
      })
      res <- list(estimates = estimates, formula = formula)
      class(res) <- "blblm"
      invisible(res)
      return(res)
    }
    estimates <- parallel::parLapply(cl, data_list, function(x) {
      lm_each_subsample(formula = formula, data = x, n = n, B = B, parallel = TRUE)
    })
    parallel::stopCluster(cl)
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}


#' split data into m parts of approximated equal sizes
#' @param data Data Frame
#' @param m Number of Desired Splits
#' @param seed Optional seed for testing
#' @return Returns a list of m subsections of data.
split_data <- function(data, m, seed = NaN) {
  if (!is.nan(seed)) {
    set.seed(seed)
  }

  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @param formula The regression formula to be applied
#' @param n Number of rows in the data
#' @param data Dataframe
#' @param B Number of replications to be used in the bootstrap
#' @param parallel Determines if the function will utilize multiple cores
#' @param seed optional seed for testing
lm_each_subsample <- function(formula, data, n, B, parallel = FALSE, seed = NaN) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  if (!is.nan(seed)) {
    replicate(B, lm1(X, y, n, seed = seed))
  }
  else {
    replicate(B, lm1(X, y, n), simplify = FALSE)
  }
}


#' compute the regression estimates for a blb dataset
#' @param X Model Matrix for Regression
#' @param y Model Response for Regression
#' @param n Number of Rows in Data
#' @param seed Optional seed for testing
lm1 <- function(X, y, n, seed = NaN) {
  if (!is.nan(seed)) {
    set.seed(seed)
  }
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit A Weighted Linear Model
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit A Weighted Linear Model
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @title Print Blblm
#' @param x blblm object
#' @param ... ...
#' @details
#' Called via the standard print() in R, it prints the formula used by the blblm object.
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", utils::capture.output(x$formula))
  cat("\n")
}


#' @export
#' @title Sigma
#' @details
#' Returns the calculated mean of all estimated sigmas from a given blblm fitted object. If confidence = TRUE,
#' returns a confidence interval at the specified level for the sigma value
#' @param object A blblm object
#' @param confidence Boolean, TRUE if confidence interval is desired, FALSE by default
#' @param level The desired level of confidence for the confidence interval (default .95)
#' @param ... ...
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @title Coef
#' @details
#' Returns the mean of all estimated coefficients for a given blblm fitted object.
#' @method coef blblm
#' @param object A blblm object
#' @param ... ...
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @title Confint blblm
#' @details
#' Returns a confidence interval for the given parameters (if no parameters are given, gives a confidence interval for all parameters in the model formula)
#' @method confint blblm
#' @param object A blblm object
#' @param parm The given parameters to be estimated
#' @param level The given level of confidence (default .95)
#' @param ... ...
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @title Predict
#' @details
#' Returns the predicted outputs for the values of the new data given. If confidence = TRUE, returns a confidence
#' interval for the predicted output value at the specified level of confidence.
#' @method predict blblm
#' @param object A blblm object
#' @param new_data A set of new data corresponding to the predictor variables of the formula (if formula is y ~ x, data is x = ...)
#' @param confidence Boolean, If TRUE will generate a confidence interval at the specified level, default FALSE
#' @param level The desired level of confidence for the confidence interval. Default .95
#' @param ... ...
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
