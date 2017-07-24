#' @export
rkan <- function(x, y, nlambda1 = 20, nlambda2 = 20, nk=p,lambda1 = NULL, lambda2 = NULL, k=NULL, g=p,
                 lambda1.min=0.05, lambda2.min=0.001, beta0 = NULL, w0 = NULL, 
                 initial = c("uniform","rkan"), intercept = TRUE, standardize = TRUE){
  n = length(y)
  p = dim(x)[2]
  ## check error
  if (class(x) != "matrix") {
    tmp <- try(x <- as.matrix(x), silent = TRUE)
    if (class(tmp)[1] == "try-error") 
      stop("x must be a matrix or able to be coerced to a matrix")
  }
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent = TRUE)
    if (class(tmp)[1] == "try-error") 
      stop("y must numeric or able to be coerced to numeric")
  }
  if (any(is.na(y)) | any(is.na(x))) 
    stop("Missing data (NA's) detected.Take actions to eliminate missing data before passing 
         X and y to pawls.")
  initial <- match.arg(initial)
  
  if (!is.null(lambda1)) 
    nlambda1 <- length(lambda1)
  if (!is.null(lambda2)) 
    nlambda2 <- length(lambda2)
  if (!is.null(k))
    nk <- length(k)

  ## set initial
  if (initial == "rkan") {
    init = rkan(x, y)
    beta0 = ifelse(abs(init$beta)<1e-7, 0.001, init$beta)
    w0 = ifelse(init$w == 1, 0.999, init$w)
  } else if (initial == "uniform") {
    if (is.null(beta0)) {
      beta0 = rep(1, p)
    }
    if (is.null(w0)){
      w0 = rep(0, n)
    }
  }
  
  ## sandardize
  if (standardize) {
    std <- .Call("Standardize", x, y)
    XX <- std[[1]]
    yy <- std[[2]]
    scale <- std[[3]]
  } else {
    XX <- x
    yy <- y
  }
  
  ## set tunning parameter
  if (is.null(lambda1)||is.null(lambda2)||is.null(k)) {
    param <- set_parameter(x=XX, y=yy, nlambda1=nlambda1, nlambda2=nlambda2, nk=nk,
                             lambda1.min=lambda1.min, lambda2.min=lambda2.min, beta0=beta0, w0=w0)
    if (is.null(lambda1)) 
      lambda1 <- param$lambda1
    if (is.null(lambda2)) 
      lambda2 <- param$lambda2
    if (is.null(k))
      k <- param$k
  }
  
  ## Fit 
  res1 <- rkan_grid(x=XX, y=yy, lambda1=lambda1, lambda2=lambda2, k=k, g=g, beta0=beta0, w0=w0)
  res2 <- BIC_grid(res1$wloss, res1$beta, res1$w,k=k, g=g)
  fit <- list(beta = res2$beta,
              w = res2$w,
              lambda1 = lambda1,
              lambda2 = lambda2,
              k = k,
              opt.lambda1 = lambda1[res2$index[1]],
              opt.lambda2 = lambda2[res2$index[2]],
              opt.k =k[res2$index[3]],
              iter = res1$iter,
              ws = res1$w,
              betas = res1$betas,
              raw.bic = res2$raw.bic,
              bic = res2$bic,
              how=res1$how,
              counter=res1$counter)
  
  ## unstandardize
  if (standardize) {
    scale = ifelse(scale == 0, 0, 1/scale)
    fit$beta = fit$beta * scale
  }
  fit
}