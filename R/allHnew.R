#' @title Estimate marginal exposure-response function on new grid of points
#'
#' @param fit An object containing the results from `rffbkmr` function.
#' @param Xrange Vector of quantiles for new grid points.
#' @param center The quantile of all exposures at which the exposure–response function equals zero.
#' @param alpha 100(1-alpha)% posterior interval (default: 0.05).
#'
#' @importFrom abind abind
#' @import stats
#' @import utils
#'
#' @return A data frame with all exposure values, posterior mean estimate, and 100(1-alpha)% posterior interval.
#' @export
#'
#' @examples
#' data(simudat)
#' fit <- rffbkmr(y = simudat$y,
#' X = simudat[,1:5],
#' Z = simudat[,6:10],
#' Z.eq = as.formula("~ z1 + as.factor(z2) + z3 + z4 + as.factor(z5)"),
#' K = 1000,
#' J = 10,
#' verbose = FALSE)
#' est_all <- allHnew(fit = fit)
allHnew = function (fit, Xrange = seq(0.01, 0.99, 0.01), center = NULL, alpha = 0.05) {

  if (sum(Xrange > 1 | Xrange < 0) != 0) {
    stop("All quantiles for new grid points should be between 0 and 1.")
  }

  K = fit$K
  Z = fit$Z
  X = fit$X

  X.new = apply(X, 2, quantile, probs = Xrange)

  if (!is.null(center)) {

    if (center < 0 | center > 1) {
      stop("The quantile for center should be between 0 and 1.")
    }

    tmp = apply(X, 2, quantile, probs = center)
    X.new = rbind(X.new, tmp)
  }

  phi.cos = cos (X.new %*% t(fit$omega.hat))
  phi.sin = sin (X.new %*% t(fit$omega.hat))

  H.new =  as.vector(cbind (phi.cos, phi.sin) %*% fit$beta.hat[-c(1:ncol(Z))])


  mutiply1 = function(omega) {
    result = X.new %*% t(omega)
    return(result)
  }

  X.omega.CI = apply(fit$omega.save[,,(K/2+1):K], 3, mutiply1, simplify = F)
  X.omega.CI = array(unlist(X.omega.CI), dim = (c(nrow(X.omega.CI[[1]]), ncol(X.omega.CI[[1]]), length(X.omega.CI))))

  phi.cos.CI = cos(X.omega.CI)
  phi.sin.CI = sin(X.omega.CI)

  phi = abind(phi.cos.CI, phi.sin.CI, along = 2)

  H.CI = matrix(rep(NA,K/2*length(H.new)), ncol = K/2)
  for (i in 1:(K/2)) {
    temp = as.vector(phi[,,i] %*% fit$beta.save[-c(1:ncol(Z)), K/2+i])
    H.CI[,i] = temp
  }
  H.CI = as.data.frame(H.CI)

  q = alpha/2

  if (!is.null(center)) {
    ngrid = length(H.new) - 1
    tmp = H.new[1:ngrid] - H.new[(ngrid+1)]
    H.new = tmp
    vcov.matrix = cov(t(H.CI))
    tmp.var = diag(vcov.matrix)[1:ngrid] + vcov.matrix[(ngrid+1),(ngrid+1)] - 2*vcov.matrix[1:ngrid,(ngrid+1)]
    H.lower = H.new + qnorm(q)*sqrt(tmp.var)
    H.upper = H.new + qnorm(1-q)*sqrt(tmp.var)
    output = as.data.frame(X.new[1:ngrid,])
  } else {
    H.lower = apply(H.CI, 1, FUN = quantile, probs = q)
    H.upper = apply(H.CI, 1, FUN = quantile, probs = 1-q)
    output = as.data.frame(X.new)
  }


  output$H.est = H.new
  output$low = H.lower
  output$up = H.upper
  output$quantile = Xrange

  return(output)

}
