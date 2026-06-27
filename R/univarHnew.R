#' @title Estimate univariate exposure-response function on new grid of points
#'
#' @param fit An object containing the results from `rffbkmr` function.
#' @param whichX Vector identifying which exposure should be selected to varying. The number of selected exposures should be one.
#' @param Xgrid Values at grid points to cover the range of selected exposure.
#' @param fixlevel Quantile that other exposures will be fixed at (default: 0.5).
#' @param center The quantile of the selected exposure at which the exposure–response function equals zero.
#' @param ngrid The number of grids points to cover the range of selected exposure (default: 50).
#' @param alpha 100(1-alpha)% posterior interval (default: 0.05).
#'
#' @return A data frame with selected exposure name, exposure value, posterior mean estimate, and 100(1-alpha)% posterior interval.
#' @export
#'
#' @importFrom abind abind
#' @import stats
#' @import utils
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
#' est_uni <- univarHnew(fit = fit, whichX = "x1")
univarHnew = function (fit, whichX, Xgrid = NULL, fixlevel = 0.5, center = NULL, ngrid = 50, alpha = 0.05) {

  if (length(whichX) != 1) {
    stop("The number of selected exposures should be one.")
  }

  if (fixlevel < 0 | fixlevel > 1) {
    stop("All other exposures should be held constant at quantile values between 0 and 1.")
  }

  K = fit$K
  Z = fit$Z
  X = fit$X

  X.fix = apply(X, 2, quantile, probs = fixlevel)

  if (is.null(Xgrid)) {
    X.new = matrix(rep(X.fix, each = ngrid), nrow = ngrid, byrow = F)
    colnames(X.new) = colnames(X)

    X.sel = unlist(as.vector(X[, whichX]))
    X.new[, whichX] = seq(min(X.sel), max(X.sel), length.out = ngrid)

  } else {
    X.new = matrix(rep(X.fix, each = length(Xgrid)), nrow = length(Xgrid), byrow = F)
    colnames(X.new) = colnames(X)

    X.sel = unlist(as.vector(X[, whichX]))
    X.new[, whichX] = Xgrid
  }

  if (!is.null(center)) {

    if (center < 0 | center > 1) {
      stop("The quantile for center should be between 0 and 1.")
    }

    tmp = X.new[nrow(X.new),]
    if(is.null(Xgrid)) {
      tmp[whichX] = quantile(X[, whichX], probs = center)
    } else {
      tmp[whichX] = Xgrid[round(center*length(Xgrid))]
    }
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
  H.lower = apply(H.CI, 1, FUN = quantile, probs = q)
  H.upper = apply(H.CI, 1, FUN = quantile, probs = 1-q)

  if (!is.null(center)) {
    ngrid = length(H.new) - 1
    tmp = H.new[1:ngrid] - H.new[(ngrid+1)]
    H.new = tmp
    vcov.matrix = cov(t(H.CI))
    tmp.var = diag(vcov.matrix)[1:ngrid] + vcov.matrix[(ngrid+1),(ngrid+1)] - 2*vcov.matrix[1:ngrid,(ngrid+1)]
    H.lower = H.new - 1.96*sqrt(tmp.var)
    H.upper = H.new + 1.96*sqrt(tmp.var)
  }


  if (is.null(center)) {
    output = data.frame(var = colnames(X)[whichX], value = X.new[, whichX],
                        H.est = H.new, low = H.lower, up = H.upper)
  } else {
    output = data.frame(var = colnames(X)[whichX], value = X.new[1:ngrid, whichX],
                        H.est = H.new, low = H.lower, up = H.upper)
  }


  return(output)

}
