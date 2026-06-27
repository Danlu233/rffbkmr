#' @title Estimate bivariate exposure-response function on new grid of points
#'
#' @param fit An object containing the results from `rffbkmr` function.
#' @param whichX Vector identifying which exposure should be selected to varying. The number of selected exposures should be two.
#' @param fixlevel Quantile that other exposures will be fixed at (default: 0.5).
#' @param center The quantile of the selected exposure at which the exposure–response function equals zero.
#' @param ngrid The number of grids points to cover the range of each selected exposures (default: 50).
#' @param alpha 100(1-alpha)% posterior interval (default: 0.05).
#'
#' @return A data frame with selected two exposure names, two exposure values, posterior mean estimate, and 100(1-alpha)% posterior interval.
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
#' est_bi <- bivarHnew(fit = fit, whichX = c("x1", "x2"))
bivarHnew = function (fit, whichX, fixlevel = 0.5, center = NULL, ngrid = 50, alpha = 0.05) {


  if (length(whichX) != 2) {
    stop("The number of selected exposures should be two")
  }

  if (fixlevel < 0 | fixlevel > 1) {
    stop("All other exposures should be held constant at quantile values between 0 and 1.")
  }

  K = fit$K
  Z = fit$Z
  X = fit$X

  X.fix = apply(X, 2, quantile, probs = fixlevel)


  X1 = X[,whichX[1]]
  X2 = X[,whichX[2]]

  X1.grid = seq(quantile(X1, 0.01), quantile(X1, 0.99), length.out = ngrid)
  X2.grid = seq(quantile(X2, 0.01), quantile(X2, 0.99), length.out = ngrid)

  grid = expand.grid(X1.grid, X2.grid)
  colnames(grid) = whichX

  if(ncol(X) == 2) {

    X.new = grid

  } else {

    med = apply(X, 2, median)

    X.median = matrix(rep(med, each = nrow(grid)), nrow = nrow(grid), byrow = F)
    colnames(X.median) = colnames(X)
    X.median = as.data.frame(X.median)[,!colnames(X.median) %in% whichX]

    X.new = cbind(grid, X.median)

  }

  X.new = as.matrix(X.new[,colnames(X)])

  if (!is.null(center)) {

    if (center < 0 | center > 1) {
      stop("The quantile for center should be between 0 and 1.")
    }

    tmp = X.new[nrow(X.new),]
    tmp[whichX] = apply(X[, whichX], 2, quantile, probs = center)

    X.new = rbind(X.new, tmp)
  }



  phi.cos = cos (as.matrix(X.new) %*% t(fit$omega.hat))
  phi.sin = sin (as.matrix(X.new) %*% t(fit$omega.hat))

  H.new =  as.vector(cbind (phi.cos, phi.sin) %*% fit$beta.hat[-c(1:ncol(Z))])

  ## CI
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
    tmp = H.new[1:(ngrid^2)] - H.new[(ngrid^2+1)]
    H.new = tmp
    vcov.matrix = cov(t(H.CI))
    tmp.var = diag(vcov.matrix)[1:(ngrid^2)] + vcov.matrix[(ngrid^2+1),(ngrid^2+1)] - 2*vcov.matrix[1:(ngrid^2),(ngrid^2+1)]
    H.lower = H.new - 1.96*sqrt(tmp.var)
    H.upper = H.new + 1.96*sqrt(tmp.var)
  }

  if (is.null(center)) {
    output = data.frame(var1 = colnames(X)[whichX[1]], var2 = colnames(X)[whichX[2]],
                        value1 = X.new[, whichX[1]], value2 = X.new[, whichX[2]],
                        H.est = H.new, low = H.lower, up = H.upper)

  } else {
    output = data.frame(var1 = colnames(X)[whichX[1]], var2 = colnames(X)[whichX[2]],
                        value1 = X.new[1:(ngrid^2), whichX[1]], value2 = X.new[1:(ngrid^2), whichX[2]],
                        H.est = H.new, low = H.lower, up = H.upper)

  }



  return(output)

}
