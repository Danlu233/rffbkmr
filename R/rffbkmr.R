
#' @title Approximate Bayesian Kernel Machine Regression via Random Fourier Features
#'
#' @param y A numeric vector for the outcome.
#' @param X A numeric matrix or data.frame for exposure variables.
#' @param Z A numeric matrix or data.frame for covariates (confounders).
#' @param Z.eq A formula specifying covariate transformation (default: ~ . -1).
#' @param K Number of MCMC iterations (default: 20000).
#' @param J Number of basis functions (default: 20).
#' @param a Shape parameter for inverse-gamma priors (default: 0.001).
#' @param b Rate parameter for inverse-gamma priors (default: 0.001).
#' @param s2_gamma Prior variance for fixed effects (default: 100).
#' @param theta Initial theta values (default: 0.5).
#' @param sig2 Initial sigma^2 (default: 1).
#' @param tau2 nitial tau^2 (default: 1).
#' @param verbose Logical; whether to print tuning info (default: TRUE).
#' @param omega_tune_rate Tuning rate for omega (default: 0.2).
#' @param beta_tune_rate Tuning rate for beta (default: 0.2).
#'
#' @return A list with posterior summaries and samples.
#' @export
#'
#' @import MASS
#' @import stats
#' @import utils
#' @examples
#' data(simudat)
#' fit <- rffbkmr(y = simudat$y,
#' X = simudat[,1:5],
#' Z = simudat[,6:10],
#' Z.eq = as.formula("~ z1 + as.factor(z2) + z3 + z4 + as.factor(z5)"),
#' K = 1000,
#' J = 10,
#' verbose = FALSE)
rffbkmr <- function(y,
                      X,
                      Z,
                      Z.eq = as.formula("~ .-1"),
                      K = 20000,
                      J = 20,
                      a = 0.001,
                      b = 0.001,
                      s2_gamma = 100, #prior for fixed effect coef
                      theta = rep(0.5, ncol(X)),
                      sig2 = 1,
                      tau2 = 1,
                      verbose = T,
                      omega_tune_rate = 0.2,
                      beta_tune_rate = 0.2) {

  start.time = Sys.time()


  X = as.matrix(X)

  #check if all exposures are continuous
  if (sum(apply(X, 2, is.numeric)) != ncol(X)){
    stop("All exposures should be continuous")
  }

  #categorical variable to dummy variable
  Z = as.data.frame(Z)
  Z = model.matrix(Z.eq, data = Z)

  ### initial setting

  M = ncol(X)  #Num of exposure
  P = ncol(Z)  #Num of confounders
  n = length(y) #Sample size

  omega = rnorm (J*M, rep(0,J*M), sqrt(rep(2*theta,each=J)))
  omega = matrix (omega, ncol = M)

  XXX = X %*% t(omega)
  phi_cos = cos (XXX)
  phi_sin = sin (XXX)

  B = cbind (Z, phi_cos, phi_sin)
  phi = cbind (phi_cos, phi_sin)

  d = c( rep(s2_gamma, P), rep(tau2/J, 2*J)) # this is the diagonal element in S

  BtB = t(B)%*%B

  S = diag (d)
  beta0 = solve(1/sig2 * BtB + solve(S)) %*% (1/sig2*t(B)%*%y)
  mu0 = B%*%beta0
  h0 =  phi %*% beta0[-c(1:P)]


  sig2 = sum((y-mu0)^2)/n


  beta = beta0

  ### derivative functions

  # HMC update 1
  U_beta = function (beta){   -1/(2*sig2)*sum( (y-B%*%beta)^2) - 1/2*sum (beta^2/d) }

  grad_U_beta = function (beta){  (1/sig2)*t(B)%*%(y-B%*%beta) - (beta/d) }

  # HMC update 2
  U_omega = function (omega_vec){
    omega = matrix (omega_vec, ncol = M)
    phi_cos = cos (X %*% t(omega))
    phi_sin = sin (X %*% t(omega))
    B = cbind (Z, phi_cos, phi_sin)
    -1/(2*sig2)*sum( (y-B%*%beta)^2) - 1/2*sum(diag(omega %*% diag(1 / theta) %*% t(omega)))
  }

  grad_U_omega = function (omega_vec){
    omega = matrix (omega_vec, ncol = M)
    phi_cos = cos (X %*% t(omega))
    phi_sin = sin (X %*% t(omega))
    B = cbind (Z, phi_cos, phi_sin)
    R = (y-B%*%beta)
    tmp = phi_sin %*% beta[(P+1):(P+J)] - phi_cos %*% beta[(P+J+1):(P+2*J)]
    grad = -1/(sig2) * t(matrix(rep(tmp, J), ncol = J, byrow = T)) %*% (matrix(rep(R,M), ncol=M) * X) - omega %*% diag(1/(2*theta))

    c(grad)
  }


  ## parameters to save
  beta.save = matrix (NA, ncol = K, nrow = 2*J+P)
  omega.save = array (NA, c(J,M,K))
  theta.save = matrix (NA, ncol = M, nrow = K)
  h.save = matrix (NA, ncol = K, nrow = n)
  sig2.save = tau2.save = rep (NA, K)
  log.pd = matrix(NA, ncol = n, nrow = K)

  # setting for leap frog
  epsilon_beta = epsilon_omega = 0.0002
  L = 5
  acc.check = 200

  leap_frog = function (theta, r, epsilon, U, grad_U){
    r_tilde = r + (epsilon/2)*grad_U(theta)
    theta_tilde = theta + epsilon*r_tilde
    r_tilde = r_tilde + (epsilon/2)*grad_U(theta_tilde)
    return (list (theta = theta_tilde, r = r_tilde))
  }

  HMC = function (U, grad_U, epsilon, L, theta)
  {
    r0 = rnorm(length(theta),0,1) # independent standard normal variates
    theta_tilde = theta
    r_tilde = r0

    for (i in 1:L){
      leap = leap_frog (theta_tilde, r_tilde, epsilon, U, grad_U)
      theta_tilde = leap$theta
      r_tilde = leap$r
    }

    # Evaluate potential and kinetic energies at start and end of trajectory
    curr_lik = U(theta)-0.5*sum(r0^2)
    prop_lik = U(theta_tilde)-0.5*sum(r_tilde^2)
    alpha = min(1,exp(prop_lik-curr_lik))
    if (runif(1) < alpha)
    {
      return (list(theta=theta_tilde, acc = TRUE)) # accept
    }
    else {
      return (list(theta=theta, acc = FALSE)) # reject
    }
  }

  counter_beta = counter_omega = 0
  for (k in 1:K){

    if (k %% 1000 == 0){print (k)}


    #Update beta
    hmc.k = HMC (U_beta, grad_U_beta, epsilon_beta, L, beta)
    counter_beta = counter_beta + hmc.k$acc
    beta=hmc.k$theta

    #Update omega
    hmc.k = HMC (U_omega, grad_U_omega, epsilon_omega, L, c(omega))
    counter_omega = counter_omega + hmc.k$acc
    omega=matrix(hmc.k$theta, ncol = M)
    phi_cos = cos (X %*% t(omega))
    phi_sin = sin (X %*% t(omega))
    B = cbind (Z, phi_cos, phi_sin)


    ##Tune during adaptive phase
    half.K = ceiling(K/2)
    if (k < half.K & k%%acc.check ==0){
      if (counter_beta/acc.check < 0.65){ epsilon_beta=epsilon_beta * (1 - beta_tune_rate)}
      if (counter_beta/acc.check > 0.85){ epsilon_beta=epsilon_beta * (1 + beta_tune_rate)}
      if (counter_omega/acc.check < 0.65){ epsilon_omega=epsilon_omega * (1 - omega_tune_rate)}
      if (counter_omega/acc.check > 0.85){ epsilon_omega=epsilon_omega * (1 + omega_tune_rate)}
      counter_beta = counter_omega = 0
      if (verbose == T) {
        print (c(epsilon_beta, epsilon_omega))
      }
    }

    ##Update Theta (no variable selection)
    theta = (1/rgamma (M, a+J/2, apply(omega^2,2,sum)/2 + b ) )/2



    #Update sigma2
    sig2 = 1/rgamma (1, a+n/2, sum( (y-B%*%beta)^2)/2 + b )

    #Update tau2
    tau2 = (1/rgamma (1, a+J, sum( (beta[-c(1:P)])^2)/2 + b ))*J

    d = c( rep(s2_gamma, P), rep(tau2/J, 2*J))

    beta.save[,k] = beta
    sig2.save[k] = sig2
    tau2.save[k] = tau2
    theta.save[k,] = theta
    omega.save[,,k] = omega
    h.save[,k] =  cbind (phi_cos, phi_sin) %*% beta[-c(1:P)]

    mu_k = h.save[,k] + Z %*% beta[(1:P)]
    log.pd[k,] = log(1/sqrt(2*pi*sig2)) - (y - mu_k)^2/(2*sig2)
  }

  beta.hat = rowMeans (beta.save[,half.K:K])
  names(beta.hat) = c(colnames(Z), paste0("a",seq(1,J)), paste0("b",seq(1,J)))

  ## computation for WAIC
  lppd = sum(log(apply(exp(log.pd), 2, mean)))
  pWAIC = sum(apply(log.pd, 2, var))
  WAIC = -2*(lppd - pWAIC)

  end.time = Sys.time()

  output = list(beta.hat = beta.hat, h.hat = rowMeans(h.save[,half.K:K]),
                sig2.hat = mean(sig2.save[half.K:K]), tau2.hat = mean(tau2.save[half.K:K]),
                theta.hat = colMeans(theta.save[half.K:K,]),
                omega.hat = rowMeans(omega.save[,,half.K:K], dims = 2),
                counter.beta = counter_beta/(half.K + acc.check), counter.omega = counter_omega/(half.K + acc.check),
                beta.save = beta.save, sig2.save = sig2.save, tau2.save = tau2.save,
                theta.save = theta.save, omega.save = omega.save, h.save = h.save,
                time = end.time - start.time, WAIC = WAIC,
                epsilon_beta = epsilon_beta, epsilon_omega = epsilon_omega,
                log.pd = log.pd, K = K,
                X = X, Z = Z)

  return(output)
}
