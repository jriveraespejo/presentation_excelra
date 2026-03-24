
#' Inverse logit function
#'
#' @param x vector
inv_logit = function(x){
  return( 1/(1+exp(-x)) )
}


#' Simulate covariates with maximum ranges
#'
#' @param n number of simulations
#' @param mr maximum ranges for covariates
#'
#' @returns matrix size [n, length(mr)] 
#' @export
#'
#' @examples
true_covariates = function( 
    mr = c(14, 100, 24*4, 1000, 30, 100, 100),
    n=10 ){
  
  # # test
  # mr = c(14, 100, 24*4, 1000, 30, 100, 100)
  # n=10
  
  # storage
  q = length(mr)
  X = matrix( NA, nrow=n, ncol=q)
  
  # simulation
  for( j in 1:q){
    X[,j] = round( runif(n=n, min=0, max=mr[j]), 1) 
  }
  
  # return object
  return(X)
}


#' Simulate true protein yields
#'
#' @param n number of simulations
#' @param X covariates of interest [ dim(x) = c(n,q) ]
#' @param mr maximum ranges for covariates [ length(mr) = q ]
#' @param cp center points [ length(cp) = q ]
#' @param twi list of specific two way interactions 
#' @param pq vector with pure quadratic terms
#' @param betas representing relationship between X and Y
#' @param sigma_e inherent variability
#'
#' @returns data.frame
#' @export
#'
#' @examples
true_protein_yield = function(
    n = 10,
    X = NULL,
    mr = c(14, 100, 24*4, 1000, 30, 100, 100),
    cp = c(12, 70, rep(0, 5) ),
    twi = list( c(1,2) ),
    pq = 1:2,
    betas = c(85, -0.4, -0.02, rep(0.01, 5), 0.1),
    sigma_e = 2 ){
  
  # # test
  # n = 10
  # X = d
  # mr = c(14, 100, 24*4, 1000, 30, 100, 100)
  # cp = c(7.2, 37, rep(0, 5) )
  # twi = list( c(1,2) )
  # pq = 1:2
  # betas = c(4.4, -0.4, -0.005, rep(0.001, 5), 0.1)
  # betas = c(85, -0.4, -0.02, rep(0.01, 5), 0.1)
  # sigma_e = 1

  # true covariates
  if( is.null(X) ){
    X = true_covariates( n=n, mr=mr )
  }
  
  # calculations
  # subtracting center points
  Xd = X - matrix( cp, nrow=nrow(X), ncol=ncol(X), byrow=T )
  
  # two way interactions
  for( k in 1:length(twi) ){
    dm = matrix( Xd[ ,twi[[k]] ], nrow=nrow(Xd), ncol=length(twi[[k]]) )
    Xd = cbind(Xd, apply( dm, 1, prod) ) # interaction  
  }
  
  # pure quadratic terms
  Xd[,pq] = Xd[,pq]^2 # quadratic terms
  
  # error
  e = rnorm( n=nrow(X), mean=0, sd=sigma_e )
  
  # outcome
  mu = betas[1] + Xd %*% betas[2:length(betas)] 
  Y = as.vector( mu + e )
   
  # return object
  return( Y )
}


#' RBF Kernel function
#'
#' @param X1 matrix of dimensions (n, d)
#' @param X2 matrix of dimensions (m, d)
#' @param l length scale
#' @param sigma_f scale parameter 
#'
#' @return matrix of dimensions (n, m)
rbf_kernel = function(X1, X2, l=1, sigma_f=1 ){
  
  # # test
  # X1=x_train[,-ncol(x_train)]
  # X2=x_train[,-ncol(x_train)]
  # l=1.0
  # sigma_f=1.0
  
  # calculating distances
  dist_sq = matrix( NA, nrow=nrow(X1), ncol=nrow(X2) )
  for( i in 1:nrow(X1) ){
    dist_sq[i,] = sapply( 
      1:nrow(X2), 
      function(j){ sum( (X1[i,] - X2[j,] )^2 ) } )
  }
  
  # return object
  return( sigma_f^2 * exp( -dist_sq /(2 * l^2) ) )
  
}


#' Posterior Gaussian Process
#'
#' @param X_predict matrix (m, d) of prediction points
#' @param X_train matrix (n, d) of training points
#' @param Y_train column vector (n, d) of training observations
#' @param l length scale
#' @param sigma_f scale parameter 
#' @param noise scalar of observation noise
#'
#' @return list of mean (mu) and covariance (sigma) for the Gaussian
posterior = function(X_train, X_predict, Y_train, l=1, sigma_f=1, noise=1e-8 ) {
  
  # # test
  # X_train = x_train
  # X_predict = x_predict
  # Y_train = y_train
  # l=1
  # sigma_f=1
  # noise=1e-8
  
  # covariance calculations
  Stt = rbf_kernel( X_train, X_train, l=l, sigma_f=sigma_f  ) # consider noise: + noise**2 * diag(dim(X_train)[[1]])
  Stp = rbf_kernel( X_train, X_predict, l=l, sigma_f=sigma_f )
  Spp = rbf_kernel( X_predict, X_predict, l=l, sigma_f=sigma_f ) # consider noise: + noise**2 * diag(dim(X_predict)[[1]])
  Stt_inv = solve(Stt)
  
  # mu and sigma
  w = t(Stp) %*% Stt_inv
  mu = as.vector( w %*% Y_train )
  sigma = Spp - ( t(Stp) %*% Stt_inv) %*% Stp
  
  # return object
  return( list( mu=mu, sigma=sigma ) )
  
}

#' Negative log-Likelihood of a Kernel
#'
#' @param X_train matrix (n, d) of training points
#' @param Y_train column vector (n, d) of training observations
#' @param noise scalar of observation noise
#'
#' @return function with kernel parameters as input and negative log likelihood
#' as output
negative_loglik = function( X_train, Y_train, noise=1e-8, ... ) {
  
  function(params) {
    n = dim(X_train)[[1]]
    K = rlang::exec( rbf_kernel, X1=X_train, X2=X_train, !!!params )
    L = chol( K ) # cholesky decomposition, consider noise: + noise**2 * diag(n)
    a = backsolve( r=L, x=forwardsolve( l=t(L), x=Y_train) )
    
    # return object
    return( 0.5*t(Y_train)%*%a + sum(log(diag(L))) + 0.5*n*log(2*pi) )
  }
}


#' Gaussian Process Regression
#'
#' @param kernel kernel function
#' @param X_train matrix (n, d) of training points
#' @param y_train column vector (n, d) of training observations
#' @param noise scalar of observation noise
#' @param ... parameters of the kernel function with initial guesses. Due to the
#' optimizer used, all parameters must be given and the order unfortunately
#' matters
#'
#' @return function that takes a matrix of prediction points as input and
#' returns the posterior predictive distribution for the output
gaussian_process_regression = function( 
    X_train, X_predict, Y_train, l=1, sigma_f=1, noise=1e-8 ) {
  
  # # test
  # X_train = x_train
  # X_predict = x_predict
  # Y_train = y_train
  # l=1.5
  # sigma_f=1
  # noise=1e-8
  
  # parameter optimization
  kernel_nll = negative_loglik( X_train, Y_train, l=l, sigma_f=sigma_f, noise=noise )
  param = c( l, sigma_f )
  opt = optim( par=rep(1, length(param)), fn=kernel_nll )
  opt_param = opt$par
  
  # posterior
  post = posterior(
    X_train = X_train,
    X_predict = X_predict,
    Y_train = Y_train,
    l=opt_param[1], 
    sigma_f=opt_param[2], 
    noise=noise)
  
  # return object
  return( list(
    mu=post$mu, sigma=diag(post$sigma),
    parameters = c(l=opt_param[1], sigma_f=opt_param[2] ) 
  ) )
  
}



#' Lower Confidence Bound (LCB)
#'
#' @param post posterior from a Gaussian process
#' @param lambda scalar, exploration/exploitation trade off
#'
#' @return EI, vector of length m
lower_confidence_bound = function( post, lambda=1 ) {
  
  # # test
  # post

  # return object
  return( post$mu  - lambda * post$sigma )
  
}



#' Expected Improvement (EI)
#'
#' @param post posterior from a Gaussian process
#' @param lambda scalar, exploration/exploitation trade off
#'
#' @return EI, vector of length m
expected_improvement = function( 
    X_train, X_predict, Y_train, l=1, sigma_f=1, noise=1e-8 ) {
  
  # # test
  # X_train = x_train
  # X_predict = x_predict
  # Y_train = y_train
  # l=1.5
  # sigma_f=1
  # noise=1e-8
  
  # calculating posteriors
  post_train = gaussian_process_regression( 
    X_train = X_train, 
    X_predict = X_train, 
    Y_train = Y_train, 
    l=l, sigma_f=l, noise=noise )
  
  post_predict = gaussian_process_regression( 
    X_train = X_train, 
    X_predict = X_predict, 
    Y_train = Y_train, 
    l=l, sigma_f=l, noise=noise )
  
  # expected improvement
  sigma = post_predict$sigma
  improvement = min(post_train$mu) - post_predict$mu
  Z = improvement / sigma
  ei = improvement * pnorm(Z) + sigma * dnorm(Z)
  ei[ sigma == 0.0 ] = 0.0
  
  # return object
  return( ei )
  
}