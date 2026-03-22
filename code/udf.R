
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
    betas = c(90, -0.4, -0.02, rep(0.01, 5), 0.1),
    sigma_e = 2 ){
  
  # # test
  # n = 10
  # X = d
  # mr = c(14, 100, 24*4, 1000, 30, 100, 100)
  # cp = c(7.2, 37, rep(0, 5) )
  # twi = list( c(1,2) )
  # pq = 1:2
  # betas = c(95, -0.6428, -0.0214, rep(0.01, 5), 0.0857),
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

