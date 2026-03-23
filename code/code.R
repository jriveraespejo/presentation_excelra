# working environment ####

# cleaning R start
rm(list=ls()); gc()

# loading libraries
libraries = c( 
  # general purpose
  'here','tidyverse','magrittr',
  # experimental design
  'FrF2','daewr','rsm'
) 
# sapply(libraries, install.packages, character.only=T)
sapply(libraries, require, character.only=T)

# sourcing used defined functions (udf)
source( file.path(here(), 'code', 'udf.R') )

seed = 12345


# Motivating example ####

# Current operation conditions
k = 7 # number of factors
d0 = data.frame( matrix( c(2, 10, 5, 100, 10, 20, 20), nrow=1, ncol=7) )
names(d0) = LETTERS[1:k]

set.seed( seed )
d0$Y = true_protein_yield( X=as.matrix(d0) )


# plot
twi = list( c(1,2), c(3,2),
            c(1,3), c(3,4) )
mr = c(14, 100, 24*4, 1000, 30, 100, 100)

par( mfrow=c(2,2) )
for( j in 1:length(twi) ){
  idx = !names(d0) %in% c( names(d0)[ twi[[j]] ], 'Y' )
  nam = paste( 
    names(d0)[ twi[[j]][1] ],'\n',
    paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )
  
  plot( d0[,twi[[j]]], 
        pch=19, col=rgb(0,0,0,0.5),
        ylim=c( 0, mr[ twi[[j]][2] ] ),
        xlim=c( 0, mr[ twi[[j]][1] ] ), 
        ylab=names(d0)[ twi[[j]][2] ],
        xlab=nam )
  text( d0[,twi[[j]]] + 1, labels =round(d0$Y,1) )
}
par( mfrow=c(1,1) )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,1:2], pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, mr[ 2 ] ), xlim=c( 0, mr[ 1 ] ), 
      ylab=names(d0)[ 2 ], xlab=nam )
text( d0[,1:2] + 0.5, labels =round(d0$Y,1) )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,1:2], pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, mr[ 2 ] ), xlim=c( 0, mr[ 1 ] ), 
      ylab=names(d0)[ 2 ], xlab=nam )
rect( 0, 0, 4, 20, lty=2 )
text( d0[,1:2] + 0.5, labels =round(d0$Y,1) )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,1:2], pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, 20 ), xlim=c( 0, 4 ), 
      ylab=names(d0)[ 2 ], xlab=nam )
text( d0[,1:2] + 0.15, labels =round(d0$Y,1) )



# Step 1: Identifying relevant factors ####

## 1. Design ####

# design CRFF
set.seed( seed )
dsg = FrF2( nfactors=k, resolution=3 )

# # design PB
# dsg = pb( nruns=2^{k-p}+4 )
# dsg = dsg[,1:k]

# checking design
colormap(dsg, mod=2)

# checking design
print('Aliases')
y = runif( nrow(dsg), 0, 1)
aliases( lm( y~(.)^2, data=dsg) )


# convert to data.frame
d = data.frame( as.matrix(dsg[,1:k]) ) # design
d = sapply(d, as.integer)
d = data.frame(d)

# create natural variables
nv = rbind( 
  d0[,-ncol(d0)], # central points
  c(3-1, 15-5, 6-4, 150-50, 15-5, 30-10, 30-10) # ranges to test
)
for( j in 1:k ){
  d[,j] = d[,j]*nv[2,j]/2 + nv[1,j] 
}


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,1:2], pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, 20 ), xlim=c( 0, 4 ), 
      ylab=names(d0)[ 2 ], xlab=nam )
points( d[,1:2], pch=19, col=rgb(0,0,1,0.5) )
legend( 'topleft', legend=c( 'CRFF design','current OC' ),
        fill=c(rgb(0,0,1,0.5), rgb(0,0,0,0.5)), bty='n' )




## 2. Measure ####
set.seed( seed )
d$Y = true_protein_yield( X = as.matrix( d[,1:k] ) )


## 3. Analyze ####

# model for all First Order (FO) effects
model1 = lm( Y ~ A + B + C + D + E + F + G, data=d)
summary(model1)

# model for all First Order (FO) effects
fullnormal( coef(model1), alpha=0.025 )

# reduced model for First Order (FO) effects
model2 = lm( Y ~ A + B + C, data=d)
summary(model2)

# reduced model for First Order (FO) effects
model3 = lm( Y ~ A + B, data=d)
summary(model3)



# Step 2: Steepest ascend ####

## 1. Design ####

# experiment block
d$block = 'screen'

# add some center points
for( i in 1:4 ){
  d = rbind( d, data.frame( d0[,-ncol(d0)], Y=NA, block='center') )
}


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d[,1:2], pch=19, col=rgb(0,0,1,0.5),
      ylim=c( 0, 20 ), xlim=c( 0, 4 ), 
      ylab=names(d0)[ 2 ], xlab=nam )
legend( 'topleft', legend=c( 'CRFF design +\ncenter points' ),
        fill=rgb(0,0,1,0.5), bty='n' )


## 2. Measure ####
set.seed( seed )
d$Y[is.na(d$Y)] = true_protein_yield( X = as.matrix( d[is.na(d$Y),1:k] ) )


## 3. Analyze ####

# model
model_rsm = rsm( Y ~ FO(A,B), data=d )
summary(model_rsm)

# direction
sa = steepest(model_rsm, dist=1, descent=F)


# new point
new = d0
check = T
set.seed( seed )
while( check ){
  new$A = new$A + 1
  if( new$A > mr[1] ){
    break
  }
  new$B = new$B + with(sa, A/B)*(new$A-2)
  new$Y = true_protein_yield( X = as.matrix( new[,-ncol(new)] ) )
  d = rbind( d, data.frame(new, block='ascend') )
  check = with(d, Y[nrow(d)] - Y[nrow(d)-1] > 0 )
}


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 0, 20 ),
         xlim=c( 0, 4 ) )
points( d[ d$block!='ascend', 1:2 ], 
        pch=19, col=rgb(0,0,0,0.5) )

# equation of a line
point1 = d[which( d$block=='center')[1], 1:2 ] 
point2 = d[which( d$block=='ascend')[1], 1:2 ] 
b = ( point2[,2] - point1[,2] ) / ( point2[,1] - point1[,1] )
a = point1[,2] - b*point1[,1]
abline( a=a, b=b, lty=2 )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 0, 100 ),
         xlim=c( 0, 14 ) )
points( d[, 1:2 ], 
        pch=19, col=rgb(0,0,0,0.5) )
text( d[ d$block=='ascend', 1:2 ] - 0.5, 
      labels = round( d$Y[ d$block=='ascend'],1) )

# equation of a line
abline( a=a, b=b, lty=2 )



# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 0, 100 ),
         xlim=c( 0, 14 ) )
points( d[, 1:2 ], 
        pch=19, col=rgb(0,0,0,0.5) )
text( d[ d$block=='ascend', 1:2 ] - 0.5, 
      labels = round( d$Y[ d$block=='ascend'],1) )

# equation of a line
abline( a=a, b=b, lty=2 )
rect( 11, 60, 14, 100, lty=2 )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 60, 100 ),
         xlim=c( 11, 14 ) )
points( d[, 1:2 ], 
        pch=19, col=rgb(0,0,0,0.5) )
text( d[ d$block=='ascend', 1:2 ] - 0.1, 
      labels = round( d$Y[ d$block=='ascend'],1) )





# Step 3: Bayesian optimization ####

## 1. Design ####

# initial training points
dBO = tail(d, 3)
x_train = dBO[, c('A','B') ]
y_train = dBO$Y

# prediction points
A = seq(11,14,by=0.2)
B = seq(60,100,by=5)
x_predict = data.frame( expand.grid( A=A, B=B ) )


# plot
par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ 1:2 ], 'Y' )
nam = paste( 
  names(d0)[ 1 ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d[, 1:2 ], 
        pch=19, col=rgb(0,0,0,0.5),
        ylim=c( 59, 115 ),
        xlim=c( 10.95, 14.05 ) )
text( d[ d$block=='ascend', 1:2 ] - 0.1, 
      labels = round( d$Y[ d$block=='ascend'],1) )
points( x_predict, pch=19, col=rgb(0,0,1,0.5) )
legend( 'topleft', legend=c( 'Training data', 'Candidate points' ),
        fill=c( rgb(0,0,0,0.5), rgb(0,0,1,0.5) ), bty='n' )

d


## 2. Measure and Analyze ####

# analyze
post = gaussian_process_regression( 
  X_train = x_train, 
  X_predict = x_predict, 
  Y_train = y_train, 
  l=1, sigma_f=1, noise=1e-8 )

# plot
plot( d[, 1:2 ], pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 59, 115 ), xlim=c( 10.95, 14.05 ) )
text( d[ d$block=='ascend', 1:2 ] - 0.1, 
      labels = round( d$Y[ d$block=='ascend'],1) )
points( x_predict, pch=19, col=rgb(0,0,1,0.5) )
legend( 'topleft', legend=c( 'Training data', 'Candidate points' ),
        fill=c( rgb(0,0,0,0.5), rgb(0,0,1,0.5) ), bty='n' )

mu = matrix( post$mu, nrow=length(A), ncol=length(B), byrow=T)
contour( x=A, y=B, z=mu, add=T, col=hcl.colors(10, "Spectral") )



## 3. Analyze ####







# The reality ####


# setting seed for plot replication
set.seed( seed )

# simulation
X = true_covariates( n=1000 )
Y = true_protein_yield( X=X )
dT = data.frame( X, Y)
m = rsm( Y ~ FO(X1,X2) + TWI(X1,X2) + PQ(X1,X2), data=dT )


# plot
contour( m, ~ X1+X2, image=T,
         xlabs=c("pH", "Temperature (ºC)") )
abline( v=12, h=70, lty=2, col=rgb(0,0,0,0.3) )


# plot
contour( m, ~ X1+X2, image=T,
         xlabs=c("pH", "Temperature (ºC)") )
abline( v=12, h=70, lty=2, col=rgb(0,0,0,0.3) )
points( d[ d$block!='ascend', 1:2 ], pch=19, col=rgb(0,0,0,0.5),
        ylim=c( 59, 115 ), xlim=c( 10.95, 14.05 ) )
points( d[ d$block=='ascend', 1:2 ], pch=19, col=rgb(0,0,1,0.5),
      ylim=c( 59, 115 ), xlim=c( 10.95, 14.05 ) )
# points( d[ d$block=='BO', 1:2 ], pch=19, col=rgb(1,0,0,0.5),
#         ylim=c( 59, 115 ), xlim=c( 10.95, 14.05 ) )
legend( 'left', legend=c( 'Screen', 'Ascend', 'Optimize' ),
        fill=c( rgb(0,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,0,0.5) ), bty='n' )

