# working environment ####

# cleaning R start
rm(list=ls()); gc()

# loading libraries
libraries = c( 
  # general purpose
  'here',
  # experimental design
  'FrF2','daewr','rsm'
) 
# sapply(libraries, install.packages, character.only=T)
sapply(libraries, require, character.only=T)

# sourcing used defined functions (udf)
source( file.path(here(), 'code', 'udf.R') )

seed = 12345


# Step 1: screening experiment ####

## Current operation conditions ####
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
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,twi[[1]]], 
      pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, mr[ twi[[1]][2] ] ),
      xlim=c( 0, mr[ twi[[1]][1] ] ), 
      ylab=names(d0)[ twi[[1]][2] ],
      xlab=nam )
text( d0[,twi[[1]]] + 0.5, labels =round(d0$Y,1) )



# plot
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,twi[[1]]], 
      pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, mr[ twi[[1]][2] ] ),
      xlim=c( 0, mr[ twi[[1]][1] ] ), 
      ylab=names(d0)[ twi[[1]][2] ],
      xlab=nam )
rect( 0, 0, 4, 20)
text( d0[,twi[[1]]] + 0.5, labels =round(d0$Y,1) )


# plot
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,twi[[1]]], 
      pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, 20 ),
      xlim=c( 0, 4 ), 
      ylab=names(d0)[ twi[[1]][2] ],
      xlab=nam )
text( d0[,twi[[1]]] + 0.15, labels =round(d0$Y,1) )



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
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

plot( d0[,twi[[1]]], 
      pch=19, col=rgb(0,0,0,0.5),
      ylim=c( 0, 20 ),
      xlim=c( 0, 4 ), 
      ylab=names(d0)[ twi[[1]][2] ],
      xlab=nam )
points( d[,1:2], pch=19, col=rgb(0,0,1,0.5) )




## 2. Measure ####
set.seed( seed )
d$Y = true_protein_yield( X = as.matrix( d[,1:k] ) )


## 3. Analyze ####

# model for all First Order (FO) effects
model1 = lm( Y ~ A + B + C + D + E + F + G, data=d)
summary(model1)

halfnorm( coef(model1) )

# reduced model for First Order (FO) effects
model2 = lm( Y ~ A + B, data=d)
summary(model2)



# Step 2: steepest ascend ####

# experiment block
d$block = 'screen'

# add some center points
for( i in 1:4 ){
  d = rbind( d, data.frame( d0[,-ncol(d0)], Y=NA, block='center') )
}

set.seed( seed )
d$Y[is.na(d$Y)] = true_protein_yield( X = as.matrix( d[is.na(d$Y),1:k] ) )


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
  new$B = new$B + with(sa, A/B)*(new$A-2)
  new$Y = true_protein_yield( X = as.matrix( new[,-ncol(new)] ) )
  d = rbind( d, data.frame(new, block='ascend') )
  check = with(d, Y[nrow(d)] - Y[nrow(d)-1] > 0 )
}


# plot
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 0, 20 ),
         xlim=c( 0, 4 ) )
points( d[ d$block!='ascend', twi[[1]] ], 
        pch=19, col=rgb(0,0,0,0.5) )

# equation of a line
point1 = d[which( d$block=='center')[1], twi[[1]] ] 
point2 = d[which( d$block=='ascend')[1], twi[[1]] ] 
b = ( point2[,2] - point1[,2] ) / ( point2[,1] - point1[,1] )
a = point1[,2] - b*point1[,1]
abline( a=a, b=b, lty=2 )


# plot
twi = list( c(1,2) )

par( mfrow=c(1,1) )
idx = !names(d0) %in% c( names(d0)[ twi[[1]] ], 'Y' )
nam = paste( 
  names(d0)[ twi[[1]][1] ],'\n',
  paste( names(d0)[idx], '=', d0[,idx], collapse=', ') )

contour( model_rsm, ~ A+B, image=T,
         ylim=c( 0, 100 ),
         xlim=c( 0, 14 ) )
points( d[, twi[[1]] ], 
        pch=19, col=rgb(0,0,0,0.5) )
text( d[ d$block=='ascend', twi[[1]] ] - 0.5, 
      labels = round( d$Y[ d$block=='ascend'],1) )

# equation of a line
abline( a=a, b=b, lty=2 )



# RSM ####



# BO ####




# final validation ####


## Unknown "true" model ####

# setting seed for plot replication
set.seed(1256)

# simulation
X = true_covariates( n=1000 )
Y = true_protein_yield( X=X )
d = data.frame( X, Y)
m = rsm( Y ~ FO(X1,X2) + TWI(X1,X2) + PQ(X1,X2), data=d )
# summary(m)

# plot
contour( m, ~ X1+X2, image=T,
         xlabs=c("pH", "Temperature (ºC)") )
abline( v=10, h=70, lty=2, col=rgb(0,0,0,0.3) )


