# R code for
#### The pseudo-stochastic EM (PSEM) algorithm for Sub-Gaussian alpha-Stable Mixture Model #######
#
# By: Shaho Zarei and Adel Mohammadpour
#
# PSEMfun(sgdata , K , N , TEK )
#
# Inputs:
#
# sgdata=  Sub-Gaussian heterogeneous data 
# K     =  Number of mixture components
# N     =  Monte Carlo sample size
# TEK   =  Number of iterations of the PSEM algorithm
# Outputs
#          
# predicted.classes= Labels of input sgdata
# mixture.weights  = Cluster weights
# The following outputs are estimated for each components
# alpha            = Tail index 
# shape.matrices   = Dispersion matrix 
# location.vectors = Location vectors 
# BIC              = Bayesian Information Criterion
#
###############################################################################################################################
# Example
###########
#
# Simulate data to test PSEMfun
#
  K <- 3                       # Number of mixture components
  d <- 2                       # Number of variables
  n <- 1500                    # Sample size
  prior   <- rep(1/K, K)       # Mixture weights
  library(stable)
#############-------values of  parameters for generated data-------################
#
alphax <- c(1.45, 1.25, 1.8)                                 # Tail Indices 
deltax <- matrix(c(0,5 ,5,0, -5,0),d,K)                      # Location parameters
sigmax <- 0.5*array(c(matrix(c(2,0.5,0.5,0.5),d,d),diag(d),
               matrix(c(2,-0.5,-0.5,0.5),d,d)), c(d,d,K))    # Dispersion Matrices
#
# Note: Sigma should multiply by 0.5 because stable package is based on the elliptical distributions and dispersion matrix
# of sub-Gaussian alpha stable distributions is half of dispersion matrix of the elliptical distributions.
#
realclass <- sample(1:K, size = n, replace = TRUE, prob = prior) # Determine sample size (component size) for  each component  
dataX <- NULL
for (i in 1:n){

  g <- realclass[i]
  x <- rmvstable( mvstable.elliptical(alphax[g], sigmax[,,g], deltax[,g]), 1) 
   dataX  <- cbind(dataX,x)

  }
 dim(dataX)
 sgdata <- t(dataX)  # Columns should show covariates

# Now we can test the function
 
PSEMfun(sgdata, K = 3, N = 2000, TEK = 30)
################
# End of Example
###############################################################################################################################

PSEMfun  <-  function(sgdata, K, N, TEK)

  { 

require(stable)       # The stable package available at http://www.robustanalysis.com
require(cluster)


n       <- nrow(sgdata)      # sample size
d       <- ncol(sgdata)      # Number of covariates

dataSIM <- data.frame(sgdata) 
prior   <- rep(1/K, K)       ## Mixture weights

########----- function for generating positive alpha-stable distribution-----################

rstab <- function (n, alpha, beta, gamma=1, delta=0, pm)
 
{
  stopifnot(0 < alpha, alpha<=  2, length(alpha)==1, -1 <= 
              beta, beta <=  1, length(beta)==1, 0 <=  gamma, length(pm)==
              1, pm %in% 0:2)
  if (pm  ==  1) {
    delta <- delta + beta * gamma * .om(gamma, alpha)
  }
  else if (pm==2) {
    delta <- delta - alpha^(-1/alpha) * gamma * stableMode(alpha,
                                                           beta)
    gamma <- alpha^(-1/alpha) * gamma
  }
  theta <- pi * (runif(n) - 1/2)
  w <- -log(runif(n))
  result <- if (alpha ==  1 & beta==0) {
    rcauchy(n)
  }
  else {
    b.tan.pa <- beta * tan(pi/2 * alpha)
    theta0 <- min(max(-pi/2, atan(b.tan.pa)/alpha), pi/2)
    c <- (1 + b.tan.pa^2)^(1/(2 * alpha))
    a.tht <- alpha* (theta + theta0)
    r <- (c * sin(a.tht)/(cos(theta))^(1/alpha)) * (cos(theta -
                                                          a.tht)/w)^((1 - alpha)/alpha)
    r - b.tan.pa
  }
  result * gamma + delta
}

.om <- function(gamma,alpha) {
  if(alpha!=  round(alpha))
    tan(pi/2*alpha)
  else if(alpha ==  1)
    (2/pi)*log(gamma)
  else 0 }

###################################################################################
 
  WEIGHT <- NULL
  ALPHAX <- NULL
  SIGMAX <- NULL
  DELTAX <- NULL
    Diff <- c()
     
    ##-------------------------------------------------------------###
       ### Determining initial values for the PSEM algorithm ###
    ##-------------------------------------------------------------###

kmedoid   <- pam(dataSIM ,k=K)
DATA      <- data.frame(cbind(dataSIM, kmedoid$clustering))

alphax.in <- c()                                     ## Tail Index
sigmax.in <- array(0,c(d,d,K))                       ## Dispersion Matrix
deltax.in <- matrix(0,d,K)                           ## Location parameter

pri <- table(kmedoid$clustering)

prior.in<- NULL
for ( h in 1:K){
prior.in <- cbind(prior.in, pri[h]/n)
}

for ( g in 1:K){

  DAg  <- DATA[DATA$kmedoid.clustering==g,]
  DAgM <- matrix(unlist(DAg), nrow(DAg), (d+1))
  fitM <- mvstable.fit.elliptical(t(DAgM[,c(1,2)]),1)
  alphax.in[g]   <- fitM$alpha
  deltax.in[,g]  <- fitM$delta
  sigmax.in[,,g] <- fitM$R

}
    
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
    ##-------------------Start of Loop-------------------------------###
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

  for (iter in 1:TEK){

##################------ Update mixture weights ----################

  comp1     <- lapply(1:K, function(j)prior.in[j]*
  dmvstable( mvstable.elliptical(alphax.in[j], sigmax.in[,,j], deltax.in[,j]), dataX))
  comp      <- sapply(comp1,cbind)
  compsum   <- apply(comp,1,sum)
  E1ijT     <- comp/compsum
  new.prior <- apply(E1ijT,2,mean)      # new mixture weights

############------ value of log-likelihood function----##############

obsloglik <- sum(log(compsum))        

#################--------- E[1/P|X] ---------#################

expectx <- function( DATAA, alfa, dispersion, loc){
 
  EijT  <- matrix(0,n,K)
 
  for (j in 1:K){
    k<-round(min(160,160/alfa[j]))        # for avoid computational error
    a0 <- alfa[j]
    bb <- seq(1,k)
    vv<-2+(exp(lgamma(a0*k/2+a0/2+1)+lgamma(a0*k/2+a0/2+d/2+1)-lgamma(a0*k/2+1)-
                 lgamma(a0*k/2+d/2+1))/((k+1)))^(2/a0)
   
    for (i in 1:n){
     
      dis <- t( DATAA[,i]-loc[,j])%*%solve(dispersion[,,j])%*%( DATAA[,i]-loc[,j])

      if (dis>vv){
        s1<-sum((-1)^(bb-1)*(dis/2)^(-a0*bb/2-d/2-1)*exp(lgamma(a0*bb/2+1)+lgamma(a0*bb/2+d/2+1)
                                                         -lgamma(bb+1))*sin(bb*pi*a0/2))
        s2<-sum((-1)^(bb-1)*(dis/2)^(-a0*bb/2-d/2)*exp(lgamma(a0*bb/2+1)+lgamma(a0*bb/2+d/2)
                                                       -lgamma(bb+1))*sin(bb*pi*a0/2))
        EijT[i,j] <- s1/s2}else{
         
          tu <- rstab(N,alfa[j]/2,1,(cos(pi*alfa[j]/4))^(2/alfa[j]),0,1)
          EijT[i,j]<-sum(na.omit(tu^(-d/2-1)*exp(-.5*dis/tu)))/
            sum(na.omit(tu^(-d/2)*exp(-.5*dis/tu)))
        }
    }
  }
  return(list(EijT = EijT))
}

ExijT <- expectx(dataX,alphax.in,sigmax.in,deltax.in)$EijT

##################------ Update location parameters ---################

dmujt1 <- lapply(1:K, function(j)t(E1ijT)[j,]%*%ExijT[,j])
dmujt  <- sapply(dmujt1,cbind)

new.deltax <- NULL
    for(j in 1:K){

        nm1   <- lapply(1:n, function(i) E1ijT[i,j]*ExijT[i,j]*dataX[,i])
        nm    <- sapply(nm1,cbind)
        nmsum <- apply(nm,1,sum)
        muj   <- nmsum/dmujt[j]
        new.deltax <- cbind(new.deltax,muj)

            }

############----- Update dispersion matrices ---------################## 

new.sigma <- array(0,dim = c(d,d,K))
num       <- apply(E1ijT,2,sum)

for(j in 1:K){

	sig1 <- lapply(1:n ,function(i) (1/num[j])*E1ijT[i,j]*ExijT[i,j]*(dataX[,i]-new.deltax[,j])
		%*%t(dataX[,i]-new.deltax[,j]))
	sig2 <- sapply(sig1,cbind)
	sig3 <- apply(sig2,1,sum)
	new.sigma[,,j] <- matrix(sig3,d,d)

            }
 
##################------------  Update tail indices  ------------####################  

alphaSEM1<-c()
calyij<-NULL
calyij<-array(0,dim = c(n,d,K))

for (j in 1:K){

calyij[,,j] <- dataX - new.deltax[,j]

}

DenPSatable<-matrix(0,n,K)

for(j in 1:K){
i=1
Vij=c()
sdstable <- c()
nj=0
t0<-alphax.in[j]
stany<-calyij[,,j]
sdstable <- dmvstable(mvstable.elliptical(alpha=t0, R=new.sigma[,,j], new.deltax[,j] ),
              t(dataSIM))

#### simulate data from conditional distribution  of (P|given x) with reject sampling ####
	
	  while(nj < n){
		vij <- rstable(1,  t0/2, beta=1, gamma=(cos(t0*pi*0.25))^(2/t0), delta=0, param=1)
            aaa=t(stany[i,])%*%solve(new.sigma[,,j])%*%stany[i,]
            maxf.y.p <- {exp(-d/2)*aaa^(-d/2)*d^(0.5*d)} / ((2*pi)^(d/2)*sqrt(det(new.sigma[,,j])))
            bond <- maxf.y.p/sdstable[i] 
		u= runif(1,0,1)
 			if(u < bond) 
                		     	{
     					Vij=cbind(Vij,vij)
            			nj=nj+1
                          i=i+1
               			}
                         }

DenPSatable[,j]<-Vij

fits=stable.fit.mle(DenPSatable[,j], param = 1)

alphaSEM1[j]<-2*fits[1]      ## the estimated value of tail index form f(p|alpha) should be multiplied by 2 for f(p|x,alpha)

}

alphaSEM <- alphaSEM1

#############-------------- save updates -------------############################

  WEIGHT <- cbind(WEIGHT, new.prior)
  ALPHAX <- cbind(ALPHAX, alphaSEM)
  SIGMAX <- cbind(SIGMAX, new.sigma)
  DELTAX <- cbind(DELTAX, new.deltax)

##############------------ Update Parameters for iterration (t+1)th----------############################

  prior.in  <- new.prior
  alphax.in <- alphaSEM
  sigmax.in <- new.sigma
  deltax.in <- new.deltax

comp2 <- lapply(1:K, function(j)prior.in[j]*
               dmvstable(mvstable.elliptical(alphax.in[j],sigmax.in[,,j],deltax.in[,j]),dataX ))

comp       <- sapply(comp2,cbind)
compsum    <- apply(comp,1,sum)
newloglik  <- sum(log(compsum))
Diff[iter] <- newloglik-obsloglik
print(iter)

                }    #### End of loop

############--- draw graph of between log-likelihood differences----##############

dplot <- plot(Diff, type = "l", xlab="iterations", ylab="Log-Likelihood differences")

###############################################################################

T0 = trunc(0.7*TEK)            ## T_{0}: An approximate value for burn-in time 

mw1 = lapply(T0:TEK, function(t)WEIGHT[,t])
mw2 = sapply(mw1,cbind)
mixweigth = apply(mw2,1,mean)

a1 = lapply(T0:TEK, function(t)ALPHAX[,t])
a2 = sapply(a1,cbind)
alphahatX = apply(a2,1,mean)

s1 = lapply(T0:TEK, function(t)SIGMAX[,t])
s2 = sapply(s1,cbind)
SIGMAX1 = apply(s2,1,mean)
sigmahatX = array(SIGMAX1,c(d,d,K))

deltahatX = matrix(0,d,K)
for (j in 1:K){
  mc = matrix(0,2,1)
  for (t in T0:(TEK-1)){
    nc = K*t+j
    mc = mc+DELTAX[,nc]
  }
  deltahatX[,j] = mc/(TEK-T0)
}

############------ pridected classes based on PSEM algorithm  ------########

prior.in1  <- mixweigth
alphax.in1 <- alphahatX
sigmax.in1 <- sigmahatX
deltax.in1 <- deltahatX

WG11<-NULL
for (j in 1:K){

  WG11<-cbind(WG11,dmvstable(mvstable.elliptical(alphax.in1[j],sigmax.in1[,,j],deltax.in1[,j],param=1),dataX))

}

compsum <- apply(WG11,1,sum)

E1ijT.hat <- WG11/compsum

SGaSMM.class <- c()   


for (i in 1:n){

  SGaSMM.class[i] <- which(E1ijT.hat[i,] == max(E1ijT.hat[i,]))   ##  determine estimated group of observations

              }
##########------ Compute Bayesian information criterion (BIC)------#############

comp3   <- lapply(1:K, function(j)mixweigth[j]*
           dmvstable(mvstable.elliptical(alphahatX[j],sigmahatX[,,j],deltahatX[,j]),dataX ))

compb    <- sapply(comp3,cbind)
compsumB <- apply(compb,1,sum)
Bloglik  <- sum(log(compsumB))
m        <- (K-1)+K+K*d+K*d*(d+1)/2
BIC      <- 2*Bloglik-m*log(n)         ## The highest BIC is preferred

##########------ The estimated values of parameters   ------#############

  mixweigth 
  alphahatX
  sigmahatX
  deltahatX

return (list(plot=dplot, predicted.classes=SGaSMM.class, mixture.weights=mixweigth, 
             alpha=alphahatX, shape.matrices=sigmahatX, 
             location.vectors=deltahatX,  BIC=BIC))

               }         

