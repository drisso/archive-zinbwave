#simulations

n=1000
J=20

#Start with known latent factor matrix:
U=cbind(runif(n)*seq(1:n)/100,runif(n)*seq(1:n)/100)
#simulate true U-coefficients for M
V=rbind(runif(J)*seq(1:20)/10,runif(J)*seq(1:20)/10)

#design matrix X.M
X.M=matrix(0,nrow=n,ncol=2)
X.M[,1]=c(rep(1,500),rep(0,500))
X.M[,2]=c(rep(0,500),rep(1,500))

#overall design matrix for M
X=cbind(X.M,U)

#coefficients for X.M
alpha.M=matrix(0,nrow=2,ncol=J)
alpha.M[1,]=10*runif(J)
alpha.M[2,]=15*runif(J)

#matrix M
M=exp(X.M%*%alpha.M+U%*%V)

#simulate true U-coefficients for Pi
W=-rbind(c(rep(1,10),rep(2,10)),c(rep(4,10),rep(3,10)))

#design matrix X.pi, take the same for the moment
X.pi=matrix(0,nrow=n,ncol=2)
X.pi[,1]=c(rep(1,500),rep(0,500))
X.pi[,2]=c(rep(0,500),rep(1,500))

#true coefficients for X.pi
alpha.pi=rbind(c(rep(5,10),rep(10,10)),c(rep(15,10),rep(8,10)))

#matrix of probabilities of dropout
Pi=binomial()$linkinv(X.pi%*%alpha.pi+U%*%W)
min(Pi)

#simulate negative binomial from M (matrix of expressions)
exprs=NULL
for(i in 1:(n*J)){
  exprs[i]=rnbinom(1,mu=as.vector(M)[i],size=log(4))
}
exprs=matrix(exprs,nrow=n)


indic=NULL
for(i in 1:(n*J)){
  indic[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
}
indic=matrix(indic,nrow=n)

Y=matrix(0,nrow=n,ncol=J)
Y[indic==T]=exprs[indic==T]
sum(Y==0)
which(apply(Y==0,2,sum)==max(apply(Y==0,2,sum)))
sum(Y[,17]==0)


#X_M design matrix for M, X_pi same for Pi
#U latent factor matrix

#matrix of covariates for M regression 
#(kx=nb of known M-factors + nb of latent factors)
X=cbind(X.M,U)

#matrix of covariates for Pi regression
#(kz=nb of known Pi-factors + nb of latent factors)
Z=cbind(X.pi,U)
#Z=X.pi

#for each gene, parms is a vector of size kx+kz+1 
#(total nb of known + nb of latent factors + theta)

offsetx=0 #no offset for the moment
offsetz=0 #no offset for the moment

#Y is matrix of counts n times J
n <- nrow(Y) #number of cells
kx <- NCOL(X) #number of M-factors
kz <- NCOL(Z) #number of Pi-factors
Y0 <- Y <= 0 #==1 if counts is 0
Y1 <- Y > 0 #==1 if count is not 0

linkstr <- "logit"
linkobj <- make.link(linkstr)
linkinv <- linkobj$linkinv

#loglikelihood
#parms should be a vector of parameters for gene j
#first all M-parameters, then all Pi-parameters, then theta

ziNegBin <- function(parms,j,X,Z,Y) {
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + 
                             offsetz))
  theta <- exp(parms[(kx + kz) + 1])
  loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                                                   size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y[,j], 
                                                     size = theta, mu = mu, log = TRUE))
  Y0 <- Y <= 0 #==1 if counts is 0
  Y1 <- Y > 0 #==1 if count is not 0
  loglik <- sum(loglik0[Y0[,j]]) + sum(loglik1[Y1[,j]])
  loglik
}

#gradient function
gradNegBin <- function(parms,j,X,Z,Y) {
  eta <- as.vector(X %*% parms[1:kx] + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
  muz <- linkinv(etaz)
  theta <- exp(parms[(kx + kz) + 1])
  Y0 <- Y <= 0 #==1 if counts is 0
  Y1 <- Y > 0 #==1 if count is not 0
  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- muz * (1 - as.numeric(Y1[,j])) + exp(log(1 - muz) + 
                                                  clogdens0)
  wres_count <- ifelse(Y1[,j], Y[,j] - mu * (Y[,j] + theta)/(mu + theta), 
                       -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) - 
                              log(mu + theta) + log(mu)))
  wres_zero <- ifelse(Y1[,j], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
  wres_theta <- theta * ifelse(Y1[,j], digamma(Y[,j] + theta) - 
                                 digamma(theta) + log(theta) - log(mu + theta) + 1 - 
                                 (Y[,j] + theta)/(mu + theta), exp(-log(dens0) + log(1 - 
                                                                                       muz) + clogdens0) * (log(theta) - log(mu + theta) + 
                                                                                                              1 - theta/(mu + theta)))

  colSums(cbind(wres_count * X, wres_zero * Z, wres_theta))
}


optim(fn=ziNegBin,gr=gradNegBin,j=2,X=X,Z=Z,Y=Y,par=c(1,1,0.5,0.1,1,1,1,1,1),control=list(fnscale=-1,trace=6),method="BFGS")
trueparam=rbind(alpha.M,V,alpha.pi,W,log(4))
trueparam[,2]
######################################################################################
#Optimization with respect to U
#log(M)=X.M%*%alpa.M+U%*%V, transpose the equation and consider X.M%*%alpha.M as known


# parms here are elements of i-th column of t(U) (i-th line of U) 
# all the other arguments are known
# i is the number of line of U to be optimized over
# V,W are from previous iteration
# offset.M = t(alpha.M)%*%t(X.M)
# offset.pi = t(alpha.pi)%*%t(X.pi)
# counts is the matrix of counts
# thetas is vector of thetas of all genes 

#loglik which will take U as parameter

ziNegBin.U <- function(parms,i,V,W,offset.M,offset.pi,counts,thetas){
  thetas=exp(thetas)
  is0counts <- counts == 0
  n0counts <- counts != 0
  mu <- as.vector(exp(t(V) %*% parms + offset.M[,i]))
  phi <- as.vector(linkinv(t(W) %*% parms + 
                             offset.pi[,i]))
  loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, 
                                          size = thetas, mu = mu, log = TRUE))))
  # counts[i,] is i-th column of t(counts)
  loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(counts[i,], 
                                                     size = thetas, mu = mu, log = TRUE))
  loglik <- sum(loglik0[is0counts[i,]]) + sum(loglik1[n0counts[i,]])
  loglik
}
# sum(sapply(seq(1000),function(i){ziNegBin.U(U[i,],i,V,W,t(alpha.M)%*%t(X.M),t(alpha.pi)%*%t(X.pi),counts,thetas)}))
# gradient function for optimization in U
# sum(sapply(seq(20),function(j){ziNegBin(trueparam[,j],j)}))
# gradient function for optimization in U

gradNegBin.U <- function(parms,i,V,W,offset.M,offset.pi,counts,thetas) {
  is0counts <- counts == 0
  n0counts <- counts != 0
  eta <- as.vector(t(V) %*% parms + offset.M[,i])
  mu <- exp(eta)
  etaz <- as.vector(t(W) %*% parms + offset.pi[,i])
  muz <- linkinv(etaz)
  thetas <- exp(thetas)
  clogdens0 <- dnbinom(0, size = thetas, mu = mu, log = TRUE)
  dens0 <- muz * (1 - as.numeric(n0counts[i,])) + exp(log(1 - muz) + 
                                                  clogdens0)
  wres_count <- ifelse(n0counts[i,],  counts[i,] - mu * (counts[i,] + thetas)/(mu + thetas), 
                     -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(thetas) - 
                            log(mu + thetas) + log(mu)))
wres_zero <- ifelse(n0counts[i,], -1/(1 - muz) * linkobj$mu.eta(etaz), 
                    (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
   colSums(wres_count * t(V) + wres_zero * t(W))
}

optim(fn=ziNegBin.U,gr=gradNegBin.U,i=448,par=c(0.7,0.5),counts=Y,
      thetas=rep(0.3266343,20),V=V,W=W,
      offset.M=t(alpha.M)%*%t(X.M),offset.pi=t(alpha.pi)%*%t(X.pi),
      control=list(fnscale=-1,trace=1),method="BFGS")

#Put everything together
#alt.number - number of alternations
#alt.number=10
alt.number=1
#initialization
p=2
U0=U+matrix(3,nrow=n,ncol=p)
V0=V+matrix(5,nrow=p,ncol=J)
alpha.M0=alpha.M*2
alpha.pi0=alpha.pi*2
W0=W*3
theta0=0.7
kxM=ncol(X.M)
kxV=p
kxPi=ncol(X.pi)
kxW=p
par.est=matrix(0,nrow=kxM+kxPi+2*p+1,ncol=J)
for (alt in 1:alt.number){
  for (gene in 1:J){
    par.est[,gene]=optim(fn=ziNegBin,Y=Y,j=gene,gr=gradNegBin,
          par=c(alpha.M0[,gene],V0[,gene],alpha.pi0[,gene],W0[,gene],theta0),
          X=cbind(X.M,U0),Z=cbind(X.pi,U0),
          control=list(fnscale=-1),method="BFGS")$par
  }
  alpha.M0=par.est[1:kxM,]
  V0=par.est[(kxM+1):(kxM+p),]
  alpha.pi0=par.est[(kxM+p+1):(kxM+p+kxPi),]
  W0=par.est[(kxM+p+kxPi+1):(kxM+2*p+kxPi),]
  thetas=par.est[kxM+2*p+kxPi+1,]
  for(cell in 1:n){
    U0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U0[cell,],counts=Y,
                    thetas=thetas,V=V0,W=W0,
                    offset.M=t(alpha.M0)%*%t(X.M),offset.pi=t(alpha.pi0)%*%t(X.pi),
                    control=list(fnscale=-1),method="BFGS")$par
  }
}























#test de optim
f=function(param,j){param[1]^2*2+j*param[2]^2}
f(c(1,1),3) #1^2*2+3*1=5 ok
gradf=function(param,j){c(4*param[1],2*j)}
gradf(c(1,1),3) #ok
#call optim with fixed j
optim(fn=f,gr=gradf,j=17,par=c(10,28),method="CG") #ok pour passer l'argument directement dans l'optim

