library(pscl)
source('~/GitHub/zinb/R/functions_svetlana.R')

# goal: simulate data K times with the same value of parameters (we take one gene only) 
# and check the variability of estimators

# sample size
n=100  #cell number
K=500  #number of repeats
J=100  #number of genes
    
# matrix to store the results of simulations
results=matrix(0,nrow=K,ncol=5)
X1=2*runif(n)
U1=3*(1+runif(n))
for (k in 1:K){

    X=cbind(X1,U1)
    
    alphaM=c(1.5,0.8)
    alphaPi=c(-0.3,-0.05)
    
    M=exp(X%*%alphaM)
    Pi=binomial()$linkinv(X%*%alphaPi)
    gene.exp=NULL
    for(i in 1:n){
        gene.exp[i]=rnbinom(1,mu=M[i],size=1)
    }
    gene.is1=NULL
    for(i in 1:n){
        gene.is1[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
    }
    gene.exp.zeroinf=gene.exp
    gene.exp.zeroinf[gene.is1==0]=0
    
    tab=cbind(gene.exp.zeroinf,X)
    colnames(tab)=c("expr","X1","X2")
    testdata=data.frame(tab)
    
#     tab2=cbind(gene.exp,X)
#     colnames(tab2)=c("expr","X1","X2")
#     testdata2=data.frame(tab2)
    
    estimate=zeroinfl(expr~X1+X2-1|X1+X2-1,data=testdata,dist="negbin",link="logit")
    results[k,1:4]=c(estimate$coefficients$count,estimate$coefficients$zero)
    results[k,5]=estimate$theta
    #glm.nb(expr~X1+X2-1,data=testdata2)
}

#plot histograms for estimators of each parameter
hist(results[,1],breaks=50,xlab="alpha_M (true = 1.5)",main="sample size n=150, 500 simulations")
hist(results[,2],breaks=50)
hist(results[,3],breaks=50,xlab="alpha_Pi (true = 0.3)",main="sample size n=150, 500 simulations")
hist(results[,4],breaks=50,xlab="alpha_Pi (true = 0.05)",main="sample size n=150, 500 simulations")
hist(results[,5],breaks=50,xlab="theta (true = 1.38)",main="sample size n=150, 500 simulations")

#########################################################################################################
#   ESTIMATE U*V and U*W without design matrices
#   data simulated from one factor model where U is a column, V and W are two lines, all of them will be supposed unknown
#########################################################################################################

U=matrix(runif(n)*3+1,ncol=1)
V=matrix(runif(J)*2,nrow=1)
W=-matrix(runif(J)*(2+2*runif(J)),nrow=1)
theta=1.1+runif(J)  #a known theta same for all genes for the moment 
theta.matr=matrix(0,nrow=n,ncol=J)
for(j in 1:length(theta)){
theta.matr[,j]=theta[j]
}
#mean of NB
M=exp(U%*%V)
#zero inflation probabilities
Pi=binomial()$linkinv(U%*%W)

gene.exp=NULL

for(i in 1:(n*J)){
    gene.exp[i]=rnbinom(1,mu=as.vector(M)[i],size=as.vector(theta.matr)[i])
}
gene.exp=matrix(gene.exp,ncol=J)

gene.is1=NULL
for(i in 1:(n*J)){
    gene.is1[i]=1-rbinom(1,size=1,prob=as.vector(Pi)[i])
}
gene.is1=matrix(gene.is1,ncol=J)

gene.exp.zeroinf=gene.exp
gene.exp.zeroinf[gene.is1==0]=0


#initialization
U.0=matrix(1,nrow=n,ncol=1)
V.0=matrix(1,nrow=1,ncol=J)
W.0=matrix(1,nrow=1,ncol=J)
theta0=rep(1,J)
#for this simulation : X.M=0, X.pi=0 alpha.M=0 alpha.pi=0

alt.number=30
for (alt in 1:alt.number){
    for (gene in 1:J){
        if(min(gene.exp.zeroinf[,gene])==0){
            estimate=zeroinfl(X1~X2-1|X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),dist="negbin",link="logit")
            V.0[,gene]=estimate$coefficients$count
            W.0[,gene]=estimate$coefficients$zero
            theta0[gene]=theta[gene]#estimate$theta
            #if no zeros, fit negative binomial (would be better to include a test of zero inflation)
        }else{
            estimate=glm.nb(X1~X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),link="log")
            V.0[,gene]=estimate$coefficients
            W.0[,gene]=0
            theta0[gene]=estimate$theta
        }
    }

    for(cell in 1:n){
        U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
                        theta=theta0,V=V.0,W=W.0,X.M=F,X.pi=F,
                        alpha.M=F,alpha.pi=F,
                        control=list(fnscale=-1),method="BFGS")$par        
    }
}

plot(as.vector(U%*%V),as.vector(U.0%*%V.0))

plot(V,V.0)
plot(U,U.0)
plot(W,W.0)
mean(V[1,]/V.0[1,])
mean(U.0[1,]/U[1,])
mean(W.0[1,]/W[1,])
plot(as.vector(U%*%V),as.vector(U.0%*%V.0))
plot(as.vector(binomial()$linkinv(U%*%W)),as.vector(binomial()$linkinv(U.0%*%W.0)))










# old version (to keep) 23 January 2016
# alt.number=30
# for (alt in 1:alt.number){
#     for (gene in 1:J){
#         #       V.0[,gene]=optim(fn=ziNegBin,gr=gradNegBin,j=gene,
#         #       par=c(V[,gene]+0.3,W[,gene]+0.4,theta+0.1),
#         #       X.M=NULL,U=U,X.pi=NULL,Y=gene.exp.zeroinf,
#         #       control=list(fnscale=-1),method="BFGS")$par
#         if(min(gene.exp.zeroinf[,gene])==0){
#             estimate=zeroinfl(X1~X2-1|X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),dist="negbin",link="logit")
#             V.0[,gene]=estimate$coefficients$count
#             W.0[,gene]=estimate$coefficients$zero
#             theta0[gene]=theta[gene]#estimate$theta
#             #if no zeros, fit negative binomial (would be better to include a test of zero inflation)
#         }else{
#             estimate=glm.nb(X1~X2-1,data=data.frame(cbind(gene.exp.zeroinf[,gene],U.0)),link="log")
#             V.0[,gene]=estimate$coefficients
#             W.0[,gene]=0
#             theta0[gene]=estimate$theta
#         }
#     }
#     #    alpha.M0=par.est[1:kxM,]
#     #    V0=par.est[(kxM+1):(kxM+p),]
#     #    alpha.pi0=par.est[(kxM+p+1):(kxM+p+kxPi),]
#     #    W0=par.est[(kxM+p+kxPi+1):(kxM+2*p+kxPi),]
#     #    thetas=par.est[kxM+2*p+kxPi+1,]
#     for(cell in 1:n){
#         #    U0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U0[cell,],counts=Y,
#         #    thetas=thetas,V=V0,W=W0,
#         #    offset.M=t(alpha.M0)%*%t(X.M),offset.pi=t(alpha.pi0)%*%t(X.pi),
#         #    control=list(fnscale=-1),method="BFGS")$par
#         U.0[cell,]=optim(fn=ziNegBin.U,gr=gradNegBin.U,i=cell,par=U.0[cell,],Y=gene.exp.zeroinf,
#                          theta=theta0,V=V.0,W=W.0,X.M=F,X.pi=F,
#                          alpha.M=F,alpha.pi=F,
#                          control=list(fnscale=-1),method="BFGS")$par
#         
#     }
# }





#test that my functions give the same result on these data
# matrix with one gene only
#Y=matrix(gene.exp.zeroinf,ncol=1)
#j=1
#same predictors for M and for Pi
#Z=X

#optim(fn=ziNegBin,gr=gradNegBin,j=j,X.M=matrix(X1,ncol=1),X.pi=matrix(X1,ncol=1),U=matrix(U1,ncol=1),Y=Y,par=c(1,1,0.5,0.1,1),control=list(fnscale=-1),method="BFGS")
#optim(fn=ziNegBin,gr=gradNegBin,j=j,X.M=NULL,X.pi=NULL,U=cbind(X1,U1),Y=Y,par=c(1,1,0.5,0.1,1),control=list(fnscale=-1),method="BFGS")


