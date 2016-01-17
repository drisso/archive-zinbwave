library(pscl)
source('~/GitHub/zinb/R/functions_svetlana.R')

# goal: simulate data K times with the same value of parameters (we take one gene only) 
# and check the variability of estimators

# sample size
n=50
K=500

# matrix to store the results of simulations
results=matrix(0,nrow=K,ncol=5)
X1=2*runif(n)
U1=3*(1+runif(n))
for (k in 1:K){

    X=cbind(X1,U1)
    
    alphaM=c(1.5,0.8)
    alphaPi=c(0.3,0.05)
    
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


#quantiles
quantile(results[,3],probs=c(0.05,0.25,0.5,0.75,0.95))

#now test that my functions give the same result on these data
# matrix with one gene only
Y=matrix(gene.exp.zeroinf,ncol=1)
j=1
#same predictors for M and for Pi
Z=X

optim(fn=ziNegBin,gr=gradNegBin,j=j,X.M=matrix(X1,ncol=1),X.pi=matrix(X1,ncol=1),U=matrix(U1,ncol=1),Y=Y,par=c(1,1,0.5,0.1,1),control=list(fnscale=-1),method="BFGS")


