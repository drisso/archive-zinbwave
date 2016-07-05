######################
# FUNCTIONS
######################
simulateZINB <- function(j = 100, s = 1, n = 10,
                         nbDisp = .1, gammaShape = 2, gammaRate = 2,
                         fc = 5, nDE = 0, upDEprop = .5, cellBiasSD = .25,
                         ZIprop = .5){
  
  stopifnot(upDEprop >= 0 | upDEprop <= 1)
  stopifnot(ZIprop >= 0 | ZIprop <= 1)
  
  #### mean NB
  # cell-specific bias for cell j
  log2theta = lapply(n, function(x) rnorm(x, 0, cellBiasSD))
  theta = lapply(log2theta, function(x) 2^x)
  
  # fold change
  lambdaBaseline = rgamma(j, shape = gammaShape, rate = gammaRate)
  fcList = lapply(1:s, function(x){
    fc = rep(1, j)
    if (nDE[x] != 0){
      de = sample(j, nDE[x])
      fc[de] = 0 # DE downregulated
      propUp = round(nDE[x] * upDEprop[x])
      if (propUp > 0){
        upDe = sample(de, propUp)
        fc[upDe] = FC[x] # DE upregulated
      }
    }
    fc
  })
  lambda = lapply(fcList, function(x) x * lambdaBaseline)
  NBmean = lapply(1:s, function(x) lambda[[x]] %*% t(theta[[x]]))
  NBmean = do.call(cbind, NBmean)
  
  #### counts
  Y = sapply(1:sum(n), function(y){
    sapply(1:j, function(x){
      rnbinom(1, mu = NBmean[x, y], size = nbDisp)  
    })
  })
  
  #### ZI 
  sumZeroBefore = sum(Y == 0)
  randomZero = sample(length(Y), round(ZIprop * prod(dim(Y))))
  Y[randomZero] = 0
  stopifnot(sumZeroBefore <= sum(Y == 0))
  Y
}


#######################
# SIMULATE
#######################
J            = 10000 # number of genes
S            = 2 # number of subpopulations
N            = rep(250, S) #number of cells for each subpopulation
NB_DISP      = .5
SHAPE        = 2
RATE         = 2
FC           = rep(5, S) #fold change
DE           = rep(1000, S) # number of DE genes
UP_DE        = rep(.2, S) # proportion of DE genes that are upregulated
CELL_BIAS_SD = .25 #gaussian
ZI           = .5 # prop of genes with technical zero

Y = simulateZINB(j = J, s = S, n = N, nbDisp = NB_DISP, gammaShape = SHAPE,
                 gammaRate = RATE, fc = FC, nDE = DE, upDEprop = UP_DE,
                 cellBiasSD = CELL_BIAS_SD, ZI = .5)


#####################
# PLOTS
#####################
color = c(rep('red', 250), rep('blue', 250))

genezero = apply(Y, 2, function(x) length(x[x == 0])/length(x))
par(lty = 0)
barplot(genezero, col = color, main = 'Perc. of genes with zero count', space=0)
par(lty = 1)

ngeneDetected = apply(Y, 1, function(x) length(x[x != 0]))
plot(table(ngeneDetected), main= 'Number of cells in which a gene is detected')
plot(density(ngeneDetected), main = 'Number of cells in which a gene is detected')
abline(v = mean(ngeneDetected), col = 'red')

aveCounts = apply(Y, 1, mean)
h = hist(aveCounts, plot = F)
h$counts = h$counts/sum(h$counts)
plot(h, xlim = c(0, 10), ylim = c(0, 1))
par(new = T)
plot(density(rgamma(J, SHAPE, RATE)), xlim = c(0,10), ylim = c(0, 1))















