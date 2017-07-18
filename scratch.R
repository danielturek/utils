
## testing Chris's new sampler_categorical for dcat nodes:
library(nimble)

code <- nimbleCode({
    x ~ dcat(p[1:5])
    ##y ~ dnorm(x, 10)
})
constants <- list(p = c(.1, .15, .2, .25, .3))
data <- list()##y=1)
inits <- list(x = 1)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes = NULL)
##conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addSampler('x', 'categorical')
conf$printSamplers()

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE)

niter <- 100000
set.seed(0); samples <- runMCMC(Cmcmc, niter)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

table(samples[,1])/niter
##      1       2       3       4       5 
##0.09990 0.15075 0.20022 0.24818 0.30095 

##      1       2       3       4       5 
##  0.100    0.15     0.2    0.25     0.3






## Nick's test of fakeDist() and mixedSizes

library(nimble)

dFakeDist <- nimbleFunction(
    run = function(x = double(0), y = double(2), log = integer(0, default = 0)) {
        returnType(double(0))
        return(0)
    })

rFakeDist <- nimbleFunction(
    run = function(n = integer(0),  y = double(2)) {
        returnType(double(0))
        if(n != 1) nimPrint('rmyexp only allows n = 1; using n = 1')
        return(0)
    }
)

registerDistributions(list(
    dFakeDist = list(
        BUGSdist = 'dFakeDist(y)',
        types = c('value = double(0)', 'y = double(2)'),
        mixedSizes = TRUE)
))

fakeCode <- nimbleCode({
    x ~ dFakeDist(y[1:3, 1:2])
})

fakeModel <- nimbleModel(code = fakeCode, data = list(x = 0), constants = list(y = matrix(nrow = 3, ncol = 2)))



## old MCMC control list defaults
##log = FALSE,
##reflective = FALSE,
##adaptive = TRUE,
##adaptScaleOnly = FALSE,
##adaptInterval = 200,
##scale = 1,
##propCov = 'identity',
##sliceWidth = 1,
##sliceMaxSteps = 100,
##sliceAdaptFactorMaxIter = 15000,  ##factorBurnIn = 15000,
##sliceAdaptFactorInterval = 1000,  ##factorAdaptInterval = 1000,
##sliceAdaptWidthMaxIter = 512,     ##sliceBurnIn = 512,
##sliceAdaptWidthTolerance = 0.1,
##scaleAdaptInterval = 200,
##sliceWidths = 'oneVec',
##pfNparticles = 1000,
##pfResample = FALSE,
##pfOptimizeNparticles = FALSE,
##pfType = 'bootstrap',
##pfLookahead = 'simulate',
##carUseConjugacy = TRUE




## test of current samples from AF_slice sampler
library(nimble)
load('~/github/hybridBlockSamplers/data/model_litters.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$addSampler(c('a','b'), 'AF_slice')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); system.time(samples <- runMCMC(Cmcmc, 20000))
if(!all(c(round(samples[1,],7) == c(14.7334414, 1.2708328, 2.8708787, 0.9404118),
          round(samples[1000,],7) == c(321.6325351, 1.5995729, 34.4197040, 0.5525827),
          round(samples[10000,2:4],6) == c(4.443445, 329.645757, 1.776704),
          round(samples[20000,],6) == c(2803.235783, 4.217422, 390.150481, 1.218532))))
    stop('something broke') else message('OK')




## test of conf$removeSamplers('a', 'b', 'c')
library(nimble)
load('~/github/hybridBlockSamplers/data/model_litters.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$removeSamplers('a', 'b')
conf$printSamplers()





## dcar_proper example from WinBUGS geoBUGS user manual:
## for running in NIMBLE
library(nimble)

code <- nimbleCode({
    ## Set up 'data' to define spatial dependence structure
    ## =====================================
    ## The vector C[] required as input into the car.proper distribution is a vector respresention
    ## of the weight matrix with elements Cij. The first J1 elements of the C[] vector contain the
    ## weights for the J1 neighbours of area i=1; the (J1+1) to J2 elements of the C[] vector contain
    ## the weights for the J2 neighbours of area i=2; etc.
    ## To set up this vector, we need to define a variable cumsum, which gives the values of J1,
    ## J2, etc.; we then set up an index matrix pick[,] with N columns corresponding to the
    ## i=1,...,N areas, and with the same number of rows as there are elements in the C[] vector 
    ## (i.e. L). The elements C[ (cumsum[i]+1):cumsum[i+1] ] correspond to 
    ## the set of weights Cij associated with area i, and so we set up ith column of the matrix pick[,]
    ## to have a 1 in all the rows k for which cumsum[i] < k <= cumsum[i+1], and 0's elsewhere. 
    ## For example, let N=4 and cumsum=c(0,3,5,6,8), so area i=1 has 3 neighbours, area i=2 has 2 
    ## neighbours, area i=3 has 1 neighbour and area i=4 has 2 neighbours. The the matrix pick[,] is:
    ##                pick                  
    ##             1, 0, 0, 0,                   
    ##             1, 0, 0, 0,                   
    ##             1, 0, 0, 0,                    
    ##             0, 1, 0, 0,                   
    ##             0, 1, 0, 0,                   
    ##             0, 0, 1, 0,                   
    ##             0, 0, 0, 1,                   
    ##             0, 0, 0, 1,                   
    ##
    ## We can then use the inner product (inprod(,)) function in WinBUGS and the kth row of pick to 
    ## select which area corresponds to the kth element in the vector C[]; likewise, we can use inprod(,) 
    ## and the ith column of pick to select the elements of C[] which correspond to area i.
    ##
    ## Note: this way of setting up the C vector is somewhat convoluted!!!! In future versions, we hope the 
    ## GeoBUGS adjacency matrix tool will be able to dump out the relevant vectors required. Alternatively, 
    ## the C vector could be created using another package (e.g. Splus) and read into WinBUGS as data.
    ##
    cumsum[1] <- 0
    for(i in 2:(N+1)) {
        cumsum[i] <- sum(num[1:(i-1)])
    }
    for(k in 1:L) {
        for(i in 1:N) {
            pick[k,i] <- step(k - cumsum[i] - epsilon)  * step(cumsum[i+1] - k)
            ##  pick[k,i] = 1    if     cumsum[i] < k <= cumsum[i=1];  otherwise, pick[k,i] = 0
        }
        C[k] <- sqrt(E[adj[k]] / inprod(E[1:N], pick[k,1:N]))    # weight for each pair of neighbours
    }
    epsilon <- 0.0001
    ## Model
    ## =====
    ## Priors:
    alpha  ~ dnorm(0, 0.0001)
    prec  ~ dgamma(0.5, 0.0005)     ## prior on precision
    v <- 1/prec                     ## variance
    sigma <- sqrt(1 / prec)         ## standard deviation
    gamma.min <- min.bound(C[1:L], adj[1:L], num[1:N], M[1:N])
    gamma.max <- max.bound(C[1:L], adj[1:L], num[1:N], M[1:N])
    gamma ~ dunif(gamma.min, gamma.max)
    S[1:N] ~ car.proper(theta[1:N], C[1:L], adj[1:L], num[1:N], M[1:N], prec, gamma)
    ## Likelihood:
    for(i in 1:N) {
        log(mu[i]) <- log(E[i]) + S[i]
        Y[i] ~ dpois(mu[i])
        RR[i] <- exp(S[i])      ## Area-specific relative risk
        theta[i] <- alpha
    }
})

N <- 56
E <- c(1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1, 2.3, 2.0, 6.6,
       4.4, 1.8, 1.1, 3.3, 7.8, 4.6, 1.1, 4.2, 5.5, 4.4,
       10.5,22.7, 8.8, 5.6,15.5,12.5, 6.0, 9.0,14.4,10.2,
       4.8, 2.9, 7.0, 8.5,12.3,10.1,12.7, 9.4, 7.2, 5.3,
       18.8,15.8, 4.3,14.6,50.7, 8.2, 5.6, 9.3,88.7,19.6,
       3.4, 3.6, 5.7, 7.0, 4.2, 1.8)
M <- 1/E
##X <- c(16,16,10,24,10,24,10, 7, 7,16, 7,16,10,24, 7,16,10,
##       7, 7,10, 7,16,10, 7, 1, 1, 7, 7,10,10,7,24,10, 7, 7,
##       0,10, 1,16, 0, 1,16,16, 0, 1, 7, 1, 1, 0, 1,1, 0, 1, 1,16,10)
num <- c(3, 2, 1, 3, 3, 0, 5, 0, 5, 4, 0, 2, 3, 3, 2, 6, 6, 6, 5, 3,
         3, 2, 4, 8, 3, 3, 4, 4, 11, 6, 7, 3, 4, 9, 4, 2, 4, 6, 3, 4, 
         5, 5, 4, 5, 4, 6, 6, 4, 9, 2, 4, 4, 4, 5, 6, 5)
adj <- c(19, 9, 5, 10, 7, 12, 28, 20, 18, 19, 12, 1, 
         17, 16, 13, 10, 2, 29, 23, 19, 17, 1, 22, 16, 7, 2, 
         5, 3, 19, 17, 7, 35, 32, 31, 29, 25, 29, 22, 21, 17,
         10, 7, 29, 19, 16, 13, 9, 7, 56, 55, 33, 28, 20, 4,
         17, 13, 9, 5, 1, 56, 18, 4, 50, 29, 16, 16, 10, 39, 34, 29, 9,
         56, 55, 48, 47, 44, 31, 30, 27, 29, 26, 15, 43, 29, 25,
         56, 32, 31, 24, 45, 33, 18, 4, 50, 43, 34, 26, 25, 23, 21,
         17, 16, 15, 9, 55, 45, 44, 42, 38, 24, 47, 46, 35, 32, 27, 24, 14, 
         31, 27, 14, 55, 45, 28, 18, 54, 52, 51, 43, 42, 40, 39, 29, 23,
         46, 37, 31, 14, 41, 37, 46, 41, 36, 35, 54, 51, 49, 44, 42, 30,
         40, 34, 23, 52, 49, 39, 34, 53, 49, 46, 37, 36, 51, 43, 38, 34, 30,
         42, 34, 29, 26, 49, 48, 38, 30, 24, 55, 33, 30, 28, 53, 47, 41,
         37, 35, 31, 53, 49, 48, 46, 31, 24, 49, 47, 44, 24, 54, 53, 52,
         48, 47, 44, 41, 40, 38, 29, 21, 54, 42, 38, 34, 54, 49, 40, 34,
         49, 47, 46, 41, 52, 51, 49, 38, 34, 56, 45, 33, 30, 24, 18,
         55, 27, 24, 20, 18)
L <- length(adj)
constants <- list(N=N, L=L, E=E, num=num, adj=adj, M=M)

Y <- c(9, 39, 11, 9, 15, 8, 26, 7, 6, 20, 13, 5, 3, 8, 17, 9, 2, 7, 9,
       7, 16, 31, 11, 7, 19, 15, 7, 10, 16, 11, 5, 3, 7, 8, 11, 9, 11,
       8, 6, 4, 10, 8, 2, 6, 19, 3, 2, 3, 28, 6, 1, 1, 1, 1, 0, 0)
data <- list(Y=Y)

inits <- list(alpha = 3, prec = 1, gamma = 0.1,
              S = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))


## dcar_proper example from WinBUGS geoBUGS user manual:
## this is the model as I ran it in WinBUGS,
## and the RESULTS FROM WINBUGS for comparison.

model
{
  for(i in 1 : N) { 
     m[i] <- 1/E[i]       # scaling factor for variance in each cell
  } 
  cumsum[1] <- 0
  for(i in 2:(N+1)) {
     cumsum[i] <- sum(num[1:(i-1)])
  }	
  for(k in 1 : sumNumNeigh) { 
 	 for(i in 1:N) {
          pick[k,i] <- step(k - cumsum[i] - epsilon)  * step(cumsum[i+1] - k)   
      }                                                       
      C[k] <- sqrt(E[adj[k]] / inprod(E[], pick[k,]))    # weight for each pair of neighbours
  }
  epsilon <- 0.0001
  for (i in 1 : N) {
      Y[i]  ~ dpois(mu[i])
      log(mu[i]) <- log(E[i]) + S[i]
      RR[i] <- exp(S[i])      # Area-specific relative risk 
      theta[i] <- alpha
  }
  S[1:N] ~ car.proper(theta[], C[], adj[], num[], m[], prec, gamma)
  alpha  ~ dnorm(0, 0.0001)  
  prec  ~ dgamma(0.5, 0.0005)     # prior on precision
   v <- 1/prec                               # variance
  sigma <- sqrt(1 / prec)               # standard deviation
  gamma.min <- min.bound(C[], adj[], num[], m[])
  gamma.max <- max.bound(C[], adj[], num[], m[])
  gamma ~ dunif(gamma.min, gamma.max)
}

## data
list(N = 56,  
     Y   = c(    9,   39,   11,    9,   15,    8,   26,    7,    6,   20,
         13,    5,    3,    8,   17,    9,    2,    7,    9,    7,
         16,   31,   11,    7,   19,   15,    7,   10,   16,   11,
         5,    3,    7,    8,   11,    9,   11,    8,    6,    4,
         10,    8,    2,    6,   19,    3,    2,    3,   28,    6,
         1,    1,    1,    1,    0,    0),
     E = c(1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1, 2.3, 2.0, 6.6, 4.4, 1.8, 1.1, 3.3, 7.8, 4.6, 1.1, 4.2, 5.5, 4.4,
         10.5,22.7, 8.8, 5.6,15.5,12.5, 6.0, 9.0,14.4,10.2, 4.8, 2.9, 7.0, 8.5,12.3,10.1,12.7, 9.4, 7.2, 5.3,
         18.8,15.8, 4.3,14.6,50.7, 8.2, 5.6, 9.3,88.7,19.6, 3.4, 3.6, 5.7, 7.0, 4.2, 1.8),
     num = c(3, 2, 1, 3, 3, 0, 5, 0, 5, 4, 0, 2, 3, 3, 2, 6, 6, 6, 5, 3, 3, 2, 4, 8, 3, 3, 4, 4, 11, 6, 7, 3,
         4, 9, 4, 2, 4, 6, 3, 4, 5, 5, 4, 5, 4, 6, 6, 4, 9, 2, 4, 4, 4, 5, 6, 5),
     adj = c(19, 9, 5, 10, 7, 12, 28, 20, 18, 19, 12, 1, 17, 16, 13, 10, 2, 
         29, 23, 19, 17, 1, 22, 16, 7, 2, 5, 3, 19, 17, 7, 35, 32, 31, 
         29, 25, 29, 22, 21, 17, 10, 7, 29, 19, 16, 13, 9, 7, 56, 55, 33, 28, 20, 4, 
         17, 13, 9, 5, 1, 56, 18, 4, 50, 29, 16, 16, 10, 39, 34, 29, 9, 56, 55, 48, 47, 44, 31, 30, 27, 
         29, 26, 15, 43, 29, 25, 56, 32, 31, 24, 45, 33, 18, 4, 50, 43, 34, 26, 25, 23, 21, 17, 16, 15, 9, 
         55, 45, 44, 42, 38, 24, 47, 46, 35, 32, 27, 24, 14, 31, 27, 14, 55, 45, 28, 18, 
         54, 52, 51, 43, 42, 40, 39, 29, 23, 46, 37, 31, 14, 41, 37, 46, 41, 36, 35, 54, 51, 49, 44, 42, 30, 
         40, 34, 23, 52, 49, 39, 34, 53, 49, 46, 37, 36, 51, 43, 38, 34, 30, 42, 34, 29, 26, 49, 48, 38, 30, 24, 
         55, 33, 30, 28, 53, 47, 41, 37, 35, 31, 53, 49, 48, 46, 31, 24, 49, 47, 44, 24, 54, 53, 52, 48, 47, 44, 41, 40, 38, 
         29, 21, 54, 42, 38, 34, 54, 49, 40, 34, 49, 47, 46, 41, 52, 51, 49, 38, 34, 56, 45, 33, 30, 24, 18, 55, 27, 24, 20, 18),
     sumNumNeigh = 234)

## inits
list(alpha=3, prec=1, S=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), gamma=0.1)

## results from WinBUGS for dcar_proper:
## dcar_proper example from WinBUGS geoBUGS user manual:
## for running in NIMBLE
##  node	 mean	 sd	2.5%	median	97.5%
##  alpha	-0.1391	0.1592	-0.4531	-0.1396	0.1732
##  gamma	0.1634	0.01875	0.1129	0.1689	0.1823
##  prec	0.3342	0.09314	0.1859	0.3226	0.5491

winbugsResults <- matrix(c(-0.1391,0.1592,-0.4531,-0.1396,0.1732,0.1634,0.01875,0.1129,0.1689,0.1823,0.3342,0.09314,0.1859,0.3226,0.5491), nrow=3, ncol=5, byrow=TRUE)
rownames(winbugsResults) <- c('alpha', 'gamma', 'prec')
colnames(winbugsResults) <- c('mean', 'sd', '2.5%', 'median', '97.5%')
winbugsResults
##          mean      sd    2.5%  median  97.5%
## alpha -0.1391 0.15920 -0.4531 -0.1396 0.1732
## gamma  0.1634 0.01875  0.1129  0.1689 0.1823
## prec   0.3342 0.09314  0.1859  0.3226 0.5491

##S[1]	1.743	0.3439	0.001097	1.016	1.762	2.368	501	99500
##S[2]	1.397	0.1632	5.225E-4	1.067	1.402	1.707	501	99500
##S[3]	1.143	0.3117	9.69E-4	0.4913	1.157	1.714	501	99500
##S[4]	1.126	0.3423	0.001084	0.4093	1.143	1.752	501	99500
##S[5]	1.129	0.2605	7.779E-4	0.59	1.14	1.611	501	99500
##S[6]	1.018	0.3665	0.001187	0.2505	1.037	1.681	501	99500
##S[7]	1.084	0.1942	5.992E-4	0.6894	1.09	1.449	501	99500
##S[8]	0.9187	0.3893	0.001243	0.09141	0.9405	1.618	501	99500
##S[9]	0.9715	0.4085	0.001309	0.1061	0.9946	1.7	501	99500
##S[10]	1.036	0.2199	6.416E-4	0.5859	1.043	1.446	501	99500
##S[11]	0.9193	0.2842	9.177E-4	0.3267	0.9312	1.442	501	99500
##S[12]	0.8776	0.4498	0.001464	-0.07956	0.9062	1.679	501	99500
##S[13]	0.8362	0.5806	0.001909	-0.4341	0.8821	1.84	501	99500
##S[14]	0.7018	0.3602	0.001163	-0.05503	0.7194	1.357	501	99500
##S[15]	0.6455	0.2399	7.549E-4	0.1524	0.6531	1.092	501	99500
##S[16]	0.6387	0.3135	9.77E-4	-0.01128	0.6519	1.214	501	99500
##S[17]	0.5477	0.6598	0.002289	-0.8985	0.606	1.674	501	99500
##S[18]	0.3313	0.3711	0.001216	-0.4439	0.3496	1.005	501	99500
##S[19]	0.4437	0.3107	9.139E-4	-0.2032	0.4586	1.013	501	99500
##S[20]	0.3253	0.364	0.001044	-0.4397	0.3424	0.9856	501	99500
##S[21]	0.2978	0.2382	7.383E-4	-0.1894	0.3042	0.7465	501	99500
##S[22]	0.245	0.1679	5.286E-4	-0.09421	0.2487	0.5618	501	99500
##S[23]	0.1375	0.2776	8.674E-4	-0.4362	0.1469	0.6536	501	99500
##S[24]	-0.09633	0.3869	0.001406	-0.9039	-0.07928	0.6078	501	99500
##S[25]	0.1488	0.2098	7.077E-4	-0.278	0.1544	0.5439	501	99500
##S[26]	0.1014	0.2378	7.583E-4	-0.3835	0.1082	0.5479	501	99500
##S[27]	0.02065	0.3505	0.001112	-0.7117	0.03649	0.6625	501	99500
##S[28]	0.007397	0.2895	9.733E-4	-0.589	0.01767	0.5463	501	99500
##S[29]	0.0935	0.2217	6.994E-4	-0.3598	0.09998	0.5113	501	99500
##S[30]	-0.1558	0.2916	9.563E-4	-0.7564	-0.1461	0.3874	501	99500
##S[31]	-0.111	0.4157	0.001364	-0.9859	-0.08949	0.6387	501	99500
##S[32]	-0.05884	0.5173	0.001618	-1.162	-0.02871	0.8635	501	99500
##S[33]	-0.155	0.3465	0.001048	-0.8779	-0.1405	0.4833	501	99500
##S[34]	-0.266	0.3345	0.001164	-0.9582	-0.2517	0.3513	501	99500
##S[35]	-0.1531	0.2626	8.8E-4	-0.6914	-0.1448	0.3414	501	99500
##S[36]	-0.1871	0.2935	9.404E-4	-0.7946	-0.1764	0.3574	501	99500
##S[37]	-0.2241	0.2653	8.216E-4	-0.7694	-0.2154	0.2753	501	99500
##S[38]	-0.4511	0.3403	0.001183	-1.152	-0.4382	0.1817	501	99500
##S[39]	-0.2213	0.3515	0.001097	-0.9473	-0.2063	0.428	501	99500
##S[40]	-0.517	0.4516	0.001434	-1.464	-0.4969	0.3053	501	99500
##S[41]	-0.641	0.2547	8.162E-4	-1.164	-0.6343	-0.1627	501	99500
##S[42]	-0.5665	0.2712	9.55E-4	-1.124	-0.557	-0.06094	501	99500
##S[43]	-0.619	0.5202	0.001687	-1.717	-0.5893	0.3157	501	99500
##S[44]	-0.7859	0.3004	9.681E-4	-1.404	-0.7759	-0.2268	501	99500
##S[45]	-0.6555	0.172	7.411E-4	-1.009	-0.6506	-0.3324	501	99500
##S[46]	-0.8712	0.4106	0.001378	-1.728	-0.8558	-0.1191	501	99500
##S[47]	-1.154	0.5385	0.001751	-2.284	-1.128	-0.1731	501	99500
##S[48]	-0.997	0.3996	0.001411	-1.83	-0.9819	-0.2598	501	99500
##S[49]	-0.8482	0.1408	6.438E-4	-1.136	-0.8441	-0.5841	501	99500
##S[50]	-0.6966	0.2663	0.001019	-1.253	-0.6847	-0.2069	501	99500
##S[51]	-1.083	0.6658	0.002157	-2.505	-1.042	0.1048	501	99500
##S[52]	-1.299	0.68	0.002363	-2.738	-1.26	-0.07973	501	99500
##S[53]	-1.406	0.566	0.002052	-2.597	-1.376	-0.3815	501	99500
##S[54]	-1.354	0.5102	0.001772	-2.43	-1.328	-0.4254	501	99500
##S[55]	-1.422	0.6701	0.002154	-2.848	-1.378	-0.2378	501	99500
##S[56]	-1.364	0.9762	0.003289	-3.505	-1.283	0.2978	501	99500




## test of compilation of dcar_proper()
## CURRENTLY FAILING
## NEED TO GET HELP FROM CJP
##
## ALSO:
## building nimble package gives this disturbing warning:
## ** testing if installed package can be loaded
## Warning: failed to assign RegisteredNativeSymbol for C_dcar_normal to C_dcar_normal since C_dcar_normal is already defined in the ‘nimble’ namespace
## * DONE (nimble)
##
library(nimble)

code <- nimbleCode({
    x[1:3] ~ dcar_proper(mu[1:3], C[1:6], adj[1:6], num[1:3], M[1:3], t, g)
})
constants <- list(mu = 1:3, C = 1:6, adj = 1:6, num = 1:3, M = 1:3, t = 1, g = 1)
data <- list()
inits <- list(x = c(1, 2, NA))

Rmodel <- nimbleModel(code, constants, data, inits, calculate=FALSE)

Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)





## problem from TR J/A 2
## deck of 52, A=1, face = 10
## drawing cards w/ replacement
## draw 2, prob(sum is even) ?
## draw 4, or draw 100?
pE <- 8/13
pO <- 5/13

N <- 100
pevensum <- numeric(N)
pevensum[1] <- pE

if(N > 1) for(i in 2:N) {
    pevensum[i] <- (1-pevensum[i-1])*pO + pevensum[i-1]*pE
    }
pevensum

pevensum[99] - pevensum[100]
plot(1:N, pevensum, type = 'l', xlim=c(1,10))

## testing as.carAdjacency() conversions
library(nimble)
as.carAdjacency
as.carAdjacency(matrix(c(0,0,0,0,0,2,0,2,0), 3, 3))

as.carAdjacency(list(numeric(0), 3, 2), list(numeric(0), 2, 2))

## testing about providing numIslands argument to dcar_normal()

library(nimble)

code <- nimbleCode({
    for(i in 1:L) {   weights[i] <- 1   }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    ##weights = c(1,1,1,1,1,1,1,1),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )

code <- nimbleCode({
    for(i in 1:L) {   weights[i] <- 1   }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1, 3)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    ##weights = c(1,1,1,1,1,1,1,1),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )

code <- nimbleCode({
    for(i in 1:L) {   weights[i] <- 1   }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1, numI)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    ##weights = c(1,1,1,1,1,1,1,1),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5), numI = 4)

code <- nimbleCode({
    for(i in 1:L) {   weights[i] <- 1   }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1, numI)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    ##weights = c(1,1,1,1,1,1,1,1),
    L = 8,
    numI = 1)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )

## 5-8

code <- nimbleCode({
    S[1:N] ~ car.normal(adj[1:L], num = num[1:N], tau = 1)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )

code <- nimbleCode({
    S[1:N] ~ car.normal(adj[1:L], num = num[1:N], tau = 1, numIslands = 3)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )

code <- nimbleCode({
    S[1:N] ~ car.normal(adj[1:L], num = num[1:N], tau = 1, numIslands = numI)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    L = 8)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5), numI = 4)

library(nimble)
code <- nimbleCode({
    ##S[1:N] ~ car.normal(adj[1:L], num = num[1:N], tau = 1, numIslands = numI)
    S[1:N] ~ car.normal(adj[1:L], num = num[1:N], tau = 1, c = numI, zero_mean = 1)
    for(i in 1:N) { y[i] ~ dnorm(S[i], 1) }
})
constants <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    L = 8,
    numI = 1)
data <- list( y = c(1,2,3,4,5) )
inits <- list( S = c(1,2,3,4,5) )
inits <- list( S = rep(NA,5) )


Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmodel$S
Cmodel$S

Rmodel$calculate()
Cmodel$calculate()
niter <- 20
set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
round(Rsamples - Csamples, 3)
Rsamples
Csamples
round(apply(Rsamples, 1, sum), 5)
round(apply(Csamples, 1, sum), 5)


debug(Rmcmc$run)

##[20,] -1.81243476 -1.10311996  0.31018093  0.69365670 1.9117171


niter <- 10000
set.seed(0); Cmcmc$run(niter)
Csamples <- as.matrix(Cmcmc$mvSamples)
samplesPlot(Csamples)
samplesPlot(Csamples, var = 1:2)



##options(warn=2)
##options(error = recover)



## testing how WinBUGS imposes the sum-to-zero constraint...
## on a per-island-basis?

library(nimble)

code <- nimbleCode({
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1)
    for(i in 1:N) {
        y[i] ~ dnorm(S[i], 1)
    }
    sum1 <- S[1] + S[2]
    sum2 <- S[3] + S[4] + S[5]
    sum3 <- S[1] + S[2] + S[3] + S[4] + S[5]
})
##
data <- list(
    N = 5,
    num = c(1,1,2,2,2),
    adj = c(2,1,   4,5,   3,5,   3,4),
    weights = c(1,1,1,1,1,1,1,1),
    L = 8,
    y = c(1,2,3,4,5)
)
##
inits <- list(
    S = c(1,2,3,4,5)
)

catCode(code, data, inits, file='~/temp/BUGS.txt')

Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['y'], inits)


conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)




## learning hose to use tapply()

require(stats)
rbinom(32, n = 5, prob = 0.4)
as.factor(rbinom(32, n = 5, prob = 0.4))
groups <- as.factor(rbinom(32, n = 5, prob = 0.4))
groups
tapply(groups, groups, length) #- is almost the same as
table(groups)

## contingency table from data.frame : array with named dimnames
str(warpbreaks)
head(warpbreaks)
warpbreaks$breaks
identical(warpbreaks[,-1], warpbreaks[-1])
str(warpbreaks[-1])

tapply(warpbreaks$breaks, warpbreaks[,-1], sum)
tapply(warpbreaks$breaks, warpbreaks[,-1], sum, simplify=FALSE)
tapply(warpbreaks$breaks, warpbreaks[,-1], function(x) x+1000)

class(tapply(warpbreaks$breaks, warpbreaks[,-1], function(x) x+1000))
typeof(tapply(warpbreaks$breaks, warpbreaks[,-1], function(x) x+1000))
tapply(warpbreaks$breaks, warpbreaks[,-1], function(x) x+1000)[1,1]

warpbreaks[, 3, drop = FALSE]
warpbreaks[, 3]
tapply(warpbreaks$breaks, warpbreaks[, 3, drop = FALSE], sum)
tapply(warpbreaks$breaks, warpbreaks[, 3], sum)

aggregate(warpbreaks$breaks, warpbreaks[, 3, drop=FALSE], sum)
aggregate(warpbreaks$breaks, warpbreaks[, -1, drop=FALSE], sum)




n <- 17
rep_len(1:3, n)
factor(rep_len(1:3, n), levels = 1:5)
fac <- factor(rep_len(1:3, n), levels = 1:5)
fac
table(fac)
tapply(1:n, fac, sum)
tapply(1:n, fac, sum, default = 0) # maybe more desirable
tapply(1:n, fac, sum, simplify = FALSE)
##tapply(1:n, fac, range, default=c(1,3))
tapply(1:n, fac, quantile)
tapply(1:n, fac, length) ## NA's
tapply(1:n, fac, length, default = 0) # == table(fac)



## example of ... argument: find quarterly means
tapply(presidents, cycle(presidents), mean, na.rm = TRUE)

ind <- list(c(1, 2, 2), c("A", "A", "B"))
ind
table(ind)
tapply(1:3, ind) #-> the split vector
tapply(1:3, ind, sum)




## checking car car_calcNumIslands()

myfun <- function() CAR_calcNumIslands(adj, num)

c_CAR_calcNumIslands <- compileNimble(CAR_calcNumIslands)
myfun <- function() c_CAR_calcNumIslands(adj, num)

num <- c(1, 1)
adj <- c(2, 1)
myfun()
## c = 1

library(nimble)
myfun <- function() {
    N <- length(num)
    L <- sum(num)
    weights <- rep(1, L)
    code <- nimbleCode({
        S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
    })
    constants <- list(adj=adj, weights=weights, num=num, N=N, L=L)
    data <- list(S = rep(0,N))
    inits <- list(tau = 1)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
}

num <- c(1, 1)
adj <- c(2, 1)
myfun()
## c = 1

num <- c(2, 1, 1)
adj <- c(2,3,   1,  1)
myfun()
## c = 1

num <- c(2, 2, 2)
adj <- c(2,3,   1,3,   1,2)
myfun()
## c = 1

num <- c(1, 0, 1)
adj <- c(3,    1)
myfun()
## c = 2

num <- c(2, 2, 2, 0)
adj <- c(2,3,   1,3,   1,2)
myfun()
## c = 2

num <- c(2, 2, 2, 0, 0)
adj <- c(2,3,   1,3,   1,2)
myfun()
## c = 3

num <- c(2, 2, 2, 0, 1, 1)
adj <- c(2,3,   1,3,   1,2,     6,    5)
myfun()
## c = 3

num <- c(2, 2, 2, 0, 1, 1, 1, 1)
adj <- c(2,3,   1,3,   1,2,     6,    5,    8,     7)
myfun()
## c = 4

num <- c(2, 2, 2, 0, 1, 1, 2, 2, 2)
adj <- c(2,3,   1,3,   1,2,     6,    5,    8,9,     7,9,   7,8)
myfun()
## c = 4

num <- c(2, 2, 2, 0, 1, 1, 2, 2, 2, 0, 0)
adj <- c(2,3,   1,3,   1,2,     6,    5,    8,9,     7,9,   7,8)
myfun()
## c = 6




## equivalance of geoBUGS / winBUGS and nimble
## CAR dcar_normal density evaluation for tau?

N <- 4
x <- 1:4
num <- c(3,2,2,1)
adj <- c(2,3,4,   1,3,    1,2,   1)
weights <- c(2,2,2,   2,3,   2,3,  2)

k <- length(x)
c <- 1
lp <- 0
count <- 1
for(i in 1:k) {
    if(num[i] == 0)   c <- c + 1
    xi <- x[i]
    for(j in 1:num[i]) {
        xj <- x[adj[count]]
        lp <- lp + weights[count] * (xi-xj)^2
        count <- count + 1
    }
}
lp

k <- length(x)
lp <- 0
count <- 1
qf <- 0
for(i in 1:k) {
    xi <- x[i]
    mu <- 0
    wPlus <- 0
    for(j in 1:num[i]) {
        mu <- mu + x[adj[count]] * weights[count]
        wPlus <- wPlus + weights[count]
        count <- count + 1
    }
    mu <- mu / wPlus
    qf <- qf + 1/2 * xi * (xi - mu) * wPlus
}
qf
qf*4


if(count != (length(adj)+1)) stop('something wrong')
lp <- lp * (-1/2) * tau / 2
lp <- lp + (k-c)/2 * log(tau/2/pi)
lp




## testing island nodes = NA values
library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
##
code <- nimbleCode({
    alpha0 ~ dflat()
    tau ~ dgamma(0.001, 0.001)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
    for(i in 1:3) {
        mu[i] <- S[i] + alpha0
        Y[i] ~ dnorm(mu[i], 2)
    }
    z[1] ~ dnorm(S[6], 1)
    z[2] ~ dnorm(S[7], 1)
    z[3] ~ dnorm(S[8] + 0*S[8]^2, 1)
    z[4] ~ dnorm(S[9] + 0*S[9]^2, 1)
    z[5] ~ dnorm(S[10] + 0*S[10]^2, 1)
})
##
data <- list(
    N = 10,
    num = c(1, 2, 1, 0, 0, 0, 0, 0, 0, 0),
    adj = c(2,   1,3,    2),
    weights = c(1,1,1,1),
    L = 4,
    Y = c(1,2,3),
    z = c(0, 0, 0, 0, 0)
)
##
inits <- list(
    alpha0 = 0,
    tau = 1,
    S = c(0, 0, 0, 3, NA, 0, NA, 0, NA, NaN)
)
##
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data[c('Y','z')], inits)

conf <- configureMCMC(Rmodel, control = list(log=TRUE))
conf$printSamplers()
conf$addMonitors('S')
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE)

niter <- 200000#500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)
dim(samples)
length(samples)
##
means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
medians <- apply(samples, 2, median)
q025 <- apply(samples, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
q975 <- apply(samples, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
res <- cbind(means, sds, medians, q025, q975)
res

               means         sds      medians        q025        q975
S[1]   -3.941210e-01   0.4893330 -0.225937468 -1.59038344   0.1910843
S[2]    1.565404e-05   0.2975876  0.000138869 -0.67216986   0.6661552
S[3]    3.941054e-01   0.4887008  0.227770323 -0.19113529   1.5807439
alpha0  2.003152e+00   0.4080147  2.003949375  1.19972110   2.8056662
tau     9.120936e+01 280.8845605  4.164618835  0.04471161 861.5163085

## winBUGS output


conf$printSamplers()
Rmodel$getLogProb('alpha0')
Rmodel$getLogProb('tau')
Rmodel$getLogProb('S')
Rmodel$getLogProb('Y')

debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)

niter <- 30
set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
sampnames <- dimnames(Rsamples)[[2]]
Rsamples[, sampnames]

Csamples[, sampnames]

Rsamples[, sampnames] - Csamples[, sampnames]


catListContents <- function(lst) {
    cat('list(')
    for(i in seq_along(lst)) {
        cat(names(lst)[i])
        cat(' = ')
        val <- lst[[i]]
        if(length(val) == 1) {
            cat(val)
        } else {
            cat('c(')
            cat(paste(val, collapse=', '))
            cat(')')
        }
        if(i < length(lst)) cat(', ')
    }
    cat(')\n\n')
}
catCode <- function(code, data, inits, file) {
    if(!missing(file)) sink(file)
    cat('\n\nmodel\n')
    print(code)
    cat('\n\n## data\n')
    catListContents(data)
    cat('\n## inits\n')
    catListContents(inits)
    if(!missing(file)) sink()
}

catCode(code, data, inits)
catCode(code, data, inits, file='~/temp/BUGS.txt')


x <- runif(100000, 0, 10000000)
hist(x)
t <- 1/x
hist(t, breaks=1000, xlim=c(0, 0.0001))



## finding BUGS tau sampling (11)
## with a gamma(0.001, 0.001) prior on tau
library(nimble)
##
code <- nimbleCode({
    tau ~ dgamma(0.001, 0.001)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 4,
    num = c(1,2,1,0),
    adj = c(2,   1,3,    2),
    weights = c(1,1,1,1),
    L = 4,
    S = c(1,0,0,0)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 1000000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:1000000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
##q025 <- apply(samples, 2, function(x) quantile(x, 0.025))
##q975 <- apply(samples, 2, function(x) quantile(x, 0.975))
mean(samples)
var(samples)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, q025, medians, q975)
res
##        means     sds       q025 medians     q975
##2.5% 2.006003 1.99855 0.05131201  1.3909 7.396893



## true posterior (for tau):
WSSD <- 1
N <- 4
c <- 2
a <- (N-c)/2 
b <- WSSD/2
a
b

## NIMBLE
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
a
b


## winBUGS output
##node  mean   sd     2.5%     median  97.5%
##tau   2.008  2.006  0.04963  1.389   7.42

m <- 2.008
s <- 2.006
a <- m^2 / s^2
b <- m / s^2
a
b







## finding BUGS tau sampling (10)
## now with a transformation on tau
## prior is uniform on variance
## this time, with two more islands!
library(nimble)
##
code <- nimbleCode({
    sigma2 ~ dunif(0, 10000)
    tau <- 1/sigma2
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(             
    N = 9,
    num = c(1,2,2,2,2,2,1,0,0),
    adj = c(2,   1,3,    2,4,   3,5,   4,6,   5,7,   6),
    weights = c(1,1,1,1,1,1,1,1,1,1,1,1),
    L = 12,
    S = c(1,0,0,0,0,0,0,0,0)
)
##
inits <- list(
    sigma2 = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addMonitors('tau')
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 1000000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:1000000,]
dim(samples)
length(samples)
##

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
medians <- apply(samples, 2, median)
q025 <- apply(samples, 2, function(x) quantile(x, 0.025))
q975 <- apply(samples, 2, function(x) quantile(x, 0.975))
##mean(samples)
##var(samples)
##means <- mean(samples)
##q025 <- quantile(samples, 0.025)
##q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, q025, medians, q975)
res
##           means      sds       q025   medians      q975
##sigma2 0.5124162 1.035113 0.08940902 0.2994832  2.081309
##tau    3.9918965 2.838581 0.48046693 3.3390858 11.184554


## true posterior (for tau):
WSSD <- 1
N <- 9
c <- 3
a <- (N-c)/2 - 1
b <- WSSD/2
a
b

## NIMBLE
tau.samp <- samples[, 'tau']
mean(tau.samp)
var(tau.samp)

## theoretical gamma, for tau:
## mean = a/b
a / b
## var = a / b^2
a / b^2

sigma2.samp <- samples[, 'sigma2']
mean(sigma2.samp)
var(sigma2.samp)

## theoretical inverse-gamma, for sigma2:
## mean = b/(a-1), for a>1   (from wikipedia)
b / (a-1)
## var = b^2 / ((a-1)^2 * (a-2)),   for a>2   (from wikipedia)
b^2 / ((a-1)^2 * (a-2))



## winBUGS output
##node     mean   sd      2.5%     median  97.5%
##sigma2   0.487  0.8645  0.09028  0.2972  2.009
##tau      4.004  2.81    0.4977   3.365   11.08


    
## finding BUGS tau sampling (9)
## now with a transformation on tau
## prior is uniform on variance
library(nimble)
##
code <- nimbleCode({
    sigma2 ~ dunif(0, 10000)
    tau <- 1/sigma2
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(             
    N = 9,
    num = c(1,2,2,2,2,2,2,2,1),
    adj = c(2,   1,3,    2,4,   3,5,   4,6,   5,7,   6,8,   7,9,   8),
    weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
    L = 16,
    S = c(1,0,0,0,0,0,0,0,0)
)
##
inits <- list(
    sigma2 = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addMonitors('tau')
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 1000000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:1000000,]
dim(samples)
length(samples)
##

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
medians <- apply(samples, 2, median)
q025 <- apply(samples, 2, function(x) quantile(x, 0.025))
q975 <- apply(samples, 2, function(x) quantile(x, 0.975))
##mean(samples)
##var(samples)
##means <- mean(samples)
##q025 <- quantile(samples, 0.025)
##q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, q025, medians, q975)
res
##           means      sds       q025  medians       q975
##sigma2 0.2515433 0.237622 0.06912374 0.188005  0.8239204
##tau    5.9778155 3.469859 1.21370951 5.319008 14.4668096

## true posterior (for tau):
WSSD <- 1
N <- 9
c <- 1
a <- (N-c)/2 - 1
b <- WSSD/2
a
b

## NIMBLE
tau.samp <- samples[, 'tau']
mean(tau.samp)
var(tau.samp)

## theoretical gamma, for tau:
## mean = a/b
a / b
## var = a / b^2
a / b^2

sigma2.samp <- samples[, 'sigma2']
mean(sigma2.samp)
var(sigma2.samp)

## theoretical inverse-gamma, for sigma2:
## mean = b/(a-1), for a>1   (from wikipedia)
b / (a-1)
## var = b^2 / ((a-1)^2 * (a-2)),   for a>2   (from wikipedia)
b^2 / ((a-1)^2 * (a-2))



## winBUGS output
##node     mean    sd      2.5%      median   97.5%
##sigma2   0.2467  0.2254  0.06937   0.186    0.7953
##tau      6.014   3.438   1.257     5.375    14.41





## finding BUGS tau sampling (8)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 10000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 4,
    num = c(1,2,1,0),
    adj = c(2,   1,3,    2),
    weights = c(1,1,1,1),
    L = 4,
    S = c(1,0,0,0)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
mean(samples)
var(samples)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res

## true posterior:
WSSD <- 1
N <- 4
c <- 2
a <- (N-c)/2 + 1
b <- WSSD/2
a
b

## NIMBLE
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
a
b
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
 node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	2.002	2.01	0.009006	0.04942	1.389	7.469	1	50000

##a = mean^2  / sd^2
##b = mean    / sd^2

m <- 2.001
s <- 2.01
a <- m^2 / s^2
b <- m / s^2
a
b

a <- 1
b <- 0.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))

WSSD <- 1
N <- 4
c <- 4 ## BUGS is using 4 islands here
a <- (N-c)/2 + 1
b <- WSSD/2
a
b





## finding BUGS tau sampling (7)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 10000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 4,
    num = c(1,2,2,1),
    adj = c(2,   1,3,    2,4,   3),
    weights = c(1,1,1,1,1,1),
    L = 6,
    S = c(1,0,0,0)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
mean(samples)
var(samples)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res

## true posterior:
WSSD <- 1
N <- 4
c <- 1
a <- (N-c)/2 + 1
b <- WSSD/2
a
b

## NIMBLE
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
a
b
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
node	 mean	sd	MC error  2.5%	 median	97.5%	start	sample
tau	3.012  	2.463	0.01127	  0.2203 2.373	9.356	1	50000

##a = mean^2  / sd^2
##b = mean    / sd^2

m <- 3.012
s <- 2.463
a <- m^2 / s^2
b <- m / s^2
a
b

a <- 1.5
b <- 0.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))

WSSD <- 1
N <- 4
c <- 3 ## BUGS is using 3 islands here
a <- (N-c)/2 + 1
b <- WSSD/2
a
b



## finding BUGS tau sampling (6)
## transformation of tau (to sd)
library(nimble)
##
code <- nimbleCode({
    ##sd ~ dgamma(0.001, 0.001)
    ##tau ~ dgamma(0.001, 0.001)
    b ~ dbern(0.5)
    tau <- b+1
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
    ##for(i in 1:N) {
    ##    mu[i] <- S[i] + alpha0
    ##    Y[i] ~ dnorm(mu[i], 1)
    ##}
})
##
data <- list(
    N = 3,
    num = c(1,2,1),
    adj = c(2,   1,3,    2),
    weights = c(1,1,1,1),
    L = 4,
    ##Y = c(1,2,3)
    S = c(0, 0, 0)
)
##
inits <- list(
    ##tau = 1
    b=0
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel, control = list(log=TRUE))
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 100000#500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)##[100001:500000,]
dim(samples)
length(samples)
##
##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res

##      means      sds  medians     q025     q975
##sd 500.2603 287.7832 500.1425 26.28444 974.5512


## winBUGS output
node   mean	sd	2.5%	  median	97.5%	
sd     7.2166	31.5099	0.276503  0.630709	2.09074	


x0 <- sum(samples==0)
x1 <- sum(samples==1)

x1/(x0+x1)


Rmodel$sd <- Cmodel$sd
Rmodel$sd
Rmodel$tau
Rmodel$calculate()
exp(Rmodel$calculate('sd'))
Rmodel$calculate('S')


Cmodel$sd
Cmodel$tau
exp(Cmodel$calculate('sd'))
Cmodel$calculate('S')

debug(Rmcmc$samplerFunctions[[1]]$run)
Rmcmc$run(10)

samples <- as.matrix(Cmcmc$mvSamples)
length(samples)

samplesPlot(samples, ind=1:1000)

## finding BUGS tau sampling (5)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 8,
    num = c(2,2,2,2,2,2,2,2),
    adj = c(2,8,    1,3,    2,4,     3,5,     4,6,     5,7,      6,8,     1,7),
    weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
    L = 16,
    S = c(0,1,0,1,0,1,0,1)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
 node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	0.874155	0.467162	0.00145191	0.214319	0.792857	1.99685	1	100000

##a = mean^2  / sd^2
##b = mean    / sd^2
## ===> WinBUGS is getting a gamma(1.5,3.5) posterior

a <- 1.5
b <- 3.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))



## finding BUGS tau sampling (4)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 5,
    num = c(1,3,3,3,2),
    adj = c(2,   1,3,4,    2,4,5,    2,3,5,    3,4),
    weights = rep(1,12),
    L = 12,
    S = c(-1, 0, 1, 2, 2)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	0.499872	0.354029	0.00117622	0.0613616	0.419798	1.39341	1	100000

##a = mean^2  / sd^2
##b = mean    / sd^2
## ===> WinBUGS is getting a gamma(1.5,3.5) posterior

a <- 1.5
b <- 3.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))



## finding BUGS tau sampling (3)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 4,
    num = c(1,3,2,2),
    adj = c(2,   1,3,4,    2,4,   2,3),
    weights = c(1,1,1,1,1,1,1,1),
    L = 8,
    S = c(-1, 0, 1, 2)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	0.429886	0.35268	0.00114715	0.031528	0.338512	1.33672	1	100000

##a = mean^2  / sd^2
##b = mean    / sd^2
## ===> WinBUGS is getting a gamma(1.5,3.5) posterior

a <- 1.5
b <- 3.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))


## finding BUGS tau sampling (2)
library(nimble)
##
code <- nimbleCode({
    tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
})
##
data <- list(
    N = 4,
    num = c(1,2,2,1),
    adj = c(2,   1,3,    2,4,   3),
    weights = c(1,1,1,1,1,1),
    L = 6,
    S = c(-1, 0, 1, 2)
)
##
inits <- list(
    tau = 1
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
q025 <- quantile(samples, 0.025)
q975 <- quantile(samples, 0.975)
res <- cbind(means, sds, medians, q025, q975)
res
## a
a <- mean(samples)^2 / var(samples)
## b
b <- mean(samples) / var(samples)
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))
a
b

## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	1.003	0.8229	0.002677	0.07357	0.7899	3.119	1	100000

##a = mean^2  / sd^2
##b = mean    / sd^2
## ===> WinBUGS is getting a gamma(1.5,1.5) posterior

a <- 1.5
b <- 1.5
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))


## finding BUGS tau sampling (1)
library(nimble)
##
code <- nimbleCode({
    ##alpha0 ~ dflat()
    ##sd ~ dunif(0, 1000)
    ##tau <- 1/(sd*sd)
    tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
    ##for(i in 1:N) {
    ##    mu[i] <- S[i] + alpha0
    ##    Y[i] ~ dnorm(mu[i], 1)
    ##}
})
##
data <- list(
    N = 3,
    num = c(1,2,1),
    adj = c(2,   1,3,    2),
    weights = c(1,1,1,1),
    L = 4,
    ##Y = c(1,2,3)
    S = c(-1, 0, 1)
)
##
inits <- list(
    tau = 1
    ##alpha0 = 0,
    ##S = c(0,0,0)
)
##
catCode(code, data, inits, file='~/temp/BUGS.txt')
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','weights','L')], data=data['S'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##
niter <- 500000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[100001:500000,]
dim(samples)
length(samples)
##

##means <- apply(samples, 2, mean)
##sds <- apply(samples, 2, sd)
##medians <- apply(samples, 2, median)
means <- mean(samples)
sds <- sd(samples)
medians <- median(samples)
res <- cbind(means, sds, medians)
res
        means      sds  medians
[1,] 2.517552 1.586795 2.194357

## a
mean(samples)^2 / var(samples)
## b
mean(samples) / var(samples)


a <- 5
b <- 4
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))


## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
tau	1.002	1.008	0.004465	0.02487	0.6956	3.749	1	47200
## ===> WinBUGS is getting a gamma(1,1) posterior
a <- 1
b <- 1
n <- 1e5
x <- rgamma(n,a,b)
mean(x)
sd(x)
median(x)
quantile(x, probs=c(0.025, 0.975))


## preproducible example of problem with weights[k] <- 1
## for Perry
## FILED ON GITHUB 5/31/17
library(nimble)
##
code <- nimbleCode({
    ##for(i in 1:L) {   weights[i] <- 1   }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 3)
    for(i in 1:N) {
	Y[i] ~ dnorm(S[i], 1)
    }
})
##
constants <- list(
    N = 3,
    num = c(1,2,1),
    adj = c(2,   1,3,    2),
    weights = rep(1,4),
    L = 4
)
##
data <- list(
    Y = c(1,2,3)
)
##
inits <- list(
    S = c(0,0,0)
)
##
Rmodel <- nimbleModel(code, constants, data, inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$samplerFunctions[[1]]$componentSamplerFunctions[[1]]$dcar$neighborWeights
Rmcmc$samplerFunctions[[1]]$componentSamplerFunctions[[2]]$dcar$neighborWeights
Rmcmc$samplerFunctions[[1]]$componentSamplerFunctions[[3]]$dcar$neighborWeights

niter <- 100
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)

samples[98:100, ]
## [98,]  0.190857546 -3.875851e-01  0.19672751
## [99,] -0.473204489  1.077377e-01  0.36546680
##[100,] -0.370861918 -8.497980e-02  0.45584172




## NOT AGREEING WITH SAMPLING FOR TAU
## CAR scalar_RW sampler
## no islands
## also sampling tau
library(nimble)
##
code <- nimbleCode({
    alpha0 ~ dflat()
    for(k in 1:L) {
        weights[k] <- 1
    }
    sd ~ dunif(0, 1000)
    tau <- 1/(sd*sd)
    ##tau ~ dunif(0, 1000)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau)
    for(i in 1:N) {
        log(mu[i]) <- alpha0 + S[i]
        Y[i] ~ dpois(mu[i])
    }
})
##
data <- list(
    N = 6,
    num = c(3,4,4,3,2,2),
    adj =     c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
    L = 18,
    Y = c(10,12,12,10,14,10)
)
##
inits <- list(
    alpha0 = 0,
    sd = 1,
    S = c(0,0,0,0,0,0)
)
##
Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','L')], data=data['Y'], inits)
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addMonitors('S','sd','tau')
Rmcmc <- buildMCMC(conf)
##
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 200000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)[50000:200000,]
##
##samplesPlot(samples, 'alpha0')
##samplesPlot(samples, 'tau')
##
means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
medians <- apply(samples, 2, median)
res <- cbind(means, sds, medians)
res


## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
sd	0.3021	0.2928	0.004865	0.01032	0.2264	1.047	501	49500





## CAR scalar_RW, scalar_conjugate, and scalar_postPred samplers
## no islands
## fixed tau = 3
library(nimble)

code <- nimbleCode({
    alpha0 ~ dflat()
    for(k in 1:L) {
        weights[k] <- 1
    }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 3)
    for(i in 1:N) {
        mu[i] <- alpha0 + S[i]
    }
    for(i in 1:2) {
        log(lambda[i]) <- mu[i]
        Y[i] ~ dpois(lambda[i])
    }
    Y[3] ~ dnorm(mu[3], 3)
    ymean4 <- 5*mu[4]
    Y[4] ~ dnorm(ymean4, 7)
    ymean5 <- 2*mu[5]
    Y[5] ~ dnorm(ymean5, 1)
})

data <- list(
    N = 6,
    num = c(3,4,4,3,2,2),
    adj =     c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
    L = 18,
    Y = c(10,12,15,20,24)
)

inits <- list(
    alpha0 = 0,
    S = c(0,0,0,0,0,0)
)

Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','L')], data=data['Y'], inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 100000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
res <- cbind(means, sds)
res
##            means       sds
##S[1]   -1.6762604 0.1988165
##S[2]   -1.3264934 0.1621055
##S[3]    1.8720905 0.2261232
##S[4]   -0.8194848 0.1493609
##S[5]    3.0231445 0.2757768
##S[6]   -1.0729963 0.3532078
##alpha0  4.8462442 0.1313957

## winBUGS output
node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
S[1]	-1.675	0.1776	0.001971	-2.027	-1.674	-1.328	1	10000
S[2]	-1.323	0.1479	0.0013	-1.616	-1.32	-1.036	1	10000
S[3]	1.867	0.2275	0.001976	1.417	1.867	2.316	1	10000
S[4]	-0.8196	0.1365	0.001492	-1.088	-0.819	-0.5559	1	10000
S[5]	3.024	0.276	0.002745	2.479	3.025	3.564	1	10000
S[6]	-1.073	0.3552	0.003718	-1.773	-1.072	-0.3791	1	10000
alpha0	4.846	0.1308	0.001459	4.589	4.847	5.098	1	10000

## agrees!!





## CAR scalar_RW sampler
## no islands
## fixed tau = 3
library(nimble)

code <- nimbleCode({
    alpha0 ~ dflat()
    for(k in 1:L) {
        weights[k] <- 1
    }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 3)
    for(i in 1:N) {
        log(mu[i]) <- alpha0 + S[i]
        Y[i] ~ dpois(mu[i])
    }
})

data <- list(
    N = 6,
    num = c(3,4,4,3,2,2),
    adj =     c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
    L = 18,
    Y = c(10,12,15,20,24,16)
)

inits <- list(
    alpha0 = 0,
    S = c(0,0,0,0,0,0)
)

Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','L')], data=data['Y'], inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 100000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
res <- cbind(means, sds)
res
##              means       sds
##S[1]   -0.265103019 0.1899968
##S[2]   -0.144243381 0.1713349
##S[3]   -0.028341197 0.1699618
##S[4]    0.125589974 0.1722319
##S[5]    0.310280206 0.1774359
##S[6]    0.001817417 0.1920584
##alpha0  2.744200920 0.1039359

## winBUGS output
 
##node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
##S[1]	-0.2619	0.1891	0.001877	-0.6427	-0.2585	0.1017	1	10000
##S[2]	-0.1415	0.1718	0.001695	-0.4899	-0.1373	0.1827	1	10000
##S[3]	-0.0275	0.1687	0.001798	-0.3593	-0.02419	0.2992	1	10000
##S[4]	0.1236	0.1738	0.001715	-0.2203	0.1239	0.4595	1	10000
##S[5]	0.3083	0.1755	0.001711	-0.04451	0.3102	0.6477	1	10000
##S[6]  -0.001004	0.1907	0.001934	-0.3914	0.004113	0.3587	1	10000
##alpha0	2.744	0.103	0.001042	2.535	2.746	2.939	1	10000

## agrees!!





## CAR scalar_RW sampler
## no islands
## fixed tau = 1
library(nimble)

code <- nimbleCode({
    alpha0 ~ dflat()
    for(k in 1:L) {
	weights[k] <- 1
    }
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], 1)
    for(i in 1:N) {
        log(mu[i]) <- alpha0 + S[i]
        Y[i] ~ dpois(mu[i])
    }
})

data <- list(
    N = 6,
    num = c(3,4,4,3,2,2),
    adj =     c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
    ##weights = c(1,1,1,   1,1,1,1,   1,1,1,1,   1,1,1,  1,1,   1,1),
    L = 18,
    Y = c(10,12,15,20,24,16)
)

inits <- list(
    alpha0 = 0,
    S = c(0,0,0,0,0,0)
)

Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','L')], data=data['Y'], inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 100000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)

library(coda)
apply(samples, 2, effectiveSize)

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
res <- cbind(means, sds)
res
##             means       sds
##S[1]   -0.36240674 0.2350176
##S[2]   -0.20235429 0.2146017
##S[3]   -0.03327436 0.2060798
##S[4]    0.19690032 0.1976192
##S[5]    0.38665926 0.1936092
##S[6]    0.01447581 0.2166168
##alpha0  2.72605237 0.1060945

## winBUGS output
##node	 mean	 sd	 MC error	2.5%	median	97.5%	start	sample
##alp0	2.726	0.1055	0.001117	2.514	2.728	2.929	1	10000
##S[1]	-0.3636	0.2341	0.002337	-0.8376	-0.3549	0.08407	1	10000
##S[2]	-0.2063	0.216	0.002338	-0.6496	-0.2007	0.1983	1	10000
##S[3]	-0.0335	0.203	0.001985	-0.4438	-0.0292	0.3531	1	10000
##S[4]	0.1971	0.1987	0.001868	-0.209	0.2003	0.5715	1	10000
##S[5]	0.3891	0.1914	0.00197	0.009187	0.3927	0.7583	1	10000
##S[6]	0.01712	0.2166	0.002037	-0.4136	0.02405	0.4305	1	10000

## agrees!!




## testing default weights in dcar_normal() distribution
library(nimble)

code <- nimbleCode({
    alpha0 ~ dflat()
    ##for(i in 1:L) {
    ##    weights[i] <- 1
    ##}
    S[1:N] ~ car.normal(adj = adj[1:L], num = num[1:N], tau = 1)
    for(i in 1:N) {
        log(mu[i]) <- alpha0 + S[i]
        Y[i] ~ dpois(mu[i])
    }
})
data <- list(
    N = 6,
    num = c(3,4,4,3,2,2),
    adj =     c(2,3,4,   1,3,5,6,   1,2,4,5,   1,3,6,  2,3,   2,4),
    ##weights = rep(1, 18),
    L = 18,
    Y = c(10,12,15,20,24,16)
)
inits <- list(
    alpha0 = 0,
    ##weights = rep(1, 18),
    S = c(0,0,0,0,0,0)
)

Rmodel <- nimbleModel(code, constants=data[c('N','num','adj','L')], data=data['Y'], inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmodel$calculate('S')
Cmodel$calculate('S')

niter <- 100000
set.seed(0); Cmcmc$run(niter)
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
res <- cbind(means, sds)
res
##             means       sds
##S[1]   -0.36240674 0.2350176
##S[2]   -0.20235429 0.2146017
##S[3]   -0.03327436 0.2060798
##S[4]    0.19690032 0.1976192
##S[5]    0.38665926 0.1936092
##S[6]    0.01447581 0.2166168
##alpha0  2.72605237 0.1060945





library(nimble)
##nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
##
## Scottish Lip Cancer example from WinBUGS/geoBUGS of dcar_normal() distribuion
##
## data copied from WinBUGS
N <- 56
E <- c(1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1, 2.3, 2.0, 6.6,
       4.4, 1.8, 1.1, 3.3, 7.8, 4.6, 1.1, 4.2, 5.5, 4.4,
       10.5,22.7, 8.8, 5.6,15.5,12.5, 6.0, 9.0,14.4,10.2,
       4.8, 2.9, 7.0, 8.5,12.3,10.1,12.7, 9.4, 7.2, 5.3,
       18.8,15.8, 4.3,14.6,50.7, 8.2, 5.6, 9.3,88.7,19.6,
       3.4, 3.6, 5.7, 7.0, 4.2, 1.8)
X <- c(16,16,10,24,10,24,10, 7, 7,16, 7,16,10,24, 7,16,10,
       7, 7,10, 7,16,10, 7, 1, 1, 7, 7,10,10,7,24,10, 7, 7,
       0,10, 1,16, 0, 1,16,16, 0, 1, 7, 1, 1, 0, 1,1, 0, 1, 1,16,10)
num <- c(3, 2, 1, 3, 3, 0, 5, 0, 5, 4, 0, 2, 3, 3, 2, 6, 6, 6, 5, 3,
         3, 2, 4, 8, 3, 3, 4, 4, 11, 6, 7, 3, 4, 9, 4, 2, 4, 6, 3, 4, 
         5, 5, 4, 5, 4, 6, 6, 4, 9, 2, 4, 4, 4, 5, 6, 5)
adj <- c(19, 9, 5, 10, 7, 12, 28, 20, 18, 19, 12, 1, 
         17, 16, 13, 10, 2, 29, 23, 19, 17, 1, 22, 16, 7, 2, 
         5, 3, 19, 17, 7, 35, 32, 31, 29, 25, 29, 22, 21, 17,
         10, 7, 29, 19, 16, 13, 9, 7, 56, 55, 33, 28, 20, 4,
         17, 13, 9, 5, 1, 56, 18, 4, 50, 29, 16, 16, 10, 39, 34, 29, 9,
         56, 55, 48, 47, 44, 31, 30, 27, 29, 26, 15, 43, 29, 25,
         56, 32, 31, 24, 45, 33, 18, 4, 50, 43, 34, 26, 25, 23, 21,
         17, 16, 15, 9, 55, 45, 44, 42, 38, 24, 47, 46, 35, 32, 27, 24, 14, 
         31, 27, 14, 55, 45, 28, 18, 54, 52, 51, 43, 42, 40, 39, 29, 23,
         46, 37, 31, 14, 41, 37, 46, 41, 36, 35, 54, 51, 49, 44, 42, 30,
         40, 34, 23, 52, 49, 39, 34, 53, 49, 46, 37, 36, 51, 43, 38, 34, 30,
         42, 34, 29, 26, 49, 48, 38, 30, 24, 55, 33, 30, 28, 53, 47, 41,
         37, 35, 31, 53, 49, 48, 46, 31, 24, 49, 47, 44, 24, 54, 53, 52,
         48, 47, 44, 41, 40, 38, 29, 21, 54, 42, 38, 34, 54, 49, 40, 34,
         49, 47, 46, 41, 52, 51, 49, 38, 34, 56, 45, 33, 30, 24, 18,
         55, 27, 24, 20, 18)
L <- length(adj)
weights <- rep(1, L)
constantsIslands <- list(N=N, L=L, E=E, X=X, num=num, adj=adj, weights=weights)
##
## alternate data with no islands, copied from Breslow and Clatyon (1993)
## no islands in this data
num <- c(4, 2, 2, 3, 5, 2, 5, 1, 6, 4, 4, 3, 4, 3, 3, 6, 6, 6, 5,  ##different!
         3, 3, 2, 6, 8, 3, 4, 4, 4, 11, 6, 7, 4, 4, 9, 5, 4, 5, 6,
         5, 5, 7, 6, 4, 5, 4, 6, 6, 4, 9, 3, 4, 4, 4, 5, 5, 6)
adj <- c(5,9,11, 19, 7,10, 6,12, 18, 20, 28, 1, 11, 12,13,19, 3,  ## different!!
         8, 2,10,13,16,17, 6, 1, 11, 17,19,23,29, 2, 7, 16, 22,
         1,5,9,12, 3, 5, 11, 5, 7,17,19, 31, 32, 35, 25, 29, 50,
         7,10,17,21, 22,29, 7,9,13,16,19,29, 4,20,28,33,55,56, 1,
         5,9,13,17, 4,18,55, 16,29,50, 10,16, 9,29,34,36,37,39, 27,
         30,31,44,47,48,55,56, 15,26,29, 25,29,42,43, 24,31,32,55,
         4,18,33,45, 9, 15, 16, 17, 21, 23, 25, 26, 34, 43, 50, 24,
         38,42,44,45,56, 14,24,27,32,35,46,47, 14,27,31,35, 18,28,
         45,56, 23, 29, 39, 40, 42, 43, 51, 52, 54, 14, 31, 32, 37,
         46, 23, 37, 39, 41, 23, 35, 36, 41, 46, 30,42,44,49,51,54,
         23,34,36,40,41, 34,39,41,49,52, 36,37,39,40,46,49,53, 26,
         30, 34, 38, 43, 51, 26, 29, 34, 42, 24, 30, 38, 48, 49, 28,
         30, 33, 56, 31, 35, 37, 41, 47, 53, 24, 31, 46, 48, 49, 53,
         24,44,47,49, 38, 40, 41, 44, 47, 48, 52, 53, 54, 15,21, 29,
         34, 38, 42, 54, 34,40,49,54, 41, 46, 47, 49, 34, 38, 49,
         51, 52, 18, 20, 24, 27, 56, 18, 24, 30, 33, 45, 55)
L <- length(adj)
weights <- rep(1, L)
constantsNoIslands <- list(N=N, L=L, E=E, X=X, num=num, adj=adj, weights=weights)
##
## observations:
Y <- c(9, 39, 11, 9, 15, 8, 26, 7, 6, 20, 13, 5, 3, 8, 17, 9, 2, 7, 9,
       7, 16, 31, 11, 7, 19, 15, 7, 10, 16, 11, 5, 3, 7, 8, 11, 9, 11,
       8, 6, 4, 10, 8, 2, 6, 19, 3, 2, 3, 28, 6, 1, 1, 1, 1, 0, 0)
data <- list(Y=Y)
##
## model code (simple)
## from geoBUGS manual
## and used in Breslow and Clatyon (1993)
codeSimple <- nimbleCode({
    alpha0  ~ dflat()  
    alpha1 ~ dnorm(0, 0.00001)
    tau ~ dgamma(0.5, 0.0005)
    sigma <- sqrt(1 / tau)
    ##for(k in 1:L)    ## need to make this work!
    ##    weights[k] <- 1
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 1)
    for(i in 1:N) {
        log(mu[i]) <- log(E[i]) + alpha0 + alpha1 * X[i]/10 + S[i]
        Y[i] ~ dpois(mu[i])
        RR[i] <- exp(alpha0 + alpha1 * X[i]/10 + S[i])
    }
})
##initsSimple <- list(tau = 1, alpha0 = 0, alpha1 = 0,    ## need to make this work!!
##                    S=c(0,0,0,0,0,NA,0,NA,0,0,
##                        NA,0,0,0,0,0,0,0,0,0,
##                        0,0,0,0,0,0,0,0,0,0,
##                        0,0,0,0,0,0,0,0,0,0,
##                        0,0,0,0,0,0,0,0,0,0,
##                        0,0,0,0,0,0))
initsSimple <- list(tau = 1, alpha0 = 0, alpha1 = 0, 
                    S=c(0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0))
##
## more complex model code
## also used in geoBUGS manual, includes results for this:
## code from geoBUGS example model
codeComplex <- nimbleCode({
    alpha0 ~ dflat()
    alpha1 ~ dnorm(0, 0.00001)
    tau ~ dgamma(0.5, 0.0005)
    tau.h ~ dgamma(0.5, 0.0005)
    S[1:N] ~ car.normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 1)
    for(i in 1:N) {
        H[i] ~ dnorm(0, tau.h)
        log(mu[i]) <- log(E[i]) + alpha0 + alpha1 * X[i]/10 + S[i] + H[i]
        Y[i] ~ dpois(mu[i])
        RR[i]  <- exp(alpha0 + alpha1 * X[i]/10 + S[i] + H[i])
        residRR[i] <- exp(S[i] + H[i])
        X.pred[i] <- alpha1 * X[i]/10
        lE[i] <- log(E[i])
    }
    rr.x <- exp(alpha1)
    sdS <- sd(S[1:N])
    sdH <- sd(H[1:N])
    sdX <- sd(X.pred[1:N])
    sdE <- sd(lE[1:N])
    sumvar <- sdS^2 + sdH^2 + sdX^2 + sdE^2
    pS <- sdS^2 / sumvar
    pH <- sdH^2 / sumvar
    pX <- sdX^2 / sumvar
    pE <- sdE^2 / sumvar
})
initsComplex <- list(alpha0=0, alpha1=0, tau=1, tau.h=1, S=rep(0,N), H=rep(0,N))

## this mimics Breslow and Clatyon (1993)
Rmodel <- nimbleModel(codeSimple, constantsNoIslands, data, initsSimple)

## this mimics geoBUGS manual
Rmodel <- nimbleModel(codeComplex, constantsIslands, data, initsComplex)

Rmodel$calculate()

conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$printMonitors()
conf$addMonitors('sigma')      ## for (B&C 1993) model
conf$addMonitors('sigma', 'S')
conf$addMonitors('sigma', 'S', 'RR')
conf$addMonitors('rr.x', 'pS', 'pH', 'pX', 'pE')        ## for geoBUGS manual model
conf$addMonitors('rr.x', 'pS', 'pH', 'pX', 'pE', 'RR')

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 300000
nburnin <- 50000
set.seed(0); system.time(Cmcmc$run(niter))   ## ~30 sec (B&C 1993 model), ~45 sec (geoBUGS manual model)
samples <- as.matrix(Cmcmc$mvSamples)[(nburnin+1):niter, ]
dim(samples)
colnames(samples)

library(coda)
apply(samples, 2, effectiveSize)

means <- apply(samples, 2, mean)
sds <- apply(samples, 2, sd)
res <- cbind(means, sds)

res
res[c('alpha0','alpha1','sigma'), ]
res[c('alpha0','alpha1','sigma','tau'), ]
res[c('pE','pH','pS','pX','rr.x'), ]

## this mimics Breslow and Clatyon (1993)
## their results:
## alpha0: -0.18 (sd = 0.12)
## alpha1: -0.35 (sd = 0.12)
## sigma:   0.73 (sd = 0.13)
## NIMBLE:
##            means       sds
##alpha0 -0.2099181 0.1201011
##alpha1  0.3582498 0.1279591
##sigma   0.7176594 0.1239661

## this mimics geoBUGS manual
## their results:
##pE:   0.6238 (sd = 0.0367)
##pH:   0.0182 (sd = 0.0308)
##pS:   0.2718 (sd = 0.0549)
##pX:   0.0862 (sd = 0.0407)
##rr.x: 1.587  (sd = 0.1907)
## NIMBLE:
##           means        sds
##pE   0.600743362 0.03484688
##pH   0.005349068 0.01306911
##pS   0.320298241 0.05209477
##pX   0.073609329 0.04107101
##rr.x 1.537977214 0.20651663


## hand-implement the "ranked" statistics to verify against BUGS book
colnames(samples)
rrnames <- paste0('RR[', 1:56, ']')
rrsamples <- samples[, rrnames]
colnames(rrsamples)
dim(rrsamples)
apply(rrsamples, 2, effectiveSize)
 
a <- rrsamples[50001:100000,]
low <- 6
high <- 51
out <- array(NA, c(dim(a)[1], 2))
 
for(i in 1:dim(a)[1]) {
    vals <- a[i,]
    sor <- sort(vals)
    out[i,1] <- sor[low]
    out[i,2] <- sor[high]
}
 
ratio <- out[,2] / out[,1]
mean(ratio)
median(ratio)
sd(ratio)







## test models for dcar_normal distribution

library(nimble)
N <- 15
num <- c(9,    rep(1,9), 0,0,0,0,0)
adj <- c(2:10, rep(1,9))
L <- length(adj)
weights <- rep(1, L)
code <- nimbleCode({
    t ~ dunif(0, 100)
    x[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau=t)
    y1 ~ dexp(x[1])
    y2 ~ dnorm(x[2], 1)
    z <- x[3] + 3
    y3 ~ dnorm(z, 1)
    y4 ~ dnorm(x[4]*x[4], 1)
    y51 ~ dnorm(x[5]+5, 1)
    y52 ~ dnorm(x[5]/2, 2)
    for(i in 1:9) {
        mu[i] <- i
    }
    mu[10] <- x[6]
    y6[1:10] ~ dmnorm(mu[1:10], C[1:10,1:10])
})
constants <- list(N = N, L = L, adj = adj, weights = weights, num = num, C = diag(10))
data <- list(y1=1, y2=1, y3=1, y4=1, y51=1, y52=1, y6=rep(0,10))
inits <- list(x = rep(1,N), t = 1)
Rmodel <- nimbleModel(code, constants, data, inits)


library(nimble)
N <- 5
num <- c(0, 0, 0, 1, 1)
adj <- c(5, 4)
L <- length(adj)
weights <- rep(1, L)
code <- nimbleCode({
    t ~ dunif(0, 100)
    x[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau = t)
    y2 ~ dnorm(x[2], 1)
    y3 ~ dexp(x[3])
})
constants <- list(N = N, L = L, adj = adj, weights = weights, num = num)
data <- list(y2 = 2, y3 = 3)
inits <- list(x = 1:N, t = 1)
Rmodel <- nimbleModel(code, constants, data, inits)



conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, control=list(adaptInterval=1000000))
##conf <- configureMCMC(Rmodel, control=list(carUseConjugacy=FALSE))
##conf$printSamplers(displayControlDefaults = TRUE)
conf$printSamplers()
##conf$addMonitors('x')
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE)

Rmodel$calculate('x')
Cmodel$calculate('x')
## -8.270447
## -1.837877  (3 islands model)
## -47.7848   (Scottish lip cancer)

niter <- 10
set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
sampnames <- dimnames(Rsamples)[[2]]
Rsamples[, sampnames]
Csamples[, sampnames]
Rsamples[, sampnames] - Csamples[, sampnames]

niter <- 100000
##niter <- 10
##set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)

samples <- as.matrix(Cmcmc$mvSamples)

samples[niter,]
##         t       x[1]       x[2]       x[3]       x[4]       x[5]       x[6] 
## 0.1287901  0.5934936 -0.4500610 -1.7883696 -1.6046201 -0.9227700  0.6587114 
##      x[7]       x[8]       x[9]      x[10] 
## 4.7227748 -0.3847551  0.7452862  2.3125083 

apply(samples, 2, mean)
apply(samples, 2, sd)

samplesPlot(samples, 1)
samplesPlot(samples, 2:6)
samplesPlot(samples, c(3,4))

library(coda)
apply(samples, 2, effectiveSize)




> apply(samples, 2, mean)
          t        x[1]        x[2]        x[3]        x[4]        x[5] 
 0.10855950  1.15413864  0.99105430 -1.73847094  0.08252756 -1.81795004 
       x[6]        x[7]        x[8]        x[9]       x[10] 
 0.08001419  1.15626026  1.14893480  1.16751225  1.13389828 

> apply(samples, 2, sd)
        t      x[1]      x[2]      x[3]      x[4]      x[5]      x[6]      x[7] 
0.1131858 0.8444364 0.9587355 0.9704785 0.9258072 0.8092274 0.9547935 5.8191558 
     x[8]      x[9]     x[10] 
5.6441207 5.6973708 5.8381872 

> apply(samples, 2, effectiveSize)
        t      x[1]      x[2]      x[3]      x[4]      x[5]      x[6]      x[7] 
 2106.414  7062.471 93635.441 44922.136 20619.570 48727.286 22121.622 87259.953 
     x[8]      x[9]     x[10] 
79181.284 64071.703 78410.282 




## example of summing 0-length vector in nimble
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)

nfDef <- nimbleFunction(
    setup = function(model, nodes) {},
    run = function() {
        a <- values(model, nodes)
        b <- sum(a)
    }
)

Rnf <- nfDef(Rmodel, character(0))

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)

Rnf$run()
Cnf$run()





## timing comparison of named vs. indexed lookup for named vectors
N <- 5
n <- 10^N
vec <- 1:n
nn <- paste0('x', 1:n)
names(vec) <- nn
##system.time(for(i in 1:n) { thisname <- paste0('x', i); out <- vec[thisname] })
system.time(for(i in 1:n) { thisname <- paste0('x', i); out <- vec[i] })   ## WAYYYYY FASTER !!!!!!!



## making conditional code block for samplerAssignmentRules

library(nimble)
nimbleOptions('MCMCuseSamplerAssignmentRules')
nimbleOptions(MCMCuseSamplerAssignmentRules = TRUE)
nimbleOptions('MCMCuseSamplerAssignmentRules')

rules <- nimbleOptions('MCMCdefaultSamplerAssignmentRules')
rules


code <- nimbleCode({
  a ~ dnorm(0, rs)
  b <- a + 1
  c ~ dnorm(b, 1)
  d <- c-1
})
constants <- list()
data <- list()
inits <- list(a = 0, c=0, rs=1)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)

conf$printSamplers()




## very quick test of compilation of dcar_normal(), for Chris P.
library(nimble)
code <- nimbleCode({
    x[1:3] ~ dcar_normal(adj[1:3], weights[1:3], num[1:3], t)
})
constants <- list(adj = 1:3, weights = 1:3, num = 1:3, t = 0)
data <- list()
inits <- list()


Rmodel <- nimbleModel(code, constants, data, inits, calculate=FALSE)

Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)




## bug in getDependencies(..., self=FALSE) !
library(nimble)
code <- nimbleCode({
    x[1:4] ~ dmnorm(mu[1:4], C[1:4,1:4])
    y ~ dexp(x[1])
})
constants <- list(mu = rep(0,4), C = diag(4))
Rmodel <- nimbleModel(code, constants)

Rmodel$getDependencies('x[1:4]')
Rmodel$getDependencies('x[1:4]', self=FALSE)
Rmodel$getDependencies('x[1:4]', returnScalarComponents = TRUE)
Rmodel$getDependencies('x[1:4]', self=FALSE, returnScalarComponents = TRUE)

Rmodel$getDependencies('x[1]')
Rmodel$getDependencies('x[1]', self=FALSE)
Rmodel$getDependencies('x[1]', returnScalarComponents = TRUE)
Rmodel$getDependencies('x[1]', self=FALSE, returnScalarComponents = TRUE)


## looking at how it always lifts chol()
nms <- ls(Rmodel$nodes)
i <- 2
ls(Rmodel$nodes[[nms[i]]])
Rmodel$nodes[[nms[i]]]$calculate
Rmodel$nodes[[nms[i]]]$simulate

## bug in model building!
library(nimble)
code <- nimbleCode({
    x[1:4] ~ dmnorm(mu[1:4], C[1:4,1:4])
    y[1] ~ dnorm(x[2], 1)
})
constants <- list(mu = rep(0,4), C = diag(4))
Rmodel <- nimbleModel(code, constants)




## testing speed comparisons of new sampler assignment rules system


library(nimble)
nimbleOptions('MCMCuseSamplerAssignmentRules')

nimbleOptions(MCMCuseSamplerAssignmentRules = FALSE)
nimbleOptions(MCMCuseSamplerAssignmentRules = TRUE)

nimbleOptions('MCMCuseSamplerAssignmentRules')

timer <- function(n) {
    N <- 10^n
    library(nimble)
    code <- nimbleCode({
        for(i in 1:N) {
            a[i] ~ dnorm(0, 1)
            b[i] ~ dexp(a[i])
        }
    })
    constants <- list(N = N)
    data <- list(b = rep(10,N))
    inits <- list(a = rep(0,N))
    Rmodel <- nimbleModel(code, constants, data, inits)
    t <- system.time(conf <- configureMCMC(Rmodel))
    ##conf$printSamplers()
    return(as.numeric(t[3]))
}

sapply(1:4, function(n) timer(n))


## old system:
old <- c(0.049, 0.255, 2.986, 21.344)
## new sampler assignment rules:
new <- c(0.061, 0.460, 5.660, 52.530)   ## original
0.078  0.505  5.552 49.853  ## ruleSelectionCodeBlock, and only 1 eval()
0.066  0.474  5.688 49.230  ##
0.164  0.430  4.555 43.934  ## with isEndNode[], and nodeDistributions[] defined
0.139  0.376  3.271 32.278  ## is isEndNode[i], isBinary[i], isMultivariate[i], isConjugate[i]
0.034  0.276  3.316 35.299
0.082  0.393  3.625 35.707
0.105  0.287  3.261 32.724  ## ruleSelectFunction()
0.033  0.282  3.384 33.004

new/old

plot(1:4, new, col='red', type='l')
lines(1:4, old)


## my attempt at profiling the time for samplerAssignmentRules
?Rprof
?summaryRprof

N <- 5000
code <- nimbleCode({
    for(i in 1:N) {
        a[i] ~ dnorm(0, 1)
        b[i] ~ dexp(a[i])
    }
})
constants <- list(N = N)
data <- list(b = rep(10,N))
inits <- list(a = rep(0,N))
Rmodel <- nimbleModel(code, constants, data, inits)

profFile <- '~/temp/prof.txt'
Rprof(profFile)

system.time(conf <- configureMCMC(Rmodel))

Rprof(NULL)

profFile <- '~/temp/prof.txt'
summaryRprof(profFile)



## table of quantiles of t-distribution

dfs <- c(1,2,5,10,100,1000)
coverages <- c(.9, .95, .99)

ps <- 1-(1-coverages)/2
dfarg <- rep(dfs, each=length(coverages))
probarg <- rep(ps, length(dfs))
tvals <- qt(probarg, dfarg)
tab <- t(matrix(tvals, nrow=length(coverages), ncol=length(dfs)))
dimnames(tab) <- list(`df` = dfs, `CI Coverage` = paste0(100*coverages,'%'))
tab





## testing that nimble MCMC doesn't put samplers on deterministic nodes

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, rs)
    b <- a + 1
    c ~ dnorm(b, 1)
    d <- c-1
})
constants <- list()
data <- list()
inits <- list(a = 0, c=0, rs=1)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$addSampler('a', 'RW')
conf$addSampler('b', 'RW')
conf$addSampler('d', 'RW')
conf$addSampler('rs', 'RW')
conf$addSampler(c('a','rs'), 'RW_block')

Rmodel$getNodeNames()

conf <- configureMCMC(Rmodel, nodes = c('a', 'b', 'c', 'rs'))
conf$printSamplers()
conf <- configureMCMC(Rmodel, nodes = Rmodel$getNodeNames(includeRHSonly=TRUE))


Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## example of Central Limit Theorem (CLT) for STAT201
unifMeans <- function(n, rep=500, breaks=100, curve=TRUE) {
    um <- replicate(rep, mean(runif(n)))
    hist(um, xlim=c(0,1), breaks=breaks, prob=TRUE, main=paste0('Uniform sample mean (n = ', n, ')'))
    if(curve) curve(dnorm(x, mean(um), sd(um)), col='red', add=TRUE)
}

n <- 4
par(mfrow=c(n,1), mar=c(2,2,2,2))
unifMeans(1, curve=FALSE)
##unifMeans(2)
##unifMeans(3)
unifMeans(4)
unifMeans(16)
unifMeans(64)





## playing with and testing new MCMC sampler assignment rules
library(nimble)
getNimbleOption('MCMCdefaultSamplerAssignmentRules')

code <- nimbleCode({
    for(i in 1:10) {
        a[i] ~ dnorm(0, 1)
    }
    for(i in 1:5) {
        b[i] ~ dexp(a[i] + 1)
    }
    x[1:10] ~ dmnorm(mu[1:10], C[1:10,1:10])
    y ~ dexp(x[1])
})
constants <- list(mu = rep(0,10), C = diag(10))
Rmodel <- nimbleModel(code, constants = constants)

conf <- configureMCMC(Rmodel)
conf$printSamplers()

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b1 ~ dnorm(0, 1)
    b2 ~ dnorm(0, 1)
    c1 ~ dnorm(0, 1)
    c2 ~ dnorm(0, 1)
    c3 ~ dnorm(0, 1)
    d ~ dnorm(0, 1)
    e1 ~ dnorm(0, 1)
    e2 ~ dnorm(0, 1)
    f ~ dnorm(0, 1)
    g ~ dnorm(0, 1)
    h ~ dnorm(0, 1)
    i ~ dnorm(0, 1)
    j ~ dnorm(0, 1)
    k1 ~ dnorm(0, 1)
    k2 ~ dnorm(0, 1)
    l1 ~ dnorm(0, 1)
    l2 ~ dnorm(0, 1)
    l3 ~ dnorm(0, 1)
    z1 ~ dnorm(0, 1)
    z2 ~ dnorm(0, 1)
})
Rmodel <- nimbleModel(code)

a_sampler <- nimbleFunction(setup = function() {}, run = function() 2)
x_sampler <- nimbleFunction(setup = function() {}, run = function() 3)

em <- samplerAssignmentRules(empty=TRUE, print=TRUE)
em$addRule(quote(grepl('^a', node)), 'a_sampler')
em$addRule(quote(grepl('^b1', node)), a_sampler)
em$addRule(quote(grepl('^b2', node)), a_sampler, name = 'sampler_for_node_b')
em$addRule(quote(grepl('^c1', node)), nimbleFunction(setup = function() {}, run = function() 4))
ll <- list()
ll[[1]] <- nimbleFunction(setup = function() {}, run = function() 5)
em$addRule(quote(grepl('^c2', node)), ll[[1]])
em$addRule(quote(grepl('^c3', node)), ll[[1]], name = 'sampler_for_c3')
em$addRule(quote(grepl('^d', node)), 'RW')
em$addRule(quote(grepl('^e1', node)), 'sampler_RW')
em$addRule(quote(grepl('^e2', node)), 'sampler_RW', name = 'sampler_for_e2')
em$addRule(quote(grepl('^f', node)), sampler_RW)
em$addRule(quote(grepl('^g', node)), quote(addSampler(target=node, type=sampler_RW)))
em$addRule(quote(grepl('^h', node)), quote(addSampler(target=node, type='RW')))
em$addRule(quote(grepl('^i', node)), quote(addSampler(target=node, type='sampler_RW')))
em$addRule(quote(grepl('^j', node)), quote(addSampler(target=node, type='x_sampler')))
em$addRule(quote(grepl('^k1', node)), quote(addSampler(target=node, type=x_sampler)))
em$addRule(quote(grepl('^k2', node)), quote(addSampler(target=node, type=x_sampler, name='sampler_for_k2')))
em$addRule(quote(grepl('^l1', node)), quote(addSampler(target=node, type=nimbleFunction(setup = function() {}, run = function() 6))))
em$addRule(
    quote(grepl('^l2', node)),
    quote(addSampler(
        target=node,
        type=nimbleFunction(setup = function() {}, run = function() 6),
        name='sampler_for_node_l2')))
em$addRule(
    quote(grepl('^l3', node)),
    quote(addSampler(
        target=node,
        type=nimbleFunction(setup = function() {}, run = function() 6))),
          name='sampler_for_node_l3')
em

conf <- configureMCMC(Rmodel, rules = em)
conf$printSamplers()

conf <- configureMCMC(Rmodel, rules = em, warnNoSamplerAssigned = FALSE)
warnings()
conf$printSamplers()

nimbleOptions(MCMCdefaultSamplerAssignmentRules = em)
getNimbleOption('MCMCdefaultSamplerAssignmentRules')

conf <- configureMCMC(Rmodel)
warnings()
conf$printSamplers()

nimbleOptions(MCMCdefaultSamplerAssignmentRules = samplerAssignmentRules())
getNimbleOption('MCMCdefaultSamplerAssignmentRules')

conf <- configureMCMC(Rmodel)
conf$printSamplers()


samplerAssignmentRules()

defaults$ruleList[[3]]
defaults$printRules(3)
defaults$printRules(c(3,5))

defaults$addRule(1, 5, print=TRUE, position=12)
defaults
defaults$orderRules(1:5, print=TRUE)

defaults$addRule(quote({a+b;TRUE;4}), 5, print=TRUE, position=3)
defaults$addRule(quote({a+b;TRUE;4}), quote({for(i in 1:10) {a + b; TRUE}}), print=TRUE)
defaults

defaults$addRule(TRUE, mean, print=TRUE)
defaults$addRule(TRUE, sampler_RW, print=TRUE)

defaults$reorder(c(1,4,5), print=TRUE)
defaults


defaults$reorder(c(1,3,2))
defaults$reorder(2:12)
defaults$reorder(c(1,4,0,3))
defaults$reorder(c(1,4,5.5,3))
defaults$addRule(quote(grepl('^a', node)), 'a_sampler', position=1, print=TRUE)


em <- samplerAssignmentRules(empty=TRUE, print=TRUE)
em$addRule(quote(grepl('^a', node)), 'a_sampler', position=1, print=TRUE)

em <- samplerAssignmentRules(empty=TRUE, print=TRUE)
em$addRule(quote(grepl('^a', node)), a_sampler, position=1, print=TRUE)

em$printRules()

conf <- configureMCMC(Rmodel, rules = em)
conf$printSamplers()


conf$printSamplers()


Rmodel$isMultivariate('a[1]')
Rmodel$isMultivariate('a[1:2]')
Rmodel$isMultivariate('x')
Rmodel$isMultivariate('x[1:5]')
Rmodel$isMultivariate(c('x[1:5]', 'x', 'a'))


## bug in copying models?

library(nimble)

code <- nimbleCode({
    for(i in 1:N) {
        x[i] ~ dnorm(0, 1)
    }
})
constants <- list(N=10)

def <- nimbleModel(code, constants, returnDef = TRUE)

Rmodel1 <- def$newModel(check = FALSE)

Rmodel2 <- def$newModel(check = FALSE)

Rmodel1$newModel(check=FALSE)
                 

## read and XML document

library(XML)

url <- 'https://www.w3schools.com/xml/plant_catalog.xml'
url

?xmlTreeParse
formals(xmlTreeParse)
args(xmlTreeParse)

xmlfile <- xmlTreeParse(url)
xmlfile <- xmlTreeParse(url, isURL=TRUE)

                                        # the xml file is now saved as an object you can easily work with in R:
class(xmlfile)
# Use the xmlRoot-function to access the top node
xmltop = xmlRoot(xmlfile)
# have a look at the XML-code of the first subnodes:
print(xmltop)[1:2]


## error in getDependencies

library(nimble)

code <- nimbleCode({
    x[1:10] ~ dmnorm(mu[1:10], C[1:10,1:10])
    for(i in 1:10) {
        y[i] ~ dnorm(x[i], 1)
    }
})

constants <- list(mu = rep(0,10), C = diag(10))
data <- list()
inits <- list(x = rep(0,10), y = rep(0,10))

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$getDependencies('x[5]')

Rmodel$calculate('x')

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)

samplesPlot(samples)

library(coda)
apply(samples, 2, effectiveSize)






## testing output of printSamplers() and getSamplers()

library(nimble)

N <- 10
code <- nimbleCode({
    mu[1:N] ~ dmnorm(ones[1:N], C[1:N,1:N])
    a[1:N] ~ dmnorm(mu[1:N], C[1:N,1:N])
    x ~ dnorm(z, 1)
    y ~ dnorm(x, 1)
    z ~ dgamma(1, 1)
    for(i in 1:N) {
        qq[i] ~ dnorm(a[i], 1)
    }
})
constants <- list(N=N, C=diag(N), ones=rep(0,N))
data <- list(y = 0)
inits <- list(a = rep(0,N), x=0, qq=rep(0,N))
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)

conf$printSamplers()
conf$getSamplerDefinition(3)

conf$addSampler('a', 'RW_block', control=list(scale=2, adaptInterval=500, propCov=10*diag(N)))
conf$addSampler('z', 'RW', control=list(scale=3))
conf$addSampler('z', 'RW', control=list(scale=1:10))
conf$addSampler('z', 'RW', control=list(scale=c(1,4,6,7)))

conf$printSamplers()
conf$printSamplers(displayConjugateDep = TRUE)
conf$printSamplers(displayControlDefaults = TRUE)
conf$printSamplers(displayNonScalars = TRUE)
conf$getSamplers()

conf$samplerConfs[[1]]
class(conf$samplerConfs[[1]])


## testing using if() statements inside BUGS code

library(nimble)

code <- nimbleCode({
    if(orderx == 1) {
        y ~ dnorm(0,1)
    } else {
        y ~ dpois(1)
    }
})

orderx <- 1

m <- nimbleModel(code, constants=list(orderx=1))


## demo of using splines
library(splines)

data <- data.frame(value = c(0.01590000, 0.04040333, 0.08487666,
                             0.36868000, 0.31937660),
                   depth = c(10, 25, 50, 100, 180),
                   t = 1:5)

value.sp <- interpSpline(value ~ t, data)
depth.sp <- interpSpline(depth ~ t, data)
class(value.sp)
value.sp
plot(value.sp)
plot(depth ~ value, data)
lines(predict(value.sp)$y, predict(depth.sp)$y)

## demo of coin flips and LLN

set.seed(0)
n <- 15
flips <- rbinom(n, 1, 0.5)
flips

cumsum(flips)

cumsum(flips) / 1:n

running_avg <- cumsum(flips) / 1:n
running_avg

par(mfrow=c(1,1))
plot(running_avg, type='l')
abline(h=0.5, col='red', lty=3)

ks <- 2:4
par(mfrow=c(length(ks),1), mar=c(2,2,2,2))
for(i in ks) {
    set.seed(0)
    n <- 10^i
    flips <- rbinom(n, 1, 0.5)
    running_avg <- cumsum(flips) / seq_along(flips)
    plot(running_avg, type='l')
    abline(h=0.5, col='red', lty=3)
}



## testing use of rep() in NIMBLE DSL code

library(nimble)

code <- nimbleCode({
    for(i in 1:5) {
        a[i] ~ dnorm(0, 1)
    }
})
constants <- list()
data <- list()
inits <- list(a = 1:5)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

nfDef <- nimbleFunction(
    setup = function(model) {
        d <- 5
        node <- 'a[1:5]'
    },
    run = function(a = double()) {
        model[[node]] <<- rep(a, d)
    }
)

Rnf <- nfDef(Rmodel)

Cnf <- compileNimble(Rnf, project=Rmodel)

Rmodel$a
Rnf$run(3)
Rmodel$a

Cmodel$a
Cnf$run(3)
Cmodel$a



## looking at Hurricane Damage data file for STAT201 takehome midterm

df <- read.csv("~/github/courses/stat201/data/HurricaneDamage.csv")
head(df)
dim(df)
df
with(df, plot(Year, Damage))
i <- which.max(df$Damage)
i
df[i,]

with(df, plot(Year[-1], Damage[-1]))


## trying to figure out what's meant by github issue for setData()

library(nimble)

code <- nimbleCode({
    for(i in 1:3) {
        a[i] ~ dnorm(0,1)
    }
})

constants <- list()
data <- list()
inits <- list(a = 1:3)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$a
Rmodel$getNodeNames(dataOnly=TRUE)

Rmodel$setData('a')

Rmodel$getNodeNames(dataOnly=TRUE)

1



## examining why the slope is *not* inverted
## when you switch x & y

set.seed(0)
x <- 1:200
y <- rnorm(x, 2*x, 100)#200)
x <- x-mean(x)
y <- y-mean(y)
lim <- c(-400, 400)
plot(x,y, pch=20, xlim=lim, ylim=lim)
points(y,x, pch=20, col='red')
abline(a=0,b=1,col='blue')
m1 <- lm(y~x)
m2 <- lm(x~y)
b1 <- m1$coef[2]
b2 <- m2$coef[2]
abline(m1)
abline(m2, col='red')
b1
1/b2
b1*b2

sd(x)
sd(y)

b1
b2




## STAT 201 class survey dataset

survey <- read.csv("~/Downloads/STAT201ClassInfo.csv")

names(survey)

attach(survey)
plot(Height, Minutes)
plot(Height, Friends)
cor(Height, Friends)
lm(Height ~ Friends)

str(survey)

survey$Siblings
table(survey$Siblings)
prop.table(table(survey$Siblings))


## warning messing in NIMBLE for forgetting type declaration

library(nimble)

## forgot type declaration.  This case is caught
nfDef <- nimbleFunction(
    setup = TRUE,
    run = function(a) { }
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)  # useful warning message:
 ## compiling... this may take a minute. Use 'showCompilerOutput = TRUE' to see C++ compiler details.
## Error: Type declaration missing for argument(s) a

## forgot type declaration, but gave a default value.
## this case gives a totally inconprehensible error
nfDef <- nimbleFunction(
    setup = TRUE,
    run = function(a = 1) { }
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)  # would love a good error message here:
compiling... this may take a minute. Use 'showCompilerOutput = TRUE' to see C++ compiler details.
Error in AT$default : $ operator is invalid for atomic vectors




## playing with calculating r^2

a <- 1
b <- 0.5
sigma <- 1
n <- 10
x <- 1:n
y <- rnorm(n, a+b*x, sigma)
plot(x,y)

r <- cor(x,y)
r^2

m <- lm(y~x)
summary(m)
r^2


yhat <- m$fitted.values
ybar <- mean(y)

1 - sum((y-yhat)^2)/sum((y-ybar)^2)
r^2




#### NIMBLE occupancy model for ZIB GLMM model of potato psyllid occupancy
#### Originally developed by Daniel Turek
## trouble-shooting his issues with NIMBLE v0.6-3
library(nimble)

dCustom <- nimbleFunction(
    run = function(x = double(), a = double(), log = integer(0, default=0)) {
        returnType(double())
        if(log) return(0) else return(1)
    }
)
rCustom <- nimbleFunction(
    run = function(n = integer(), a = double()) {
        returnType(double())
        return(0)
    }
)
registerDistributions(list(
    dCustom = list(
        BUGSdist = "dCustom(a)",
        discrete = TRUE
    )
))

code <- nimbleCode({
    a ~ dnorm(0, 1)
    y ~ dCustom(0.5)
})
constants <- list()
data <- list()
inits <- list(y = 0, a=0)
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)




library(nimble)
setwd("~/Downloads")

inputData <- readRDS('output/data_nimble_zib.rds')
source('R_functions/nimble_definitions.R')



code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:9) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[4]*aet[i] + beta[5]*tmn[i] + beta[6]*tmx[i] + beta[7]*year[i] + beta[8]*month[i] + beta[9]*month2[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*year_list_length[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       aet=aet, tmn=tmn, tmx=tmx, 
                       year=year, 
                       month=month,
                       month2=month2,
                       list_length=list_length, 
                       year_list_length=year_list_length, 
                       siteID=siteID))

data <- with(inputData, list(y=y))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,9))#, betaseason=rep(0,4))

modelInfo_month <- list(code=code, constants=constants, data=data, inits=inits, name='month_model')

Rmodel <- nimbleModel(modelInfo_month$code,
                      modelInfo_month$constants,
                      modelInfo_month$data,
                      modelInfo_month$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:9]')
spec$addSampler('beta[1:3]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[4:9]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors(c('p_occ')) # add a monitor to get p_occ in output
#spec$addMonitors(c('p_obs')) # add a monitor to get p_obs in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 15000
burnin <- 5000

ti <- Sys.time()
samplesList <- lapply(3, mcmcClusterFunction)
tf <- Sys.time()

# The time it took to run MCMC
tf-ti




## helping to debug Dao's autoAdapt algorithm.
## NOTE: need to build and install package from his codess branch
## main function is MCMC_autoAdapt.R / autoAdapt() .... (I think??)

library(nimble)
setwd('~/Downloads')
load('./data/model_spatial.RData')
Rmodel <- nimbleModel(code, constants, data, inits)

CandidateSamplerList <- list(
#  sampler_conjugate= list(type = "sampler_conjugate", name = "conjugate"),
  sampler_RW = list(type = "sampler_RW", control =  list(adaptive = TRUE,adaptInterval = 20), name = "adaptive"),
  sampler_RWlog = list(type = "sampler_RW", control =  list(adaptive = TRUE,adaptInterval = 20,  log =TRUE), name = "adaptivelog"),
  sampler_slice = list(type = "sampler_slice", control =  list(adaptive = TRUE,adaptInterval = 20), name = "Myslice"),
  sampler_RW_block = list(type = "sampler_RW_block", control =  (list(adaptive = TRUE, adaptInterval = 20)), name = "block_RW"),
  sampler_RW_block = list(type = "sampler_AF_slice", control =  (list(adaptive = TRUE, adaptInterval = 20)), name = "block_AFS"),
  sampler_RW_block = list(type = "RW_rotated_block", control =  (list(adaptive = TRUE, adaptInterval = 20)), name = "block_Rotated")
)
##
Candidates<-list(2,3,4,5,6,7)  
##
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$getSamplers()
##
monitor <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
monitor
##
n =length(monitor)
print(n)            
##
DefaultSamplerList <- vector(mode="list", length=n)
DefaultSamplerList
##
for(i in 1: n){
  names(DefaultSamplerList)[i]<-'RW'
  DefaultSamplerList[[i]]$name <-monitor[i]
  DefaultSamplerList[[i]]$target <-monitor[i]
  DefaultSamplerList[[i]]$control <-list()
  DefaultSamplerList[[i]]$type <-'sampler_RW'
  DefaultSamplerList[[i]]$oldtype <-DefaultSamplerList[[i]]$type
}
##
DefaultSamplerList[1:5]
##
DefaultSamplerList
Candidates
monitor


debug(autoAdapt)
debug(autoAdaptClass)
debug(ab$run)

ab <- autoAdapt(Rmodel, niter=1000, DefaultSamplerList=DefaultSamplerList, Candidates = Candidates, monitor=monitor, iteration =5)




## tobit model for Chris, Ann Raiho <ann.raiho@gmail.com>, Mike Dietze <dietze@bu.edu>
## experimenting with "turning off samplers" in the MCMC

## Daniel, Perry, the idea is that in the following code, 'y.ind' is an indicator vector that says whether each element of a data vector is positive or zero.  Then in the BUGS code, y.censored is either the actual data value, or if the data value is 0, y.censored[i] is NA.  In that case, we want to do MCMC sampling on that element of y.censored using a univariate sampler (we recognize that y.censored has a dmnorm distribution -- in this case we only want to update the values for which the corresponding data values are 0).
## When building the model, y.ind is flagged as 'data', y.censored is not flagged as 'data'. The sampler on y.censored is removed, and then a univariate sampler is assigned to each element of y.censored that contains an NA (which corresponds to elements of y.ind that are 0).
## This same model is used multiple times during a data assimilation, with the input data vector (and therefore the y.ind and y.censored vectors) changing at each time. Goal is to using the same compiled MCMC at each time, but dynamically controlling which components of y.censored are sampled.


library(nimble)

sampler_toggle <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        type <- control$type
        nested_sampler_name <- paste0('sampler_', type)
        control_new <- nimbleOptions('MCMCcontrolDefaultList')
        control_new[[names(control)]] <- control
        nested_sampler_list <- nimbleFunctionList(sampler_BASE)
        nested_sampler_list[[1]] <- do.call(nested_sampler_name, list(model, mvSaved, target, control_new))
        toggle <- 1
    },
    run = function() {
        if(toggle == 1)
            nested_sampler_list[[1]]$run()
    },
    methods = list(
        reset = function()
            nested_sampler_list[[1]]$reset()
    )
)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a*a + 1, 2)
    c ~ dnorm(a + b*b, 5)
})
constants <- list()
data <- list()
inits <- list(a = 0, b=0, c=0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
conf$addMonitors('a','b','c')
conf$printMonitors()
conf$printSamplers()

conf$addSampler('a', 'RW')
conf$addSampler('b', 'RW')
conf$addSampler('c', 'RW')

conf$addSampler('c', 'slice')

conf$addSampler('a', 'toggle', control = list(type='RW'))
conf$addSampler('b', 'toggle', control = list(type='RW'))
conf$addSampler('c', 'toggle', control = list(type='RW'))

conf$addSampler('c', 'toggle', control = list(type='slice'))

conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$samplerFunctions$contentsList[[3]]$toggle
valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[3]], 'toggle')
Rmcmc$samplerFunctions$contentsList[[3]]$toggle <- 0
valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[3]], 'toggle', 0)
Rmcmc$samplerFunctions$contentsList[[3]]$toggle
valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[3]], 'toggle')

niter <- 50
set.seed(0); Rsamples <- runMCMC(Rmcmc, niter)
set.seed(0); Csamples <- runMCMC(Cmcmc, niter)
Rsamples
Csamples
Rsamples - Csamples

saveSamples <- Rsamples
Rsamples - saveSamples

##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)
samplesPlot(samples)


code <- nimbleCode({ 
    ##q[1:N,1:N]  ~ dwish(R = aq[1:N,1:N], df = bq) ## aq and bq are estimated over time
    ##Q[1:N,1:N] <- inverse(q[1:N,1:N])
    ##X.mod[1:N] ~ dmnorm(muf[1:N],prec = pf[1:N,1:N]) ## Model Forecast ##muf and pf are assigned from ensembles
    #### add process error
    X[1:N]  ~ dmnorm(X.mod[1:N],prec = q[1:N,1:N])
    ## Analysis
    y.censored[1:N] ~ dmnorm(X[1:N], prec = r[1:N,1:N]) 
    for(i in 1:N){
        y.ind[i] ~ dconstraint(y.censored[i] > 0)
    }
})
N <- 5
constants <- list(N=N, X.mod=rep(10,N), q=diag(N), r=diag(N))
## value of y.ind and y.censored are just placeholders...
## BUT, important, y.censored *must* be specified as data at this point.
## it's a long story why, has to do with the initializeModel routine
## at the beginning of MCMC execution
data <- list(y.ind=rep(1,N), y.censored=rep(10,N))  
inits <- list(X=rep(10,N))

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, print=TRUE)

## important!
## this is needed for correct indexing later
samplerNumberOffset <- length(conf$getSamplers())

for(i in 1:N) {
    node <- paste0('y.censored[',i,']')
    conf$addSampler(node, 'toggle', control=list(type='RW'))
    ## could instead use slice samplers, or any combination thereof, e.g.:
    ##conf$addSampler(node, 'toggle', control=list(type='slice'))
}

conf$printSamplers()

conf$printMonitors()
## could monitor y.censored, if you wish, to verify correct behaviour
conf$addMonitors('y.censored')

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## new dataset (1)
y.ind <- c(1, 1, 1, 0, 0)
## important!
## change: rather than providing NA for the non-data values (those to be sampled),
## you'll have to provide some values here.
## that's because we effective disabled the model initialization routine earlier
y.censored <- c(9, 9, 11, 20, 20)

Cmodel$y.ind <- y.ind
Cmodel$y.censored <- y.censored

for(i in 1:N) {
    ## ironically, here we have to "toggle" the value of y.ind[i]
    ## this specifies that when y.ind[i] = 1,
    ## indicator variable is set to 0, which specifies *not* to sample
    valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle', 1-y.ind[i])
}

## check the values of indicator variables, if you wish:
##for(i in 1:N) {
##    print(valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle'))
##}

niter <- 1000
set.seed(0)
samples <- runMCMC(Cmcmc, niter)

head(samples)
tail(samples)


## new dataset (2)
y.ind <- c(1, 1, 1, 1, 1)   ## everything is data
## important!
## change: rather than providing NA for the non-data values (those to be sampled),
## you'll have to provide some values here.
## that's because we effective disabled the model initialization routine earlier
y.censored <- c(9, 9, 11, 11, 12)

Cmodel$y.ind <- y.ind
Cmodel$y.censored <- y.censored

for(i in 1:N) {
    ## ironically, here we have to "toggle" the value of y.ind[i]
    ## this specifies that when y.ind[i] = 1,
    ## indicator variable is set to 0, which specifies *not* to sample
    valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle', 1-y.ind[i])
}

## check the values of indicator variables, if you wish:
##for(i in 1:N) {
##    print(valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle'))
##}

niter <- 1000
set.seed(0)
samples <- runMCMC(Cmcmc, niter)

head(samples)
tail(samples)



## new dataset (3)
y.ind <- c(0, 0, 0, 0, 0)   ## nothing is data
## important!
## change: rather than providing NA for the non-data values (those to be sampled),
## you'll have to provide some values here.
## that's because we effective disabled the model initialization routine earlier
y.censored <- c(20, 20, 20, 20, 20)

Cmodel$y.ind <- y.ind
Cmodel$y.censored <- y.censored

for(i in 1:N) {
    ## ironically, here we have to "toggle" the value of y.ind[i]
    ## this specifies that when y.ind[i] = 1,
    ## indicator variable is set to 0, which specifies *not* to sample
    valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle', 1-y.ind[i])
}

## check the values of indicator variables, if you wish:
##for(i in 1:N) {
##    print(valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[samplerNumberOffset+i]], 'toggle'))
##}

niter <- 1000
set.seed(0)
samples <- runMCMC(Cmcmc, niter)

head(samples)
tail(samples)






## prep for STAT 201 lecture

df <- read.csv('~/Downloads/TVhours.csv')
df <- read.delim('~/github/courses/stat201/data/Global.Temperature.txt')
df <- read.csv('~/Downloads/Global.Temperature.csv')
df
str(df)
barplot(df$animal)
df$animal
df

hist(df$tvhours)
hist(df$tvhours, breaks=30)
mean(df$tvhours)
sd(df$tvhours)


## testing speed of new RW sampler, which checks prior logProb first,
## to see if it's a valid proposal

library(nimble)
sampler_RW_new <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        depNodes       <- model$getDependencies(target, self = FALSE)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        ##scaleHistory  <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
        } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        model[[target]] <<- propValue
############### revisons starting here
        ##logMHR <- calculateDiff(model, calcNodes) + propLogScale
        ##jump <- decide(logMHR)
        priorLP0 <- getLogProb(model, target)
        priorLP1 <- calculate(model, target)
        if(is.nan(priorLP1) | priorLP1 == -Inf ) { jump <- FALSE
                               ##print('automatic rejection')
                           } else { logMHR <- calculateDiff(model, depNodes) + priorLP1 - priorLP0 + propLogScale
                                    ##print(logMHR)
                                    jump <- decide(logMHR) }
############### same after here
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                ##setSize(scaleHistory, timesAdapted)         ## scaleHistory
                ##scaleHistory[timesAdapted] <<- scale        ## scaleHistory
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            gamma1 <<- 0
        }
    )
)
##
N <- 100
code <- nimbleCode({
    a ~ dnorm(1, 1)
    ##a ~ T(dnorm(0, 1), -1, 1)
    ##for(i in 1:N) {
    ##    b[i] ~ dnorm(a, 0.001)
    ##}
    ##b ~ dgamma(a, 1)
})
constants <- list()
data <- list()##b = rnorm(N))
inits <- list(a = 0)
##
Rmodel1 <- nimbleModel(code, constants, data, inits)
Cmodel1 <- compileNimble(Rmodel1)
Rmodel2 <- nimbleModel(code, constants, data, inits)
Cmodel2 <- compileNimble(Rmodel2)
##
conf1 <- configureMCMC(Rmodel1, nodes=NULL)
conf1$addSampler('a', 'RW')
conf1$printSamplers()
Rmcmc1 <- buildMCMC(conf1)
Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
##
conf2 <- configureMCMC(Rmodel2, nodes=NULL)
conf2$addSampler('a', 'RW_new')
conf2$printSamplers()
Rmcmc2 <- buildMCMC(conf2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Rmodel2)


niter <- 100

set.seed(0); Rmcmc1$run(niter); Rsamples1 <- as.matrix(Rmcmc1$mvSamples)
set.seed(0); Cmcmc1$run(niter); Csamples1 <- as.matrix(Cmcmc1$mvSamples)
set.seed(0); Rmcmc2$run(niter); Rsamples2 <- as.matrix(Rmcmc2$mvSamples)
set.seed(0); Cmcmc2$run(niter); Csamples2 <- as.matrix(Cmcmc2$mvSamples)


Rsamples1 - Csamples1
Rsamples2 - Csamples2
Rsamples1 - Rsamples2  ## *should be* different when there's truncation



niter <- 1e7

t1 <- system.time(Cmcmc1$run(niter))[1]
t2 <- system.time(Cmcmc2$run(niter))[1]

t1
t2

(t2-t1)/t1*100




library(nimble)
code <- nimbleCode({
    a ~ dbeta(1,1)
    b ~ dgamma(a, 3)
})
constants <- list()
data <- list()
inits <- list(a = -1, b=3)
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)


Rmodel$a
Cmodel$a
calculate(Rmodel, 'a')
calculate(Cmodel, 'a')

Rmodel$b
Cmodel$b
calculate(Rmodel, 'b')
calculate(Cmodel, 'b')

    



## demo of how to use autoBlock on spatial model
## for Dao

library(nimble)
load('~/github/automated-blocking-examples/data/model_spatial.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel, autoBlock=TRUE)
conf$printSamplers()



## sampler to record scale and acceptanceRate history, for Dao

sampler_RW_record <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes  <- model$getDependencies(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory  <- c(0, 0)   ## scaleHistory
        acceptanceRateHistory <- c(0, 0)   ## scaleHistory
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { propLogScale <- rnorm(1, mean = 0, sd = scale)
                       propValue <- currentValue * exp(propLogScale)
                   } else         propValue <- rnorm(1, mean = currentValue,  sd = scale)
        if(reflective) {
            lower <- model$getBound(target, 'lower')
            upper <- model$getBound(target, 'upper')
            while(propValue < lower | propValue > upper) {
                if(propValue < lower) propValue <- 2*lower - propValue
                if(propValue > upper) propValue <- 2*upper - propValue
            }
        }
        model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes) + propLogScale
        jump <- decide(logMHR)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                setSize(scaleHistory,          timesAdapted)         ## scaleHistory
                setSize(acceptanceRateHistory, timesAdapted)         ## scaleHistory
                scaleHistory[timesAdapted]          <<- scale        ## scaleHistory
                acceptanceRateHistory[timesAdapted] <<- acceptanceRate        ## scaleHistory
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        getScaleHistory = function()          { returnType(double(1)); return(scaleHistory) },          ## scaleHistory
        getAcceptanceRateHistory = function() { returnType(double(1)); return(acceptanceRateHistory) },          ## scaleHistory
        ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
        ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
        ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
        ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
        ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
        ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            scaleHistory           <<- scaleHistory * 0    ## scaleHistory
            acceptanceRateHistory  <<- acceptanceRateHistory * 0    ## scaleHistory
            gamma1 <<- 0
        }
    ), where = getLoadingNamespace()
)





## Perry's demo of making R and C external calls in NIMBLE
## This is a demo for beta (or gamma?) users of external call features in nimble
## I have actually built TWO kinds of external calls:
## 1. Call arbitrary (up to limited argument conventions) compiled code as long as you can provide a .h file and a .o file.  Chris suggests maybe we should also (instead) allow the user (that's you) to provide a DLL (.so in Linx / OS X, .dll in Windows).
## 2. Call arbitrary R functions from within NIMBLE (including compiled NIMBLE), as long as argument and return types are guaranteed in advance.  Obviously if you have a need to call external code that cannot be handled by mechanism 1, you could use mechanism 2 and have your R function call your external code in whatever complicated way is needed.  For the case of deSolve, one could call a function that calls deSolve from R.  There is some overhead to using mechanism 2 since it requires copies of arguments and return value to be made, and of course the R steps of execution will be at the speed of R.

## You will need to install from the "call_external_c_function" branch (now mis-named because it also includes calling R functions)
library(devtools)
install_github("nimble-dev/nimble", ref = "call_external_c_function", subdir = "packages/nimble")

## Throughout this demo, ignore warning messages about functions that may not be defined.  I haven't updated the list of keywords used in that error-trapping.

library(nimble)

## Calling external C:

## Say we want a C function that adds 1.5 to a vector of values:
## We can only give non-scalars results by pointer argument.
## Any NIMBLE non-scalar can become a double* of contiguously allocated memory.
## Like in LAPACK or similar low-level libraries, you need a separate argument to say how long the allocated memory is.
sink('add1p5.h')
cat('
extern "C" {
  void my_internal_function(double *p, double*ans, int n);
}
')
sink()

sink('add1p5.cpp') 
cat('
#include "stdio.h"
#include "add1p5.h"

void my_internal_function(double *p, double *ans, int n) {
  printf("In my_internal_function\\n"); /* cat reduces the double slash to single slash */ 
  for(int i = 0; i < n; i++) 
    ans[i] = p[i] + 1.5;
}
')
sink()

## A major limitation is that we cannot have compiled code return anything except a scalar.  So returned non-scalars will have to be via arguments.  That meansto use this in BUGS code we'll have to wrap it in another nimbleFunction.  It will make sense below.

## compile that to .o
## I am not an expert on all compilation twists, but on my system I need to do it with g++ instead of gcc (even though it is pure C) in order for the linker to be later happy linking it to compiled C++ from NIMBLE.  Using g++ on C code gives me a clang warning about deprecated behavior, but it's ok.
system('g++ add1p5.cpp -c -o add1p5.o')

## now to create a nimbleFunction that interfaces to add1p5:

Radd1p5 <- nimbleExternalCall(function(x = double(1), ans = double(1), n = integer()){}, Cfun = 'my_internal_function', headerFile = 'add1p5.h', oFile = 'add1p5.o')
## The first argument uses the format of a nimbleFunction with argument type declarations, but the body of the function can be empty ('{}')
## You can choose any names for the arguments in R. They don't have to match the C code.
## Ignore the warning here and later warnings
Radd1p5 ## this is a nimbleFunction with some internal special sauce, but from here on out it should behave like a nimbleFunction
## We'll wait to use it in BUGS code

## Radd1p5 doesn't return anything and would only be able to return a scalar, so we can wrap it like this:
wrappedRadd1p5 <- nimbleFunction(
    run = function(x = double(1)) {
        ans <- numeric(length(x))
        Radd1p5(x, ans, length(x))
        return(ans)
        returnType(double(1))
    })
## Chris, if you are reading this, it looks like the DSLcode check is weird here.

## calling an R function
## Say we want an R function that adds 2.5 to every value in a vector
add2p5 <- function(x) {
    x + 2.5 ## This can be pure R
}

## now create a nimbleFunction that interfaces to add2p5, not necessary when running uncompiled but necessary when running compiled
Radd2p5 <- nimbleRcall(function(x = double(1)){}, Rfun = 'add2p5', returnType = double(1), envir = .GlobalEnv)
## Similar to above.  The function prototype and the returnType represent a promise that add2p5 will always take and return these types.  Also an environment is needed and it defaults to R's global environment but I'm showing it explicitly.
## Again, ignore the more extensive warnings
Radd2p5 ## this is also a nimbleFunction with more special sauce

## Now let's use these in a model

demoCode <- nimbleCode({
    for(i in 1:4) {x[i] ~ dnorm(0,1)} ## just to get a vector
    y[1:4] <- wrappedRadd1p5(x[1:4])
    z[1:4] <- Radd2p5(x[1:4])
    })

demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)))
## Again ignore the error during checking.  We'll have to trap and handle that, but right now I'm focused on core functionality.
## The model will not work uncompiled!

nimbleOptions(showCompilerOutput = TRUE) ## so you can see what is happening
CdemoModel <- compileNimble(demoModel, dirName = '.') ## last arg puts the C++ code in your working directory so you can look at it if you like

CdemoModel$x
CdemoModel$calculate()
CdemoModel$y - CdemoModel$x
CdemoModel$z - CdemoModel$x

## The correct calculations happened.




library(nimble)
nimbleOptions(MCMCprogressBar = FALSE)


code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


progressBarOption = TRUE
progressBarOption = FALSE

Cmcmc$run(30, progressBar = progressBarOption)
Cmcmc$run(30000, progressBar = progressBarOption)





## testing new sampler function: RW_dirichlet
## for ddirch or ddirich Dirichlet nodes
library(nimble)
##nimbleOptions(showCompilerOutput=TRUE)
n <- 100
alpha <- c(10, 30, 15)#, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09)#, .54, .01)
y <- rmulti(1, n, p)
code <- quote({
    p[1:K]  ~ ddirch(alpha[1:K])
    ##p2[1:K] ~ ddirch(alpha[1:K])
    y[1:K]  ~ dmulti(p[1:K], n)
    ##y2[1:K] ~ dmulti(p2[1:K], n)
    ##x ~ dnorm(p[1], 1)
    ##x2 ~ dnorm(p2[1],1)
    ##for(i in 1:K) {
    ##    alpha[i] ~ dgamma(.001, .001);
    ##}
})
inits <- list(p = rep(1/K, K), alpha = alpha)#, p2 = rep(1/K,K))##, x=0)
constants <- list(n=n, K=K)
data <- list(y = y)#, y2 = y)
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

##Rmodel$alpha
##Rmodel$alpha2
##Rmodel$getParam('p2', 'alpha')
##Rmodel$getParam('p[5:6]', 'alpha')

conf <- configureMCMC(Rmodel, nodes=NULL)
conf$addSampler('p[1:3]',  'RW_dirichlet', print=TRUE)
conf$addSampler('p2[1:3]', 'conjugate',    print=TRUE)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Rmodel$alpha
Rmodel$p
Rmodel$y

debug(Rmcmc$run)
Rmcmc$run(1001)

debug(samplerFunctions[[1]]$reset)
debug(samplerFunctions[[1]]$run)


Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 100000

set.seed(0)
Rsamples <- runMCMC(Rmcmc, niter)

set.seed(0)
Csamples <- runMCMC(Cmcmc, niter)

colnames(Csamples)
apply(Csamples, 2, mean)

for(i in 1:5) {
    samplesPlot(Csamples, paste0(c('p[', 'p2['), i, c(']')))
}


Rmodel$getLogProb()
Cmodel$getLogProb()
Rmodel$calculate()
Cmodel$calculate()


Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

dimnames(Rsamples)
dimnames(Csamples)

Rsamples - Csamples
Csamples[1:200,]

mean(diff(Csamples[,1])==0)
(1-0.44)^3



## reproducible example of problem with model$values(...) <- ...

library(nimble)
code <- nimbleCode({
    for(i in 1:5) {
        a[i] ~ dnorm(0, 1)
    }
})
constants <- list()
data <- list()
inits <- list(a = rep(0,5))
Rmodel <- nimbleModel(code, constants, data, inits)

nfDef <- nimbleFunction(
    setup = function(model) {
        newValues <- rep(1, 5)
        nodes <- 'a'
    },
    run = function() {
        ##values(model, nodes) <<- newValues    ## this one works fine
        model$values(nodes) <<- newValues    ## this causes failure
    }
)

##Error in makeAssgnFcn(e[[1]]) : 
##  'model$values' is not a valid function in complex assignments

Rnf <- nfDef(Rmodel)

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)

Rmodel$a
Rnf$run()
Rmodel$a

Cmodel$a
Cnf$run()
Cmodel$a





## trying to use nodes[i] functionality,
## works for: model[[ multivariateNodeName ]][i] <- scalar
## does not work for:
## model[[ nodes[i] ]] <- x, or
## values(model, nodes[i]) <- x

library(nimble)

code <- nimbleCode({
    for(i in 1:3) {
        a[i] ~ dnorm(0, 1)
    }
})
constants <- list()
data <- list()
inits <- list(a = rep(0,3))
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$a

nfDef <- nimbleFunction(
    setup = function(model) {
        nodes <- 'a'
        d <- length(model$expandNodeNames(nodes))
    },
    run = function() {
        for(i in 1:d) {
            model[[nodes]][i] <<- i*10
        }
    }
)

Rnf <- nfDef(Rmodel)

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)  ## compilation fails

Rmodel$a
Rnf$run()
Rmodel$a

Cmodel$a
Cnf$run()
Cmodel$a


## looking into nimble gthub issue 250, possible infinite
## recursion in getDependencyPaths in checkConjugacy()
## raised by Chris P.

library(nimble)

code <- nimbleCode({
  for(r in 1:R){ 
    for(i in 1:I){
      beta[r,i] ~ dnorm(0,.04)
      beta.pine[r,i] ~ dnorm(0,.04)
    }  
  }
  # if loop over 2nd index, number of dependencies would be greatly reduced
  for(i in 1:I) {
  phi.first[1:J,i] <- Z[1:J,1:R]%*%beta[1:R,i]
  pine.phi[1:J,i] <- Z[1:J,1:R]%*%beta.pine[1:R,i]
  exp.phi[1:J, i] <- exp(phi.first[1:J, i])
  exp.pine.phi[1:J, i] <- exp(pine.phi[1:J, i])
  }
    for(j in 1:J){
    	p.true[j,1] ~ dbeta(exp.phi[j,1],exp.pine.phi[j,1])
    for(i in 2:(I-1)){
        p.rel[j,i]  ~ dbeta(exp.phi[j,i],exp.pine.phi[j,i]) 
        p.true[j,i] <-  p.rel[j,i] * (1 - sum(p.true[j,1:(i-1)]))
    }	
       p.true[j,21] <- 1 - sum(p.true[j,1:20])
    }  
  for(j in 1:J){
    Y[j,1:I] ~ dmulti(size = n[j], prob = p.true[j,1:I])
  }
})

J <- 1
I <- 21
R <- 5

data = list(Y = matrix(0, nrow = J, ncol = I))

constants = list(n = rep(500, J), R = R, I = I, J = J,  Z =  matrix(0, nrow = J, ncol = R))
# try with Z as variable

model <- nimbleModel(code, constants = constants, data = data)

## infinite loop?:
model$checkConjugacy('p.rel[1, 2]')




## generating correlated chains of samples for acf acfplot question
## in stat365 final exam

library(coda)
set.seed(0)
n <- 10000
r1 <- 0.45
r2 <- 0.85
x <- y <- numeric(n)
x[1] <- y[1] <- 0
for(i in 2:n) {
    x[i] <- r1*x[i-1] + rnorm(1)
    y[i] <- r2*y[i-1] + rnorm(1)
}
samp <- cbind(alpha=x, beta=y)
effectiveSize(samp)
setwd('~/github/private/courses/stat365/exams')
pdf('acfPlot.pdf',width=6,height=4,paper='special')
acfplot(as.mcmc(samp), ylim=c(-0.05, 1), lag.max=35)
dev.off()




## looking at interval-censored Weibull regression model
## for Intekhab

library(nimble)
file <- read.csv('~/Downloads/365Project_data.csv')

##dim(file)
##colnames(file)
##head(file)
##str(file)
##file$t

code <- nimbleCode({
    b0 ~ dnorm(0, sd = 1000)
    b_T ~ dnorm(0, sd = 1000)
    b_karnof~ dnorm(0, sd = 1000)
    b_cd4 ~ dnorm(0,sd=1000)
    b_age ~ dnorm(0, sd=1000)
    sigma_N ~ dunif(0,100)
    for (i in 1:N) {
        error[i] ~ dnorm(0,sd=sigma_N)  #random effect
        eta[i] <- b0 + b_T*Tx[i] + b_karnof*karnof[i]+ b_cd4*cd4[i] + b_age*age[i]+error[i]
        mu[i] <- exp(eta[i])
        censored[i] ~ dinterval(t[i], c[i])
        t[i] ~ dweib(1, mu[i]) 
    }
})

file[c(1,14), c('censor','c','t')]

constants <- list(N = nrow(file), Tx=file$tx, karnof=file$karnof, cd4=file$cd4, age=file$age,c=file$c)
data <- list(censored =file$censor,t=file$t)
inits <- list(b0 = 0, b_T = 0, b_karnof=0, b_cd4=0, b_age=0, sigma_N=5, error=rnorm(nrow(file),0,5))

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$simulate('t')
Rmodel$calculate()
Rmodel$calculate('censored')

i <- 1
Rmodel$c[i]
Rmodel$t[i]
Rmodel$censored[i]
Rmodel$mu[i]
Rmodel$logProb_t[i]
Rmodel$logProb_censored[i]


constants
data
inits

Rmodel$c
range(Rmodel$mu)
i <- 1
Rmodel$censored[1]
Rmodel$logProb_censored[14]
range(Rmodel$logProb_error)


Rmodel$getNodeNames(dataOnly=TRUE)

Rmodel$calculate()
Rmodel$simulate()
Rmodel$calculate()


conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 500)

samples
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)
samplesPlot(samples)


## facebook use time comparison example for STAT201


df <- read.csv('~/Downloads/timeFB_FRvsSO.csv')

head(df)

hist(df$FBtime)

dim(df)

t.test(df$FBtime)
t.test(df$FBtime, mu=64)
t.test(df$FBtime, mu=64, alternative='greater')

fr <- subset(df$FBtime, df$class=='FR')
so <- subset(df$FBtime, df$class=='SO')
length(fr)
length(so)

t.test(fr, so)


prop.test(x=6, n=10, correct=FALSE)
prop.test(x=60, n=100, correct=FALSE)  ## by default, p0 = 0.5
prop.test(x=600, n=1000, correct=FALSE, p=.6)




head(df)

## trying different samplers for
## "seizures" example for STAT365

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

library(nimble)
library(coda)
load('~/Downloads/seizures.RData')
N <- dim(seizures$Counts)[1]


code <- nimbleCode({
    b0 ~ dnorm(0, sd=1000)
    bbase ~ dnorm(0, sd=1000)
    bage ~ dnorm(0, sd=1000)
    btrt ~ dnorm(0, sd=1000)
    sigma ~ dunif(0, 1000)
    sigma_patient ~ dunif(0, 1000)
    for(i in 1:N) {
        g[i] ~ dnorm(0, sd=sigma_patient)
        for(j in 1:4) {
            eps[i, j] ~ dnorm(0, sd=sigma)
            log(lambda[i,j]) <- b0 + bbase * log(baseline[i]) + bage * age[i] + btrt * treatment[i] + eps[i,j] + g[i]
            y[i,j] ~ dpois(lambda[i,j])
        }
    }
})
constants <- list(N = dim(seizures$Counts)[1], baseline=seizures$Baseline, age=seizures$Age, treatment=seizures$Treatment)
data <- list(y=seizures$Counts)
inits <- list(b0=1, bbase=0, bage=0, btrt=0, sigma=1, sigma_patient=1, g=rep(0,N), eps = array(0,c(N,4)))
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf <- configureMCMC(Rmodel, onlySlice=TRUE)
conf <- configureMCMC(Rmodel, autoBlock=TRUE)
conf <- configureMCMC(Rmodel, onlySlice=TRUE, control=list(adaptInterval=100))

conf$addSampler(c('b0','bage','bbase'), 'RW_block')

conf$printSamplers()
##conf$removeSamplers('b0')
##conf$addSampler('b0', 'slice')
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)

system.time(samples <- runMCMC(Cmcmc, 10000, nburnin=2000))

samplesPlot(samples)


effectiveSize(samples)
## default RW, took 2 seconds
##           b0          bage         bbase          btrt         sigma 
##     9.639895     10.662738     13.488713     59.645279    553.986654 
##sigma_patient 
##    24.535715 




## fitting classification model of oak tree heights


load('~/Downloads/oak_heights.RData')
ls()

hist(oak_heights, breaks=20)

library(nimble)

N <- length(oak_heights)

code <- nimbleCode({
    mu0 ~ dnorm(0, sd=10000)
    mu1 ~ dnorm(0, sd=10000)
    sigma0 ~ dunif(0, 10000)
    sigma1 ~ dunif(0, 10000)
    con ~ dconstraint(mu0 < mu1)
    for(i in 1:N) {
        z[i] ~ dbern(0.5)
        means[i]  <- equals(z[i],0)*mu0    + equals(z[i],1)*mu1
        sigmas[i] <- equals(z[i],0)*sigma0 + equals(z[i],1)*sigma1
        y[i] ~ dnorm(means[i], sd = sigmas[i])
    }
})

constants <- list(N=N)
data <- list(y = oak_heights, con=1)
inits <- list(mu0=0, mu1=1, sigma0=1, sigma1=1, z=rep(0,N))

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$calculate()




## experimenting with t.test and prop.test()
## for STAT201

x1 <- 347
n1 <- 11535
x2 <- 327
n2 <- 14035
p1 <- x1/n1
p2 <- x2/n2
p1
p2

se <- sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
se
ci <- p1-p2 + c(-1,1) * 1.96*se
ci

p_pooled <- (x1+x2)/(n1+n2)
p_pooled
se_pooled = sqrt(p_pooled*(1-p_pooled)*(1/n1+1/n2))
se_pooled

z <- (p1-p2) / (se_pooled)
z
2*(1-pnorm(z))


prop.test(c(x1,x2), c(n1,n2), correct=FALSE)

prop.test(220,400,correct=FALSE)
prop.test(219,400,correct=FALSE)
prop.test(c(219,220),c(400,400),correct=FALSE)

xbar1 <- 33
xbar2 <- 19.9
n1 <- 476
n2 <- 496
s1 <- 21.9
s2 <- 14.6

se <- sqrt(s1^2/n1 + s2^2/n2)
se
ci <- xbar1-xbar2 + c(-1,1)*1.96*se
ci
z99 <- qnorm(0.995)
z99
ci99 <- xbar1-xbar2 + c(-1,1)*z99*se
ci99


f <- function(..., a=3) {
    browser()
    x <- list(...)
    print(x)
    print(a)
}

f(a=3, 5)


## experimentation with dconstraint and truncation

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    con1 ~ dconstraint(a > 0)
    b ~ T(dnorm(.5, 1),a,)
    c ~ dnorm(0, 1)
    d ~ dnorm(c, 1)
    e ~ dexp(c+d)
    f ~ dexp(e^2)
})
constants <- list()
data <- list(con1=1)
inits <- list(a=1, b=2.5, c=0, d=0, e=1, f=1)
Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$getDependencies('a')

Rmodel$a
Rmodel$con1
exp(Rmodel$calculate('a'))
exp(Rmodel$calculate('con1'))

Rmodel$a <- -1
Rmodel$a
Rmodel$calculate()
exp(Rmodel$calculate('a'))
exp(Rmodel$calculate('con1'))

Rmodel$a <- 1
Rmodel$a
Rmodel$con1 <- 0
Rmodel$con1
Rmodel$calculate()
exp(Rmodel$calculate('a'))
exp(Rmodel$calculate('con1'))

Rmodel$a <- -10
Rmodel$a
Rmodel$con1 <- 0
Rmodel$con1
Rmodel$calculate()
exp(Rmodel$calculate('a'))
exp(Rmodel$calculate('con1'))

conf <- configureMCMC(Rmodel)

conf$getMonitors()
conf$getMonitors2()
conf$printMonitors()

conf$resetMonitors()
conf$addMonitors(c('b'), 'e', c('f', 'a'), print=FALSE)
conf$addMonitors2(c('d', 'e', 'f'), print=FALSE)
conf$addMonitors()
conf$getMonitors()


Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, 500000)


colnames(samples)
##apply(samples, 2, mean)
##samplesPlot(samples)
par(mfrow=c(2,1), mar=c(2,1,1,1))
hist(samples[,'a'], breaks=200, xlim=c(-3,6))
hist(samples[,'b'], breaks=200, xlim=c(-3,6))


## testing new initializeModel() function

library(nimble)
code <- nimbleCode({
    p1 ~ dnorm(a, b)
    p2 <- p1 + 1
    p3 <- p2 + p2
    x[1] ~ dnorm(0,10)
    for(i in 2:10) {
        mu[i] <- x[i-1] + p3
        mu2[i] <- mu[i] + p2
        x[i] ~ dnorm(mu2[i], 1)
    }
    muy <- x[10] +100
    y ~ dnorm(muy, 1)
})
constants <- list()
data <- list(y = 100)
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
RinitModel <- initializeModel(Rmodel)
CinitModel <- compileNimble(RinitModel, project=Rmodel)

model <- Rmodel
ini <- RinitModel

model <- Cmodel
ini <- CinitModel

model$a <- 0
model$b <- 3

set.seed(0)
ini$run()
model$calculate()

model$a
model$b
model$lifted_d1_over_sqrt_oPb_cP
##model$simulate('p1')
model$p1
model$p2
model$p3
model$x
model$mu
model$mu2
model$muy
model$y

model$getDependencies('a')
model$getDependencies('b')





## group membership / classification problem for STAT365
## generating data

mu1 <- 24
sd1 <- 4
mu2 <- 32
sd2 <- .5
n1 <- 300
n2 <- 80
set.seed(0)
y1 <- rnorm(n1, mu1, sd1)
y2 <- rnorm(n2, mu2, sd2)
y <- c(y1, y2)
y <- sort(y)
hist(y, breaks=20)

oak_heights <- y
save(oak_heights, file='~/Downloads/oak_heights.RData')

## fit classification model

library(nimble)

load('~/Downloads/oak_heights.RData')
N <- length(oak_heights)
y <- oak_heights

code <- nimbleCode({
    mu0 ~ dnorm(0, sd=10000)
    mu1 ~ dnorm(0, sd=10000)
    constraint ~ dconstraint(mu0 < mu1)
    sigma0 ~ dunif(0, 10000)
    sigma1 ~ dunif(0, 10000)
    for(i in 1:N) {
        x[i] ~ dbinom(size=1, prob=0.5)
        mu[i] <- equals(x[i],0)*mu0 + equals(x[i],1)*mu1
        sigma[i] <- equals(x[i],0)*sigma0 + equals(x[i],1)*sigma1
        y[i] ~ dnorm(mu[i], sd=sigma[i])
    }
})
constants <- list(N=N)
data <- list(y=oak_heights, constraint=1)
inits <- list(mu0=1, mu1=1, sigma0=1, sigma1=1, x=rep(0,N))

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$getMonitors()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 50000, nburnin=10000)
dim(samples)

topLevel <- c('mu0','mu1','sigma0','sigma1')
apply(samples[,topLevel], 2, mean)
samplesPlot(samples, topLevel)

hist(y, breaks=20)

c(y[253], y[260], y[275])
apply(samples[, c('x[253]','x[260]','x[275]')], 2, mean)


samplesPlot(samples, c('x[253]', 'x[260]', 'x[275]'), ind=1001:1050)



## doing continuous-valued state-space model
## for STAT365 problem set 12

## PHASE 1
## CREATE YOUR DATA

n <- 100
a <- 6
b <- 0.8
sigmaPN <- 2
sigmaOE <- 4

set.seed(0)
x <- numeric(n)
y <- numeric(n)
x[1] <- 1
y[1] <- rnorm(1, x[1], sigmaOE)

for(i in 2:n) {
    x[i] <- rnorm(1, a+b*x[i-1], sigmaPN)
    y[i] <- rnorm(1, x[i], sigmaOE)
}

tail(y)

plot(1:n, y, type='l')
lines(1:n, x, col='blue', lwd=2)

## NOW PHASE TWO
## ALL WE HAVE ARE THE Y's (the data)
## ANALYZE


library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, sd=10000)
    b ~ dnorm(0, sd=10000)
    sigmaPN ~ dunif(0, 10000)
    sigmaOE ~ dunif(0, 10000)
    x[1] ~ dnorm(0, sd=10000)   ## x1 technically a top-level param.
    y[1] ~ dnorm(x[1], sd=sigmaOE)
    for(t in 2:n) {
        x[t] ~ dnorm(a+b*x[t-1], sd=sigmaPN)
        y[t] ~ dnorm(x[t], sd=sigmaOE)
    }
})

constants <- list(n = length(y))
data <- list(y=y)
inits <- list(a=0, b=0, sigmaOE=1, sigmaPN=1, x=y)

Rmodel <- nimbleModel(code, constants, data, inits)




## error "Negative Probability Given to Rank Sample" sent to users list
## by Colin <clewisbeck@gmail.com>
## never did figure this one out...

library(nimble)
## Simulate Data#
set.seed(10)
t<-1:(365)
v<-rnorm(365,sd=0.5)
seasonal<-sin(2*pi*t/365) + cos(pi*2*t/365)
y<- seasonal + v
##plot(y,type="p")
##Create Constant W, G and F matrices to pass to NIMBLE
W<-diag(2)
F<-matrix(c(1,0), nrow = 1, ncol = 2)
##G matrix; can generalize this later
G<-matrix(c(cos(2*pi/365),-sin(2*pi/365),sin(2*pi/365),cos(2*pi/365)),nrow=2,ncol=2)
smosCode <- nimbleCode({
    ##Initial Values for State: x0 ~ Normal(c(0,0),diag(1000))
    m0[1]<-0
    m0[2]<-0
    C0[1:2,1:2]<-1000*W[1:2,1:2]
    ##x0[1:2]~dmnorm(m0[1:2],C0[1:2,1:2])
    x0[1:2]~dmnorm(m0[1:2],cov=C0[1:2,1:2])              ## using covariance
    ##initial values
    m1[1:2,1]<-G[1:2,1:2]%*%x0[1:2]
    C1[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]
    ##x[1:2,1]~dmnorm(m1[1:2,1],C1[1:2,1:2])
    x[1:2,1]~dmnorm(m1[1:2,1],cov=C1[1:2,1:2])           ## using covariance
    y[1]~dnorm(inprod(F[1,1:2],x[1:2,1]),sigmaSquaredInvv)
    ##Model
    for (t in 2:T){
        mu[1:2,t]<-G[1:2,1:2]%*%x[1:2,t-1]
        sigs[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]
        ##x[1:2,t] ~ dmnorm(mu[1:2,t],sigs[1:2,1:2])
        x[1:2,t] ~ dmnorm(mu[1:2,t],cov=sigs[1:2,1:2])   ## using covariance
        y[t] ~ dnorm(inprod(F[1,1:2],x[1:2,t]),sigmaSquaredInvv)
    }
    ##Priors 
    sigmaSquaredInvv~dgamma(5,20)
    sigmaSquaredInvw~dgamma(5,20)
})
smosModel<-nimbleModel(code=smosCode,constants=list(T=365,pi=pi,G=G,W=W,F=F),data=list(y=y), inits=list(sigmaSquaredInvv=1,sigmaSquaredInvw=1))
smosLiuWestFilter<-buildLiuWestFilter(model=smosModel,nodes='x',params=c('sigmaSquaredInvv','sigmaSquaredInvw'))

csmosLiuWestFilter <- compileNimble(smosModel, smosLiuWestFilter)

set.seed(0); smosLiuWestFilter$run(100)
head(as.matrix(smosLiuWestFilter$mvWSamples))[,c(3,4,5,1,2)]
##     sigmaSquaredInvv[1] sigmaSquaredInvw[1]      wts[1]     x[1]      x[2]
##[1,]         0.001055908           0.3644654 -0.04219564 35.55898 -15.77173
##[2,]         0.001056058           0.3644703 -0.04545424 35.64636 -15.62196
##[3,]         0.001056334           0.3644903 -0.04317495 35.58381 -15.96475
##[4,]         0.001056401           0.3644907  0.02075697 35.27854 -14.09417
##[5,]         0.001055928           0.3644734  0.02530712 35.15605 -15.57586
##[6,]         0.001056553           0.3645051  0.03630814 34.85197 -14.91214
tail(as.matrix(smosLiuWestFilter$mvWSamples))[,c(3,4,5,1,2)]
##       sigmaSquaredInvv[1] sigmaSquaredInvw[1]       wts[1]     x[1]      x[2]
##[95,]          0.001065339           0.3649603 -0.025964321 35.14237 -15.52160
##[96,]          0.001065207           0.3649535 -0.015406686 34.85541 -15.56082
##[97,]          0.001065683           0.3649811 -0.060543006 36.06590 -15.44988
##[98,]          0.001065355           0.3649594 -0.048862101 36.16558 -13.87091
##[99,]          0.001065157           0.3649503  0.003280856 34.76580 -14.95481
##[100,]         0.001065097           0.3649547  0.009111265 34.93088 -15.21369

set.seed(0); csmosLiuWestFilter[[2]]$run(10000)
head(as.matrix(csmosLiuWestFilter[[2]]$mvWSamples))
##     sigmaSquaredInvv[1] sigmaSquaredInvw[1]      wts[1]     x[1]      x[2]
##[1,]         0.001055908           0.3644654 -0.04219564 35.55898 -15.77173
##[2,]         0.001056058           0.3644703 -0.04545424 35.64636 -15.62196
##[3,]         0.001056334           0.3644903 -0.04317495 35.58381 -15.96475
##[4,]         0.001056401           0.3644907  0.02075697 35.27854 -14.09417
##[5,]         0.001055928           0.3644734  0.02530712 35.15605 -15.57586
##[6,]         0.001056553           0.3645051  0.03630814 34.85197 -14.91214
tail(as.matrix(csmosLiuWestFilter[[2]]$mvWSamples))
##       sigmaSquaredInvv[1] sigmaSquaredInvw[1]       wts[1]     x[1]      x[2]
##[95,]          0.001065339           0.3649603 -0.025964321 35.14237 -15.52160
##[96,]          0.001065207           0.3649535 -0.015406686 34.85541 -15.56082
##[97,]          0.001065683           0.3649811 -0.060543006 36.06590 -15.44988
##[98,]          0.001065355           0.3649594 -0.048862101 36.16558 -13.87091
##[99,]          0.001065157           0.3649503  0.003280856 34.76580 -14.95481
##[100,]         0.001065097           0.3649547  0.009111265 34.93088 -15.21369



## playing with NIMBLE model dinterval, dconstraint

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, sd=10)
    d ~ dinterval(a, int)
    c ~ dconstraint(a > 0 & int==1)
})

Rmodel <- nimbleModel(code)
Rmodel$getNodeNames()
Rmodel$getNodeNames(includeRHSonly=TRUE)
Rmodel$getDependencies('a')
ls(Rmodel$nodes)
Rmodel$nodes$lifted_a_gt_0_and_int__1
Rmodel$getDependencies('d')
Rmodel$getDependencies('c')

Rmodel$a <- 0
Rmodel$a
exp(Rmodel$calculate('a'))

Rmodel$d <- 0
Rmodel$int <- -2
Rmodel$d
Rmodel$int
exp(Rmodel$calculate('d'))

Rmodel$simulate('d'); Rmodel$d
exp(Rmodel$calculate('d'))

Rmodel$int <- 1
Rmodel$c <- 1
Rmodel$a <- .1
Rmodel$a
Rmodel$int
Rmodel$c
invisible(Rmodel$calculate()); exp(Rmodel$calculate('c'))




## doing wwwhours, web usage problem for STAT201

##df <- read.csv('~/Downloads/wwwhours.csv')
df <- read.csv('~/github/courses/stat201/data/wwwhours.csv')
head(df)
summary(df)
dim(df)
hours <- df$hours
n <- length(hours)
n
hist(hours)
mean(hours)
se <- sd(hours)/sqrt(n)
z <- qnorm(0.975)
ci <- mean(hours) + c(-1,1) * z * se
ci

iter <- 10000
means <- numeric(iter)
for(i in 1:iter) {
    ind <- sample(1:n, replace=TRUE)
    bootstrapSample <- hours[ind]
    means[i] <- mean(bootstrapSample)
}
hist(means)
abline(v=ci, lwd=2, col='red')


bootstrapCI <- quantile(means, c(0.025, 0.975))
bootstrapCI
ci
abline(v=bootstrapCI, lwd=2, col='blue')


## starting to think about BUGS example mice, which is
## a Weibull regression survival analysis with censoring, censored data

l <- list(t = structure(.Data = 
                            c(12, 1, 21, 25, 11, 26, 27, 30, 13, 12, 21, 20, 23, 25, 23, 29, 35, NA, 31, 36, 
                              32, 27, 23, 12, 18, NA, NA, 38, 29, 30, NA, 32, NA, NA, NA, NA, 25, 30, 37, 27, 
                              22, 26, NA, 28, 19, 15, 12, 35, 35, 10, 22, 18, NA, 12, NA, NA, 31, 24, 37, 29, 
                              27, 18, 22, 13, 18, 29, 28, NA, 16, 22, 26, 19, NA, NA, 17, 28, 26, 12, 17, 26),
              .Dim = c(4, 20)),
          t.cen = structure(.Data = 
                                c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 
                                  0, 0, 0, 0, 0, 40, 40, 0, 0, 0, 40, 0, 40, 40, 40, 40, 0, 0, 0, 0, 
                                  0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 0, 40, 40, 0, 0, 0, 0, 
                                  0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 29, 10, 0, 0, 0, 0, 0, 0),
              .Dim = c(4, 20)),
          M = 4, N = 20)



t <- l$t
t.censored <- l$t.cen



t
t.censored
censored <- t.censored==0


library(nimble)

code <- nimbleCode({
    for(i in 1:4) {
        for(j in 1:20) {
            t[i,j] ~ dweib(r, mu[i])
            t[i,j] 
        }
    }
})

constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)
samplesPlot(samples)



## generating Dolphin dolphin capture-recapture data for STAT365


phi <- .8
p <- 0.5
k <- 10
n <- 145
y <- array(NA, c(n,k))
for(i in 1:n) {
    y[i,1] <- 1
    x <- 1
    for(t in 2:k) {
        x <- rbinom(1, prob=x*phi, size=1)
        y[i,t] <- rbinom(1, prob=p*x, size=1)
    }
}

sightings <- y
save(sightings, file='~/Downloads/dolphins.RData')


load('~/github/courses/stat365/data/dolphins.RData')
load('~/downloads/dolphins.RData')



## generate time series data from
## first-order autoregressive process

## simulation parameters
n <- 100
a <- 6
b <- .8
sigmaOE <- 4
sigmaPN <- 2

## generate time-series data: y1, y2, ..., y100
set.seed(0)
x <- numeric(n)
y <- numeric(n)
x[1] <- 1
y[1] <- rnorm(1, x[1], sigmaOE)
for (i in 2:n) {
    x[i] <- rnorm(1, a + b*x[i-1], sigmaPN)
    y[i] <- rnorm(1, x[i], sigmaOE)
}

## visualise the data
plot(1:n, y, type='l')
lines(1:n, x, col='blue', lwd=2)

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 0.0001)
    b ~ dnorm(0, 0.0001)
    sigmaOE ~ dunif(0, 10000)
    sigmaPN ~ dunif(0, 10000)
    x[1] ~ dnorm(0, 0.00001)
    y[1] ~ dnorm(x[1], sd=sigmaOE)
    a[1] ~ dnorm(0, 0.00001)
    for(i in 2:n) {
        a[i] <- p*i
        x[i] ~ dnorm(a[i] + b * x[i-1], sd=sigmaPN)
        y[i] ~ dnorm(x[i], sd=sigmaOE)
    }
})
constants <- list(n=n)
data <- list(y=y)
inits <- list(c=.1, b=0, sigmaOE=1, sigmaPN=1, x=rep(0, n))

Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$resetMonitors()
conf$addMonitors(c('p', 'b', 'sigmaOE', 'sigmaPN'))
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 100000, nburnin=20000)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)
samplesPlot(samples, burnin=1000)



## modelling dolphin capture-recapture data for STAT 365

load('~/Downloads/dolphins.RData')

library(nimble)

dim(sightings)
nind <- dim(sightings)[1]
nocc <- dim(sightings)[2]

code <- nimbleCode({
    p   ~ dunif(0,1)
    phi ~ dunif(0,1)
    for(ind in 1:nind) {
        x[ind,1] <- 1
        y[ind,1] <- 1
        for(t in 2:nocc) {
            x[ind,t] ~ dbinom(size=1, prob=phi*x[ind,t-1])
            y[ind,t] ~ dbinom(size=1, prob=p*x[ind,t])
        }
    }
})
constants <- list(nind=nind, nocc=nocc)
data <- list(y=sightings)
inits <- list(x=array(1,c(nind,nocc)), p=0.5, phi=0.5)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()

conf <- configureMCMC(Rmodel)

conf$printSamplers('p')
conf$printSamplers('phi')
conf$printSamplers()
conf$getMonitors()





## testing new addition to NIMBLE: conf$addSampler('node', 'conjugate')
library(nimble)

code <- nimbleCode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)

conf$printSamplers()
##[1]  conjugate_dgamma_dnorm_dpois_dgamma sampler: x,  dep_dnorm: a, a2,  dep_dpois: b, b2,  dep_dgamma: c
##[2]  conjugate_dnorm_dnorm_dlnorm sampler: a,  dep_dnorm: jNorm[2], jNorm[3],  dep_dlnorm: kLogNorm[2], kLogNorm[3]
##[3]  posterior_predictive sampler: a2
##[4]  posterior_predictive sampler: b
##[5]  posterior_predictive sampler: b2
##[6]  RW sampler: c
##[7]  posterior_predictive sampler: kLogNorm[2]
##[8]  posterior_predictive sampler: kLogNorm[3]
##[9]  posterior_predictive sampler: jNorm[2]
##[10] posterior_predictive sampler: jNorm[3]

pr <- TRUE
pr <- FALSE
nd <- 'x'
nd <- 'a'
nd <- 'c'
nd <- 'a2'
nd <- 'kLogNorm[3]'
conf$addSampler(nd, 'RW',        print=pr)
conf$addSampler(nd, 'slice',     print=pr)
conf$addSampler(nd, 'conjugate', print=pr)

conf$printSamplers()

conf$addSampler(nd, print=TRUE)
conf$addSampler(nd, print=TRUE)
conf$addSampler(nd, print=TRUE)





par(mfrow=c(1,3))
ns <- c(50, 500, 5000)
for(n in ns) {
    x <- rbinom(n=10000, size=n, prob=0.8) / n
    hist(x, xlim=c(0,1))
}




library(nimble)
df <- read.csv('~/Downloads/UsedCars.csv')
code <- nimbleCode({
    b0  ~ dnorm(0, sd=10000)
    bage ~ dnorm(0, sd=10000)
    bhp ~ dnorm(0, sd=10000)
    btype ~ dnorm(0, sd=10000)
    sigma ~ dunif(0, 50000)
    for(i in 1:N) {
        y[i] ~ dnorm(mu[i], sd = sigma)
        mu[i] <- b0 + bage*age[i] + bhp*hp[i] + btype*type[i]
    }
})
constants <- list(N = dim(df)[1], age=df$Age, hp=df$HP, type=df$Type)
data <- list(y = df$Price)
inits <- list(b0=0, bage=0, bhp=0, btype=0, sigma=1)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


samples <- runMCMC(Cmcmc, 20000, nburnin=10000)
dim(samples)

samplesPlot(samples)
samplesPlot(samples, 'btype')

apply(samples, 2, effectiveSize)






1

a <- acf(samples)
a
a$acf
plot(a)
acfplot(as.mcmc(samples), thin=10)

library(coda)
geweke.diag(codaMCMC[[1]][-(1:1000),])
geweke.plot(codaMCMC[[1]])
samplesPlot(codaMCMC[[1]])





?acf
acfplot

colnames(samples)
quantile(samples[, 'bage'], c(0.025, 0.975))

dim(samples)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var='btype', burnin=2000)
samplesPlot(samples, var='btype', ind=1:200)
samplesPlot(samples, var='btype', ind=2001:10000)
samplesPlot(samples, var=c('btype', 'bage'), ind=2001:10000)
samplesPlot(samples, var=c('btype', 'bage'), ind=2001:10000, densityplots=FALSE)

samplesPlot(samples, burnin=2000)
samplesPlot(samples, var=1:3, burnin=2000)
samplesPlot(samples, var=4, burnin=5000)
samplesPlot(samples, var=1, burnin=2000)

quantile(samples[-(1:2000), 1], c(0.025, 0.975))








## "seizures" example for STAT365

library(nimble)
load('~/Downloads/seizures.RData')
N <- dim(seizures$Counts)[1]


code <- nimbleCode({
    b0 ~ dnorm(0, sd=1000)
    bbase ~ dnorm(0, sd=1000)
    bage ~ dnorm(0, sd=1000)
    btrt ~ dnorm(0, sd=1000)
    sigma ~ dunif(0, 1000)
    sigma_patient ~ dunif(0, 1000)
    for(i in 1:N) {
        g[i] ~ dnorm(0, sd=sigma_patient)
        for(j in 1:4) {
            eps[i, j] ~ dnorm(0, sd=sigma)
            log(lambda[i,j]) <- b0 + bbase * log(baseline[i]) + bage * age[i] + btrt * treatment[i] + eps[i,j] + g[i]
            y[i,j] ~ dpois(lambda[i,j])
        }
    }
})
constants <- list(N = dim(seizures$Counts)[1], baseline=seizures$Baseline, age=seizures$Age, treatment=seizures$Treatment)
data <- list(y=seizures$Counts)
inits <- list(b0=1, bbase=0, bage=0, btrt=0, sigma=1, sigma_patient=1, g=rep(0,N), eps = array(0,c(N,4)))
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, onlySlice=TRUE)
conf$printSamplers()
##conf$removeSamplers('b0')
##conf$addSampler('b0', 'slice')
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samplesList <- runMCMC(Cmcmc, 500000, nburnin=10000, nchains=3, returnCodaMCMC=TRUE)

samples <- samplesList[[1]]
samplesPlot(samples)
samplesPlot(samples, 'btrt')
samplesPlot(samples, 'b0')
##samplesPlot(samples, 'bbase')
library(coda)
apply(samples, 2, effectiveSize)
##dim(samples)
gelman.diag(samplesList)

apply(samples, 2, effectiveSize)
apply(samples, 2, mean)
sqrt(apply(samples, 2, var)) / sqrt(apply(samples, 2, effectiveSize))



codalist <- runMCMC(Cmcmc, 100000, nchains=3, returnCodaMCMC=TRUE)

gelman.diag(codalist)
geweke.diag(codalist)
geweke.plot(codalist)

apply(codalist[[1]], 2, effectiveSize)
apply(codalist[[1]], 2, var) / apply(codalist[[1]], 2, effectiveSize)


## "seeds" example for STAT365
write.csv(df, '~/Downloads/Seeds.csv')


df <- read.csv('~/Downloads/Seeds.csv')
df

cucumber <- as.numeric(df$plant) - 1
fertB <- as.numeric(df$fertilizer) - 1
y <- df$germinations
n <- df$seeds
N <- dim(df)[1]

library(nimble)


sds <- c(.1, .5, 1, 5, 10)
samps <- array(0, c(10000, length(sds)))

for(i in seq_along(sds)) {
    sdC <- sds[i]
    code <- nimbleCode({
        sigma ~ dunif(0, 100)
        b0 ~ dnorm(0, sd=1000)
        bCuc ~ dnorm(0, sd=sdC)
        bFertB ~ dnorm(0, sd=1000)
        for(i in 1:N) {
            mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
            logit(p[i]) ~ dnorm(mu[i], sd=sigma)
            y[i] ~ dbinom(size=n[i], prob=p[i])
        }
    })
    constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n, sdC=sdC)
    data <- list(y=y)
    inits <- list(sigma=1, b0=1, bCuc=0, bFertB=0)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$mu
    calculate(Rmodel)
    conf <- configureMCMC(Rmodel)
    conf$printSamplers()
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    set.seed(0)
    samples <- runMCMC(Cmcmc, 100000)
    ##Cmcmc$run(10000)
    ##samples <- as.matrix(Cmcmc$mvSamples)
    samps[, i] <- samples[ 50001:60000, 'bCuc']
}

dim(samps)
dimnames(samps)
colnames(samps) <- paste0('sd', as.character(sds))
dimnames(samps)

samplesPlot(samples)

## prior sensitity analysis for different
## parametrizations of the 'sigma' term
niter <- 100000
samps <- array(0, c(niter, 3))
colnames(samps) <- c('sd', 'var', 'tau')

code <- nimbleCode({
    sigma ~ dunif(0, 1000)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], sd=sigma)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(sigma=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 1] <- samples[, 'sigma']
code <- nimbleCode({
    v ~ dunif(0, 1000)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], var=v)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(v=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 2] <- samples[, 'v']
code <- nimbleCode({
    tau ~ dgamma(0.001, 0.001)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], tau=tau)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(tau=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 3] <- samples[, 'tau']

samps[, 'var'] <- sqrt(samps[, 'var'])
samps[, 'tau'] <- 1/sqrt(samps[, 'tau'])


head(samps)

dim(samps)
dimnames(samps)

samplesPlot(samps, burnin=50000)








## "pairs"
## studying estimation of normal variance / covariance,
## from paired observations from a common variance
## (and possibly different mean)


## this shows that for pairs x1,x2 ~ N(mean=[changing], var= constant v)
## using the average of (x2-x1)^2 / 2
## gives an unbiased estimate of the common variance v
estimateV <- function(n, v) {
    pairs <- t(replicate(n, rnorm(2, mean=runif(1,0,100), sd=sqrt(v))))
    vest <- apply(pairs, 1, function(x) ((x[2]-x[1])^2)/2)
    mean(vest)
}
n <- 1000
v <- 2
n.average.over=1000
estimates <- replicate(n.average.over, estimateV(n=n, v=v))
mean(estimates)

## this shows that the above (vest = (x2-x1)^2 / 2)
## is *not* the same as a traditional sd(x), if the x are used to
## generate a vector of n consecutive "pairs"
n <- 10
v <- 2
xs <- rnorm(n, sd=sqrt(v))
sd(xs)^2
pairs <- array(NA, c(n,2))
for(i in 1:(n-1)) pairs[i,] <- xs[i:(i+1)]
pairs[n,] <- xs[c(n,1)]
xs
pairs
vest <- apply(pairs, 1, function(x) ((x[2]-x[1])^2)/2)
sum(vest) / (n-1)
sd(xs)^2


## this shows that for multivariate-normal pairs
## x1,x2 ~ MVN(mean=[changing], Sigma = constant V)
## using the average of dif=x2-x1, dif %*% t(dif) / 2
## gives an unbiased estimate of the common covariance matrix V
estimateV <- function(n, V) {
    pairs <- lapply(1:n, function(x) {
        mus <- runif(2,0,100)
        list(rmvnorm(1, mus, sigma=V)[1,],
             rmvnorm(1, mus, sigma=V)[1,])
    })
    vest <- lapply(pairs, function(x) {
        dif <- x[[2]] - x[[1]]
        (dif %*% t(dif)) / 2
    })
    apply(array(unlist(vest), dim=c(2,2,n)), c(1,2), sum) / n
}
## this new way does the "better" "newer" way of estimating
## empirical covariance.  I'm really just implementing it,
## because if it works, will be easier to modify the block sampler.
estimateV2 <- function(n, V) {
    pairs <- lapply(1:n, function(x) {
        mus <- runif(2,0,100)
        list(rmvnorm(1, mus, sigma=V)[1,],
             rmvnorm(1, mus, sigma=V)[1,])
    })
    ## original:
    ##vest <- lapply(pairs, function(x) {
    ##    dif <- x[[2]] - x[[1]]
    ##    (dif %*% t(dif)) / 2
    ##})
    ##apply(array(unlist(vest), dim=c(2,2,n)), c(1,2), sum) / n
    ## following is, essentially, vest2:
    difs <- lapply(pairs, function(x) { x[[2]] - x[[1]] })
    ## keep track of array of the after-before differences:
    empirSamp <- matrix(unlist(difs), nrow=n, ncol=2, byrow=TRUE)
    ## estimate covariance as:
    empirCov <- t(empirSamp) %*% empirSamp / (2*n)  ## this works!
    empirCov
}
library(mvtnorm)
n <- 200
v1 <- 2
v2 <- 5
rho <- 0.6
V <- array(c(v1^2, rho*v1*v2, rho*v1*v2, v2^2), c(2,2))
n.average.over <- 200
set.seed(0); estimates  <- replicate(n.average.over, estimateV( n=n, V=V))
set.seed(0); estimates2 <- replicate(n.average.over, estimateV2(n=n, V=V))
V
## two methods are identical:
apply(estimates,  c(1,2), mean)
apply(estimates2, c(1,2), mean)

## still working on "pairs"
## now, testing that the new block sampler, RW_block_new
## appears to be correctly saving before/after samples, etc..

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a/2, 1)
    ##c ~ dgamma(a^2, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0, b=0)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
conf$printSamplers()
conf$addSampler(c('a', 'b'), 'RW_block_new', control=list(adaptInterval=10))
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

iter <- 5
set.seed(0); samples <- runMCMC(Rmcmc, iter)
set.seed(0); samples <- runMCMC(Cmcmc, iter)



## still continuing "pairs"
## now, running correlated state space model, and seeing if
## new RW_block_new sampler adapts towards a new, different covariance?
## and whether it achieves good sampling quicker???
## re-creating the deep-dive deep dive into correlated state-space state space model
## to figure out why the block sampler is taking so long to adapt (from Dao).
## where 'a' and 'b' don't mix until after 150,000 iterations
## this time trying the RW_block_new sampler

library(nimble)
library(coda)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
#####nimbleOptions(showCompilerOutput=TRUE)
##
code <- nimbleCode({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(1e-04, 1)
    sigOE ~ dunif(1e-04, 1)
    x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for (i in 2:t) {
        x[i] ~ dnorm(x[i - 1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
constants <- list(t = 100)
data <- list(y = c(20.24405,20.57693,20.49357,20.34159,20.45759,20.43326,20.20554,20.12860,20.14756,20.20781,20.23022,20.26766,20.22984,20.37703,20.13641,20.05309,19.95709,20.19303,20.30562,20.54443,20.91010,20.70580,20.42344,20.19795,20.28816,20.31894,20.76939,20.77023,20.83486,20.29335,20.40990,20.19601,20.04083,19.76056,19.80810,19.83129,19.69174,19.90069,19.87623,19.63371,19.62360,19.72630,19.64450,19.86779,20.17104,20.34797,20.32968,20.48027,20.46694,20.47006,20.51676,20.40695,20.18715,19.97552,19.88331,19.67831,19.74702,19.47502,19.24408,19.37179,19.38277,19.15034,19.08723,19.37051,19.14274,19.46433,19.62459,19.77971,19.54194,19.39081,19.61621,19.51307,19.34745,19.17019,19.26829,19.58943,19.77143,19.83582,19.71198,19.67746,19.75053,20.40197,20.49363,20.37079,20.19005,20.55862,20.48523,20.33071,19.97069,19.79758,19.83811,19.79728,19.86277,19.86836,19.92481,19.88095,20.24899,20.55165,20.22707,20.11235))
inits <- list(a = 0.95, b=1, sigPN = 0.2, sigOE=0.05, x = c(20.26036,20.51331,20.57057,20.35633,20.33736,20.47321,20.22002,20.14917,20.19216,20.26969,20.21135,20.22745,20.20466,20.41158,20.13408,20.08023,19.98956,20.13543,20.32709,20.55840,20.88206,20.74740,20.47671,20.14012,20.29953,20.33778,20.80916,20.75773,20.84349,20.35654,20.41045,20.20180,20.02872,19.74226,19.80483,19.81842,19.69770,19.84564,19.88211,19.70559,19.56090,19.73728,19.66545,19.88158,20.13870,20.39163,20.37372,20.47429,20.39414,20.42024,20.55560,20.40462,20.15831,19.89425,19.79939,19.72692,19.74565,19.42233,19.22730,19.36489,19.37289,19.19050,19.00823,19.35738,19.14293,19.48812,19.67329,19.82750,19.58979,19.43634,19.61278,19.56739,19.38584,19.19260,19.32732,19.65500,19.65295,19.84843,19.68285,19.69620,19.77497,20.31795,20.45797,20.32650,20.24045,20.60507,20.51597,20.30076,19.98100,19.86709,19.85965,19.74822,19.86730,19.90523,19.86970,19.87286,20.28417,20.46212,20.22618,20.13689))
##
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## [1] 183.3436
Rmodel2 <- Rmodel$newModel(replicate=TRUE)
Rmodel2$calculate()   ## [1] 183.3436
Rmodel3 <- Rmodel$newModel(replicate=TRUE)
Rmodel3$calculate()   ## [1] 183.3436
##
conf <- configureMCMC(Rmodel, nodes = NULL)
conf$addSampler(c('a', 'b'), 'RW_block')
conf$printSamplers()
conf$addSampler('sigOE', 'RW')
conf$addSampler('sigPN', 'RW')
for(node in Rmodel$expandNodeNames('x'))
    conf$addSampler(node, 'RW')
conf$resetMonitors()
conf$addMonitors(c('a', 'b', 'sigOE', 'sigPN'))
conf$getMonitors()
##
conf2 <- configureMCMC(Rmodel2, nodes = NULL)
conf2$addSampler(c('a', 'b'), 'RW_block_new', control=list(useAcceptedOnly=FALSE))
##conf2$addSampler(c('a', 'b'), 'AF_slice', control=list(sliceWidths=c(1,1)))
conf2$printSamplers()
conf2$addSampler('sigOE', 'RW')
conf2$addSampler('sigPN', 'RW')
for(node in Rmodel$expandNodeNames('x'))
    conf2$addSampler(node, 'RW')
conf2$resetMonitors()
conf2$addMonitors(c('a', 'b', 'sigOE', 'sigPN'))
conf2$getMonitors()
##
conf3 <- configureMCMC(Rmodel3, nodes = NULL)
conf3$addSampler(c('a', 'b'), 'RW_block_new', control=list(useAcceptedOnly=TRUE))
##conf3$addSampler(c('a', 'b'), 'AF_slice', control=list(sliceWidths=c(1,1)))
conf3$printSamplers()
conf3$addSampler('sigOE', 'RW')
conf3$addSampler('sigPN', 'RW')
for(node in Rmodel$expandNodeNames('x'))
    conf3$addSampler(node, 'RW')
conf3$resetMonitors()
conf3$addMonitors(c('a', 'b', 'sigOE', 'sigPN'))
conf3$getMonitors()

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()   ## [1] 183.3436
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmodel2$calculate()   ## [1] 183.3436
Cmcmc2 <- compileNimble(Rmcmc2, project = Rmodel2)
##
Rmcmc3 <- buildMCMC(conf3)
Cmodel3 <- compileNimble(Rmodel3)
Cmodel3$calculate()   ## [1] 183.3436
Cmcmc3 <- compileNimble(Rmcmc3, project = Rmodel3)

##niter <- 100000
##niter <- 300000
##niter <- 500000
niter <- 1000000

set.seed(0); system.time(Cmcmc$run(niter))
set.seed(0); system.time(Cmcmc2$run(niter))
set.seed(0); system.time(Cmcmc3$run(niter))

samples  <- as.matrix(Cmcmc$mvSamples)
samples2 <- as.matrix(Cmcmc2$mvSamples)
samples3 <- as.matrix(Cmcmc3$mvSamples)

#### Plot1
dev.new(width=8, height=6)
par(mfrow=c(3,2), mar=c(1,1,1,1))
indToPlot <- 1:700000
plot( samples[indToPlot,'a'], type='l', ylab='a')
plot( samples[indToPlot,'b'], type='l', ylab='b')
plot(samples2[indToPlot,'a'], type='l', ylab='a')
plot(samples2[indToPlot,'b'], type='l', ylab='b')
plot(samples3[indToPlot,'a'], type='l', ylab='a')
plot(samples3[indToPlot,'b'], type='l', ylab='b')

plot( samples[,'a'], type='l', ylab='a')
plot( samples[,'b'], type='l', ylab='b')
plot(samples2[,'a'], type='l', ylab='a')
plot(samples2[,'b'], type='l', ylab='b')
plot(samples3[,'a'], type='l', ylab='a')
plot(samples3[,'b'], type='l', ylab='b')

## plots of sigma term sampling look just fine, in both cases
##samplesPlot(samples, var=c('sigOE','sigPN'))
##samplesPlot(samples2, var=c('sigOE','sigPN'))

pdf('~/downloads/RWblock.pdf')
par(mfrow=c(2,1), mar=c(2,1,2,1))
indToPlot <- 1:200000
plot( samples[indToPlot,'a'], type='l', ylab='a', main = 'MCMC Samples of Parameter 1')
plot( samples[indToPlot,'b'], type='l', ylab='b', main = 'MCMC Samples of Parameter 2')
dev.off()


block_scales <- Cmcmc$samplerFunctions$contentsList[[1]]$getScaleHistory()
length(block_scales)
block_propCovHistory <- Cmcmc$samplerFunctions$contentsList[[1]]$getPropCovHistory()
## create block_propCovScale
block_propCovScale <- block_propCovHistory
for(i in 1:length(block_scales))   block_propCovScale[i,,] <- block_scales[i] * block_propCovHistory[i,,]
block_scale_a <- apply(block_propCovScale, 1, function(x) sqrt(x[1,1]))
block_scale_b <- apply(block_propCovScale, 1, function(x) sqrt(x[2,2]))
block_cors <- apply(block_propCovHistory, 1, function(x) cov2cor(x)[1,2])
ar <- cbind(block_scales, block_scale_a, block_scale_b, block_cors)
colnames(ar) <- c('scale', 'sig_a', 'sig_b', 'cor')
samplesPlot(ar)

block_scales2 <- Cmcmc2$samplerFunctions$contentsList[[1]]$getScaleHistory()
length(block_scales2)
block_propCovHistory2 <- Cmcmc2$samplerFunctions$contentsList[[1]]$getPropCovHistory()
## create block_propCovScale
block_propCovScale2 <- block_propCovHistory2
for(i in 1:length(block_scales2))   block_propCovScale2[i,,] <- block_scales2[i] * block_propCovHistory2[i,,]
block_scale_a2 <- apply(block_propCovScale2, 1, function(x) sqrt(x[1,1]))
block_scale_b2 <- apply(block_propCovScale2, 1, function(x) sqrt(x[2,2]))
block_cors2 <- apply(block_propCovHistory2, 1, function(x) cov2cor(x)[1,2])
ar2 <- cbind(block_scales2, block_scale_a2, block_scale_b2, block_cors2)
colnames(ar2) <- c('scale', 'sig_a', 'sig_b', 'cor')
samplesPlot(ar2)

block_scales3 <- Cmcmc3$samplerFunctions$contentsList[[1]]$getScaleHistory()
length(block_scales3)
block_propCovHistory3 <- Cmcmc3$samplerFunctions$contentsList[[1]]$getPropCovHistory()
## create block_propCovScale
block_propCovScale3 <- block_propCovHistory3
for(i in 1:length(block_scales3))   block_propCovScale3[i,,] <- block_scales3[i] * block_propCovHistory3[i,,]
block_scale_a3 <- apply(block_propCovScale3, 1, function(x) sqrt(x[1,1]))
block_scale_b3 <- apply(block_propCovScale3, 1, function(x) sqrt(x[2,2]))
block_cors3 <- apply(block_propCovHistory3, 1, function(x) cov2cor(x)[1,2])
ar3 <- cbind(block_scales3, block_scale_a3, block_scale_b3, block_cors3)
colnames(ar3) <- c('scale', 'sig_a', 'sig_b', 'cor')
samplesPlot(ar3)


## here's propCov that RW_block adapts towards:
block_propCovScale[dim(block_propCovScale)[1],,]
##            [,1]        [,2]
##[1,]  0.00293348 -0.05870338
##[2,] -0.05870338  1.17509114
##(from values of scale):
block_scales[length(block_scales)]
##[1] 1.686427
##(and propCovHistory):
block_propCovHistory[dim(block_propCovHistory)[1],,]
##             [,1]        [,2]
##[1,]  0.001739464 -0.03480932
##[2,] -0.034809317  0.69679328

## empirical covariance between a&b (what RW_block adapts towards)
ind <- 800001:1000000
cov(samples[ind, c('a','b')])
##             a           b
##a  0.001805492 -0.03613168
##b -0.036131683  0.72328234
cov(samples2[ind, c('a','b')])  ## same!
##             a          b
##a  0.001802017 -0.0360598
##b -0.036059804  0.7217960

## here's propCov that RW_block_NEW adapts towards (conditional cov):
block_propCovScale2[dim(block_propCovScale2)[1],,]
##             [,1]        [,2]
##[1,]  0.001016699 -0.02033899
##[2,] -0.020338986  0.40709698
block_scales2[length(block_scales2)]
##[1] 3.65329
##(and propCovHistory):
block_propCovHistory2[dim(block_propCovHistory2)[1],,]
##              [,1]         [,2]
##[1,]  0.0002782968 -0.005567306
##[2,] -0.0055673059  0.111432961


## what we're calling the "conditional covariance" of a&b:
ind <- 800001:1000000
difpairs <- diff(samples[ind,c('a','b')])
(t(difpairs) %*% difpairs) / (2*dim(difpairs)[1])
##             a            b
##a  0.000400004 -0.008005672
##b -0.008005672  0.160272494
difpairs2 <- diff(samples2[ind,c('a','b')])
(t(difpairs2) %*% difpairs2) / (2*dim(difpairs2)[1])
##              a            b
##a  0.0003152922 -0.006307686
##b -0.0063076856  0.126247785


## using these FINAL ADAPTED scale and propCov values, how about ESS??
dim(samples)
dimnames(samples)
ind <- 800001:1000000
## for original RW_block:
ess <- effectiveSize(samples[ind,])
round(ess)
##    a     b sigOE sigPN 
##22954 22969   251 16831 
## for NEW RW_block_NEW (adapting to conditional covariance):
ess2 <- effectiveSize(samples2[ind,])
round(ess2)
##    a     b sigOE sigPN 
##18452 18443   185 15305 
ess3 <- effectiveSize(samples3[ind,])
round(ess3)
##    a     b sigOE sigPN 
##21923 21937   147 17625

##
## conclusion: mixing is WORSE using RW_block_NEW, quite reasonably so.....
## conclusion2: with using accepted proposals only, it isn't quite as bad.


## 
#### final constant that scale approaches:
##block_scales[length(block_scales)]
####propCov adapts very nicely to true covariance between 'a' and 'b'
##cov(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
##block_propCovHistory[dim(ar)[1],,]
#### final adapted (and scaled) proposal corrleation is very accurate:
##cor(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
##cov2cor(block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])
#### final adapted (and scaled) proposal standard deviations for 'a' and 'b':
##sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[1,1])
##sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[2,2])
## 
#### expand cor, sig_a, and sig_b by adaptInterval:
##length(block_scale_a)
##aI <- 200
##block_scale_a_ex <- rep(block_scale_a, each=aI)
##block_scale_b_ex <- rep(block_scale_b, each=aI)
##block_cors_ex    <- rep(block_cors,    each=aI)
##block_scales_ex  <- rep(block_scales,  each=aI)
##samples_block_info <- cbind(samples[,'a'], samples[,'b'], block_scale_a_ex, block_scale_b_ex, block_cors_ex, block_scales_ex)
##dimnames(samples_block_info)[[2]] <- c('a', 'b', 'sig_a', 'sig_b', 'cor', 'scale')
## 
####samplesPlot(samples_block_info, ind=1:300000, var=c('b', 'a'))
####samplesPlot(samples_block_info, ind=1:300000, var=c('sig_a', 'sig_b', 'cor', 'scale'))
## 
#### this one is best:
#### artificially trim 'b' samples:
##samples_block_info_trim <- samples_block_info
##samples_block_info_trim[,'b'] <- pmin(samples_block_info_trim[,'b'], 2)
##samplesPlot(samples_block_info_trim, ind=1:300000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'))
##dev.copy2pdf(file='~/Downloads/plot2.pdf')
## 
#### same thing, on the "early" time scale
##samplesPlot(samples_block_info_trim, ind=1:2000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
##dev.copy2pdf(file='~/Downloads/plot3.pdf')
## 
##samplesPlot(samples_block_info_trim, ind=1:20000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
##dev.copy2pdf(file='~/Downloads/plot4.pdf')





## testing the new adaptive properties for covariance in RW_block sampler,
## not adapting until acceptance rate >= 0.15 at least once
## uses scaleHistory and propCovHistory

library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

code <- nimbleCode({
    a ~ dnorm(0, sd=100)
    b ~ dnorm(0, sd=100)
    c ~ dnorm(a/2, sd=.1)
    for(i in 1:n1) {
        y1[i] ~ dnorm(a, sd=0.1)
    }
    for(i in 1:n2) {
        y2[i] ~ dnorm(b, sd=0.1)
    }
})
n1 <- 10
n2 <- n1
y1 <- rnorm(n1, 3, 0.1)
y2 <- rnorm(n2, y1+5, 0.01)
constants <- list(n1=n1, n2=n2)
data <- list(y1=y1, y2=y2)
inits <- list(a = 0, b=0, c=0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
##conf$addSampler('a', 'RW')
conf$addSampler(c('a', 'b', 'c'), 'RW_block')
conf$printSamplers()
conf$addMonitors(c('a', 'b', 'c'))
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(4000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

samplesPlot(samples)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])


samplesPlot(samples)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])








## getting model running for Colin Lewis-Beck at Iowa State

## Create Constant W, G and F matrices to pass to NIMBLE
library(nimble)
y1 <- rep(0,2190)
W<-diag(2)
F1<-matrix(c(1,0), nrow = 1, ncol = 2)
G<-matrix(c(cos(2*pi/365),-sin(2*pi/365),sin(2*pi/365),cos(2*pi/365)),nrow=2,ncol=2)
smosCode <- nimbleCode({
    ## Initial Values for States
    c0[1]<-0
    c0[2]<-0
    m0[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]   ## CHANGE 1: add full indexing to m0[...]
    x0[1:2]~dmnorm(c0[1:2],m0[1:2,1:2])        ## CHANGE 2: add full indexing to x0[...]
    ## initial values
    m[1:2,1]<-G[1:2,1:2] %*% x0[1:2]           ## CHANGE 3: add indexing to m[...], inprod() only returns a scalar (not a vector or matrix) so use %*%
    var0[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2] ## CHANGE 4: add full indexing to var0[...]
    x[1:2,1]~dmnorm(m[1:2,1],var0[1:2,1:2])
    y[1]~dnorm(inprod(F1[1,1:2],x[1:2,1]),sigmaSquaredInvv)
    ## Model
    for (t in 2:T){
        mu[1:2, t] <- G[1:2,1:2] %*% x[1:2,t-1]     ## CHANGE 5: something was funny with your mu[] declaration.  In addition to not having
                                                    ## the required indexing, it should should be indexed by t
        sigs[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]  ## CHANGE 6: add full indexing to sigs[...]
        x[1:2,t] ~ dmnorm(mu[1:2,t],sigs[1:2,1:2])  ## CHANGE 7: added appropriate 't' indexing to mu[...]
        y[t] ~ dnorm(inprod(F1[1,1:2],x[1:2,t]),sigmaSquaredInvv)
    }
    ## Priors
    sigmaSquaredInvv~dgamma(5,20)
    sigmaSquaredInvw~dgamma(5,200)
})
smosModel<-nimbleModel(code=smosCode,name='2harm',constants=list(T=2190,pi=pi,W=W,G=G,F1=F1),data=list(y=y1), inits=list(sigmaSquaredInvv=1,sigmaSquaredInvw=1))
Cmodel <- compileNimble(smosModel)




## doing the Surgeries example for STAT 365 homework

library(nimble)
df <- read.csv('~/Downloads/Surgeries.csv')
df

code <- nimbleCode({
    for(i in 1:N) {
        p[i] ~ dbeta(1, 1)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(p = rep(0.5, dim(df)[1]))

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=1)

mean(samples[,1])
median(samples[,1])
sd(samples[,1])
quantile(samples[,1], probs=c(0.025, 0.975))



code <- nimbleCode({
    mu ~ dnorm(0, 0.0001)
    sigma ~ dunif(0, 1000)
    for(i in 1:N) {
        logit(p[i]) ~ dnorm(mu, sd=sigma)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(mu=0, sigma=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=1)

mean(samples[,1])
quantile(samples[,1], probs=c(0.025, 0.975))

expit(mean(samples[,1]))
expit(quantile(samples[,1], probs=c(0.025, 0.975)))

sum(df$Mortalities)/sum(df$Surgeries)


code <- nimbleCode({
    b0 ~ dnorm(0, 0.0001)
    b1 ~ dnorm(0, 0.0001)
    sigma ~ dunif(0, 1000)
    for(i in 1:N) {
        logit(p[i]) ~ dnorm(b0 + b1*x[i], sd=sigma)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
    old_procedure <- expit(b0)
    new_procedure <- expit(b0 + b1)
    mortality_difference <- new_procedure - old_procedure
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(b0=0, b1=0, sigma=1, x=df$Procedure)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addMonitors(c('old_procedure', 'new_procedure', 'mortality_difference'))
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=3:5, burnin=500)

mean(samples[,3])
quantile(samples[,3], probs=c(0.025, 0.975))

expit(mean(samples[,1]))
expit(quantile(samples[,1], probs=c(0.025, 0.975)))

sum(df$Mortalities)/sum(df$Surgeries)


## spatial NIMBLE model example from Abhirup Datta
library(nimble)
library(mvtnorm)
nimbleOptions(showCompilerOutput=TRUE)

## model code ###
gpCode <- nimbleCode({
    for(i in 1:n){ y[i] ~ dnorm(w[i],sd=tau)}
    Vw[,] <- sigma*chol(G[,])
    w[] ~ dmnorm(muw[], cholesky=Vw[,], prec_param=0)
    sigma ~ dunif(0, 100)
    tau ~ dunif(0, 10)
})

## data ###
set.seed(1)
n=100
sigs=4
w=as.vector(rmvnorm(1,rep(0,n),sigs*0.5^as.matrix(dist(1:n))))
Y=w+sqrt(0.5*sigs)*rnorm(n)

### setting up nimble model ####
gpModel <- nimbleModel(gpCode, constants = list(n = n, G=0.5^as.matrix(dist(1:n)), muw=rep(0,n)), 
    dimensions=list(y=n, Vw=c(n, n), w=n, Vy=c(n, n), G = c(n, n)), check=FALSE)

gpModel$setData(list(y = Y))

gpModel$sigma <- 1
gpModel$tau <- 0.1
gpModel$w <- rep(0,n)

CgpModel <- compileNimble(gpModel)

### MCMC ####

gpconf <- configureMCMC(CgpModel, print=TRUE)


gpconf$printSamplers()
gpconf$removeSamplers('w')
gpconf$printSamplers()
for(node in gpModel$expandNodeNames('w', returnScalarComponents=TRUE)) {
    gpconf$addSampler(node, 'RW')
}
gpconf$printSamplers()



gpconf$addMonitors(c('w'))	

gpMCMC <- buildMCMC(gpconf)

CgpMCMC <- compileNimble(gpMCMC, project = CgpModel)

##Nmcmc=100000
Nmcmc=1000
set.seed(1)
CgpMCMC$run(Nmcmc)

gpMCMCsamples <- as.matrix(CgpMCMC$mvSamples)

plot(gpMCMCsamples[,"w[1]"])

dimnames(gpMCMCsamples)

samplesPlot(gpMCMCsamples, var=3:10)


## implementing RJMCMC for a simple 2 models, logistic regression
## yi ~ bernoulli(logit(p) = a + bx)
set.seed(0)
n <- 50
x <- runif(n, -10, 10)
a <- .3
b <- 0.1
p <- 1/(1 + exp(-1*(a + b*x)))
y <- rbinom(n, prob=p, size=1)
##plot(x,y)
m1 <- glm(y~1, family=binomial())
m1
m2 <- glm(y~x, family=binomial())
m2
aic1 <- AIC(m1)  ## for this AIC function, "lower values indicate better fit"
aic2 <- AIC(m2)
aics <- c(aic1, aic2)
aic_min <- min(aics)
daics <- aics - aic_min
w1_aic <- exp(-daics[1]/2) / sum(exp(-daics/2))
w2_aic <- exp(-daics[2]/2) / sum(exp(-daics/2))
bic1 <- BIC(m1)  ## for this AIC function, "lower values indicate better fit"
bic2 <- BIC(m2)
bics <- c(bic1, bic2)
bic_min <- min(bics)
dbics <- bics - bic_min
w1_bic <- exp(-dbics[1]/2) / sum(exp(-dbics/2))
w2_bic <- exp(-dbics[2]/2) / sum(exp(-dbics/2))
w1_aic
w1_bic

## NEED TO KEEP WORKING ON LOGISTIC EXAMPLE FROM HERE
prior_sd <- 10
rj_proposal_sd <- 10
prior_m1 <- 0.5
{
    decide <- function(lMHR) {
        if(is.nan(lMHR)) return(FALSE)
        if(log(runif(1,0,1)) < lMHR) return(TRUE) else return(FALSE)
    }
    prior <- function(x) return(dnorm(x, 0, sd=prior_sd, log=TRUE))
    ll <- function(m, a, b, y, x) {
        if(m==1) return(sum(dnorm(y, a, 1, log=TRUE)))
        if(m==2) return(sum(dnorm(y, a+b*x, 1, log=TRUE)))
    }
    updateM <- function(mab, y, x) {
        ##browser()
        if(mab[1] == 1) {
            bprop <- rnorm(1, 0, rj_proposal_sd)   ## proposal for new slope 'b'
            lMHR <- prior(bprop) + ll(2,mab[2],bprop,y,x)  -  (ll(1,mab[2],0,y,x) + dnorm(bprop,0,rj_proposal_sd,log=TRUE)) + log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(2,mab[2],bprop)) else return(mab)
        }
        if(mab[1] == 2) {
            lMHR <- ll(1,mab[2],0,y,x) + dnorm(mab[3],0,rj_proposal_sd,log=TRUE)  -  (prior(mab[3]) + ll(2,mab[2],mab[3],y,x)) - log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(1,mab[2],0)) else return(mab)
        }
    }
    updateAB <- function(mab, y, x) {
        ## update a
        sd.aprop <- 0.1
        aprop <- rnorm(1,mab[2],sd.aprop)
        lMHR <- prior(aprop) + ll(mab[1],aprop,mab[3],y,x) - prior(mab[2]) - ll(mab[1],mab[2],mab[3],y,x)
        if(decide(lMHR)) mab <- c(mab[1],aprop,mab[3])
        ## update b
        if(mab[1] == 2) {
            sd.bprop <- 0.1
            bprop <- rnorm(1,mab[3],sd.bprop)
            lMHR <- prior(bprop) + ll(mab[1],mab[2],bprop,y,x) - prior(mab[3]) - ll(mab[1],mab[2],mab[3],y,x)
            if(decide(lMHR)) mab <- c(mab[1],mab[2],bprop) else mab
        }
        return(mab)
    }
}

set.seed(0)
iter <- 100000
mab <- c(1, 0, 0)  ## c(1, 0, 0)
##samp <- data.frame(a=rep(NA,iter), b=NA, c=NA)
samp <- array(NA, c(iter,3))
colnames(samp) <- c('m', 'a', 'b')
for(i in 1:iter) {
    mab <- updateM(mab, y, x)
    mab <- updateAB(mab, y, x)
    samp[i,] <- mab
}
df <- as.data.frame(samp)

mean(df$m==1)
##head(df,10)

c(mean(df$a[df$m==1]), m1$coef[1])
c(mean(df$a[df$m==2]), m2$coef[1])
c(mean(df$b[df$m==2]), m2$coef[2])

##plot(df$a, type='l')
##plot(1:iter, df$b)

prod(dnorm(y, 0, sd=sqrt(prior_sd^2 + 1^2)/10))   ## wrong for p(y|m1), don't know why

pym1 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd), 1))))
pym2 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd) + rnorm(1,0,prior_sd)*x, 1))))
pym1
pym2
ppm1 <- prior_m1*pym1 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm2 <- (1-prior_m1)*pym2 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm1
ppm2

w1_aic
w1_bic
mean(df$m==1)
ppm1

library(BMA)
bma <- bic.glm(x=data.frame(x=x), y=y, glm.family='gaussian')
##class(bma)
##str(bma)
bma$postprob[1]
bma$deviance
bma$label
bma$size
bma$which
bma$probne0
bma$postmean
bma$condpostmean
bma$mle





## implementing RJMCMC for a simple 2 models
## yi ~ N(a + bx, 1)
set.seed(0)
n <- 10
x <- runif(n, -1, 1)
a <- .3
b <- 0.1
y <- rnorm(n, a+b*x, 1)
##plot(x,y)
m1 <- lm(y~1)
m1
m2 <- lm(y~x)
m2
aic1 <- AIC(m1)  ## for this AIC function, "lower values indicate better fit"
aic2 <- AIC(m2)
aics <- c(aic1, aic2)
aic_min <- min(aics)
daics <- aics - aic_min
w1_aic <- exp(-daics[1]/2) / sum(exp(-daics/2))
w2_aic <- exp(-daics[2]/2) / sum(exp(-daics/2))
bic1 <- BIC(m1)  ## for this AIC function, "lower values indicate better fit"
bic2 <- BIC(m2)
bics <- c(bic1, bic2)
bic_min <- min(bics)
dbics <- bics - bic_min
w1_bic <- exp(-dbics[1]/2) / sum(exp(-dbics/2))
w2_bic <- exp(-dbics[2]/2) / sum(exp(-dbics/2))
w1_aic
w1_bic

prior_sd <- 10
rj_proposal_sd <- 10
prior_m1 <- 0.5
{
    decide <- function(lMHR) {
        if(is.nan(lMHR)) return(FALSE)
        if(log(runif(1,0,1)) < lMHR) return(TRUE) else return(FALSE)
    }
    prior <- function(x) return(dnorm(x, 0, sd=prior_sd, log=TRUE))
    ll <- function(m, a, b, y, x) {
        if(m==1) return(sum(dnorm(y, a, 1, log=TRUE)))
        if(m==2) return(sum(dnorm(y, a+b*x, 1, log=TRUE)))
    }
    updateM <- function(mab, y, x) {
        ##browser()
        if(mab[1] == 1) {
            bprop <- rnorm(1, 0, rj_proposal_sd)   ## proposal for new slope 'b'
            lMHR <- prior(bprop) + ll(2,mab[2],bprop,y,x)  -  (ll(1,mab[2],0,y,x) + dnorm(bprop,0,rj_proposal_sd,log=TRUE)) + log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(2,mab[2],bprop)) else return(mab)
        }
        if(mab[1] == 2) {
            lMHR <- ll(1,mab[2],0,y,x) + dnorm(mab[3],0,rj_proposal_sd,log=TRUE)  -  (prior(mab[3]) + ll(2,mab[2],mab[3],y,x)) - log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(1,mab[2],0)) else return(mab)
        }
    }
    updateAB <- function(mab, y, x) {
        ## update a
        sd.aprop <- 0.1
        aprop <- rnorm(1,mab[2],sd.aprop)
        lMHR <- prior(aprop) + ll(mab[1],aprop,mab[3],y,x) - prior(mab[2]) - ll(mab[1],mab[2],mab[3],y,x)
        if(decide(lMHR)) mab <- c(mab[1],aprop,mab[3])
        ## update b
        if(mab[1] == 2) {
            sd.bprop <- 0.1
            bprop <- rnorm(1,mab[3],sd.bprop)
            lMHR <- prior(bprop) + ll(mab[1],mab[2],bprop,y,x) - prior(mab[3]) - ll(mab[1],mab[2],mab[3],y,x)
            if(decide(lMHR)) mab <- c(mab[1],mab[2],bprop) else mab
        }
        return(mab)
    }
}

set.seed(0)
iter <- 100000
mab <- c(1, 0, 0)  ## c(1, 0, 0)
##samp <- data.frame(a=rep(NA,iter), b=NA, c=NA)
samp <- array(NA, c(iter,3))
colnames(samp) <- c('m', 'a', 'b')
for(i in 1:iter) {
    mab <- updateM(mab, y, x)
    mab <- updateAB(mab, y, x)
    samp[i,] <- mab
}
df <- as.data.frame(samp)

mean(df$m==1)
##head(df,10)

c(mean(df$a[df$m==1]), m1$coef[1])
c(mean(df$a[df$m==2]), m2$coef[1])
c(mean(df$b[df$m==2]), m2$coef[2])

##plot(df$a, type='l')
##plot(1:iter, df$b)

prod(dnorm(y, 0, sd=sqrt(prior_sd^2 + 1^2)/10))   ## wrong for p(y|m1), don't know why

pym1 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd), 1))))
pym2 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd) + rnorm(1,0,prior_sd)*x, 1))))
pym1
pym2
ppm1 <- prior_m1*pym1 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm2 <- (1-prior_m1)*pym2 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm1
ppm2

w1_aic
w1_bic
mean(df$m==1)
ppm1

library(BMA)
bma <- bic.glm(x=data.frame(x=x), y=y, glm.family='gaussian')
##class(bma)
##str(bma)
bma$postprob[1]
bma$deviance
bma$label
bma$size
bma$which
bma$probne0
bma$postmean
bma$condpostmean
bma$mle





## testing samplers used by JAGS for multivariate normal (dmnorm) nodes
## also, modified code in JAGS MNormal.cc (block sampler)
## to see what's going on

## (see http://danielturek.github.io/public/jags_block_adaptation/jags_block_adaptation.html)
## to build jags from source:
##   $ cd ~/Downloads/JAGS-4.2.0
##   $ ./configure; make -j 4; sudo make install
## then restart R

##remove.packages('rjags')
##install.packages('rjags')
library(rjags)
code <- quote({
    C[1,1] <- 1
    C[1,2] <- 0
    C[2,1] <- 0
    C[2,2] <- 1
    mu[1] <- 0
    mu[2] <- 0
    y[1:2] ~ dmnorm(mu[1:2], C[1:2, 1:2])
    d[1] <- exp(y[1]) + 4 * pow(y[1], 2)
    d[2] <- exp(y[2]) + 10 * pow(y[1], 2) + 3
    dat1 ~ dexp(d[1])
    dat2 ~ dexp(d[2])
})
constants <- list()
data <- list(dat1=3, dat2=6)
inits <- list(y = c(0,0))
##Rmodel <- nimbleModel(code, constants, data, inits)
##conf <- configureMCMC(Rmodel)
##conf$printSamplers()
niter <- 100
monitorVars <- c('y')
constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)
set.seed(0); jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE, n.adapt=40)

set.seed(0); jags_out <- rjags::coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1)

list.samplers(jags_mod)



## doing the midterm question about normal mean hypothesis test
## "weights of sheep raised on a farm"
## now also with predictive distribution

library(nimble)

y <- c(78, 81, 77, 76, 75, 74, 78, 75, 77, 75)
n <- length(y)
np <- 5

code <- nimbleCode({
    mu ~ dnorm(75, sd=10)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, sd=3)
    }
    for(i in 1:np) {
        p[i] ~ dnorm(mu, sd=3)
    }
    pmean <- mean(p[1:np])
})
constants <- list(n=n, np=np)
data <- list(y=y)
inits <- list(mu = 75, p = rep(0,np))

Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$addMonitors(c('mu', 'p', 'pmean'))
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, resetFunctions=TRUE)

set.seed(0)
Cmcmc$run(1000000)
samples <- as.matrix(Cmcmc$mvSamples)
colnames(samples)
dim(samples)
samples[1:20,]
apply(samples, 2, mean)
apply(samples, 2, sd)
apply(samples, 2, var)
1/apply(samples, 2, var)

tau0 <- 1/100
mu0 <- 75
ybar <- mean(y)
tau <- 1/9
tau_n <- tau0 + n*tau
mu_n <- (mu0*tau0 + ybar*n*tau)/(tau0+n*tau)
mu_n
tau_n
1/tau_n
var(samples[,'mu'])

v <- 1/tau
v_n <- 1/tau_n
v_p1 <- v_n + v
v_p1
apply(samples, 2, var)

v_p1/5
var(samples[,'pmean'])
v_n/5 + v
v_n + v/5

var(apply(matrix(sample(as.numeric(samples[,2:6])), ncol=5), 1, mean))




## working on used cars regression example for Bayes course 365
df <- read.csv('~/Downloads/UsedCars.csv')
str(df)
hist(df$Age)
hist(df$HP)
dim(df)
head(df)

m <- lm(Price ~ Age + HP | Type, df=df)
summary(m)

summary(lm(Price ~ Age + HP, df=subset(df, Type==0)))
summary(lm(Price ~ Age + HP, df=subset(df, Type==1)))
summary(lm(Price ~ Age + HP + Type, df=df))


sd(residuals(lm(Price ~ Age + HP, df=subset(df, Type==0))))
sd(residuals(lm(Price ~ Age + HP, df=subset(df, Type==1))))




## deep-dive deep dive into correlated state-space state space model
## to figure out why the block sampler is taking so long to adapt.
## recreating and fixing Dao's issue with the correlated SSM,
## where 'a' and 'b' don't mix until after 150,000 iterations
library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
##
code <- nimbleCode({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(1e-04, 1)
    sigOE ~ dunif(1e-04, 1)
    x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for (i in 2:t) {
        x[i] ~ dnorm(x[i - 1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
constants <- list(t = 100)
data <- list(y = c(20.24405,20.57693,20.49357,20.34159,20.45759,20.43326,20.20554,20.12860,20.14756,20.20781,20.23022,20.26766,20.22984,20.37703,20.13641,20.05309,19.95709,20.19303,20.30562,20.54443,20.91010,20.70580,20.42344,20.19795,20.28816,20.31894,20.76939,20.77023,20.83486,20.29335,20.40990,20.19601,20.04083,19.76056,19.80810,19.83129,19.69174,19.90069,19.87623,19.63371,19.62360,19.72630,19.64450,19.86779,20.17104,20.34797,20.32968,20.48027,20.46694,20.47006,20.51676,20.40695,20.18715,19.97552,19.88331,19.67831,19.74702,19.47502,19.24408,19.37179,19.38277,19.15034,19.08723,19.37051,19.14274,19.46433,19.62459,19.77971,19.54194,19.39081,19.61621,19.51307,19.34745,19.17019,19.26829,19.58943,19.77143,19.83582,19.71198,19.67746,19.75053,20.40197,20.49363,20.37079,20.19005,20.55862,20.48523,20.33071,19.97069,19.79758,19.83811,19.79728,19.86277,19.86836,19.92481,19.88095,20.24899,20.55165,20.22707,20.11235))
inits <- list(a = 0.95, b=1, sigPN = 0.2, sigOE=0.05, x = c(20.26036,20.51331,20.57057,20.35633,20.33736,20.47321,20.22002,20.14917,20.19216,20.26969,20.21135,20.22745,20.20466,20.41158,20.13408,20.08023,19.98956,20.13543,20.32709,20.55840,20.88206,20.74740,20.47671,20.14012,20.29953,20.33778,20.80916,20.75773,20.84349,20.35654,20.41045,20.20180,20.02872,19.74226,19.80483,19.81842,19.69770,19.84564,19.88211,19.70559,19.56090,19.73728,19.66545,19.88158,20.13870,20.39163,20.37372,20.47429,20.39414,20.42024,20.55560,20.40462,20.15831,19.89425,19.79939,19.72692,19.74565,19.42233,19.22730,19.36489,19.37289,19.19050,19.00823,19.35738,19.14293,19.48812,19.67329,19.82750,19.58979,19.43634,19.61278,19.56739,19.38584,19.19260,19.32732,19.65500,19.65295,19.84843,19.68285,19.69620,19.77497,20.31795,20.45797,20.32650,20.24045,20.60507,20.51597,20.30076,19.98100,19.86709,19.85965,19.74822,19.86730,19.90523,19.86970,19.87286,20.28417,20.46212,20.22618,20.13689))
##
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## [1] 183.3436
##

conf <- configureMCMC(Rmodel, nodes = NULL)

conf$addSampler(c('a', 'b'), 'RW_block')
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(adaptInterval=100))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(1,-.99,-0.99,1), c(2,2))))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(1,-.99,-0.99,1), c(2,2)), scale=0.01))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(0.001709168, -0.0341986, -0.0341986, 0.6844844), c(2,2))))


conf$printSamplers(c('a','b'))

conf$addSampler('sigOE', 'RW')
conf$addSampler('sigPN', 'RW')
for(node in Rmodel$expandNodeNames('x'))
    conf$addSampler(node, 'RW')
conf$resetMonitors()
conf$addMonitors(c('a', 'b', 'sigOE', 'sigPN'))
##conf$getMonitors()
##
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()   ## [1] 183.3436
nimbleOptions(showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 500000
system.time(Cmcmc$run(niter))
##
samples <- as.matrix(Cmcmc$mvSamples)
dim(samples)
dimnames(samples)
##

## Plot1
dev.new(width=8, height=5)
par(mfrow=c(1,2))
plot(samples[1:300000,'a'], type='l', ylab='a')
plot(samples[1:300000,'b'], type='l', ylab='b')
par(mfrow=c(1,1))
getwd()
dev.copy2pdf(file='~/Downloads/plot1.pdf')

samplesPlot(samples, var=c('b','a'), ind=1:300000)

##samplesPlot(samples, var=c('sigOE','sigPN'))

##xs <- Rmodel$expandNodeNames('x')
##samplesPlot(samples, var=xs[ 1:10])
##samplesPlot(samples, var=xs[11:20])
##samplesPlot(samples, var=xs[21:30])
##samplesPlot(samples, var=xs[31:40])

##i <-  4:13   ## x[ 1]...x[10]
##i <- 14:23   ## x[11]...x[20]
##i <-  1:8    ## (a,b), sigOE, sigPN, x[1]...x[5]
##i <-  2:3    ## sigOE, sigPN

## all the latent state scales follow sigOE scale
##i <-  c(2,4:8)    ## sigOE, x[1]...x[5]
##nms <- sapply(conf$samplerConfs[i], function(x) x$target)
##nms
##ar <- do.call(cbind, lapply(Cmcmc$samplerFunctions$contentsList[i], function(x) x$scaleHistory))
##colnames(ar) <- nms
##samplesPlot(ar)
##dim(ar)
##samplesPlot(ar, burnin=500)

block_scales <- Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
block_propCovHistory <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
## create block_propCovScale
block_propCovScale <- block_propCovHistory
for(i in 1:length(block_scales))   block_propCovScale[i,,] <- block_scales[i] * block_propCovHistory[i,,]
##dim(block_propCovScale)
block_scale_a <- apply(block_propCovScale, 1, function(x) sqrt(x[1,1]))
block_scale_b <- apply(block_propCovScale, 1, function(x) sqrt(x[2,2]))
block_cors <- apply(block_propCovHistory, 1, function(x) cov2cor(x)[1,2])
ar <- cbind(block_scales, block_scale_a, block_scale_b, block_cors)
colnames(ar) <- c('scale', 'sig_a', 'sig_b', 'cor')
samplesPlot(ar)

## final constant that scale approaches:
block_scales[length(block_scales)]
##propCov adapts very nicely to true covariance between 'a' and 'b'
cov(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
block_propCovHistory[dim(ar)[1],,]
## final adapted (and scaled) proposal corrleation is very accurate:
cor(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
cov2cor(block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])
## final adapted (and scaled) proposal standard deviations for 'a' and 'b':
sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[1,1])
sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[2,2])

## expand cor, sig_a, and sig_b by adaptInterval:
length(block_scale_a)
aI <- 200
block_scale_a_ex <- rep(block_scale_a, each=aI)
block_scale_b_ex <- rep(block_scale_b, each=aI)
block_cors_ex    <- rep(block_cors,    each=aI)
block_scales_ex  <- rep(block_scales,  each=aI)
samples_block_info <- cbind(samples[,'a'], samples[,'b'], block_scale_a_ex, block_scale_b_ex, block_cors_ex, block_scales_ex)
dimnames(samples_block_info)[[2]] <- c('a', 'b', 'sig_a', 'sig_b', 'cor', 'scale')

##samplesPlot(samples_block_info, ind=1:300000, var=c('b', 'a'))
##samplesPlot(samples_block_info, ind=1:300000, var=c('sig_a', 'sig_b', 'cor', 'scale'))

## this one is best:
## artificially trim 'b' samples:
samples_block_info_trim <- samples_block_info
samples_block_info_trim[,'b'] <- pmin(samples_block_info_trim[,'b'], 2)
samplesPlot(samples_block_info_trim, ind=1:300000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'))
dev.copy2pdf(file='~/Downloads/plot2.pdf')

## same thing, on the "early" time scale
samplesPlot(samples_block_info_trim, ind=1:2000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
dev.copy2pdf(file='~/Downloads/plot3.pdf')

samplesPlot(samples_block_info_trim, ind=1:20000, var=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
dev.copy2pdf(file='~/Downloads/plot4.pdf')







## testing addition of scaleHistory and propCovHistory
## back into RW and RW_block samplers


library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0, b=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
conf$addSampler('a', 'RW')
conf$addSampler(c('a', 'b'), 'RW_block')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])

set.seed(0)
Cmcmc$run(1000)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
a <- Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])







## playing with colorRamp and colorRampPalette
## to get gradients of colors in plots

FLcrime <- read.csv("~/github/courses/stat201/data/fl_crime.csv")
head(FLcrime[, 1:4])
FLcrime <- FLcrime[, 2:4]
colnames(FLcrime) <- c("Crime", "Education", "Urbanization")
head(FLcrime)
pairs(FLcrime, panel=panel.smooth)
pairs(FLcrime)
plot(FLcrime$Education, FLcrime$Crime)
hist(FLcrime$Urbanization)
boxplot(FLcrime$Urbanization)
quantile(FLcrime$Urbanization)

plot(FLcrime$Education, FLcrime$Crime, pch=19)
with(subset(FLcrime, Urbanization<20), points(Education, Crime, col='green', pch=19))
with(subset(FLcrime, Urbanization>80), points(Education, Crime, col='red', pch=19))


pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(nrow(FLcrime))[cut(FLcrime$Urbanization, nrow(FLcrime))]
plot(FLcrime$Education, FLcrime$Crime, col=cols, pch=19)
with(subset(FLcrime, Urbanization<20), points(Education, Crime, col='green', pch=19))
with(subset(FLcrime, Urbanization>80), points(Education, Crime, col='red', pch=19))


x <- runif(100)
dat <- data.frame(x = x,y = x^2 + 1)

?colorRampPalette
rbPal <- colorRampPalette(c('red','blue'))

cc <- rbPal(100)[cut(dat$x, 100)]
plot(dat$x, dat$y, col = cc)

cc <- colorRampPalette(1:2)
plot(1:10, 1:10, col=cc(1:10/10))


plot(FLcrime$Education, FLcrime$Crime)

## practice problem 6 for STAT 365

## likelihood: yi ~ Normal(mu, sigma)
n <- 4
y <- c(100, 110, 112, 118)
sigma <- 10

## prior: mu ~ Normal(mu0, sigma0)
mu0 <- 100
sigma0 <- 10

## derive precisions
tau0 <- sigma0^-2
tau <- sigma^-2

## posterior distribution
tau_n <- tau0 + n*tau
mu_n <- (tau0*mu0 + n*tau*mean(y)) / tau_n

sigma_n <- 1/sqrt(tau_n)

mu_n
tau_n
sigma_n

## prior, likelihood, posterior plots

## 95% BCI

## one-sided hypothesis test
## H0: mu|y <= 100
## Ha: mu|y >  100

## two-sided hypothesis test
## H0: mu|y  = 100
## Ha: mu|y != 100

## samples from posterior

## 95% BCI

## one-sided hypothesis test
## H0: mu|y <= 100
## Ha: mu|y >  100


## two-sided hypothesis test
## H0: mu|y  = 100
## Ha: mu|y != 100


.





file <- 'HurricaneDamage.csv'
data <- read.csv(paste0('~/github/courses/stat201/data/', file))
Year <- data$Year
Damage <- data$Damage
Year2 <- Year[-1]
Damage2 <- Damage[-1]
m2 <- lm(Damage2~Year2)

res <- residuals(m2)

par(mfrow=c(1,1))
hist(res, breaks=12)

par(mfrow=c(2,1))
plot(Year2, res)
plot(fitted(m2), res)






## modifying hurricanes.csv data file into Hurricane_Damage.csv,
## for use in the STAT201 miderm
file <- 'hurricanes.csv'
data <- read.csv(paste0('~/Downloads/', file))
str(data)
names(data)
names(data)[1] <- 'Rank'
names(data)[2] <- 'Tropical.Cyclone'
names(data)[3] <- 'Year'
names(data)
dim(data)
data <- data[-(2:5),]
dim(data)
mtemp <- lm(data$Damage ~ data$Year)
coef(mtemp)
mtemp <- lm(data$Damage[-1] ~ data$Year[-1])
coef(mtemp)
fitted(mtemp)
newY <- c(data$Damage[1], 0.4*data$Damage[-1] + 0.6*fitted(mtemp))
data$Damage <- newY
write.csv(data, file='~/Downloads/HurricaneDamage.csv', row.names=FALSE)


file <- 'HurricaneDamage.csv'
data <- read.csv(paste0('~/Downloads/', file))
str(data)
x <- data$Year
y <- data$Damage
m <- lm(y~x)
summary(m)
plot(x, y)
abline(m, col='red')
coef(m)

sort(y)
which(y>40000)
i <- which(y>40000)
x2 <- x[-i]
y2 <- y[-i]
m2 <- lm(y2~x2)
summary(m2)
plot(x2, y2)
abline(m2, col='red')
coef(m2)
cor(x2,y2)
cor(x2,y2)^2





## doing LAX flight departures problem
## from the STAT201 midterm

x <- 3648
y <- 25843
w <- 3407
z <- 26134
N <- x+y+w+z

year <- c(rep(2014, x+y), rep(2015, w+z))
delay <- c(rep('ayes',x), rep('no',y), rep('ayes',w), rep('no',z))
tab <- table(year, delay)
tab

## (a) display contingency table
rbind(cbind(tab, margin.table(tab,1)), c(margin.table(tab,2), N))
##     ayes    no      
##2014 3648 25843 29491
##2015 3183 26284 29467
##     6831 52127 58958

## (b) percentage in 2014?
(x+y)/N * 100
prop.table(margin.table(tab, 1))[1]
## 50.02035 %

## (c) percentage of delayed departures in 2015?
w/(w+x) * 100
## 46.5964 %

## (d) diff. in prop, delays in 2015 relative to 2014?
ptab <- prop.table(tab, 1)
ptab[2,1] - ptab[1,1]
## -0.01567962 = -1.567962 %

## (e) interpret this diff. in prop.
## The fraction of LAX December 2015 departures that were delayed was 1.57 *percentage points* lower than the fraction of LAX December 2016 departures that were delayed.

iter <- 1000
data <- data.frame(year=year, delay=delay)
nr <- nrow(data)
ratio <- numeric(iter)
for(i in 1:iter) {
    ind <- sample(1:nr)
    tab <- table(data$year, data$delay[ind])
    cond.tab <- prop.table(tab, 1)
    ratio[i] <- cond.tab[2,1] - cond.tab[1,1]
}

hist(ratio, breaks=15)
tab <- table(data$year, data$delay)
cond.tab <- prop.table(tab, 1)
observed_ratio <- cond.tab[2,1] - cond.tab[1,1]
abline(v = observed_ratio, col = 'red', lwd = 2)







## testing seq_along in NIMBLE run code
library(nimble)

nfDef <- nimbleFunction(
    setup = function() {
        a <- 1:10
    },
    run = function() {
        for(i in seq_along(a)) {
            print(i)
        }
    }
)

Rnf <- nfDef()

Rnf$run()

Cnf <- compileNimble(Rnf)

Cnf$run()




## trying out MCMC with two Gibbs samplers for Normal

y <- c(4.3, 2.5, 3.2, 3.8, 2.9, 3.1, 4.2, 4.0)
n <- length(y)
n
y
mean(y)
sd(y)

## mu ~ dnorm(mu0=0, tau0=0.001)
## tau ~ dgamma(r0=0.001, v0=0.001)
## yi ~ dnorm(mu, tau)

mu0 <- 0
tau0 <- 0.001
r0 <- 0.001
v0 <- 0.001

## iter sampling iterations to run
iter <- 100000
samp <- cbind(mu = rep(NA,iter), tau = rep(NA,iter))

## inits: mu=0, tau=1
mu <- 0
tau <- 1

set.seed(0)
for(i in 1:iter) {
    mu <- rnorm(1, (tau0*mu0+tau*sum(y))/(tau0+n*tau), (tau0+n*tau)^-0.5)
    tau <- rgamma(1, r0+n/2, v0+0.5*sum((y-mu)^2))
    samp[i, 1] <- mu
    samp[i, 2] <- tau
}


1/2 * sum((y-mu)^2)
n/2 * (mean(y)-mu)^2


head(samp)

samp <- cbind(samp, sd = samp[,'tau']^-0.5)

apply(samp, 2, mean)
apply(samp, 2, median)
mean(y)
sd(y)
sd(y)^-2

samplesPlot(samp)

library(nimble)

code <- nimbleCode({
    mu ~ dnorm(0, 0.001)
    tau ~ dgamma(0.001, 0.001)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, tau)
    }
})
constants <- list(n=n)
data <- list(y=y)
inits <- list(mu=0, tau=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

head(Rsamples, 10)
head(Csamples, 10)
head(samp, 10)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)




## fixing bug in conjugacy

library(nimble)
set.seed(0)
n <- 10
##y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=1))
y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=2), rnorm(n,20,sd=3))

code <- nimbleCode({
    for(i in 1:3) {
        ##sig[i] ~ dunif(0, 100)
        tau[i] ~ dgamma(0.001, 0.001)   ## using TAU
    }
    mu[1] ~ dnorm(0, sd=1000)
    mu[2] <- mu[1] + delta1
    delta1 ~ T(dnorm(0, sd=1000), 0, 1000)
    mu[3] <- mu[2] + delta2
    delta2 ~ T(dnorm(0, sd=1000), 0, 1000)
    for(i in 1:N) {
        ##means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2]
        means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2] + equals(z[i],3)*mu[3]
        ##sigmas[i] <- equals(z[i],1)*sig[1] + equals(z[i],2)*sig[2] + equals(z[i],3)*sig[3]
        taus[i] <- equals(z[i],1)*tau[1] + equals(z[i],2)*tau[2] + equals(z[i],3)*tau[3]   ## using TAU
        ##y[i] ~ dnorm(means[i], sd=sigmas[i])
        y[i] ~ dnorm(means[i], taus[i])   ## using TAU
        z[i] ~ dcat(pi[1:3])
    }
    for(i in 1:3) {
        pi0[i] ~ dgamma(1, 1)
        pi[i] <- pi0[i] / (pi0[1] + pi0[2] + pi0[3])
    }
})

N <- length(y)
constants <- list(N=N)
data <- list(y=y)
##inits <- list(mu=c(1,2,3), delta1=1, delta2=1, pi0=c(1,1,1), sig=c(1,1,1), z=rep(1:3, each=n))
inits <- list(mu=c(1,2,3), delta1=1, delta2=1, pi0=c(1,1,1), tau=c(1,1,1), z=rep(1:3, each=n))  ## using TAU
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$pi0
Rmodel$pi
Rmodel$z
Rmodel$mu
Rmodel$means
Rmodel$sigmas
Rmodel$taus

##undebug(Rmodel$checkConjugacy)
##Rmodel$checkConjugacy('mu[1]')
## 
##undebug(Rmodel$checkConjugacy2)
##Rmodel$checkConjugacy2('mu[1]')

conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, useConjugacy=FALSE)
conf$resetMonitors()
##conf$addMonitors(c('mu', 'z', 'sig'))
conf$addMonitors(c('mu', 'z', 'tau'))
conf$printSamplers('mu')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
calculate(Rmodel)
calculate(Cmodel)

##iter <- 20
## 
##set.seed(0)
##Rmcmc$run(iter)
## 
##set.seed(0)
##Cmcmc$run(iter)
## 
##Rsamples <- as.matrix(Rmcmc$mvSamples)
##Csamples <- as.matrix(Cmcmc$mvSamples)
## 
##sampNames <- colnames(Rsamples)
## 
##Rsamples[, sampNames]
##Csamples[, sampNames]
##Rsamples[, sampNames] - Csamples[, sampNames]

set.seed(0)
Cmcmc$run(50000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

mus <- Rmodel$expandNodeNames(c('mu'))
sigs <- Rmodel$expandNodeNames(c('sig'))
taus <- Rmodel$expandNodeNames(c('tau'))
zs <- Rmodel$expandNodeNames(c('z'))

ind <- 1:dim(samples)[1]
ind <- 1:5000
samplesPlot(samples, col = mus, ind = ind)
samplesPlot(samples, col = sigs, ind = ind)
samplesPlot(samples, col = taus, ind = ind)
samplesPlot(samples, col = c('mu[3]','sig[3]'), ind = ind)
samplesPlot(samples, col = c('mu[3]','tau[3]'), ind = ind)
samplesPlot(samples, col = c(        'tau[3]'), ind = ind)

samplesPlot(samples, col = zs, ind = ind)

samplez <- samples[, zs]
samplezsum <- cbind(samplez, z1=apply(samplez, 1, function(x) sum(x==1)), z2=apply(samplez, 1, function(x) sum(x==2)), z3=apply(samplez, 1, function(x) sum(x==3)))
dim(samplez)
dim(samplezsum)

samplesPlot(samplezsum, col = c('z1','z2','z3'), ind=1:5000)
samplesPlot(samplezsum, col = c('z1','z2','z3'), ind=1:15000)

mean(y[1:n])
mean(y[(n+1):(2*n)])
mean(y[(2*n+1):(3*n)])



## doing the midterm question about normal mean hypothesis test
## "weights of sheep raised on a farm"

tau0 <- 1/100
mu0 <- 75

y <- c(78, 81, 77, 76, 75, 74, 78, 75, 77, 75)
n <- length(y)
ybar <- mean(y)
ybar
tau <- 1/9

taup <- tau0 + n*tau
mup <- (mu0*tau0 + ybar*n*tau)/(tau0+n*tau)

mup
taup
pnorm(75, mup, 1/sqrt(taup))


library(nimble)

code <- nimbleCode({
    mu ~ dnorm(75, sd=10)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, sd=3)
    }
})
constants <- list(n=n)
data <- list(y=y)
inits <- list(mu = 75)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(500000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
mup
mean(samples[,1] < 75)
pnorm(75, mup, 1/sqrt(taup))






## inplementing a finite mixture model in NIMBLE

library(nimble)

set.seed(0)
n <- 500
y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=2), rnorm(n,20,sd=3))
##hist(y, breaks=30)

code <- nimbleCode({
    for(i in 1:3) {
        sig[i] ~ dunif(0, 100)
    }
    mu[1] ~ dnorm(0, sd=1000)
    mu[2] <- mu[1] + delta1
    delta1 ~ T(dnorm(0, sd=1000), 0, 1000)
    mu[3] <- mu[2] + delta2
    delta2 ~ T(dnorm(0, sd=1000), 0, 1000)
    for(i in 1:N) {
        means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2] + equals(z[i],3)*mu[3]
        sigmas[i] <- equals(z[i],1)*sig[1] + equals(z[i],2)*sig[2] + equals(z[i],3)*sig[3]
        y[i] ~ dnorm(means[i], sd=sigmas[i])
        z[i] ~ dcat(pi[1:3])
    }
    ##pi[1:3] ~ ddirch(alpha[1:3])
    for(i in 1:3) {
        pi0[i] ~ dgamma(1, 1)
        pi[i] <- pi0[i] / (pi0[1] + pi0[2] + pi0[3])
    }
})

N <- length(y)
constants <- list(N=N)
data <- list(y=y)
inits <- list(sig=rep(1,3), mu=c(1,2,3), delta1=1, delta2=1, z=rep(1,N), pi0=c(1,1,1))

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$pi0
Rmodel$pi
Rmodel$z
Rmodel$mu
Rmodel$means
Rmodel$sigmas

conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, useConjugacy = FALSE)

conf$printSamplers()
conf$printSamplers('mu')
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

calculate(Rmodel)
calculate(Cmodel)

##set.seed(0)
##Rmcmc$run(2)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


                   

## bootstrapping the slope 

data <- read.csv('~/Downloads/Animals\ (1).csv')
names(data)
head(data)

plot(data$gestation, data$longevity)

m <- lm(data$longevity ~ data$gestation)
m
summary(m)

s <- numeric(10000)
n <- dim(data)[1]
for(i in 1:10000) {
    ind <- sample(1:n, replace = TRUE)
    dat <- data[ind,]
    m <- lm(dat$longevity ~ dat$gestation)
    s[i] <- unname(m$coef[2])
}

hist(s)
m <- lm(data$longevity ~ data$gestation)
abline(v = unname(m$coef[2]), col='red')

              



## doing the bigmac problem

data <- read.delim('~/Downloads/bigmac.txt')
names(data)
head(data)
with(data, plot(EngSal, BigMac))
rownames(data)
colnames(data)
dimnames(data)

data <- read.delim('~/Downloads/bigmac.txt')
with(data, plot(EngSal, BigMac))
with(data, text(EngSal, BigMac, labels=rownames(data), cex=0.5, pos=4))

i <- -c(17,21)

m <- lm(data$EngSal[i] ~ data$BigMac[i])
coef(m)
summary(m)

plot(data$EngSal[i], residuals(m))

dim(data)
length(data$EngSal)
length(unname(residuals(m)))
names(m)
unname(residuals(m))
data$EngSal
data$BigMac

library(MASS)
boxcox
boxcox(m)

plot(data$EngSal[i], (data$BigMac[i])^.01)
plot(data$EngSal[i], log(data$BigMac[i]))

## likelihood: y ~ Normal(mu, sigma)
## what prior for sigma?
## sigma ~ Uniform(0, 1000)

sigma <- runif(100000, 0, 1000)

hist(sigma)



## verifying the change of variable formula

a <- runif(100000, 0, 1000)

par(mfrow=c(2,1))
plot(density(a))
plot(density(sqrt(a)))
curve(2*x/1000, col='red', add=TRUE)



## STAT 365 Monte Carlo Exercise 9.1
## comparing Bayesian and Frequentist estimators of pi

n <- 10
msef <- numeric()
mseb <- numeric()
biasf <- numeric()
biasb <- numeric()
varf <- numeric()
varb <- numeric()
pis <- c(1:5/100, 1:9/10, 95:99/100)
for(i in 1:length(pis)) {
    pi <- pis[i]
    samp <- rbinom(10000, size=n, prob=pi)
    pihatf <- samp / n
    pihatb <- (samp+1)/(n+2)
    biaspif <- mean(pihatf) - pi
    biaspib <- mean(pihatb) - pi
    varpif <- var(pihatf)
    varpib <- var(pihatb)
    msef[i] <- biaspif^2 + varpif
    mseb[i] <- biaspib^2 + varpib
    biasf[i] <- biaspif
    biasb[i] <- biaspib
    varf[i] <- varpif
    varb[i] <- varpib
}
par(mfrow = c(3,1))
plot(pis, biasf, main = 'Bias', col='red', type = 'b', ylim = range(c(biasf,biasb)))
lines(pis, biasb, col='blue', type = 'b')
plot(pis, varf, main = 'Variance', col='red', type = 'b', ylim = range(c(varf,varb)))
lines(pis, varb, col='blue', type = 'b')
plot(pis, msef, main = 'MSE', col='red', type = 'b', ylim = range(c(msef,mseb)))
lines(pis, mseb, col='blue', type = 'b')



## trying HW for STAT201

data <- read.delim('~/Downloads/Saratoga.txt')
names(data)
str(data)
attach(data)
boxplot(Price)
hist(Price)
boxplot(Living.Area)
hist(Living.Area)
plot(Living.Area, Price)

data <- read.delim('~/Downloads/Mauna-Loa-and-DJIA.txt')
attach(data)
ind <- Year>=1982
plot(DJIA, CO2.Avg)
plot(DJIA[ind], CO2.Avg[ind])
cor(DJIA[ind], CO2.Avg[ind])
str(data)

pairs(data, panel = panel.smooth)
?pairs

data <- read.delim('~/Downloads/Kentucky_Derby_2014.txt')
names(data)
str(data)
attach(data)
plot(data$Year, data$Speed..mph.)


## pulling together pieces of not having to recompile models and MCMCs, for Dao

library(nimble)
nimbleOptions(showCompilerOutput = TRUE) ### DELETE THIS later
code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dnorm(0, 1)
    c ~ dnorm(0, 1)
})
Rmodel <- nimbleModel(code, inits = list(a=0, b=0, c=0))
conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

model_orig <- Cmodel

## recover the R (uncompiled) model object
if(inherits(model_orig, 'CmodelBaseClass'))
    model_orig <- model_orig$Rmodel


##md <<- Rmodel_orig$modelDef
Rmodel <- model_orig$newModel(replicate = TRUE, check = FALSE)

conf_initial <- configureMCMC(Rmodel)

monitorsVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
conf_initial$addMonitors(monitorsVector, print=FALSE)

scalarNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
discreteInd <- sapply(scalarNodeVector, function(n) Rmodel$isDiscrete(n), USE.NAMES=FALSE)
scalarNodeVectorContinuous <<- scalarNodeVector[!discreteInd]
scalarNodeVectorContinuous

firstScalarNode <- scalarNodeVectorContinuous[1]
firstScalarNode

conf_initial$printSamplers()

samplersWeMightUse <- c('RW', 'slice', 'RW_block')
for(sampler in samplersWeMightUse)
    conf_initial$addSampler(target = firstScalarNode, type = sampler)

conf_initial$printSamplers()

Rmcmc_initial <- buildMCMC(conf_initial)
Cmodel <- compileNimble(Rmodel)
Cmcmc_initial <- compileNimble(Rmcmc_initial, project = Rmodel)

conf_new <- configureMCMC(oldConf = conf_initial)

conf_new$setSamplers()  ## remove all samplers
conf_new$printSamplers()

nodes <- c('a', 'b', 'c')
for(node in nodes)
    conf_new$addSampler(target = node, type = 'slice')
conf_new$printSamplers()

Rmcmc_new <- buildMCMC(conf_new)
Cmcmc_new <- compileNimble(Rmcmc_new, project = Rmodel)


nimCopy(from = model_orig, to = Cmodel, logProb = TRUE)
calculate(Cmodel)


## doing the random drug abuse survey using NIMBLE MCMC

library(nimble)

N <- 29
y <- 14


code <- nimbleCode({
    theta ~ dunif(0, 1)
    y ~ dbinom(size = N, prob = 0.25 + 0.5*theta)
})

constants <- list(N = N)
data <- list(y = y)
inits <- list(theta = 0.5)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)

samplesPlot(samples, col='theta')
theta <- samples[,'theta']
mean(theta)
median(theta)
quantile(theta, c(0.25, 0.5, 0.75))


## testing random drug abuse survey strategy, the variance of the results

f <- function(n, iter) {
    res <- rbinom(n=iter, size=n, prob=0.25)
    var(res)
}

## expect variance is 3*n/16
iter <- 100000

n <- 100
f(n, iter) - 3*n/16


g <- function(n, iter) {
    res <- numeric(iter)
    for(i in 1:iter) {
        a <- rbinom(n=1, size=n, prob=1/2)
        res[i] <- rbinom(n=1, size=a, prob=1/2)
    }
    var(res)
}

g(n, iter) - 3*n/16



test <- function(iter, n, theta) {
    res <- numeric()
    c <- numeric()
    d <- numeric()
    for(i in 1:iter) {
        a <- rbinom(1, size=n, prob=1/2)
        b <- n - a
        c[i] <- rbinom(1, size=a, prob=1/2)
        d[i] <- rbinom(1, size=b, prob=theta)
        res[i] <- c[i] + d[i]
    }
    var(d)
}

iter <- 10000
n <- 100
theta <- .2
test(iter, n, theta)

3*n/16

theta*n*(2-theta)/4



## solving 5 statisticians problem

one <- function() {
    cur <- 1  ## 1, 2, 3, 4, 5
    while(TRUE) {
        n <- runif(1)  # stay, +1, -1
        if(n < 1/3) {
            return(cur)
        } else if(n < 2/3) {
            cur <- ((cur - 1) - 1) %% 5 + 1
        } else cur <- ((cur + 1) - 1) %% 5 + 1
    }
}

many <- function(n) {
    ret <- numeric(n)
    for(i in 1:n)
        ret[i] <- one()
    return(ret)
}

prop.table(table(many(1000000))) * 11

##        1        2        3        4        5 
## 5.002679 2.000020 0.996589 0.999779 2.000933 
##     5/11,    2/11,    1/11,    1/11,    2/11


## playihg with Rmd

setwd('~/temp/lecTEMP')
getwd()
library(methods)
library(knitr)
library(rmarkdown)

list.files()

render('Lecture1Slides.rmd')

## testing Nick Michaud's question about nimbleFunctionLists
library(nimble)

bigFunction <- nimbleFunction(
    setup = function(N) {
        functions <- nimbleFunctionList(littleFunction_virtual)  ## CHANGE
        for(n in 1:N)
            functions[[n]] <-littleFunction(n)
    },
    run = function() {
        returnType(integer(0))
        sum_N <- 0
        for(n in 1:N) 
            sum_N <- sum_N + functions[[n]]$run()
        return(sum_N)
    }
)

## NEW
littleFunction_virtual <- nimbleFunctionVirtual(
    run = function() {
        returnType(integer(0))
    }
)

littleFunction <- nimbleFunction(
    contains = littleFunction_virtual,
    setup = function(n){},
    run = function(){
        returnType(integer(0))
        return(n)
    }
)

testFunction <- bigFunction(5)
testFunction$run()
CtestFunction <- compileNimble(testFunction)
CtestFunction$run()


## make configureMCMC respect dconstraint()
library(nimble)

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[j,i] ~ dconstraint(w[j,i] > 0)
            w[j,i] ~ dnorm(theta[j], 1)
        }
        theta[j] ~ dnorm(mu, itau2)
    }
    itau2 ~ dgamma(a, b)
    mu ~ dnorm(0, .00001)
})

J <- 3
n <- rep(2, J)

y <- matrix( sample(c(0,1), sum(n), replace = TRUE), nrow = J)
m <- nimbleModel(code, constants = list(n = n, J = J), data = list(y = y))

conf <- configureMCMC(m)
conf$printSamplers()


## nested sampler function wrapper for Dao
library(nimble)

sampler_record_wrapperNEW <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control){
      numSamples <- 0
      before <- c(0, 0)
      after <- c(0, 0)
      samplerFunctionList <- nimbleFunctionList(sampler_BASE)
    ###### make sure to provide *named* arguments to this function
    ###### shouldn't require anything in control$control, if you don't want
    controlListForNestedSampler <- mcmc_generateControlListArgument(samplerFunction = control$sampler_function, control = control$control)
    samplerFunctionList[[1]] <- eval(call( control$sampler_function, model = model, mvSaved = mvSaved, target = target, control =  controlListForNestedSampler))}, 
    run = function() {
      ## these lines are new:
      numSamples <<- numSamples + 1
      setSize(before, numSamples)
      setSize(after, numSamples)
      before[numSamples] <<- model[[target]]
      ## back to the original sampler function code:
      samplerFunctionList[[1]]$run()
      ## this line new:
      after[numSamples] <<- model[[target]]
    },
    methods = list(
        reset = function() {samplerFunctionList[[1]]$reset()}
    ))

code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)

conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$removeSamplers('sigma')
conf$printSamplers()

## SEE CHANGES HERE
conf$addSampler(target = 'sigma', type = sampler_record_wrapperNEW, control = list(sampler_function = 'sampler_slice', control=list()))
conf$printSamplers()

## SEE CHANGES HERE
conf$addSampler(target = 'sigma', type = 'sampler_record_wrapperNEW', control = list(sampler_function = 'sampler_RW', control = list()))
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



## trying out a recursive nimble function
library(nimble)

Rnf <- nimbleFunction(
    run = function(x = double()) {
        if(x == 0 || x == 1) {
            return(1)
        } else {
            a <- 1
        }
    })



## testing the new NIMBLE function: runMCMC()
library(nimble)
code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)
Rmodel$setData(list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3)))
conf <- configureMCMC(Rmodel)
conf$getMonitors()
conf$setThin(10)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


## testing new functions for addSampler(), 'name', 'libraryTag', etc...
library(nimble)
code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)
conf <- configureMCMC(Rmodel)

debug(conf$addSampler)
undebug(conf$addSampler)
conf$printSamplers()
conf$removeSamplers('sigma')
conf$printSamplers()
conf$addSampler(target = 'sigma', type = sampler_slice, name='slice1')
conf$addSampler(target = 'sigma', type = 'slice', name='slice2')
conf$addSampler(target = 'sigma', type = sampler_slice)
conf$addSampler(target = 'sigma', type = 'slice')
conf$addSampler(target = 'sigma', type = sampler_RW, name='slice1')
conf$addSampler(target = 'sigma', type = 'RW', name='slice2')
conf$addSampler(target = 'sigma', type = sampler_RW)
conf$addSampler(target = 'sigma', type = 'RW')
conf$printSamplers()
conf$controlNamesLibrary



##debug(buildMCMC)
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



            

mc <- Cmcmc
mc <- Rmcmc

pb <- TRUE
pb <- FALSE

si <- TRUE
si <- FALSE

ni <- 10
ni <- 100

nb <- 5
nb <- 9



inits <- function() list(mu = rnorm(1,0,1000), sigma = runif(1,0,10))
##inits <- function() list(mu = 1:2, sigma = runif(1,0,10), x = 3)

initsList <- list(inits(), inits(), inits())
initsList <- list(inits(), inits())
initsList <- list(inits())

debug(runMCMC)
undebug(runMCMC)

debug(mcmc$run)

runMCMC(Rmcmc, niter = 3, nchains = 3, inits = inits())
runMCMC(Cmcmc, niter = 3, nchains = 3, inits = inits())

runMCMC(Rmcmc, niter = 3, nchains = 3, inits = initsList)
runMCMC(Cmcmc, niter = 3, nchains = 3, inits = initsList)

runMCMC(Rmcmc, niter = 300, nchains = 3, inits = inits())
runMCMC(Cmcmc, niter = 300, nchains = 3, inits = ii, nburnin=10)
a <- runMCMC(Rmcmc, niter = 300, nburnin=10, nchains = 4)

runMCMC(Rmcmc, niter = 300, nchains = 3, inits = inits(), nburnin=295)
runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29995)

runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE, returnCodaMCMC = TRUE)

runMCMC(Cmcmc, niter = 30000, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE, returnCodaMCMC = TRUE)

runMCMC(Cmcmc, niter = 30000, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE)

runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29999, progressBar=TRUE, silent=TRUE, setSeed=TRUE)

runMCMC(Rmcmc, niter = 30, nchains = 3, inits = inits, nburnin=29, silent=TRUE, setSeed=TRUE)

runMCMC(Cmcmc, niter = 30000, nchains = 3, nburnin=29999, progressBar=TRUE, silent=TRUE, setSeed=TRUE)

runMCMC(Rmcmc, niter = 30, nchains = 3, nburnin=29, silent=TRUE, setSeed=TRUE)



## testing inconsistancy in dgamma() between R and C

library(nimble)

Rnf <- nimbleFunction(
    run = function(x = double(), a = double(), b = double()) {
        lp <- dgamma(x, a, b, log = 1)
        returnType(double())
        return(lp)
    }
)

Cnf <- compileNimble(Rnf)

x <- 6e-100
a <- 0.001
b <- 1.0
Rnf(x, a, b)
Cnf(x, a, b)




## making utility function for MCMC sample traceplots and density histograms

library(nimble)   ## get samples from birats2 model
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(20000)
samples <- as.matrix(Cmcmc$mvSamples)

library(nimble)   ## get samples from Dave Pleydel's multinomial test model
codeTest <- nimbleCode ({
    X[1:nGroups] ~ dmultinom(size=N, prob=pVecX[1:nGroups])
    Y[1:nGroups] ~ dmultinom(size=N, prob=pVecY[1:nGroups])
    for (ii in 1:nGroups)
        Z[ii] ~ dbeta(1 + X[ii], 1 + Y[ii]) })
nGroups   <- 5
N         <- 1E6
pVecX     <- rdirch(1, rep(1, nGroups))
pVecY     <- rdirch(1, rep(1, nGroups))
X         <- rmultinom(1, N, pVecX)[,1]
Y         <- rmultinom(1, N, pVecY)[,1]
Z         <- rbeta(nGroups, 1+X, 1+Y)
Xini      <- rmultinom(1, N, sample(pVecX))[,1]
Yini      <- rmultinom(1, N, sample(pVecY))[,1]
Constants <- list(nGroups=nGroups)
Inits     <- list(X=Xini, Y=Yini, pVecX=pVecX, pVecY=pVecY, N=N)
Data      <- list(Z=Z)
modelTest <- nimbleModel(codeTest, constants=Constants, inits=Inits, data=Data)
mcmcTest  <- buildMCMC(modelTest) 
cModelTest <- compileNimble(modelTest)
cMcmcTest <- compileNimble(mcmcTest, project=modelTest)
cModelTest$N     <- N <- 1E3
cModelTest$pVecX <- sort(rdirch(1, rep(1, nGroups)))
cModelTest$pVecY <- sort(rdirch(1, rep(1, nGroups)))
simulate(cModelTest, c('X','Y','Z'), includeData=TRUE)
niter  <- 1E4
cMcmcTest$run(niter)
samples <- as.matrix(cMcmcTest$mvSamples)

samplesPlot <- function(samples, ind=1:ncol(samples), burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright') {
    ## device window and plotting parameters
    dev.new(height=height, width=width)
    par(mfrow=c(1,2), cex=0.7, cex.main=1.5, lab=c(3,3,7), mgp=c(0,0.6,0), mar=c(2,1,2,1), oma=c(0,0,0,0), tcl=-0.3, yaxt='n', bty='l')
    ## process samples
    samples <- samples[, ind, drop=FALSE]
    if(!is.null(burnin))
        samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
    nparam <- ncol(samples)
    rng <- range(samples)
    ## traceplots
    plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
    for(i in 1:nparam)
        lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
    ## posterior densities
    xMin <- xMax <- yMax <- NULL
    for(i in 1:nparam) {
        d <- density(samples[,i])
        xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
    plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='Posterior Densities', xlab='', ylab='')
    alpha_density <- 0.2
    for(i in 1:nparam)
        polygon(density(samples[,i]), col=rainbow(nparam, alpha=alpha_density)[i], border=rainbow(nparam, alpha=alpha_density)[i])
    if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
        legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
}


dim(samples)
dimnames(samples)
samplesPlot(samples)

apply(samples, 2, mean)

samplesPlot(samples, ind=c(1,2,3,7), burnin=1000, legend.location='topleft')
samplesPlot(samples, ind=c(4,6), burnin=1000)
samplesPlot(samples, ind=c(5), burnin=1000)

## better to just use the plotting functions in coda package!!!
library(coda)
mcmcSamples <- as.mcmc(samples)
acfplot(mcmcSamples)
plot(mcmcSamples)



## playing with plot(density(x))
x <- rnorm(10000)
plot(density(x))
d <- density(x)
class(d)
ls(d)
length(d$x)
length(d$y)
plot(d$x, d$y, type='l')  ## this creates the standard plot(density(x)) plot

      plot(prior ~ param.x, ylim = yLims, type = "l", xlim = range(param.x), 
            xlab = "", ylab = "", main = "", axes = FALSE, ...)
        polygon(param.x, prior, col = "red")
        box()
        r = legend("topleft", legend = "Prior", lty = 1, bty = "n", 
            plot = FALSE)$text
        text(r$x, r$y, "Prior", adj = 0)
        plot(likelihood ~ param.x, type = "l", xlab = "", ylab = "", 
            main = "", axes = FALSE, ...)
        polygon(param.x, likelihood, col = "green")
        box()
        r = legend("topleft", legend = "Prior", lty = 1, bty = "n", 
            plot = FALSE)$text
        text(r$x, r$y, "Likelihood", adj = 0)
        plot(posterior ~ param.x, ylim = yLims, type = "l", xlab = "", 
            ylab = "", main = "", axes = F, ...)
        polygon(param.x, posterior, col = "blue")




## testing funny behavior of DSL round() and nimRound()
library(nimble)

Rnf <- nimbleFunction(
    run = function() {
        for(i in 0:50) {
            x <- i/10
            xRound <- round(x)
            print('x: ', x, ',   round(x): ', xRound)
        }
    }
)

Cnf <- compileNimble(Rnf)

Rnf()
Cnf()



## testing making my own progress bar in NIMBLE DSL nimbleFunction

library(nimble)

rfun <- nimbleFunction(
    run = function(pb = logical(default=TRUE)) {
        ##print('|')
        ##for(i in 1:20)   { a <- i %% 7; print(a) }
        ##print('|')
        ##cat('|')
        ##for(i in 1:3)   cat('-', i)
        ##cat('|')
        a <- 1
        pb <- pb & 0
        print(pb)
    }
)

##rfun()
cfun <- compileNimble(rfun)
##cfun()
cfun()


library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


library(nimble)
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(20000)
Cmcmc$run(30007, reset = FALSE)         ## continue previous run
Cmcmc$run(10000, progressBar = FALSE)   ## turn off progress bar
Cmcmc$run(100003, reset = FALSE)        ## sort of slow run
Cmcmc$run(5000)     ## faster
Cmcmc$run(2000)     ## faster still
Cmcmc$run(1000)     ## ...
Cmcmc$run(500)
Cmcmc$run(100)
Cmcmc$run(40)  ## no bar when too few iterations



## playing with R progress bars: txtProgressBar

f <- function(n = 1e6, frac = 0.01) {
    pb <- txtProgressBar(style = 3, char = '-')
    nupdate <- floor(frac * n) 
    for(i in 1:n) {
        a <- rnorm(1)
        if(i %% nupdate == 0) {
            setTxtProgressBar(pb, i/n)
        }
    }
    setTxtProgressBar(pb, 1)
    close(pb)
}

f(n=2e6, frac = .01)



library(Bolstad)

## datasets in library(Bolstad):
## bears (exercise in chapter 3)
## slug  (exercise in chapter 14)

## functions that I might possibly use:
## decomp: makes plots of prior, likelihood, posterior, but only for class=Bolstad.
##         maybe steal a bunch of the code, make one that works for general samples?
##         BUT WAIT -- this only works for sorted (x,y) pairs -- not for samples.


## getting MCMC samples for logProbs of variables

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$getMonitors()
conf$addMonitors('logProb_a')
conf$getMonitors()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##paste0('logProb_', letters)

set.seed(0)
Cmcmc$run(10)
samples <- as.matrix(Cmcmc$mvSamples)
samples
apply(samples, 2, mean)


## tracking down error in conjugacy checking system,
## this model was submitted by: Eduardo Martins <egmartins@gmail.com>
## during the ISEC 2016 NIMBLE workshop
library(nimble)
load('~/Downloads/data_inits.RData')
fec_code <- nimbleCode({
    ## Likelihood
    ## Fecundity
    for (i in 1:ny.fec) {
        y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec) # observation model
        z.fec[i] ~ dnegbin(p.fec[i], size.fec) # sampling model
        p.fec[i] <- size.fec / (size.fec + mu.fec[i])
        log(mu.fec[i]) <- a.fec.raw + b.fec * z.sl[i] + nu.fec.k.raw[idx.k[i]] + nu.fec.l.raw[idx.l[i]]
    }
    a.fec <- a.fec.raw + mean(nu.fec.k.raw[1:nk]) + mean(nu.fec.l.raw[1:nl])
    for (k in 1:nk){
        nu.fec.k.raw[k] ~ dnorm(mu.nu.fec.k, tau.nu.fec.k)
        nu.fec.k[k] <- nu.fec.k.raw[k] - mean(nu.fec.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.fec.l.raw[l] ~ dnorm(mu.nu.fec.l, tau.nu.fec.l)
        nu.fec.l[l] <- nu.fec.l.raw[l] - mean(nu.fec.l.raw[1:nl])
    }
    ## Length
    for(i in 1:ny.fec) {
        z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
        mu.sl[i] <- a.sl.raw + b.sl * z.ocr[i] + nu.sl.k.raw[idx.k[i]] + nu.sl.l.raw[idx.l[i]]
    }
    a.sl <- a.sl.raw + mean(nu.sl.k.raw[1:nk]) + mean(nu.sl.l.raw[1:nl])
    for (k in 1:nk){
        nu.sl.k.raw[k] ~ dnorm(mu.nu.sl.k, tau.nu.sl.k)
        nu.sl.k[k] <- nu.sl.k.raw[k] - mean(nu.sl.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.sl.l.raw[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
        mu.nu.sl.l[l] <- b.nu.sl.l[1] * sst6.std[l] + b.nu.sl.l[2] * pdo6.std[l] + b.nu.sl.l[3] * nsal.std[l]
        nu.sl.l[l] <- nu.sl.l.raw[l] - mean(nu.sl.l.raw[1:nl])
    }
    ## Years of ocean residence
    for(i in 1:ny.fec){
        z.ocr[i] ~ dbern(p.z.ocr[i])
        logit(p.z.ocr[i]) <- a.ocr.raw + nu.ocr.k.raw[idx.k[i]] + nu.ocr.l.raw[idx.l[i]]
    }
    a.ocr <- a.ocr.raw + mean(nu.ocr.k.raw[1:nk]) + mean(nu.ocr.l.raw[1:nl])
    for (k in 1:nk){
        nu.ocr.k.raw[k] ~ dnorm(mu.nu.ocr.k, tau.nu.ocr.k)
        nu.ocr.k[k] <- nu.ocr.k.raw[k] - mean(nu.ocr.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.ocr.l.raw[l] ~ dnorm(mu.nu.ocr.l, tau.nu.ocr.l)
        nu.ocr.l[l] <- nu.ocr.l.raw[l] - mean(nu.ocr.l.raw[1:nl])
    }
    ## Priors
    ## Fecundity
    tau.y.fec <- pow(sig.y.fec, -2)
    sig.y.fec ~ dunif(0, 100)
    num.size.fec ~ dnorm(0, 0.0016)
    den.size.fec ~ dnorm(0, 1)
    size.fec <- abs(num.size.fec / den.size.fec)
    mu.nu.fec.k ~ dnorm(0, 0.0001)
    tau.nu.fec.k <- pow(sig.nu.fec.k, -2)
    sig.nu.fec.k ~ dunif(0, 100)
    mu.nu.fec.l ~ dnorm(0, 0.0001)
    tau.nu.fec.l <- pow(sig.nu.fec.l, -2)
    sig.nu.fec.l ~ dunif(0, 100)
    a.fec.raw ~ dnorm(mean.z.fec, tau.z.fec)
    b.fec ~ dnorm(0, 0.0001)
    ## Length
    tau.sl.i <- pow(sig.sl.i, -2)
    sig.sl.i ~ dunif(0, 100)
    mu.nu.sl.k ~ dnorm(0, 0.0001)
    tau.nu.sl.k <- pow(sig.nu.sl.k, -2)
    sig.nu.sl.k ~ dunif(0, 100)
    tau.nu.sl.l <- pow(sig.nu.sl.l, -2)
    sig.nu.sl.l ~ dunif(0, 100)
    a.sl.raw ~ dnorm(0, 0.0001)
    b.sl ~ dnorm(0, 0.0001)
    for (b in 1:nb.nu.sl.l) {
        b.nu.sl.l[b] ~ dnorm(0, 0.0001)
    }
    ## Years of ocean residence
    mu.nu.ocr.k ~ dnorm(0, 0.0001)     
    tau.nu.ocr.k <- pow(sig.nu.ocr.k, -2)
    sig.nu.ocr.k ~ dunif(0, 100)
    mu.nu.ocr.l ~ dnorm(0, 0.0001)
    tau.nu.ocr.l <- pow(sig.nu.ocr.l, -2)
    sig.nu.ocr.l ~ dunif(0, 100)
    a.ocr.raw ~ dnorm(0, 0.0001)
    ## Residuals and metrics for model assessment
    ## Fecundity
    for(i in 1:ny.fec){
        res.y.fec[i] <- y.fec[i] - z.fec[i]
        new.y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec)
    }
    for(i in 1:ny.fec){
        pres.z.fec[i] <- (z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        new.z.fec[i] ~ dnegbin(p.fec[i], size.fec)
        pres.new.z.fec[i] <- (new.z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        d.z.fec[i] <- pow(pres.z.fec[i], 2)
        d.new.z.fec[i] <- pow(pres.new.z.fec[i], 2)
    }
    fit.z.fec <- sum(d.z.fec[1:ny.fec])
    fit.new.z.fec <- sum(d.new.z.fec[1:ny.fec])
    ## Length
    for (i in 1:ny.fec) {
        res.z.sl[i] <- z.sl[i] - mu.sl[i]
        new.z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
    }
    for (l in 1:nl){
        res.nu.sl.l[l] <- nu.sl.l[l] - mu.nu.sl.l[l]
        new.nu.sl.l[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
    }
})
fec_mod <- nimbleModel(fec_code, constants = dat, inits = inits())
fec_cmod <- compileNimble(fec_mod)

##options(error = recover)

## next line errors out in conjugacy check:
## Error in x[[1]] : subscript out of bounds
mcmcConf <- configureMCMC(fec_mod)

mcmcConf$printSamplers()
mcmcConf$getMonitors()
fecMCMC <- buildMCMC(mcmcConf)
cfecMCMC <- compileNimble(fecMCMC, project = fec_mod)

cfecMCMC$run(10)

## testing of new setData() functionality:
## taking a character vector of variable names

library(nimble)
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')

Rmodel$getNodeNames(includeData = FALSE)

modelLs <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R', returnModelComponents = TRUE)
modelLs$model

Rmodel$beta
Rmodel$tau.c
Rmodel$r

Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('Y'))





## testing new MCMC runtime argument: time
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    c ~ dnorm(b + b^2, 1)
})
constants <- list()
data <- list(b=0, c=0)
inits <- list(a=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000, time = TRUE)
Cmcmc$getTimes()
Cmcmc$run(100000, time = TRUE, reset = FALSE)
Cmcmc$getTimes()



## error from Perry about not finding nimArray
## but I can't reproduce this error

library(nimble)

Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
##spec <- configureMCMC(Rmodel)
##spec$printSamplers()
##spec$getSamplerDefinition(1)
##spec$getSamplerDefinition(4)
##mcmc <- buildMCMC(spec)
mcmc <- buildMCMC(Rmodel)
mcmc$run(1)
##Error in samplerFunctions[[i]]$run() : could not find function "nimArray"

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R', data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)




## Soledad's compilation error

library(nimble)
x <- c(2, 2, 3)

## works
test1 <- nimbleFunction ( run = function (indx = double(1) ) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        aa <- indx[1] 
    } else {
        newindx <- indx
        aa         <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

## doesn't work
test1 <- nimbleFunction(run = function(indx = double(1)) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        newindx <- indx[1]
        aa          <- newindx     ## THIS IS THE MAIN CHANGE IN THIS FUNCTION
    } else {
        newindx <- indx
        aa          <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

Ctest1 <- compileNimble(test1)


test1(x)
Ctest1(x)





## using a custom distribution to have flexible length data
## for floriane plard population models

library(nimble)


dxxx <- nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1), length = integer(), log.p = double()) {
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dnorm(x[i], mu[i], sd=sigma[i], log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)

rxxx <- nimbleFunction(
    run = function(n = integer(), mu = double(1), sigma = double(1), length = integer()) {
        print('this should never run')
        x <- numeric(length)
        return(x)
    }
)

registerDistributions(list(
    dxxx = list(
        BUGSdist = 'dxxx(mu, sigma, length)',
        types    = c('value = double(1)', 'mu = double(1)', 'sigma = double(1)', 'length = integer()')
    )
))

code <- nimbleCode({
    y[1:N] ~ dxxx(mu[1:N], sigma[1:N], length)
})
constants <- list(N = 10)
data <- list()
inits <- list(length = 10)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$length
Rmodel$mu
Rmodel$sigma
Rmodel$y
Rmodel$getNodeNames(dataOnly = TRUE)

length <- 10
length
mu <- ((1:10)/2)[1:length]
mu
sigma <- (1+(1:10)/10)[1:length]
sigma
y <- c(1,0,1,3,2,5,4,6,7,5)[1:length]
y
yaug <- c(y, rep(as.numeric(NA), 10-length))
yaug

Rmodel$length <- length
Rmodel$length
Rmodel$mu[1:length] <- mu
Rmodel$mu
Rmodel$sigma[1:length] <- sigma
Rmodel$sigma
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(list(y=yaug))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$y

sum(dnorm(y, mu, sd=sigma, log=TRUE))
calculate(Rmodel)


debug(exprClasses_setSizes)



## testing MCEM for pump model
library(nimble)

pumpCode <- nimbleCode({ 
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta) 
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i]) 
  }
  alpha ~ dexp(1.0) 
  beta ~ dgamma(0.1,1.0) 
}) 
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                       31.4, 1.05, 1.05, 2.1, 10.5)) 
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)) 
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N)) 
newPump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits) 

pumpMCEM <- buildMCEM(model = newPump,
                      latentNodes = 'theta',
                      burnIn = 100,
                      mcmcControl = list(adaptInterval = 20),
                      boxConstraints = list( list( c('alpha', 'beta'), 
                          limits = c(0, Inf) ) ), 
                      buffer = 1e-6)

pumpMCEM(maxit = 20, m1 = 250, m2 = 500)
pumpMCEM(maxit = 50, m1 = 1000, m2 = 5000)


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
Rfun <- nimbleFunction(run = function() {
    ##ans <- numeric(10, value = 2)
    ##
    #x <- 100 
    #ans <- integer(x) 
    #for(i in 1:x) {
    #    ans[i] <- i 
    #}
    ##
    ##ans <- matrix(1, nrow = 10, ncol = 1)
    ##y <- 3
    ##ans <- numeric(10, value = y) 
    ##ans <- array(y, dim = 10) 
    ##ans <- array(y, dim = c(10))
    ##z <- numeric(10)
    ##z[5] <- 3
    ##x <- 20
    ##y <- 30
    ##ans <- integer(z[5], value = x + y) 
    ##ans <- array(x+y, dim = z[5], type = 'integer')
    ##ans <- array(x+y, dim = c(z[5]), type = 'integer')
    ##
    ##x <- array(0, c(4,5))
    ##ans <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2]) 
    ##ans <- array(0, dim = c(dim(x)[1], dim(x)[2]))
    ##
    x <- 1
    y <- 2
    z <- 3
    ans <- array(0, dim = c(x, y, z))
    ##returnType(double(1))
    ##returnType(integer(1))
    ##returnType(double(2))
    returnType(double(3))
    return(ans)
})
Cfun <- compileNimble(Rfun)

Rfun()
Cfun()
class(Rfun()[1])
class(Cfun()[1])
class(Rfun()[1,1])
class(Cfun()[1,1])
class(Rfun()[1,1,1])
class(Cfun()[1,1,1])

## creates a length-100 integer vector, containing 1, 2, ..., 100 

## creates a 10x1 ones-matrix 

## the following three lines are equivalent 
## each creates a length-10 vector, with all elements equal to y 

## the following two lines are equivalent 
## each creates an integer vector of length z[5], with all elements equal to x+y 

## the following two lines are equivalent 
## each one creates a matrix of 0's of the same size as matrix x 

## the following creates a 3-dimensional array of 0's 



## testing adding numeric(), integer(), array(), matrix()
library(nimble)
library(testthat)
##source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

expected <- numeric(10)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(10)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))


expected <- rep(3, length = 2)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(value = 3, length = 2)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- rep(9, 3)
Rfun <- nimbleFunction(run = function() {
    x <- numeric(10, value = 3)
    ans <- integer(x[2], x[2]+x[3]*2)
    returnType(integer(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1]), 'integer'))

expected <- array(4, c(10,11))
Rfun <- nimbleFunction(run = function() {
    ans <- array(4, c(10,11))
    returnType(double(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('integer', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- matrix(as.integer(0), nrow=4, ncol=5)
Rfun <- nimbleFunction(run = function() {
    x <- 4
    y <- 5
    ans <- matrix(init=FALSE, nrow=x, ncol=y, type='integer')
    returnType(integer(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1,1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1,1]), 'integer'))




expected
Rfun()
Cfun()


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        x <- 2
        bbb <- array(3, x+10, type='integer')
        ##print(bbb)
        bbb2 <- array(value=3, dim=c(3,4), type='double')
        ##print(bbb2)
        ccc <- array(value=3, dim=c(3,4), type='double', init=FALSE)
        ##print(ccc)
        x <- 5
        y <- matrix(2, ncol = 3, nrow = x)
        ##print(y)
        zxc <- array(2, c(x, dim(y)[2], y[2,2]))
        zxc[1,1,1] <- 98
        ##returnType(double(3));  return(zxc)
        zz <- numeric(4)
        zz[1] <- 99
        ##print(zz)
        z <- numeric(value=.5, 4)
        ##print(z)
        zz[2:4] <- z[1:3] + numeric(3, 10)
        ##print(zz)
        mat <- matrix(0, 5, 5)
        mat[1, 1] <- 97
        ##print(mat)
        arr <- array(4, type = 'integer', init = FALSE, dim=c(3,4))
        ##print(arr)
        ar2 <- array(dim = dim(z)[1])
        ##print(ar2)
        ar3 <- array(dim = c(dim(z)[1]))
        ##print(ar3)
        ar4 <- array(dim = c(dim(arr)[1], dim(arr)[2]), type = 'integer')
        ##print(ar4)
    }
)
Rnf <- nfDef()
##Rnf$run
Cnf <- compileNimble(Rnf, dirName)

Rnf$run()

Cnf$run()



tempdir()

##Cnf$run()



## ???

library(nimble)
sampler_dt <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        browser()
        1
        2
        3
        ###  node list generation  ###
        calcNodes  <- model$getDependencies(target)
    },
    run = function() {
        simulate(model, target)
        calculate(model, calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)
spec$getMonitors()
spec$addMonitors(c('theta', 'alpha2'))
spec$getMonitors()
spec$printSamplers()
spec$addSampler(type = 'dt', target = 'alpha')
spec$printSamplers()

Rmcmc <- buildMCMC(spec)

## testing new function in modelBaseClass: model$getNodeFunctions()
library(nimble)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)
Cmodel <- compileNimble(Rmodel)

nodes <- c('theta[3]')
nodes <- c('theta[3]', 'x[10]')
nodes <- c('theta[3]', 'x[9:10]')
Rmodel$getNodeFunctions(nodes)$calculate
Cmodel$getNodeFunctions(nodes)[[1]]$calculate

Rmodel$getNodeNames()
modelDef <- Rmodel$getModelDef()
gids <- mdef$nodeName2GraphIDs(nodes)
gids
mdef$graphIDs2indexedNodeInfo(gids)
Rmodel$nodeFunctions[[2]]$calculate


## making dmnorm conjugate sampler work with dependent
## nodes of different size from target node
library(nimble)
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

code <- nimbleCode({
    x[1:3] ~ dmnorm(mu0[1:3], prec = ident[1:3,1:3])
    mu_y2[1:2] <- asCol(a[1:2]) + B[1:2,1:3] %*% asCol(x[1:3])
    mu_y3[1:3] <- asCol(a[1:3]) + B[1:3,1:3] %*% asCol(x[1:3])
    mu_y5[1:5] <- asCol(a[1:5]) + B[1:5,1:3] %*% asCol(x[1:3])
    y2[1:2] ~ dmnorm(mu_y2[1:2], prec = prec_y[1:2,1:2])
    y3[1:3] ~ dmnorm(mu_y3[1:3], prec = prec_y[1:3,1:3])
    y5[1:5] ~ dmnorm(mu_y5[1:5], prec = prec_y[1:5,1:5])
})

mu0 <- rep(0,3)
ident <- diag(3)
a <- 11:15
B <- matrix(1:15, nrow=5, ncol=3, byrow=TRUE)
prec_y <- diag(1:5)

constants <- list(mu0=mu0, ident=ident, a=a, B=B, prec_y=prec_y)
data <- list(y2=1:2, y3=1:3, y5=1:5)
inits <- list(x=rep(0,3))

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
##spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##         x[1]       x[2]      x[3]
##[1,] 4.966869 -1.1252473 -3.922769
##[2,] 4.972376 -0.5288582 -4.511282


##node <- quote(mu0[1:3])
##node <- quote(prec_y[1:2, 1:2])
##node <- quote(prec_y[1:3, 1:3])
##node <- quote(prec_y[1:5, 1:5])
##nimble:::cc_expandDetermNodesInExpr(Rmodel, node)
##debug(nimble:::cc_expandDetermNodesInExpr)
##undebug(nimble:::cc_expandDetermNodesInExpr)
##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
##a <- spec$getSamplerDefinition(1)
##createNamedObjectsFromList(a, writeToFile = '~/temp/del.R')







## ???
library(nimble)

d <- read.csv("http://personal.bgsu.edu/~albert/data/gateway.csv")
library(dplyr)
d$month <- as.numeric(d$Month)
d$year <- d$Year - 2002

worshipCode <- nimbleCode({ 
  for (i in 1:N){
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- mu + epsilon[i] + alpha[month[i]]
  epsilon[i] ~ dnorm(0, sd=sigma.y)
  }
  for (j in 1:J){
  alpha[j] ~ dnorm(0, sd=sigma.month)
}
mu ~ dnorm(0, sd=1000)
sigma.y ~ dunif(0, 100)
sigma.month ~ dunif(0, 100)
})

worshipConsts <- list(N=484,
                      J=12,
                      K=10,
                      month=d$month)
worshipData <- list(y=d$Count)
worshipInits <- list(mu=5, sigma.y=.4, sigma.month=.4,
                     alpha=rep(0, 12))

worship <- nimbleModel(code = worshipCode, name = 'worship', 
                    constants = worshipConsts,
                    data = worshipData,
                    inits=worshipInits)

Rmodel <- worship

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## error produced in model$check() when using equals()
library(nimble)
code <- nimbleCode({
    for (i in 1:5) {
        x[i] ~ dbern(0.5)
        ##y[i] <- equals(x[i], 0)
    }
    sss <- sum(x[1:5])
    none <- equals(sss, 0)
})
constants <- list()
data <- list()
inits <- list(x = c(0,0,1,1,1))

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$x
Rmodel$sss
Rmodel$none

## testing new log=TRUE option of RW sampler

library(nimble)
code <- nimbleCode({
    sigma ~ dunif(0,10)
    a ~ dnorm(0, sd=sigma)
    b ~ dnorm(sigma, tau=3)
})
constants <- list()
data <- list(a=1, b=4)
inits <- list(sigma=1)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes=NULL)

spec$addSampler('sigma', 'RW_log')
### --- or ---
spec$addSampler('sigma', 'RW', control=list(log=TRUE))

Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)

##         sigma
## [1,] 3.535852
## [2,] 3.535852
## [3,] 3.535852
## [4,] 5.352671
## [5,] 5.352671
## [6,] 3.986347
## [7,] 3.963423
## [8,] 3.963423
## [9,] 3.963423
##[10,] 3.963423




## Github issue #107
## initializeModel issues false warning about RHSonly variable not initialized
## reduced case
library(nimble)

code <- nimbleCode({
    for (n in 1:N){
        a0[n,1:3] ~ dmnorm(mu[1:3], Sigma[1:3,1:3])
        ##lambda[n,1:3] <- exp(a0[n,1:3]) ## alternative, which shifts the problem to lambda
        for (k in 1:K){
            lambda[n,k] <- exp(a0[n,k])
            y[n,k]~dpois(lambda[n,k])
        }
    }
})

constants <- list(N = 3, K = 3)
data <- list(y = matrix(rpois(9, 3), nrow=3))
inits <- list(mu = rep(0,3), Sigma = diag(3))

Rmodel <- nimbleModel(code, constants, data, inits)

Cmodel <- compileNimble(Rmodel)

Rmcmc <- buildMCMC(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$run(1)
Cmcmc$run(1)



## testing new binary sampler
library(nimble)

code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dbern(0.6)
    c ~ dbern(0.05)
    d ~ dbin(prob=0.2, size=1)
    e ~ dbinom(prob=0.9, size=1)
    f ~ dbern(0.5)
    g ~ dbern(0.5)
    h ~ dbern(0.5)
    for(i in 1:10)
        yf[i] ~ dnorm(f, sd = 1)
    for(i in 1:10)
        yg[i] ~ dnorm(g, sd = 1)
    for(i in 1:10)
        yh[i] ~ dnorm(h, sd = 1)
    ##x ~ dnorm(0,1)
    ##y ~ dnorm(x*x, 1)
    ##z ~ dnorm(x*y, 1)
    ##zz ~ dnorm(z + z, 1)
})
constants <- list()
data <- list(yf = c(rep(0,2), rep(1,8)), yg = c(rep(0,8), rep(1,2)), yh = c(rep(0,5), rep(1,5)))
inits <- list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
Rmodel <- nimbleModel(code, constants, data, inits)

##spec <- configureMCMC(Rmodel, autoBlock = TRUE)
##spec$printSamplers()

Rmodel$isBinary('a')
Rmodel$isBinary('b')
Rmodel$isBinary('c')
Rmodel$isBinary('d')
Rmodel$isBinary('e')
Rmodel$isBinary('f')
Rmodel$isBinary('g')
Rmodel$isBinary('h')

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'binary', print=FALSE)
spec$addSampler('b', 'binary', print=FALSE)
spec$addSampler('c', 'binary', print=FALSE)
spec$addSampler('d', 'binary', print=FALSE)
spec$addSampler('e', 'binary', print=FALSE)
spec$addSampler('f', 'binary', print=FALSE)
spec$addSampler('g', 'binary', print=FALSE)
spec$addSampler('h', 'binary', print=FALSE)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)
##Rmcmc$run(3)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)
means


# Slice sampler name conflict with library(bbmle) 'slice' function
library(bbmle)
library(nimble)  # No messages about slice are generated

pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})

pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(2000000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## testing passing a compiled model as a setup argument to a NF

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        simulate(model, node)
    }
)

Rnf <- nfDef(Rmodel, 'a')
Rnf <- nfDef(Cmodel, 'a')

Cnf <- compileNimble(Rnf, project = Rmodel)

set.seed(0)
Rmodel$a
Rnf$run()
Rmodel$a

set.seed(0)
Cmodel$a
Cnf$run()
Cmodel$a



## testing passing a compiled model object to configureMCMC(), or buildMCMC()

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

###
Rmcmc <- buildMCMC(Cmodel)

###
spec <- configureMCMC(Cmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)


## testing problem with covariance matrix in BUGS model causing a chol() error at time of model checking

library(nimble)

code <- nimbleCode({
    C[1,1] <- 1
    C[1,2] <- 0
    C[2,1] <- 0
    C[2,2] <- 1
    mu[1] <- 0
    mu[2] <- 0
    y[1:2] ~ dmnorm(mu[1:2], prec = C[1:2, 1:2])
})
constants <- list()
data <- list(y = c(0,0))
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)





##If you're interested, here's the ecological/scientific background to the project: the potato psyllid spreads a disease called Zebra Chip Disease that is damaging potato and tomato crops in California and western US and Mexico. The disease is new, having only been described about 20 years ago, but the potato psyllid is native to the area. So what's changed to cause the disease outbreaks? I'm using data from where and when potato psyllid museum specimens were collected in California to test if psyllid populations have increased over the last century and if climate change may play a role. If so, this could explain the disease outbreaks. It's complicated though because museum data is messy, specifically, the museum specimens give me presence-only information. So I'm using lists of related species collected at the same time as the psyllids to infer psyllid absences, based on an assumed correlation between number of species collected and collecting effort. This is the "list_length" covariate that's in the model.

##BLOCK SAMPLE PARAMS THAT GO IN THE SAME LINEAR PREDICTOR,
##i.e. all the betapN terms in this one

{
    ## Priors
    ## For the site random effect
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu.alpha, tau.alpha) 
    }
    mu.alpha ~ dnorm(0, 0.001)
    tau.alpha <- 1 / (sigma.alpha * sigma.alpha)
    sigma.alpha ~ dunif(0, 5)
    
    ## Grand mean
    muq ~ dnorm(0, 0.001)
    
    ## For the fixed effect coefficients
    betaq ~ dnorm(0, 0.001)
    betap1 ~ dnorm(0, 0.001)
    betap2 ~ dnorm(0, 0.001)
    betap3 ~ dnorm(0, 0.001)
    betap4 ~ dnorm(0, 0.001)
    betap5 ~ dnorm(0, 0.001)
    betap6 ~ dnorm(0, 0.001)
    betap7 ~ dnorm(0, 0.001)
    betap8 ~ dnorm(0, 0.001)

    ## Likelihood
    for (i in 1:nlist){ # i = events (year-months)
        for(j in 1:nsite) { # j = sites

            logit(p[i,j]) <- betap1*year[i,j] + betap2*pow(year[i,j],2) + # Year quadratic effects
                betap2*month[i,j] + betap3*pow(month[i,j],2) + # month quadratic effects
                    betap4*aet[i,j] + betap5*cwd[i,j] + betap6*tmn[i,j] + betap7*tmx[i,j] + # Climate effects
                        alpha[j] # Random effects

            Y[i,j] ~ dbern(p[i,j]) # Occupancy probability

            logit(q[i,j]) <- Y[i,j] + muq + betaq*list_length[i,j] # Logistic regression for detection

            detectionMatrix[i,j] ~ dbern(q[i,j])  # Distribution for random part; observed presences relating to detection probability

        } #j
    } #i
}



## bug with forwardsolve in nimbleFunctions

library(nimble)
A <- matrix(c(2,1,0,3), nrow = 2)
b <- c(1,5)

nf <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        declare(x, double(1))
        setSize(x, 2)
        x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
        ##x <- forwardsolve(A[1:2,1:2], b[1:2])
        x[1:2,1:2] <- chol(A[1:2,1:2])
        returnType(double(1))
        ##declare(x, double(2))
        ##setSize(x, 2, 2)
        ##returnType(double(2))
        return(x)
    })

cnf <- compileNimble(nf, dirName = '.')

nf(A,b)
cnf(A,b)



## bug with forwardsolve in BUGS models

library(nimble)

code <- nimbleCode({
    A[1,1] <- 2
    A[2,1] <- 1
    A[1,2] <- 0
    A[2,2] <- 3
    b[1] <- 1
    b[2] <- 5
    x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
})

Rmodel <- nimbleModel(code)
Cmodel <- compileNimble(Rmodel)

Cmodel$A
Cmodel$b
Cmodel$x
calculate(Cmodel)
## Process R floating point exception: 8 at Thu Apr  7 14:38:02 2016



## testing using forwardsolve and backsolve in a user-defined function
## in a BUGS model

library(nimble)
A <- array(NA, c(2,2))
A[1,1] <- 2
A[2,1] <- 1
A[1,2] <- 0
A[2,2] <- 3
b <- c(1,5)

myFS <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        returnType(double(1))
        ret <- forwardsolve(A, b)
        return(ret)
    }
)

myChol <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- chol(A)
        return(ret)
    }
)

myInverse <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- inverse(A)
        return(ret)
    }
)

code <- nimbleCode({
    x1[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])  ## forwardsolve() directly here does NOT work
    x2[1:2] <- myFS(A = A[1:2, 1:2], b = b[1:2])
    C1[1:2,1:2] <- chol(A[1:2, 1:2])
    C2[1:2,1:2] <- myChol(A[1:2, 1:2])
    I1[1:2,1:2] <- inverse(A[1:2, 1:2])
    I2[1:2,1:2] <- myInverse(A[1:2, 1:2])
    xx1[1:2,1:2] <- chol(A[1:2, 1:2]) + 2 * inverse(A[1:2, 1:2])
    xx2[1:2,1:2] <- myInverse(A[1:2, 1:2]) * 2 + myChol(A[1:2, 1:2])
})

Rmodel <- nimbleModel(code, constants = list(A=A, b=b))
Cmodel <- compileNimble(Rmodel)

m <- Rmodel
m <- Cmodel

calculate(m)
m$x1; m$x2
m$C1; m$C2
m$I1; m$I2
m$xx1; m$xx2

## another example of MVN conjugate sampler, for test-mcmc.R
## using both cov and prec parametrizaions of MVN,
## and various linear links

library(nimble)

set.seed(0)
mu0 <- rep(0,5)
ident <- diag(5)
a <- array(rnorm(20), c(4,5))
B <- array(NA, c(4, 5, 5))
for(i in c(1,3)) {
    A <- array(rnorm(25,i), c(5,5))
    A <- A + t(A) + 5*i*diag(5)
    B[i,,] <- A
}
B[2,,] <- diag(5)
B[4,,] <- diag(5)
M_y <- array(NA, c(4, 5, 5))
for(i in 1:4)
    M_y[i,,] <- i * diag(5)
x <- array(0, c(4,5))
y <- array(0, c(4,5))

code <- nimbleCode({
    for(i in 1:2)
        x[i,1:5] ~ dmnorm(mu0[1:5], prec = ident[1:5,1:5])
    for(i in 3:4)
        x[i,1:5] ~ dmnorm(mu0[1:5], cov  = ident[1:5,1:5])
    for(i in 1:4)
        mu_y[i,1:5] <- asCol(a[i,1:5]) + B[i,1:5,1:5] %*% asCol(x[i,1:5])
    y[1,1:5] ~ dmnorm(mu_y[1,1:5], prec = M_y[1,1:5,1:5])
    y[2,1:5] ~ dmnorm(mu_y[2,1:5], cov  = M_y[2,1:5,1:5])
    y[3,1:5] ~ dmnorm(mu_y[3,1:5], prec = M_y[3,1:5,1:5])
    y[4,1:5] ~ dmnorm(mu_y[4,1:5], cov  = M_y[4,1:5,1:5])
})
constants <- list(mu0=mu0, ident=ident, a=a, B=B, M_y=M_y)
data <- list(y=y)
inits <- list(x=x)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##        x[1, 1]     x[2, 1]     x[3, 1]    x[4, 1]     x[1, 2]   x[2, 2]
##[1,] -0.7098220 -0.92579244  0.02859492 -0.9271968 -0.03440439 0.1133790
##[2,] -0.6642612  0.09915454 -0.13079585  0.6298495 -0.08432850 0.8806305
##        x[3, 2]    x[4, 2]     x[1, 3]   x[2, 3]     x[3, 3]    x[4, 3]
##[1,] 0.09213596  0.3257091  0.27465408 -1.410726 -0.12937422 -0.2311789
##[2,] 0.01777368 -0.6051182 -0.09322287 -1.723685  0.07759651  0.4119755
##       x[1, 4]     x[2, 4]     x[3, 4]     x[4, 4]   x[1, 5]  x[2, 5]
##[1,] 0.5071505 0.140652780  0.01778143  0.51383003 0.3639871 2.146307
##[2,] 0.4781541 0.002361055 -0.02517465 -0.06995162 0.4202336 1.225317
##         x[3, 5]    x[4, 5]
##[1,] -0.02913630 -0.2563026
##[2,]  0.04303507  0.2716980




## Richard McElreath Statistical Rethinking package,
## might have lots of good data and hierarchical (Bayesian) modeling examples
install.packages('rethinking')
library(rethinking)



## tesing MVN conjugate sampler

library(nimble)

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

Rmodel <- nimbleModel(code, data= data)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
##spec$getSamplerDefinition(1)
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples
##         mu[1]    mu[2]    mu[3]
## [1,] 3.877574 2.201065 1.084173
## [2,] 2.537805 1.916830 2.084021
## [3,] 4.426193 2.326982 1.567054
## [4,] 2.456652 1.921560 1.892777
## [5,] 2.959890 2.133226 1.506502
## [6,] 3.084360 1.560040 1.941621
## [7,] 3.178667 2.072212 2.611881
## [8,] 2.996492 2.210513 2.795391
## [9,] 2.494759 1.539610 2.118232
##[10,] 2.880285 1.826356 1.805385

Rsamples - Csamples
##              mu[1] mu[2] mu[3]
## [1,] -8.881784e-16     0     0
## [2,]  0.000000e+00     0     0
## [3,] -8.881784e-16     0     0
## [4,] -4.440892e-16     0     0
## [5,] -4.440892e-16     0     0
## [6,]  0.000000e+00     0     0
## [7,] -4.440892e-16     0     0
## [8,] -4.440892e-16     0     0
## [9,]  0.000000e+00     0     0
##[10,] -4.440892e-16     0     0


## inconsistent behavior between R and C chol()

n <- 4
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
mat <- A + L

sym <- function(mat, lower) {
    ind <- if(lower) lower.tri(mat) else upper.tri(mat)
    ans <- mat
    ans[!ind] <- 0
    ans <- ans + t(ans) + diag(diag(mat))
    ans
}

## R execution
chol(mat)
chol(sym(mat, FALSE))
chol(mat) - chol(sym(mat, FALSE)) ## R takes the UPPER part of the matrix

## Rnf execution
Rnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Rnf also UPPER part of matrix

## Cnf execution
Cnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Cnf now also uses UPPER part of matrix





## test case for forwardsolve, backsolve in the DSL
## A.triangularView<Eigen::Upper>().solve(b);
## could make a file called nimble_macros.h:
##define ......

library(nimble)
## inputs
n <- 4
tests <- c(
    'chol(A)',
    'forwardsolve(A, b)', 'forwardsolve(A, B)', 'forwardsolve(L, b)', 'forwardsolve(L, B)',
       'backsolve(A, b)',    'backsolve(A, B)',    'backsolve(R, b)',    'backsolve(R, B)',
           'solve(A, b)',        'solve(A, B)',        'solve(R, b)',        'solve(R, B)'
)
## derived objects
testNames <- as.character(sapply(tests, nimble:::Rname2CppName))
testDims <- sapply(grepl('_b_', testNames), function(x) if(x) 1 else 2)
methods <- mapply(function(test, dim) eval(substitute(function() { ans <- TEST; returnType(double(DIM)); return(ans) }, list(TEST = parse(text=test)[[1]], DIM = dim))), tests, testDims)
names(methods) <- testNames
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
nfDef <- nimbleFunction(
    setup = function(A, B, b) {
        R <- chol(A)
        L <- t(R)
    },
    run = function() {},
    methods = methods
)
Rnf <- nfDef(A, B, b)

Cnf <- compileNimble(Rnf)

testOneCase <- function(test, testName, Rnf, Cnf) {
    Rans <- eval(parse(text=test)[[1]])
    Rnfans <- eval(substitute(Rnf$TEST(), list(TEST=as.name(testName))))
    Cnfans <- eval(substitute(Cnf$TEST(), list(TEST=as.name(testName))))
    dif <- max(abs(Rans - Rnfans))
    if(dif > 1E-15) return(1)
    dif <- max(abs(Rans - Cnfans))
    if(dif > 1E-15) return(1)
    return(0)
}
for(i in seq_along(tests)) {
    test <- tests[i]
    testName <- testNames[i]
    testResult <- testOneCase(test, testName, Rnf, Cnf)
    if(testResult != 0) message('failed test: ', test)
}




library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ##ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
        ans4 <- chol(A + A)
    }
)

Cnf <- compileNimble(Rnf)   ## this works




nfDef <- nimbleFunction(
    setup = function(A, B, b) {},
    run = function() {
        ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
    }
)

Rnf <- nfDef(A, B, b)
Cnf <- compileNimble(Rnf)   ## FAIL


library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ans1 <- chol(A)
        print('chol(A)')
        print(ans1)
        ##
        ##ans2 <- forwardsolve(A, b)
        ##print('forwardsolve(A, b)')
        ##print(ans2)
        ####
        ##ans3 <- forwardsolve(A, B)
        ##print('forwardsolve(A, B)')
        ##print(ans3)
        ####
        ##ans4 <- backsolve(A, b)
        ##print('backsolve(A, b)')
        ##print(ans4)
        ####
        ans5 <- backsolve(A, B)
        print('backsolve(A, B)')
        print(ans5)
    }
)

Cnf <- compileNimble(Rnf)


set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
B <- array(1:9, c(3,3))
b <- 1:3

chol(A)
forwardsolve(A, b)
forwardsolve(A, B)
backsolve(A, b)
backsolve(A, B)

Rnf(A, b, B)
Cnf(A, b, B)

##debug(exprClasses_labelForEigenization)
##undebug(exprClasses_labelForEigenization)
##debug(exprClasses_eigenize)
##undebug(exprClasses_eigenize)
##writeCode(nimDeparse(code))
##writeCode(nimDeparse(compileInfo$nimExpr))


## using compareMCMCs, and MCMC comparisons

library(nimble)

code1 <- nimbleCode({ ## really toy model
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    sigma ~ dunif(0.5, 1)
    for(i in 1:3) y[i] ~ dnorm(b, sd = sigma)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5)) ## can also provide constants
results1 <- compareMCMCs(input1, MCMCs = c('nimble')) ## run a single case

results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')) ## run two more cases
results2[[1]]$efficiency ## inspect
results2[[1]]$summary
results2[[1]]$timing
make_MCMC_comparison_pages(results2, 'model1b') ## could generate a comparison of these

## can combine results as follows
## If the same name (“nimble” in this case) was used multiple times one will have to be renamed:
results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]])
results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results')
make_MCMC_comparison_pages(results3, 'model1c')



## testing calc_dmnormConjugacyContributions

nimble:::identityMatrix(3)

calc_dmnormConjugacyContributions(diag(3), diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 7*diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 8*diag(3), 2)
calc_dmnormConjugacyContributions(3*diag(3), 1+diag(3), 1)
calc_dmnormConjugacyContributions(3*diag(3), 2+diag(3), 2)

Cnf <- compileNimble(calc_dmnormConjugacyContributions)

Cnf(diag(3), diag(3), 1)
Cnf(diag(3), 7*diag(3), 1)
Cnf(diag(3), 8*diag(3), 2)
Cnf(3*diag(3), 1+diag(3), 1)
Cnf(3*diag(3), 2+diag(3), 2)



## compiler lesson 101
library(nimble)

nf <- nimbleFunction(
    setup = function() {
        a <- 1
    },
    run = function(arg = double()) {
        b <- a + arg
        returnType(double())
        return(b)
    })

Rnf <- nf()

## nfProc object is created every time we compile a new NF generator
## nfCompileInfo data structure for compilation info RCfunctionCompileClass"
## RCfunProcessing class to manage compilation info for a *single* run/member function
proj <- nimble:::nimbleProjectClass()
nfProc <- nimble:::nfProcessing(Rnf, project = proj)
nfProc$process(control = list(debug = TRUE))

## place to define 1st arg switching at top of genCpp_sizeProcessing
## switch name to nimArr_xxx in genCpp_processSpecificCalls






## testing of new parsed posterior text in conjugacy system
posteriorText = '{a<-1;
b
c
dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                            sd   = (prior_tau + contribution_tau)^(-0.5))
}'

posteriorText = '{ R <- chol(prior_prec + contribution_prec)
                        A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
                        mu <- backsolve(R, forwardsolve(t(R), A))[,1]
                        dmnorm_chol(mean = mu, cholesky = R, prec_param = 1) }'

parsedTotalPosterior <- parse(text = posteriorText)[[1]]
parsedTotalPosterior
if(parsedTotalPosterior[[1]] != '{') parsedTotalPosterior <- substitute({POST}, list(POST = parsedTotalPosterior))
parsedTotalPosterior
prePosteriorCodeBlock <- parsedTotalPosterior[-length(parsedTotalPosterior)]
prePosteriorCodeBlock
posteriorExpr <- parsedTotalPosterior[[length(parsedTotalPosterior)]]
posteriorExpr
rDistribution <- cc_makeRDistributionName(as.character(posteriorExpr[[1]]))
rDistribution
dDistribution <- as.character(posteriorExpr[[1]])
dDistribution
argumentExprs <- as.list(posteriorExpr)[-1]
argumentExprs
argumentNames <- names(argumentExprs)
argumentNames
rCallExpr <- as.call(c(as.name(rDistribution), 1, argumentExprs))
rCallExpr
dCallExpr <- as.call(c(as.name(dDistribution), quote(VALUE), argumentExprs, log = 1))
dCallExpr
posteriorVars <- all.vars(parsedTotalPosterior)
posteriorVars
neededPriorParams <- gsub('^prior_', '', posteriorVars[grepl('^prior_', posteriorVars)])
neededPriorParams
neededContributionNames <- posteriorVars[grepl('^contribution_', posteriorVars)]
neededContributionNames
neededContributionDims <- inferContributionTermDimensions(prior)


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

code <- nimbleCode({ x ~ dnorm(mu, var = v) })
inits <- list(x=0, mu=0, v=1)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
code <- nimbleCode({ x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3]) })
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

md <- nimbleModel(code, inits = inits, returnDef = TRUE, debug=TRUE)

Rmodel <- nimbleModel(code, inits = inits)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node
Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_var
Rmodel$nodes[[nn]]$get_sd
Rmodel$nodes[[nn]]$get_tau
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_cov


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)

## parametrize in terms of prec matrix (Q):
code <- nimbleCode({
    x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3])
})
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

Rmodel <- nimbleModel(code, inits = inits)
Cmodel <- compileNimble(Rmodel)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel[[lifted]]

Q <- A
ch <- chol(Q)
V <- solve(Q)
I <- diag(3)
Q
ch
V
t(ch) %*% ch
Q
backsolve(ch, forwardsolve(t(ch), I))
V

Rmodel$nodes[[nn]]$get_mean
Rmodel$nodes[[nn]]$get_mean()
Rmodel$getParam(nn, 'mean')

Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
Rmodel$getParam(nn, 'cholesky')
ch

Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_prec()
Rmodel$getParam(nn, 'prec')
Q

Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
Rmodel$getParam(nn, 'cov')
V

xx <- forwardsolve(t(ch), I)
t(xx) %*% xx

xx <- backsolve(ch, I)  ## would this work too???
xx %*% t(xx)

t(ch) %*% ch %*% V
t(ch) %*% ch
ch %*% V
Q


nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        print('chol:')
        chh <- model$getParam(node, 'cholesky')
        print(chh)
        print('prec:')
        q <- model$getParam(node, 'prec')
        print(q)
        print('cov:')
        vv <- model$getParam(node, 'cov')
        print(vv)
    }
)

Rnf <- nfDef(Rmodel, nn)
Cnf <- compileNimble(Rnf, project=Rmodel)


ch
Q
V
Rnf$run()
Cnf$run()

## store something into model variable 'Q'
Rmodel$Q <- Q

## calculate deterministic dependents
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node

Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
ch

## request 'prec'
Rmodel$nodes[[nn]]$get_prec  ## direct: model$Q
Rmodel$nodes[[nn]]$get_prec()
Q

## request 'cov'
## here's what we want to get cov: backsolve(ch,  forwardsolve(t(ch), I)  )
Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
V







    


## figuring out R / nimble functions:
## inverse, solve, forwardsolve, backsolve

library(nimble)
set.seed(0)
A <- array(rnorm(9), c(3,3))
A <- t(A) %*% A

nfDef <- nimbleFunction(
    setup <- function(A) {
        D <- diag(3)
        b <- 1:3
    },
    methods = list(
        inv = function() {
            Ainv <- inverse(A)
            print(Ainv)
        },
        ch = function() {
            R <- chol(A)
            print(R)
        },
        fs = function() {
            R <- chol(A)
            L <- t(R)
            Linv <- forwardsolve(L, D)
            print(Linv)
        },
        bs = function() {
            R <- chol(A)
        }
    )
)

Rnf <- nfDef(A)
Cnf <- compileNimble(Rnf)

print('inverse of A')
solve(A)
Rnf$inv()
Cnf$inv()

print('chol of A')
chol(A)
Rnf$ch()
Cnf$ch()

R <- chol(A)
L <- t(R)

print('forwardward solve')
forwardsolve(L, diag(3))

print('backward solve')
backsolve(R, diag(3))



## example of using nf$getDefinition() and nf_getDefinition()
library(nimble)
nfDef <- nimbleFunction(
    setup = function(a) { },
    run = function(b = double()) {
        a <<- a + b
        returnType(double())
        return(a)
    },
    methods = list(
        nothing = function() {
            print(a)
        })
)
Rnf <- nfDef(3)
Cnf <- compileNimble(Rnf)
getDefinition(nfDef)
Rnf$getDefinition()
Cnf$getDefinition()

## example of spec$getSamplerDefinition(..)
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    for(i in 1:140)
        c[i] ~ dnorm(0,1)
})

Rmodel <- nimbleModel(code, check = FALSE)

spec <- configureMCMC(Rmodel, print = TRUE)
spec$printSamplers()
spec$getSamplerDefinition(1)
spec$getSamplerDefinition('b')

spec$getSamplers()
spec$printSamplers()
spec$addSampler(type='RW_block', target=c('a', 'b'))
spec$getSamplers()
spec$printSamplers()


## performance comparison of of new dynamic conjugate samplers, rats example
library(nimble)
rats <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
Rmodel <- nimbleModel(rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits)
spec <- configureMCMC(Rmodel, print = TRUE)
out <- MCMCsuite(
    rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits,
    niter = 100000,
    MCMCs = c('staticConj', 'dynamicConj'),
    MCMCdefs = list(
        staticConj  = quote({ nimbleOptions(useDynamicConjugacy = FALSE)
                              configureMCMC(Rmodel) }),
        dynamicConj = quote({ nimbleOptions(useDynamicConjugacy = TRUE)
                              configureMCMC(Rmodel) })),
    makePlot = FALSE
)
sampdiff <- out$samples['staticConj',,] - out$samples['dynamicConj',,]
all(sampdiff == 0)
## [1] TRUE                      ## samples from each MCMC are identical
out$timing
## staticConj    dynamicConj
##     26.309         17.801     ## ~ 33% decrease in conjugate sampler runtime



## testing proper dmnorm conjugate calculations for Chris

## code = nimbleCode({
##    y[1:3]  ~ dmnorm(mu [1:3], cov = covy [1:3,1:3])
##    mu[1:3] ~ dmnorm(mu0[1:3], cov = covmu[1:3,1:3])
## })
## Concerned about some explicit inverses in conjugate update for mu,
## compared to use of already existing Cholesky factors.  In the R sampler,
## 1) we have: inverse(prior_prec + contribution_prec)
## 2) we use get_prec, which when we specify a dmnorm with a covariance
## has code like the following (in calc_dmnormAltParams()):
##  tmp <- t(cholesky) %*% cholesky
##  return(inverse(tmp))

library(nimble)

nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions()

code = nimbleCode({
   mu[1:3] ~ dmnorm(mu0[1:3], cov=covmu[1:3,1:3])
   y[1:3] ~ dmnorm(mu[1:3], cov = covy[1:3,1:3])
   a ~ dnorm(y[1], y[1])
   b ~ dgamma(a, y[1])
})

cov <- diag(3)
constants <- list(covy = cov, covmu = cov, mu0 = rep(0,3))
data <- list(y = 1:3)
inits <- list(mu = 4:6, a = 1, b = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$printSamplers()
spec$getSamplerDefinition('mu')
Rmcmc <- buildMCMC(spec)


Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$samplerFunctions$contentsList[[2]]$scale
Cmcmc$samplerFunctions$contentsList[[2]]
Cmcmc$samplerFunctions$contentsList[[2]][['scale']]
Cmcmc$samplerFunctions$contentsList[[2]]$scale



set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## testing what works and what doesn't with nimbleFunctionLists, getParam, etc..

library(nimble)

getNimbleOption('checkModel')
nimbleOptions(checkModel = FALSE)
getNimbleOption('checkModel')

code <- quote({
    x ~ dgamma(1, 1)
    y[1] ~ dnorm(1, tau = 1*x)
    y[2] ~ dnorm(2, tau = 2*x)
    y[3] ~ dnorm(3, tau = 3*x)
})
Rmodel <- nimbleModel(code, data=list(), inits=list(x=1, y=11:13))

nfDef <- nimbleFunction(
    setup = function(model) {
        nl <- model$expandNodeNames('y')
        N <- length(nl)
        nfl <- nimbleFunctionList(node_stoch_dnorm)
        for (i in 1:N) nfl[[i]] <- model$nodeFunctions[[nl[i]]]
    },
    run = function() {
        for(i in 1:N) {
            print('i = ', i)
            val <- nfl[[i]]$get_mean()
            print('mean = ', val)
            val <- nfl[[i]]$get_tau()
            print('tau = ', val)
            val <- nfl[[i]]$get_var()
            print('var = ', val)
            val <- nfl[[i]]$get_value()
            print('value = ', val)
        }
        returnType(double(0))
        return(1)
    }
)
Rnf <- nfDef(Rmodel)
Rnf$run()
Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)
Cnf$run()


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
    z ~ dbeta(10, 10)
    zz[1] ~ dbin(z, 10)
    zz[2] ~ dbin(z, 20)
    aaa ~ dbeta(10, 10)
    bbb[1] ~ dbin(aaa, 10)
    bbb[2] ~ dbin(aaa, 20)
    bbb[3] ~ dbin(aaa, 20)
    d ~ dnorm(0, 1)
    dd ~ dnorm(0, 1)
    ddd ~ dnorm(0, 1)
    e ~ dnorm(d + dd + ddd, 1)
    f ~ dlnorm(d + dd + ddd, 1)
    ff ~ dlnorm(d + dd, 1)
    g ~ dnorm(0, 1)
    gg ~ dnorm(5+g, 1)
})
constants <- list()
inits <- list(x = 0.5, z=0.6, aaa=0.3, d=0, dd=0, ddd=0, g=0)
data = list(y = c(3,4), zz=c(5,5), bbb=c(3,3,3), ff=1, e=0, f=1, gg=0)
Rmodel <- nimbleModel(code, constants, data, inits)


spec <- configureMCMC(Rmodel)

spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Rmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(3)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]
##     x z aaa d dd ddd g
##[1,] 0 0   0 0  0   0 0
##[2,] 0 0   0 0  0   0 0
##[3,] 0 0   0 0  0   0 0
Rsamples[, Rnames]
##             x         z       aaa         d         dd        ddd         g
##[1,] 0.3254630 0.3746455 0.3658422 0.6362147 -0.2698403 -1.1333402 -3.156596
##[2,] 0.1975665 0.3995488 0.3210974 0.3695457 -0.2843177 -0.2239394 -2.711577
##[3,] 0.1900000 0.4199454 0.2194685 0.5430496 -0.9140867  0.1178770 -2.233141


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
})
data = list(y = c(3,4))
constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples
as.numeric(Rsamples) - c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)



## testing of new dynamic conjugate samplers
### checkConjugacy_demo3_run.R - various conjugacies
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})
data = list()
constants <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('x', 'c'), control = list(scale = 0.01))
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]

Rsamples[, Rnames] - cbind(c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977), c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))


## testing of new dynamic conjugate samplers
### Dirichlet-multinomial conjugacy
# as of v0.4, exact numerical results here have changed because
# ddirch now sometimes returns NaN rather than -Inf (when an
# alpha is proposed to be negative) -- this changes the RNG
# sequence because NaN values result in no runif() call in decide()
# single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- quote({
    y[1:K] ~ dmulti(p[1:K], n);
    p[1:K] ~ ddirch(alpha[1:K]);
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001);
    }
               })
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
constants <- list(n=n, K=K)
data <- list(y = y)
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)[6:10]
p
apply(samples, 2, mean)[6:10] - p



## testing of new dynamic conjugate samplers
## rats example
library(nimble)
lst <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
code <- lst$model
data <- lst$data[c('Y')]
constants <- lst$data[c('N', 'T', 'x')]
inits <- lst$inits

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Cmcmc$run(5000)
samples <- as.matrix(Cmcmc$mvSamples)

a <- samples
b <- samples
a
b
a-b
any(a-b != 0)




## phase transition regression example for Mark.

## generate data
sigma <- 5
a <- 2
b <- 0.5
t0 <- 30
x <- 1:100
n <- length(x)
y <- rnorm(n, a + (x>t0) * (b*(x-t0)), sd=sigma)
plot(x, y)

## fit NIMBLE model
library(nimble)

code <- nimbleCode({
    a ~ dunif(-100, 100)
    b ~ dunif(-100, 100)
    t0 ~ dunif(0, 100)
    sigma ~ dunif(0, 100)
    for(i in 1:n) {
        phase[i] <- step(x[i] - t0)
        mu[i] <- a + phase[i] * b * (x[i]-t0)
        y[i] ~ dnorm(mu[i], sd = sigma)
    }
})
constants <- list(n=n, x=x)
data <- list(y=y)
inits <- list(sigma=1, a=0, b=1, t0=50)

## MCMC Suite will let us do a few MCMCs quickly, and look at results easily
out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_slice', 'bugs'),
    niter = 20000,
    burnin = 2000,
    makePlot = TRUE,
    savePlot = FALSE,
    calculateEfficiency = TRUE
)

out$summary
out$efficiency
dim(out$samples)
out$timing





## nimble problem with using 'sd' name for 'tau' parameter, e.g.
library(nimble)

code <- nimbleCode({
    sd ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd)
    shape1 ~ dnorm(0, 1)
    f ~ dbeta(mean = y, sd = shape1)
})

m <- nimbleModel(code, inits = list(sd=.1))

m$sd

m$getModelDef()$printDI()

code <- nimbleCode({
    sd2 ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd2)
})
m <- nimbleModel(code, debug=TRUE)

m$sd2






     
## Plot the fake MCMC output
denoverplot(fakemcmc, fakemcmc2)
     denoverplot(fakemcmc, fakemcmc2, style="plain",
                 col=mcmcplotsPalette(3, type="grayscale"),
                 ci=0.95, greek=TRUE)
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data")
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data", greek=TRUE)

setwd('/Users/dturek/GitHub/legacy/autoBlock/')
getwd()

source("autoBlock.R")
load(file.path("data", "model_test.RData"))
saveSamples <- TRUE
niter <- 10000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dftest <- ab$summary
save(dftest, file = file.path("results_samples", "results_test.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_test_samples.RData"))
}

class(burnedSamplesList)
names(burnedSamplesList)

i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])

dimnames(ab$samples[[3]])


ret$summary

dftest <- autoBlock(code, constants, data, inits, 10000, runList)$summary
save(dftest, file = file.path("results_hclust_complete2", "results_test.RData"))


rm(list=ls())
getwd()
load('results_samples/results_test.RData')
ls()
dftest

rm(list=ls())
getwd()
load('results_samples/results_test_samples.RData')
ls()
niter
class(burnedSamplesList)
names(burnedSamplesList)
i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])







## renaming stupid data frame from the autoBlock results
setwd('~/GitHub/legacy/autoBlock')

rm(list=ls())
file <- 'results/results_samplingEfficiency.RData'
load(file)
ls()
dfsamplingEfficiency <- dfSamplingEfficiency
rm(dfSamplingEfficiency)
ls()
save(list='dfsamplingEfficiency', file = file)

rm(list=ls())
file <- 'results/results_computationalRequirement.RData'
load(file)
ls()
dfcomputationalRequirement <- dfComputationalRequirement
rm(dfComputationalRequirement)
ls()
save(list='dfcomputationalRequirement', file = file)





## working on the spatial capture-recapture (SCR) models

tr<-seq(15,85, length=10)

X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. traps

plot(X, xlim=c(0,100), ylim=c(0,100), pch=3, cex=0.75)

set.seed(10)
xlim <- c(0,100); ylim <- c(0,100)      # Area 100*100=1e4
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/10000
mu <- 50                 # Density
N <- rpois(1, mu*A) ;N   # Generate population


s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(s, pch=16, col=2)

sigma <- 5
lambda0 <- 0.4
J <- nrow(X)
K <- 5
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
    dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2]-s[,2])^2)
    lambda <- lambda0*exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        yy[,j,k] <- rpois(N, lambda)
    }
}

n <- apply(yy, c(2,3), sum)

# Plot capture events
tot<-apply(n, 1,sum)
symbols(X, circles=tot, inches=F, bg="#00000022", add=T)
points(X, pch=3, cex=0.75); points(s, pch=16, col=2)


library(nimble)

## define the model
code <- nimbleCode({
    sigma ~ dunif(0,10)
    lam0 ~ dunif(0,5)
    psi ~ dbeta(1,1)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        s[i,1] ~ dunif(xlim[1], xlim[2])
        s[i,2] ~ dunif(ylim[1], ylim[2])
        for(j in 1:J) {# Number of traps
            dist[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
            lam[i,j] <- lam0*exp(-dist[i,j]/(2*sigma^2))*z[i]
        }
    }
    for(j in 1:J){
        bigLambda[j] <- sum(lam[1:M,j])
        for(k in 1:K) {
            n[j,k] ~ dpois(bigLambda[j])
        }
    }
    N <- sum(z[1:M])
})


M<-200

constants <- list(M = M, K=K, J=J)
n1<-apply(n,1,sum)
data<-list(n=n, X=X, xlim=xlim, ylim=ylim)
s<-cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
z<-rep(1,M)
inits <- list (sigma=0.5, lam0=0.1, s=s, z=z)

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE) ## check=FALSE is faster

mcmcspec<-configureMCMC(Rmodel, print=TRUE, monitors = c("N", "lam0", "psi", "sigma"))
scrMCMC <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel) 
CscrMCMC <- compileNimble(scrMCMC, project = Rmodel)

## It's pretty slow to run so I just want an execution time first for a small sample
t1_100 <- system.time(CscrMCMC$run(100))

## And now I want a decent sample
t1_20k <- system.time(CscrMCMC$run(20000))
t1_20k_samples <- as.matrix(CscrMCMC$mvSamples)
save(t1_100, t1_20k, t1_20k_samples, file = "case1results.Rdata")


## checking dcat distribution for Chris

library(nimble)

code <- nimbleCode({
    y ~ dcat((p[1:2,1:2] %*% ones[1:2])[1:2,1])
})
constants <- list(ones = c(1,1))
data <- list()
inits <- list(y = 2, p = diag(c(.5,.5)))

Rmodel <- nimbleModel(code, constants, data, inits)




library(nimble)
Rmodel <- readBUGSmodel('rats', dir = '~/GitHub/nimble/nimble/packages/nimble/inst/classic-bugs/vol1/rats/')

customSpec <- configureMCMC(Rmodel)
customSpec$getSamplers()
customSpec$removeSamplers(1:65)
customSpec$getSamplers()
customSpec$addSampler(target = 'alpha.c', type = 'slice')
customSpec$addSampler(target = 'beta.c', type = 'slice')
customSpec$addSampler(target = 'tau.c', type = 'slice')
customSpec$addSampler(target = 'tau.alpha', type = 'slice')
customSpec$addSampler(target = 'tau.beta', type = 'slice')
customSpec$getSamplers()
alphaNames <- Rmodel$getDependencies('tau.alpha', self = FALSE, stochOnly = TRUE)
alphaNames
length(alphaNames)  ## 30
for(aN in alphaNames) customSpec$addSampler(target = aN, type = 'slice')
customSpec$getSamplers()
betaNames <- Rmodel$getDependencies('tau.beta', self = FALSE, stochOnly = TRUE)
betaNames
length(betaNames)  ## 30
for(bN in betaNames) customSpec$addSampler(target = bN, type = 'slice')
customSpec$getSamplers()
customSpec


## making a simple MCMC traceplot, and density plot

n <- 1000
x <- rnorm(n, 3, 1)
dev.new(width=4, height=3)
plot(1:n, x, type='l')
dev.copy2pdf(file = '~/Downloads/traceplot.pdf')
dev.new(width=4, height=3)
plot(density(x))
dev.copy2pdf(file = '~/Downloads/densityplot.pdf')


## make sure I get the empirical covariance right

x <- matrix(c(1,2,1,2,1,4,6,5,7,4), 5, 2)
timesRan <- 5
statSums <- apply(x, 2, sum)
statSums <- t(statSums)
statProds <- matrix(0, 2, 2)
for(i in 1:timesRan) statProds <- statProds + t(t(x[i,])) %*% t(x[i,])
statSums
statProds
(statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
cov(x)
colMeans <- apply(x, 2, mean)
xCentered <- x - rep(colMeans, each=timesRan)
xCentered
(t(xCentered) %*% xCentered) / (timesRan-1)

set.seed(0)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 4
mu <- 1:d
ch <- array(rnorm(d^2), c(d,d))
Sigma <- ch %*% t(ch)
print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$addSampler('x', 'RW_block_NEW', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(1000)
samples <- as.matrix(Cmcmc$mvSamples)
cov(samples)
           x[1]       x[2]     x[3]       x[4]
x[1]  2.0934846 -0.1465907 1.107342  1.1252451
x[2] -0.1465907  4.1696671 1.644427 -0.5950559
x[3]  1.1073420  1.6444267 2.581730  1.5097074
x[4]  1.1252451 -0.5950559 1.509707  1.8741393
samples[c(1,250,500,750,1000),]
           x[1]      x[2]       x[3]     x[4]
[1,]  1.0000000  1.000000  1.0000000 1.000000
[2,]  0.2162332  2.260594  1.9804354 2.913190
[3,]  1.8118014  3.966497  2.7563534 2.552825
[4,] -1.6015689 -1.004725 -0.9864974 1.065231
[5,]  0.8037007  2.496745  1.2659404 2.208456






## looking into RStudio hanging problem (once again... Nov 2015)

library(nimble)
code <- nimbleCode({
  a ~ dnorm(0, 1)
  sd ~ dunif(0, 100)
  for(i in 1:N) {
    b[i] ~ dnorm(a, sd = sd)
  }
})
N <- 100
constants <- list(N = N)
data <- list()
inits <- list(a = 0, sd = 1, b = rep(0, N))


Rmodel <- nimbleModel(code, constants, data, inits)


## NEED TO FIND THE PROBLEM IN MODEL CHECKING

## EXAMPLE 1:
rm(list = ls())
library(nimble)
classicDyesModel <- readBUGSmodel('dyes', dir = getBUGSexampleDir('dyes'))

## EXAMPLE 2:
set.seed(0)
library(nimble)
n_obs <- 50
### Parameters
mu <- 2
rho_a <- 1
rho_b <- 0.1
monitored_param <- c('mu', 'rho_a', 'rho_b', 'ngis')
SV =  rgamma(n=n_obs, 1, 1)
### Data generation
lambda <-  exp(mu)*SV
ngis  <-  rpois(n=n_obs, lambda=lambda )
Y <- rgamma(n=n_obs, ngis*rho_a, rho_b)
### Remove zero values from ngis
ngis <- ngis[which(ngis>0)]
abs <- which(Y==0)
n_abs <- length(abs)
pres <- which(Y>0)
n_pres <- length(pres)

# Install from source for OS X, Linux, or Windows:
#install.packages("nimble", repos = "http://r-nimble.org", type = "source")
M1_CPG_nimbleCode <- nimbleCode(
    {
        ## Covariates
        for (s in 1:n_obs){
            lambda[s] <- exp(mu)*SV[s]
        }
        ## Observation model
        ## Strictly positive Observation
        for ( s in 1:n_pres){
            ngis[s] ~ T(dpois(lambda[pres[s]]), 1, )  
            Y[pres[s]] ~ dgamma(rho_a * ngis[s], rho_b)
        }
        ## Zero Observation
        for ( s in 1:n_abs){
            proba[s] <- 1 - exp(-lambda[abs[s]])
            Y[abs[s]] ~ dbern(proba[s])
        }
        ## Prior
        rho_a ~ dgamma(0.01, 0.01)
        rho_b ~ dgamma(0.01, 0.01)
        mu ~ dnorm(0, 0.01)
    })

constants_list <- list(n_obs = length(Y), n_abs=n_abs, n_pres=n_pres, abs=abs, pres=pres)
nimble_data <- list(Y = Y, SV=SV)
## Initial values
inits_fc <- function(chain_id = 1){
  mu <- rnorm(1,0,1)
  rho_a <- runif(1,0,10)
  rho_b <- runif(1,0,10)
  ngis <- round(runif(n_pres, 2, 5), digits=0)
  list(mu=mu, rho_a=rho_a, rho_b=rho_b ,ngis=ngis)
}
init_ll <- lapply(1, function(id) inits_fc(chain_id = id))

Rmodel <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll[[1]], check = TRUE)

calculate(Rmodel, 'ngis[1]')

Rmodel$ngis

check <- TRUE
check <- FALSE

md <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll, returnDef = TRUE)

M1_CPG_nimble <- md$newModel(data =nimble_data, inits = init_ll, check = check)

debug(md$newModel)
md$newModel(data =nimble_data, inits = init_ll, check = TRUE)

debug(model$check)





## random stuff for CR paper....

rm(list=ls())
library(nimble)
##source('~/GitHub/legacy/dipper/dipperCode.R')

Rmodel <- nimbleModel(code_dipper, constants, data, inits)
lnodes <- Rmodel$expandNodeNames('x')
length(lnodes) - nind

rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)

con <- orchidDHMM$constants
con$f
con$nind
#### actually 250 individuals in original dataset!!!

con$f
12 - con$f
sum(12 - con$f)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
or <- results$df[results$df$model == 'orchid', ]
or
or[or$mcmc %in% c('jags', 'nimbleDHMM2'), ]
or1 <- or[or$param == 's[1]', ]
data.frame(mcmc = or1$mcmc, runtime_minutes = or1$timing/60)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
go <- results$df[results$df$model == 'goose', ]
go
go1 <- go[go$param == 'p[1]', ]
data.frame(mcmc = go1$mcmc, runtime_minutes = go1$timing/60, runtime_hours = go1$timing/60/60)


rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(gooseDHMM$code, gooseDHMM$constants, gooseDHMM$data, gooseDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)



## testing the MVN conjugacy test from test-mcmc.R
library(nimble)
set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

##### not currently working - see Perry's email of ~ 10/6/14
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3]), Q[1:3, 1:3])
## })

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})

## Simplest version of model w/o 'a' and 'B'
## a = rep(0,3)
## B = diag(rep(1,3))
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(mu[1:3], Q[1:3, 1:3])
## })


mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)

muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
set.seed(0)
Rmcmc$run(5)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0)
Cmcmc$run(5)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples

Rsamples - Csamples

Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
muMeanTrue

apply(samples, 2, mean) - muMeanTrue


source('~/GitHub/nimble/nimble/packages/nimble/inst/tests/test_utils.R')
##test_mcmc
test_mcmc(model = code, name = 'two-level multivariate normal', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))



## working the basketball "streaks" problem emailed out by Aaron Strauss
p <- 0.1
ns <- 1
run <- function(p, ns) {
    out <- numeric()
    n <- 0
    while(n < ns) {
        nmade <- 1
        while(rbinom(1, 1, p) == 1) nmade <- nmade + 1
        out <- c(out, nmade)
        n <- length(out)
    }
    return(out)
}

pOfStreaks <- function(p, streaks, niter=1000) {
    pMissing <- missing(p)
    ns <- length(streaks)
    suc <- 0
    for(i in 1:niter) {
        if(pMissing) p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) suc <- suc+1
    }
    prob <- suc/niter
    return(prob)
}

pOfStreaks(p=0.2, streaks=4, niter=100000)
0.2^3*(.8)

p <- 0.7
n <- 4
m <- 7
pOfStreaks(p=p, streaks=c(n,m), niter=1000000)
p^(n+m-2)*(1-p)^2

n <- 2
m <- 1
pOfStreaks(streaks=c(n,m), niter=1000000)
1/12

pOfPlessthanPROB <- function(streaks, PROB, niter=1000) {
    ns <- length(streaks)
    suc <- 0
    ranSoFar <- 0
    while(ranSoFar < niter) {
        p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) {
            if(p < PROB) suc <- suc+1
            ranSoFar <- ranSoFar+1
        }
    }
    prob <- suc/niter
    return(prob)
}

pOfPlessthanPROB(streaks=c(2,1), PROB=0.5, niter=20000)
11/16

## continusing the "sreaks" basketball problem, now in NIMBLE

library(nimble)

code <- nimbleCode({
    ##p ~ dunif(0, 1)
    p ~ dbeta(1, 1)
    y[1] ~ dnegbin(prob=p, size=1)
    y[2] ~ dnegbin(prob=p, size=1)
})
constants <- list()
data <- list(y = 0:1)  ## the 'streak' history of {2,1} in the original problem
inits <- list(p = 0.5) ## translates into c(1,0) in terms of the negative-binomial distribution

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
#Rmcmc$samplerFunctions$contentsList[[1]]$run
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
debug(Rmcmc$run)
debug(samplerFunctions[[1]]$run)
Rmcmc$run(1)
set.seed(0)
Cmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
Rsamples
Csamples
Rsamples - Csamples

apply(samples, 2, mean)

mean(samples[,1] > 0.5)
5/16

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_noConj'),
    summaryStats = c('mean', 'median', 'sd', 'function(x) mean(x>0.5)'),
    calculateEfficiency = TRUE,
    makePlot = TRUE,
    savePlot = FALSE,
    niter = 100000
)

out$timing
out$summary
11/16


## creating a really small, blank figure, for use in captions in CR-MCMC paper
dev.new(width=3, height=2)
plot(1, 1)
dev.copy2pdf(file='~/GitHub/nimble/nimblePapers/CR-MCMC/blankFigure.pdf')

## testing the problem with naming a variable 'ans' in a nimbleFunction
library(nimble)

rcFunction1 <- nimbleFunction(
    run = function() {
        declare(ans, double(1, 5))
        for(i in 1:5) ans[i] <- i
        print(ans)
    })

rcFunction1()

rcFunction2 <- nimbleFunction(
    run = function() {
        declare(functionAsList, double(1, 5))
        for(i in 1:5) functionAsList[i] <- i
        print(functionAsList)
    })

rcFunction2()

code <- nimbleCode({
    mu ~ dnorm(0, 1)
    mu2 <- timesTwo(mu)
})

Rmodel <- nimbleModel(code, inits = list(mu=1), check=FALSE)

##undebug(calculate)
##debug(nimble:::rCalcNodes)
##Rmodel$nodes$mu2$calculate
##Rmodel$nodes$mu2$calculate()
##Rmodel$nodes$mu2$simulate
##Rmodel$nodes$mu2$simulate()

calculate(Rmodel)
Cmodel <- compileNimble(Rmodel)

Rmodel$mu
Rmodel$mu2
Cmodel$mu
Cmodel$mu2

Rmodel$mu <- 2
calculate(Rmodel)
Rmodel$mu2


## finding the posterior pairwise correlations of
## the blocked parameters for the orchid model
## (multistate CR-MCMC, userDistMCMC)

library(nimble)
rm(list=ls())
source('~/GitHub/userDistMCMC/defs.R')
load('~/GitHub/userDistMCMC/models.RData')

Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$removeSamplers(numeric(0))
spec$removeSamplers(c('b[1:2]','c[1:2]'))
spec$addSampler(c('b[1]','b[2]'), 'RW_block')
spec$addSampler(c('c[1]','c[2]'), 'RW_block')
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
## n <- 10000
n <- 50000
system.time(Cmcmc$run(n)) / 60
samples <- as.matrix(Cmcmc$mvSamples)
samples2 <- samples[(n/2+1):n, ]   ## truncate to second half of samples
apply(samples2, 2, mean)
Cov <- cov(samples2)
Cov
Cor <- cov2cor(Cov)
Cor
sort(as.numeric(abs(Cor)))
Cor['b[1]', 'b[2]']
## [1] 0.9866091
Cor['c[1]', 'c[2]']
## [1] 0.9675666


## printing the numeric performance results for CR-MCMC paper, or userDistMCMC paper
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')
source('~/GitHub/userDistMCMC/defs.R')
options(digits = 3)
df <- results$df
for(mod in unique(df$model)) {
    dfmod <- df[df$model==mod, ]
    message('******************************************')
    message(mod, ' model')
    message('******************************************')
    for(mc in unique(dfmod$mcmc)) {
        dfmodmcmc <- dfmod[dfmod$mcmc==mc, ]
        message(mc, ' MCMC:')
        message('minimum: ',  min(dfmodmcmc$Efficiency))
        message('mean: ',    mean(dfmodmcmc$Efficiency))
    }
    message('******************************************')
}

df[df$model=='orchid',]
df[df$model=='goose',]

results$check()
results$quickplot()


## problem with multivariate conjugate samplers... ?
## fixed.  problem was in vector demoting in the conjugacy definition
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
mu <- rep(0, d)
Sigma <- diag(d)
constants <- list(d=d, mu=mu, Sigma=Sigma)
Rmodel <- nimbleModel(code, constants, data=list(y=rep(0,d)), inits=list(x=rep(0, d)))
Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
solve(cov(samples))

library(nimble)
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dnorm
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dbeta
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dmnorm



## testing different priors for rho in package gpmanagement
library(gpmanagement)
library(ggplot2)
example_data <- function(p, x0, Tobs){
    f <- function (x, h){
        sapply(x, function(x) {
            x <- pmax(0, x - h)
            x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
        })
    }
    sigma_g <- 0.1
    pdfn <- function(x, mu, sigma = sigma_g){
        dlnorm(x, log(mu), sdlog = sigma)
    }
    z_g <- function() rlnorm(1, 0, sigma_g)
    x <- numeric(Tobs)
    x[1] <- x0
    for(t in 1:(Tobs-1))
        x[t+1] = z_g() * f(x[t], h=0)
    obs <- data.frame(x = c(0, 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                      y = c(0, 
                          x[2:Tobs]))
    obs
}
myfit <- function(xObs, yObs, priors = NULL, niter = 100000){    ## default value for priors
  xPred <- seq(0, 1.1 * max(xObs), length = 50)
  fit <- gp_setup(xObs, yObs, xPred, priors)     ## included priors argument
  Cmcmc <- fit$Cmcmc 
  Cpred <- fit$Cpred
  Cmodel <- fit$Cmodel
  system.time(Cmcmc$run(niter))
  samples <- as.matrix(Cmcmc$mvSamples)
  system.time(Cpred$run(samples))
  E <- Cpred$getE()
  V <- sqrt(diag(Cpred$getC()))
  list(samples = tidyr::gather(as.data.frame(samples)),
       pred = data.frame(x = xPred, y = E, ymin = E - V, ymax = E + V))
}

densPlot <- function(fit, paramName, paramPrior) {
        dev.new()
        samp <- fit$samples[fit$samples$key==paramName, 'value']
        plot(density(samp), main=paste0('black posterior, red prior: ', paramName, ' ~ ', deparse(paramPrior)))
        xs <- seq(range(samp)[1], range(samp)[2], length.out = 500)
        yCall <- as.call(c(list(as.name(paramPrior[[1]])), list(quote(xs)), as.list(paramPrior[-1])))
        ys <- eval(yCall, envir = environment())
        lines(xs, ys, col='red')
}

## Carrying capcacity at 10, critical point at 5, x0=6
## 100, 50, x0=60
plotPred <- function(ccap, crit, x0, Tobs=40, plotData=FALSE, rhoPrior, plotRhoDensity=TRUE, niter=100000, plotSigGPDensity=FALSE, plotSigOEDensity=FALSE, plotPredictionIntervals=TRUE) {
    set.seed(0)
    obs <- example_data(p=c(2,ccap,crit), x0=x0, Tobs=Tobs)
    if(plotData) {dev.new(); qplot(seq_along(obs$x[-1]), obs$x[-1]) + geom_line()}
    priors <- expression({
        rho ~ X
        sigGP ~ dunif(0, 1e5)
        sigOE ~ dunif(0, 1e5)
    })
    if(missing(rhoPrior)) stop('must specify a prior for rho')
    rhoPrior[-1] <- lapply(rhoPrior[-1], function(param)
        eval(eval(substitute(substitute(PARAM, list(xObs=obs$x[-1])), list(PARAM=param)))))
    print(rhoPrior)
    priors[[1]][[2]][[3]] <- rhoPrior
    fit <- myfit(obs$x, obs$y, priors = priors, niter=niter)
    if(plotRhoDensity) densPlot(fit, 'rho', rhoPrior)
    if(plotSigGPDensity) densPlot(fit, 'sigGP', priors[[1]][[3]][[3]])
    if(plotSigOEDensity) densPlot(fit, 'sigOE', priors[[1]][[4]][[3]])
    if(plotPredictionIntervals) {
        dev.new()
        print(ggplot2::ggplot(fit$pred) + 
                  geom_ribbon(aes(x = x, y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
                      geom_line(aes(x = x, y = y), size = 1) + 
                          geom_point(data = obs, aes(x, y)) +
                              coord_cartesian(xlim = range(c(obs$x, fit$pred$x)), 
                                              ylim = range(c(obs$y, fit$pred$y))) +
                                                  ggtitle(paste0('rho ~ ', deparse(rhoPrior))))
    }
}

rp <- quote(dunif(0, diff(range(xObs))/sqrt(6)/2))

plotPred(rhoPrior = rp, ccap= 10, crit= 5, x0= 6)
plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60)

plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60, plotSigGPDensity=TRUE, plotSigOEDensity=TRUE)

plotPred(rhoPrior = quote(dunif(0,200)), ccap=10, crit=5, x0=6)

for(uMax in 5:6) {
    rp <- substitute(dunif(0, UMAX), list(UMAX = as.numeric(uMax)))
    plotPred(rhoPrior = rp, ccap=10, crit=5, x0=6, niter=10000)
}



## testing ggplot() inside a loop??
library(ggplot2)
df <- data.frame(a=1:2, b=3:4)
for(i in 1:3) {
    dev.new()
    ## need explict print(...) to make ggplot() plot, when inside a loop!!!
    print(ggplot(df, aes(a, b)) + geom_line() + ggtitle(i))
}



## testing use of substitute...
## which I used to be more comfortable with!

f <- function(a) {
    b <- substitute(a)
    print(b)
    print(class(b))
    vv <- eval(b)
    print(vv)
    print(class(vv))
}

f(dnorm(0,1))
f(4+5)

x <- f(quote(dnorm(0,1)))
x
class(x)


## testing NIMBLE compiler to handle pi
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        a <- pi
        print(a)
        b <- 2 * pi
        print(b)
    }
)
Rnf <- nfDef()
Rnf$run()
Cnf <- compileNimble(Rnf)    ## ERROR
Cnf$run()


## example of trying to use NIMBLE functions in another package
## this uses my minimalist nimtest package; same thing I sent to Duncan
library(devtools)
install_github('danielturek/nimtest')
library(nimtest)
library(nimble)
## First Test:
## a nimbleFunction defined as a package function
nfDefinition
Rnf <- nfDefinition()
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
## Second Test:
## A package function which internally defines and returns a nimbleFunction definition
## This one gives ugly looking warnings
nextFunction
def <- nextFunction()        ## WARNINGS
Rnf <- def()
Cnf <- compileNimble(Rnf)    ## WARNINGS
Rnf$run()
Cnf$run()

## alternate.  from shell:
## R CMD BUILD nimtest
## R CMD install nimtest_0.0.0.9000.tar.gz
library(nimtest)
library(nimble)
nfDefinition    ## finds it -- it created the NIMBLE function at build time
Rnf <- nfDefinition()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
nextFunction
nextFunction()
def <- nextFunction()
Rnf <- def()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()


## testing out the new elliptical slice sampler (ess)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
##mu <- rep(0, d)
mu <- 1:d
##Sigma <- diag(d);  print(Sigma)
ch <- array(c(1, 0.7, 0, 2), c(2,2));  Sigma <- ch %*% t(ch);  print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print=FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0); Cmcmc$run(20000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
Sigma


## conjugate MVN-MVN
set.seed(1)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
})
d <- 3
mu_x <- rnorm(d)
temp <- array(rnorm(d^2), c(d,d))
prec_x <- solve(temp %*% t(temp))
temp <- array(rnorm(d^2), c(d,d))
prec_y <- solve(temp %*% t(temp))
y <- rnorm(d)
constants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
data <- list(y = y)
inits <- list(x = rep(0, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

## block sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)

## ess sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)



## testing elliptical slice sampler on 'spatial' model
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_spatial.RData')

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('block', 'ess'),
    MCMCdefs = list(
        block = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'RW_block', print=FALSE)
            spec
        }),
        ess = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'ess', print=FALSE)
            spec
        })),
    niter = 10000,
    burnin = 2000,
    monitors = c('mu', 'sigma', 'rho', 'g'),
    makePlot = FALSE
)

out$timing
out$summary[,,1:5]    ## WRONG again...
cov(t(out$samples['block', -(1:3), ])) - cov(t(out$samples['ess', -(1:3), ]))






## Taken from Carl's lab notebook:
## Gpdp Via MDPtoolbox Cont; 30 Jul 2015
## http://www.carlboettiger.info/2015/07/30/gpdp-via-mdptoolbox-cont.html
## (Markov Decision Processes (MDP) toolbox)
##knitr::opts_chunk$set(eval=FALSE)
##devtools::install_github("cboettig/gpmanagement@b3b765cbceb51c9b0b8cb2724e395353ec365df9")
library("MDPtoolbox")
library("gpmanagement")
library("tidyr")
library("dplyr")
library("ggplot2")

## True model
p <- c(2, 100, 50)
f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

sigma_g <- 0.1
pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

z_g <- function() rlnorm(1, 0, sigma_g)
set.seed(0)
Tobs <- 40

x <- numeric(Tobs)
x[1] <- 60
for(t in 1:(Tobs-1)) x[t+1] = z_g() * f(x[t], h=0)
obs <- data.frame(x = c(0, pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                  y = c(0, x[2:Tobs]))
xObs <- obs$x
yObs <- obs$y
xPred <- seq(0, 1.1 * max(xObs), length = 50)
qplot(seq_along(x), x) + geom_line()

## Now the GP estimation from NIMBLE. Let’s emphasize shorter length-scales with the prior to compare:

## having this work with output from `nimbleCode` might be more natural than `expression`
priors <- expression({
  rho ~ dgamma(1, 1)
  sigGP ~ dunif(0, 1e5)
  sigOE ~ dunif(0, 1e5)
})

fit <- gp_setup(xObs, yObs, xPred)

Cmcmc <- fit$Cmcmc 
Cpred <- fit$Cpred
Cmodel <- fit$Cmodel
system.time(Cmcmc$run(100000))
samples <- as.matrix(Cmcmc$mvSamples)
## basic sanity check
testthat::expect_identical(Cmodel$getNodeNames(topOnly = TRUE), colnames(samples))

## predict from GP model using posterior MCMC samples
system.time(Cpred$run(samples))

## Posteriors
samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
df <- tidyr::gather(samples)

ggplot(df) + 
  geom_density(aes(value)) + 
  facet_wrap(~key, scale='free')

## extract predictions: E and C
E <- Cpred$getE()
C <- Cpred$getC()

obs <- data.frame(x = xObs, y = yObs)
pred <- data.frame(x = xPred, y = E, ymin = E - sqrt(diag(C)), ymax = E + sqrt(diag(C)))

ggplot2::ggplot(pred) + 
  geom_ribbon(aes(x = x,y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
  geom_line(aes(x = x, y = y), size=1) + 
  geom_point(data = obs, aes(x,y)) +
  coord_cartesian(xlim = range(c(xObs, xPred)), ylim = range(c(yObs,E))) +
  theme_bw()

## Decision theory
states <- xPred # Vector of all possible states
actions <- states # Vector of actions: harvest
## Let’s consider a slight variation of the most trivial utility function:
## one which explicitly adds a cost to completely exhausting the stock
## (or reducing the stock by more than, say 95% in this case.)
## This should be somewhat similar to the impact of no discount rate.

## Utility function
discount = 0.99

#get_utility <- function(x,h) pmin(x,h)
#R <- outer(states, actions, get_utility)

R <- sapply(actions, function(h){
    sapply(states, function(x){
        if(h < x) h else - 1 * max(states)
    })
})

## Implementing policy
z <- function() rlnorm(1, meanlog = 0, sdlog = sigma_g)

simulate_policy <- function(states, actions, policy, f, z, s0, steps = 50, utility = function(s,a) NA, discount = 1){
    s <- numeric(steps)
    a <- numeric(steps)
    u <- numeric(steps)
    s[1] <- s0
    for(t in 1:(steps-1)){
        a[t] <- actions[policy[which.min(abs(states - s[t]))]]
        s[t+1] <- z() * f(s[t], a[t])
        u[t] <- utility(s[t], a[t]) * discount ^ t
    }
    ## Final action determined but not implemented
    a[steps] <- actions[policy[which.min(abs(states - s[t]))]]
    data.frame(time = 1:steps, state = s, action = a, utility = u)
}

## GP model
gp_matrix <- function(states, actions, E, C){
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    sigmas <- sqrt(diag(C))
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            nextpop <- E[k] - actions[i]
            if(nextpop <= 0) {
                transition[k, , i] <- c(1, rep(0, K - 1))
            } else {
                transition[k, , i] <- dnorm(states, nextpop, sigmas[i]) / sum(dnorm(states, nextpop, sigmas[i]))
            }
        }
    }
    transition
}

P_gp <- gp_matrix(states, actions, E, C)
mdp_check(P = P_gp, R = R)

gp <- mdp_value_iteration(P_gp, R, discount = discount, epsilon = 0.00001, max_iter = 5e3, V0 = numeric(length(states)))

plot(states, states - actions[gp$policy],  xlab="Population size", ylab="Escapement")

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")
## With this amount of data, the gp solution is too cautious, and avoids any exploitation.

## Simulate under the true model

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
## (Average utility is approximate here since it does not include penalty;
## since a function and not a matrix is requred by this function at this time.)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")

P <- transition_matrix(states, actions, f, pdfn)
mdp_check(P = P, R = R)

mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Note that the altered award structure has almost no effect on the optimal policy
## given the true model, other than to avoid harvesting directly to zero even when
## the stock cannot persist, due to the explicit penalty for doing so.

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, mdp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")





## Taken from Carl's lab notebook:
## MDPtoolbox Allen Model; 29 Jul 2015
## http://www.carlboettiger.info/2015/07/29/mdptoolbox-allen-model.html
## (Markov Decision Processes (MDP) toolbox)
library("MDPtoolbox", quietly = TRUE)
library("ggplot2", quietly = TRUE)
K <- 10 # state space limit
states <- 0:K # Vector of all possible states
actions <- states # Vector of actions: harvest

sigma_g = 0.1
p <- c(2, 15, 5)

f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

## Utility function
discount = 0.95
get_utility <- function(x,h) {
    pmin(x,h)
}

R <- outer(states, actions, get_utility)

transition_matrix <- function(states, actions, f, pdfn){
    ## Initialize
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            ## Calculate the transition state at the next step, given the 
            ## current state k and action i (harvest H[i])
            nextpop <- f(states[k], actions[i])
            ## Population always extinct if this is negative.
            ## since multiplicitive shock z_t * f(n) < 0 for all f(n) < 0
            if(nextpop <= 0)
                transition[k, , i] <- c(1, rep(0, length(states) - 1))
            ## Implement demographic stochasticity 
            else {
                ## Cts distributions need long-tailed denominator as normalizing factor:
                fine_states <- seq(min(states), 10 * max(states), by = states[2] - states[1])
                N <- sum(pdfn(fine_states, nextpop))  
                transition[k, , i] <-pdfn(states, nextpop) / N
                ## We need to correct this density for the final capping state ("Pile on boundary")
                ## (discrete or cts case)
                ## this can be a tiny but negative value due to floating-point errors.
                ## so we take max(v,0) to avoid
                transition[k, K, i] <- max(1 - sum(transition[k, -K, i]), 0)
            }
        } 
    }
    transition
}

P <- transition_matrix(states, actions, f, pdfn)
apply(P, c(1,3), sum)  ## double-check

## Using toolbox
mdp_check(P = P, R = R)
mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Compare to Reed
## From Reed (1979) we know that the optimal solution is a constant-escapement rule
## when the growth function in convex. Note that this condition is violated by the
## growth function with alternative stable states (Allen/Ricker-Allee model),
## resulting in a very different optimal policy:
## f^prime(s^star) = 1/discount
## For growth-rate function f, where discount is the discount factor and s∗ the stock size
## for the constant escapement. Analytic solutions are clearly possible for certain
## growth functions, but here I’ve just implemented a generic numerical solution.
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)
D <- actions[mdp$policy]
# The difference between Bellman and the analytical solution is small:
plot(states, states - D,  xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)




## Carl's examples using MDPtoolbox
## MDPtoolbox Ex 2; 14 Jul 2015
## (Markov Decision Processes (MDP) toolbox)
## http://www.carlboettiger.info/2015/07/14/mdptoolbox-ex-2.html
## Adapted from Marescot et al. appendix 5, to Reed optimal control problem, including direct comparison against (semi) analytic optimum.

## step 1: define objectives
## This is a conceptual step which does not require coding
## step 2: define states
##K <- 150 # state space limit
K <- 10 # state space limit
states <- 0:K # Vector of all possible states

## step 3: define control actions
## Vector of actions: harvest
H <- states

## step 4: define dynamic model (with demographic parameters)
##p <- c(6,0.05)
p <- c(3, 0.2)
f <- function(x, h){
  A <- p[1] 
  B <- p[2] 
  s <- pmax(x-h, 0)
  A * s/(1 + B * s)   #### Bellman equation
}
sigma_g = 0.1

## step 5: define utility
## Utility function
get_utility <- function(x,h) {
    pmin(x,h)
}

## step 6: solve bellman equation with value iteration
## Initialize transition matrix
transition <- array(0, dim = c(length(states), length(states), length(H)))
## Initialize utility matrix
utility <- array(0, dim = c(length(states), length(H)))
## Fill in the transition and utility matrix
## Loop on all states
for (k in 0:K) {
    ## Loop on all actions
    for (i in 1:length(H)) {
        ## Calculate the transition state at the next step, given the 
        ## current state k and the harvest H[i]
        nextpop <- f(k, H[i])
        if(nextpop <= 0)
            transition[k+1, , i] <- c(1, rep(0, length(states) - 1))
        ## Implement demographic stochasticity by drawing probability from a density function
        else {
            ## We need to correct this density for the final capping state ("Pile on boundary").
            ## For discrete probability distribution,
            ## this is easy if `states` includes all possible
            ## discrete states below the capping state.
            ## (e.g. all non-negative integers less than K).  
            ## For a continuous distribution, this is more problematic
            ## as we have to first normalize the densities.
            ## EDIT: this can be negative, due to floating-point errors.
            ## so we take max(v,0) to avoid.
            ## Get long-tailed denominator as normalizing factor (continuous distributions only):
            fine_states <- seq(min(states), 10 * max(states), by = states[2]-states[1])
            N <- sum(dlnorm(fine_states, log(nextpop), sdlog = sigma_g))
            transition[k+1, , i] <- dlnorm(states, log(nextpop), sdlog = sigma_g) / N
            ## We need to correct this density for the final capping state ("Pile on boundary")
            transition[k+1, K+1, i] <- max(1 - sum(transition[k+1, -(K+1), i]), 0)
        }
        ## Compute utility
        utility[k+1, i] <- get_utility(k, H[i])
    } # end of action loop
} # end of state loop

## Solution calculated explicitly:
## The backward iteration consists in storing action values in the vector
## Vt which is the maximum of utility plus the future action values for all
## possible next states. Knowing the final action values, we can then backwardly
## reset the next action value Vtplus to the new value Vt.  We start The backward
## iteration at time T-1 since we already defined the action value at Tmax.

## Discount factor
discount <- 0.95
## Action value vector at tmax
Vtmax <- numeric(length(states))
## Action value vector at t and t+1
Vt <- numeric(length(states))
Vtplus <- numeric(length(states))
## Optimal policy vector
D <- numeric(length(states))
## Time horizon
##Tmax <- 150
Tmax <- 4

for (t in (Tmax - 1):1) {
    ## We define a matrix Q that stores the updated action values for 
    ## all states (rows)
    ## actions (columns)
    Q <- array(0, dim = c(length(states), length(H)))
    for (i in 1:length(H)) {
        ## For each harvest rate we fill for all states values (row) 
        ## the ith column (Action) of matrix Q
        ## The utility of the ith action recorded for all states is 
        ## added to the product of the transition matrix of the ith 
        ## action by the action value of all states 
        Q[,i] <- utility[, i] + discount * (transition[,,i] %*% Vtplus)
    } # end of the harvest loop
    ## Find the optimal action value at time t is the maximum of Q
    Vt <- apply(Q, 1, max)
    ## After filling vector Vt of the action values at all states, we 
    ## update the vector Vt+1 to Vt and we go to the next step standing 
    ## for previous time t-1, since we iterate backward
    Vtplus <- Vt
} # end of the time loop

## Find optimal action for each state
for (k in 0:K) {
    ## We look for each state which column of Q corresponds to the 
    ## maximum of the last updated value 
    ## of Vt (the one at time t + 1). If the index vector is longer than 1 
    ## (if there is more than one optimal value we chose the minimum 
    ## harvest rate)
    D[k + 1] <- H[(min(which(Q[k + 1, ] == Vt[k + 1])))]
}

## plot solution
plot(states, states - D, xlab="Population size", ylab="Escapement")

## proof of optimality: compare with analytical solution
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)

## The difference between Bellman equation solution and the analytical solution is small:
plot(states, states - D, xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)

## Using MDPtoolbox
library('MDPtoolbox')
mdp_check(P = transition, R = utility)
out <- mdp_value_iteration(transition, utility, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = Vtmax)
plot(states, states - D, xlab="Population size", ylab="Escapement")
lines(states, states - H[out$policy], col="red", lty=2)




## testing using new NIMBLE option: verbose
nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleOptions('verbose')

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)


## testing the new includeEfficiency option for MCMCsuite
## calculates N, ESS, and Eff

library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
load('~/GitHub/legacy/autoBlock/data/model_ice.RData')

out <- MCMCsuite(code, constants, data, inits,
##                 MCMCs = c('nimble', 'nimble_slice', 'autoBlock'),
                 MCMCs = c('nimble', 'nimble_slice'),
                 makePlot = FALSE,
                 calculateEfficiency = TRUE
                 )

out$timing

out$summary

apply(out$summary[, 'efficiency', ], 1, min)
apply(out$summary[, 'efficiency', ], 1, mean)

out$efficiency


## running OpenBUGS installation from OSX (using wine)

##library(R2OpenBUGS)
library(R2WinBUGS)
library(BRugs)

## schools data in the R2OpenBUGS library
data(schools)
schools

## define the model
nummodel <- function() {
    for (j in 1:J) {
        y[j] ~ dnorm(theta[j], tau.y[j])
        theta[j] ~ dnorm(mu.theta, tau.theta)
        tau.y[j] <- pow(sigma.y[j], -2)
    }
    mu.theta ~ dnorm(0, 1E-6)
    tau.theta <- pow(sigma.theta, -2)
    sigma.theta ~ dunif(0, 1000)
}

## write the model code out to a file
write.model(nummodel, 'nummodel.txt')
model.file1 = paste(getwd(), 'nummodel.txt', sep='/')
## and let's take a look:
file.show('nummodel.txt')

## prepare the data for input into OpenBUGS
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list(J = J, y = y, sigma.y = sigma.y)

## initialization of variables
inits <- function() {
    list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100), sigma.theta = runif(1, 0, 100))
}

## set the WINE working directory and the directory to OpenBUGS
## change the OpenBUGS.exe location as necessary
WINE <- '/opt/local/bin/wine'
WINEPATH <- '/opt/local/bin/winepath'
##OpenBUGS.pgm <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe'
bugs.directory <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323'

## these are the parameters to save
parameters = c('theta', 'mu.theta', 'sigma.theta')

## run the model
schools.sim <- bugs(data, inits, model.file = model.file1, parameters=parameters,
                    n.chains = 3, n.iter = 1000,
                    ##OpenBUGS.pgm = OpenBUGS.pgm,   ## syntax for R2OpenBUGS::bugs
                    program = 'OpenBUGS',
                    bugs.directory = bugs.directory,
                    WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE)

## R will pause
## When model is complete a prompt will reappear
print(schools.sim)


## debugging commands in R
## n: next
## c: continue
## s: step into function call
## f: step out of current function call



library(nimble)
modelName <- 'rats'
rats <- readBUGSmodel(model = modelName, dir = getBUGSexampleDir(modelName), returnModelComponentsOnly = TRUE)

code <- rats$model
constants <- rats$data[c('N','T','x')]
data <- rats$data[c('Y')]
inits <- rats$inits

Rmodel <- nimbleModel(code, constants, data, inits, dimensions = rats$dims)

spec <- configureMCMC(Rmodel)
spec$getSamplers()




## creating dendogram figure of autoBlocking for JSM presentation
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
Rmcmc <- buildMCMC(Rmodel, autoBlock = TRUE, makePlots = TRUE)




## working on the conjugacy system, to infer dimensions of 'contributions'
library(nimble)
getDistributionsInfo('dnorm')



## let's manually run all NIMBLE testing
## (to find what's crashing, easier than using travis)
library(nimble)
library(testthat)
setwd('~/GitHub/nimble/nimble/packages/nimble/inst/tests/')
source('test_utils.R')

testFiles <- list.files()
print(testFiles)

##### this listing is out of date!
test_package('nimble', 'copy')      ## pass
test_package('nimble', 'dsl_dists') ## pass
test_package('nimble', 'math')      ## pass
test_package('nimble', 'mcmc')      ## errors/warnings copied below
test_package('nimble', 'meta')      ## pass
test_package('nimble', 'models')    ## 
test_package('nimble', 'trunc')     ## 
test_package('nimble', 'user')      ## 



## fully rebuild NIMBLE user manual

setwd('~/GitHub/nimble/nimble-docs/UserManual/')
library(knitr) 
knit2pdf('NimbleUserManual.Rnw') 
system('open NimbleUserManual.pdf')


## investigating the bug in R M-H sampler

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'RW')
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
niter <- 1000

set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples

d <- Rsamples - Csamples
all(d == 0)



Rmodel$nodes[['a']]$calculateDiff

Rmodel$logProb_a
Cmodel$logProb_a

calculateDiff(Rmodel, 'a')
calculateDiff(Cmodel, 'a')


nimble:::rCalcDiffNodes(Rmodel, 'a')
nimble:::rCalcDiffNodes(Cmodel, 'a')


## testing nimbleOptions() system
## and also the new nimble option: MCMCcontrolDefaultList

library(nimble)

getNimbleOption('MCMCcontrolDefaultList')

environment(sampler_conjugate_dnorm)$nfRefClassDef

getNimbleOption('verifyConjugatePosteriors')
nimble:::setNimbleOption('verifyConjugatePosteriors', TRUE)
getNimbleOption('verifyConjugatePosteriors')
environment(sampler_conjugate_dnorm)$nfRefClassDef
buildConjugateSamplerFunctions()
environment(sampler_conjugate_dnorm)$nfRefClassDef


code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 2)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b = .5, c = 4, d = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

lst <- getNimbleOption('MCMCcontrolDefaultList')
lst
lst$sliceWidth <- 99
lst$scale <- 9999
lst
nimbleOptions(MCMCcontrolDefaultList = lst)
getNimbleOption('MCMCcontrolDefaultList')

spec <- configureMCMC(Rmodel)
spec$getSamplers()


## complete example of using RW_llFunction sampler,
## for the User Manual

library(nimble)

code <- nimbleCode({
    p ~ dunif(0, 1)
    y ~ dbin(p, n)
})

Rmodel <- nimbleModel(code, data = list(y = 3),
                      inits = list(p = 0.5, n = 10))

llFun <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        y <- model$y
        p <- model$p
        n <- model$n
        ll <- lfactorial(n) - lfactorial(y) - lfactorial(n-y)
                 + y*log(p) + (n-y)*log(1-p)
        returnType(double())
        return(ll)
    }
)

RllFun <- llFun(Rmodel)

mcmcspec <- configureMCMC(Rmodel, nodes = NULL)

mcmcspec$addSampler(target = 'p', type = 'RW_llFunction',
    control = list(llFunction = RllFun, includesTarget = FALSE))

Rmcmc <- buildMCMC(spec)





## testing new nimble options system

library(nimble)
nimbleOptions()
ls(nimble:::.nimbleOptions)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b=.5, c=1, d=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## working on dipper models

rm(list=ls())
library(nimble)
trunc <- TRUE
load('~/GitHub/userDistMCMC/dipperData.RData')
## optionally truncate data:
last <- apply(y, 1, function(hist) max(which(hist==1)))
yDHMM <- 2 - y
onesMatrix = matrix(1,nind,k)  ## 'ones' matrix for survival
onesVector = rep(1,nind)       ## 'ones' vector for chi
if(trunc) { ind <- 1:3;   nind<-length(ind);   first<-first[ind];   last<-last[ind];   y<-y[ind,,drop=FALSE];   yDHMM<-yDHMM[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE];   onesMatrix<-onesMatrix[ind,,drop=FALSE];   onesVector<-onesVector[ind] }
code <- quote({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:nind) {
        for(t in (first[i] + 1):last[i]) {
            y[i,t] ~ dbin(p, 1)
            ##onesMatrix[i,t] ~ dbin(phi, 1)
        }
        ##onesVector[i] ~ dbin(chi[last[i]], 1)
    }
    ##chi[k] <- 1
    ##for(t in 1:(k-1)) {
    ##chi[k-t] <- (1-phi) + phi * (1-p) * chi[k-t+1]
    ##}
})
constants <- list(k=k, nind=nind, first=first, last=last)
data      <- list(y=y)##, onesMatrix=onesMatrix, onesVector=onesVector)
inits     <- list(phi=0.6, p=0.9)

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)


Rmodel <- nimbleModel(code, constants, data, inits)



## reproducible example of problem in new igraph package
## emailed this example to Gabor


remove.packages('igraph')
install.packages('igraph')

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
})
constants <- list()
data <- list()
inits <- list(a=0, b=0)

Rmodel <- nimbleModel(code, constants, data, inits)


library(igraph)
graph <- graph.empty()
graph <- add.vertices(graph, 2, name = c('a', 'b'))
graph <- add.edges(graph, c(1, 2))
toposortReturn <- topological.sort(graph, mode = 'out')
print(toposortReturn)
class(toposortReturn)



oldGraphID_2_newGraphID <- sort(newGraphID_2_oldGraphID, index = TRUE)$ix
graph <<- permute.vertices(graph, oldGraphID_2_newGraphID)  # re-label vertices in the graph



## writing stan model file for RSBS
##source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
##install_rstan()

load('~/GitHub/legacy/automated-blocking-examples/data/model_redblue.RData')
constantsAndData <- c(constants, data)

library(rstan)

stan_mod <- stan_model(file = '~/temp/redblue.stan')

stan_out <- sampling(stan_mod, data=constantsAndData, chains=1, iter=10000, thin=1, init=list(inits))

tempArray <- extract(stan_out, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
dimnames(tempArray)[[2]] <- gsub('_', '.', dimnames(tempArray)[[2]])
if(!all(monitorNodesBUGS %in% dimnames(tempArray)[[2]])) {
    missingNames <- setdiff(monitorNodesBUGS, dimnames(tempArray)[[2]])
    warning(paste0('Stan output is missing values for: ', paste0(missingNames,collapse=', ')))
}
samplesArray <- array(0, dim = c(nkeep, length(monitorNodesBUGS)))
dimnames(samplesArray)[[2]] <- monitorNodesBUGS
monitorsWeHave <- intersect(monitorNodesBUGS, dimnames(tempArray)[[2]])
samplesArray[, monitorsWeHave] <- tempArray[(burnin+1):floor(niter/thin), monitorsWeHave, drop=FALSE]
addToOutput('stan', samplesArray, timeResult)

## playing around to build, install, run, and test gpmanagement package

q('no')
library(devtools)
setwd('~/GitHub/forks/gpmanagement')
remove.packages('gpmanagement')
install_github('danielturek/gpmanagement')
##document('gpmanagement')
system('R CMD build gpmanagement')
system('R CMD INSTALL .')

library(nimble)
library(testthat)
library(gpmanagement)

check('.')  ## this gives me the errors 

gp_setup

test_package('gpmanagement')

source('~/GitHub/forks/gpmanagement/tests/testthat/test_gp.R')

##setwd('~/GitHub/forks/gpmanagement/tests/')
##test_check('gpmanagement')

nf <- gpmanagement:::nfd(3)
nf$run()

set.seed(0)
x <- 1:100
y <- sin(x/5) + rnorm(100, 0.1)
ind <- sort(sample(1:100, 40))
xObs <- x[ind]                ## input
yObs <- y[ind]                ## input
xPred <- c(x, 101:120)        ## input

fit <- gp_setup(xObs, yObs, xPred)


## testing of 3-D array arguments, and their automatic demotion to smaller dim arrays

library(nimble)

nf1 <- nimbleFunction(
    run = function(a = double(2), b = double(3)) {
        print('a11: ', a[1,1])
        print('b111: ', b[1,1,1])
        print(dim(a)[1])
        print(dim(a)[2])
        print(dim(b)[1])
        print(dim(b)[2])
        print(dim(b)[3])
    }
)

nf <- nimbleFunction(
    run = function() {
        a <- nimArray(0, 3, 3)
        declare(b, double(3, c(3,3,3)))
        for(i in 1:3)
            for(j in 1:3) {
                a[i,j] <- 10*i + j
                for(k in 1:3)
                    b[i,j,k] <- 100*i+10*j+k
            }
        print(a)
        ##print(b)
        nf1(a[2,1:2], b[,,])
    }
)

nf()
Cnf <- compileNimble(nf)
Cnf()


## testing distribution dDHMM on a multistate example

CdDHMM <- compileNimble(dDHMM)
CrDHMM <- compileNimble(rDHMM)

checkD <- function(ind) {
    y <- yDHMM
    print(y[ind,])
    x <- y[ind, first[ind]:k]
    length <- length(x)
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(dDHMM(x, length, pi, Z, T, condition, 0))
    print(CdDHMM(x, length, pi, Z, T, condition, 0))
    print(exp(dDHMM(x, length, pi, Z, T, condition, 1)))
    print(exp(CdDHMM(x, length, pi, Z, T, condition, 1)))
}
checkD(225)   ## 0.2484
checkD(202)   ## 0.017496
checkD(167)   ## 0.1246882

checkR <- function(length) {
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(rDHMM(1, length, pi, Z, T, condition))
    print(CrDHMM(1, length, pi, Z, T, condition))
}
checkR(1)
checkR(2)
checkR(3)
checkR(10)


Z <- array(c(.9,.05,.05,.1,.8,.1,0,1/2,1/2), c(3,3))
Z
T <- array(c(0,0,1,0,0,1,0,0,1), c(3,3))
T
pi <- c(1, 1, 1)
condition <- c(1, 1, 1)

y <- c(1,1)
length <- length(y)
dDHMM( y, length, pi, Z, T, condition, 0)
CdDHMM(y, length, pi, Z, T, condition, 0)


## testing new DSL functions:
## nimVector()
## nimArray()

library(nimble)

Rnf <- nimbleFunction(
    run = function(val = double(), len = double(), r = double(), c = double()) {
        onesVec <- nimVector(1, 5)
        print('onesVec:');   print(onesVec)
        zerosVec <- nimVector(0, 5)
        print('zerosVec:');   print(zerosVec)
        otherVec <- nimVector(val, len)
        print('otherVec:');   print(otherVec)
        sumVec <- zerosVec + onesVec + otherVec[1:5]
        print('sumVec:');   print(sumVec)
        arr <- nimArray(.4, r, c)
        print('arr:');   print(arr)
        print('sum of otherVec: ', sum(otherVec))
        print('sum of arrL: ', sum(arr))
    }
)

Cnf <- compileNimble(Rnf)

Rnf(10, 7, 3, 4)
Cnf(10, 7, 3, 4)

## trying nimVector() and nimArray() in BUGS code ?!?
## it works!!
code <- nimbleCode({
    a[1:4] <- nimVector(2, 4)
    phi ~ dnorm(0, 1)
    b[2:10] <- nimVector(10, 9)
    b[1] <- phi
})

Rmodel <- nimbleModel(code)
simulate(Rmodel, 'phi')
calculate(Rmodel)

Cmodel <- compileNimble(Rmodel)

model <- Rmodel
model <- Cmodel

model$a
model$phi
model$b

simulate(model)
calculate(model)

## re-creating NA's data example error
## 11 June 2015:
## Perry said he might get to this by July
## running tests/test-copy.R also gives this same error

library(nimble)

code <- nimbleCode({
    y[2:4] ~ dmnorm(mu3[1:3], cov = cov3[1:3, 1:3])
})
constants <- list(mu3=rep(0,3), cov3=diag(3))
data <- list(y = c(NA,0,0,0))

Rmodel <- nimbleModel(code, constants, data)

calculate(Rmodel)
simulate(Rmodel)

Cmodel <- compileNimble(Rmodel)
calculate(Rmodel)
simulate(Rmodel)


## same error from test-copy.R
library(testthat)
library(nimble)
test_package('nimble', 'copy')




## Zhenglei's example, trying MCMC suite on it

library(lme4)
library(ggplot2)
library(scales)
library(nimble)
library(coda)

I <- 25
J <- 5
beta <- c(-0.2, 0.1)
theta <- 1
sigma <-0.5

simfun <- function(substance=letters[1:10],experiment=factor(1:5),log10HQ=seq(-1,3,by=0.2),theta=1,sigma=0.1,beta=c(-0.2, 0.5)){
    expdat <- expand.grid(experiment = experiment,log10HQ=log10HQ,substance=substance)
    expdat$HQ <- 10^expdat$log10HQ
    expdat$obs <- factor(seq(nrow(expdat)))
    n <- nrow(expdat)
    nsub <- length(substance)
    re <- rnorm(nsub,0,theta)
    names(re) <- substance
    fixed <- beta[1]+beta[2]*expdat$log10HQ
    expdat$y <- re[expdat$substance]+rnorm(n,0,sigma)+fixed
    expdat$effect <- gtools::inv.logit(expdat$y)
    return(expdat)
}

expdat <- simfun(substance=factor(paste0('S',1:I)),experiment = factor(1:J),log10HQ=seq(-1,3,by=0.2),theta=0.5,beta=c(-0.2,0.5),sigma=0.3)

code <- nimbleCode({
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm( 0, 0.001 )
    tau.y ~ dgamma(0.001, 0.001)  
    tau.re ~ dgamma(0.001, 0.001)
    for(j in 1:M){
        ranef.v[j]~dnorm(0, tau.re)
    }
    for (i in 1:N){
        y[i] ~ dnorm ( mu.y[i], tau.y )
        mu.y[i] <- alpha+beta*x[i] + ranef.v[unit[i]]
    }
})

constants <- list(N = nrow(expdat),M=nlevels(expdat$substance),x = expdat$log10HQ,unit=as.numeric(expdat$substance))
##data <- data.frame(y = expdat$y)
data <- list(y = expdat$y)
inits <- list(alpha = 0, beta = 0.5, tau.y=1,tau.re=1)

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'jags', 'noConj'),
    niter = 20000,
    makePlot = TRUE,
    savePlot = TRUE,
    summaryStats = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
)

out$summary[, 'effectiveSize', ] / out$timing[c('nimble', 'jags', 'noConj')]


## helping Nick get started with using node$get_mean() function for particle filter

library(nimble)

code <- nimbleCode({
    x ~ dnorm(1, 1)
    y ~ dnorm(x, 1)
})

Rmodel <- nimbleModel(code)
simulate(Rmodel)
Rmodel$x
Rmodel$y
Rmodel$nodeFunctions[['x']]$get_mean()
Rmodel$nodeFunctions[['y']]$get_mean()

## only necessary for the 'higherLevel' example
virtual_NF <- nimbleFunctionVirtual(
    ## don't need to supply prototpye for run(),
    ## since the default is no args, and returnType void
    methods = list(  ## definitely need prototype for return_mean() function
        return_mean = function() {
            returnType(double())
        }
    )
)

nfDef <- nimbleFunction(
    contains = virtual_NF,  ## only necessary for the 'higherLevel' example
    setup = function(model, node) {
        nfList <- nimbleFunctionList(node_stoch_dnorm)
        nfList[[1]] <- model$nodeFunctions[[node]]
    },
    run = function() {
        print('node value is: ', model[[node]])
        print('(assumed to be dnorm) node mean is: ', nfList[[1]]$get_mean())
    },
    methods = list(                        ## only necessary for the 'higherLevel' example
        return_mean = function() {         ##
            returnType(double())           ##
            return(nfList[[1]]$get_mean()) ##
        }                                  ##
    )                                      ##
)

nfX <- nfDef(Rmodel, 'x')
nfY <- nfDef(Rmodel, 'y')
nfX$run()
nfY$run()

compiledList <- compileNimble(list(Rmodel, nfX, nfY))
Cmodel <- compiledList[[1]]
CnfX <- compiledList[[2]]
CnfY <- compiledList[[3]]

CnfX$run()
CnfY$run()

CnfX$return_mean()
CnfY$return_mean()

## this is the 'higherLevel' example
nfDefHigherLevel <- nimbleFunction(
    setup = function(model, nodes) {
        nfList <- nimbleFunctionList(virtual_NF)
        for(i in seq_along(nodes))
            nfList[[i]] <- nfDef(model, nodes[i])
    },
    run = function() {
        declare(meansVector, double(1, length(nfList)))  ## declares a vector of doubles, with length equal to the length of nfList
        for(i in seq_along(nfList)) {
            nfList[[i]]$run()
            meansVector[i] <- nfList[[i]]$return_mean()
        }
        print('vector of the means is: ', meansVector)
    }
)

nfHigherLevel <- nfDefHigherLevel(Rmodel, c('x', 'y'))

nfHigherLevel$run()

CnfHigherLevel <- compileNimble(nfHigherLevel, project = Rmodel)

CnfHigherLevel$run()


## fixing values(model, nodes) <- value reference class '<<-' warning


library(nimble)

nfdef <- nimbleFunction(
    setup = function(model) {},
    run = function(val = double(1)) {
        values(model, 'xxx') <<- val
    }
)

code <- nimbleCode({
    xxx ~ dnorm(0, 1)
})

Rmodel <- nimbleModel(code)

nf <- nfdef(Rmodel)

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(nf, project = Rmodel)

Rmodel$xxx
nf$run(3)
Rmodel$xxx

Cmodel$xxx
Cnf$run(4)
Cmodel$xxx


## testing / trying to find some error in conjugacy system(??)  May 2015

library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')

## run entire MCMC testing system
source('~/GitHub/nimble/packages/nimble/inst/tests/test-mcmc.R')

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)
code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list(mu0 = mu0, Q0 = Q0, Q = Q, B = B, a=a)
data <- list(y=y)
m <- nimbleModel(code, constants = constants, data=data)

m <- nimbleModel(code, constants = data)

spec <- configureMCMC(m)
Rmcmc <- buildMCMC(spec)
spec$getSamplers()

options(error = recover)

test_mcmc(model = code, data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))


## testing distribution of particle filter likelihood estimates

ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))

plot(muvec, varll)
varest <- 1/2*sigx^4/sigy^4 + muvec^2*sigx^2/sigy^4
lines(muvec, varest, col='red')



## library(nimble)
## Rpf <- buildPF(Rmodel, 'x')
## Cmodel <- compileNimble(Rmodel)
## Cpf <- compileNimble(Rpf, project = Rmodel)
## Rpf$run(10)
## Cpf$run(10)
## ll <- numeric()
## for(i in seq_along(muvec)) {
##     mu <- muvec[i]
##     Cmodel$mu <- mu
##     ll[i] <- Cpf$run(10000)
## }
## plot(muvec, ll)


library(nimble)
source('~/GitHub/pfLL/pfLL.R')
m <- 1000
rep <- 1000
n <- m*rep
muvec <- seq(-5, 5, by=0.2)
sigx <- 1
sigy <- 1
taux <- 1/sigx^2; tauy <- 1/sigy^2
y <- 0   ## observation
constants <- list()
inits <- list(mu = y, sigx = sigx, sigy = sigy)
data <- list(y = y)
code <- quote({
    mu ~ dnorm(0, 0.00001)  ## these should not be necessary
    sigx ~ dunif(0, 1000)  ## these should not be necessary
    sigy ~ dunif(0, 1000)  ## these should not be necessary
    x ~ dnorm(mu, sigx)
    y ~ dnorm(x, sigy)
})
Rmodel <- nimbleModel(code, constants, data, inits)

out <- pfLL(Rmodel, 'x', param = data.frame(mu=muvec), m=m, rep=rep, makePlot=FALSE)

ll <- apply(out$ll, 1, mean)
##ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))
plot(muvec, ll)
ll.pred <- -1/2 * muvec^2 / (sigx^2+sigy^2) + log(1/sqrt(2*pi*(sigx^2+sigy^2)))
lines(muvec, ll.pred, col='red')


varll <- apply(out$ll, 1, var)
##varll <- varll - min(varll)  ### TEMPORARY
plot(muvec, varll)
var.pred <- exp(-8.8 + 1.2592*abs(muvec))
lines(muvec, var.pred, col='red')

1



n <- 100000
sigx <- 1
muvec <- seq(-3, 3, by=0.1)
nmu <- length(muvec)
X <- array(NA, c(n, nmu))
for(i in 1:nmu) X[,i] <- rnorm(n, muvec[i], sigx)

Y <- X^2
var <- apply(Y, 2, var)
pred <- muvec^2 + sigx^2
plot(muvec, var, type='p')
lines(muvec, pred, col='red')



## figuring out 3D plotting in R
rm(list=ls())
df <- expand.grid(list(x=1:10, y=1:10))
rep <- 1
ll <- array(NA, c(dim(df)[1], rep))
for(i in 1:dim(df)[1]) {
    for(j in 1:rep) {
        x <- df[i, 'x']
        y <- df[i, 'y']
        dist <- sqrt((x-4)^2 + (y-7)^2)
        ll[i,j] <- rnorm(1, 100-dist^2, sd = dist/10)
    }
}

library(plot3D)
?scatter3D

x <- rep(df$x, rep)
y <- rep(df$y, rep)
z <- as.numeric(ll)

scatter3D(x, y, z)


## minimally reproducible example of getDependencies error for Perry
library(nimble)
code <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
    y[1] ~ dnorm(x[1], 1)
    y[2] ~ dnorm(x[2], 1)
})

Rmodel <- nimbleModel(code)

Rmodel$getDependencies('x[1]')   ## should be x[1], y[1]
## [1] 'x[1]' 'y[1]' 'y[2]'

Rmodel$getDependencies('x[2]')   ## should be x[2], y[2]
## [1] 'x[2]'



## trying to figure out conjugate sampling in 'equiv' BUGS example model
library(nimble)
code <- nimbleCode({
    tau[1] ~ dgamma(0.001, 0.001)
    tau[2] ~ dgamma(0.001, 0.001)
    pi ~ dnorm(0, 1.0E-06)
    phi ~ dnorm(0, 1.0E-06)
    mu ~ dnorm(0, 1.0E-06)
    for(i in 1:N) { 
        d[i] ~ dnorm(0,tau[2])  # Subject random effect
        for (k in 1:2) {
            Treat[i,k] <- group[i]*(k-1.5) + 1.5  # treatment given
            m[i,k] <- mu + pow(-1, Treat[i,k]-1)* phi /2 + pow(-1, k-1)* pi /2 + d[i]
            Y[i,k] ~ dnorm(m[i,k], tau[1])
        }
    }
    ##theta <- exp(phi)
    ##equivalence <- step(theta - 0.8) - step(theta - 1.2)
})
Y <- structure(c(1.4, 1.64, 1.44, 1.36, 1.65, 1.08, 1.09, 1.25, 1.25, 1.3, 1.65, 1.57, 1.58, 1.68, 1.69, 1.31, 1.43, 1.44, 1.39, 1.52), .Dim = c(10, 2))
group <- c(1, 1, -1, -1, -1, 1, 1, 1, -1, -1)
N <- 10
##T <- 2
constants <- list(N=N, group=group)
data <- list(Y=Y)
inits <- list(mu=0, phi=0, pi=0, tau=c(1,1))
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 10

set.seed(0); Rmcmc$run(niter)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0); Cmcmc$run(niter)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples <- Csamples[, dimnames(Rsamples)[[2]]]
Csamples

Rsamples-Csamples


## test equiv model
library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')
test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)


1

## trying Kalman Filter and calculating likelihood for a linear SSM

# Say state process is
# X(t+1) = r * X(t) + b + nu(t)
# nu(t) ~ dnorm(0, sig2nu)
# 
# Obs process is
# Y(t) = X(t) + eps(t)
# eps(t) ~ dnorm(0, sig2eps)
# 
# Time series goes from t = 1:T
# 
# This is the simplest possible model for us.
# 
# Direct mvn approach:
#   1. Stationary variance of states:
#   V[x] = r^2 V[x] + sig2nu
# V[x] = sig2nu / (1-r^2)
# 
# 2. Deriving entries in covariance matrix of Y (from 1:T)
# Cov[X(t), X(t+1)] = r V[x]
# Cov[X(t), X(t+d)] = r^d V[x]
# Therefore Cov[Y(t), Y(t+d)] = r^d V[x]
# 
# And on the diagonal:
#   V[y] = V[x] + sig2eps
# 
# And mu(t) = E[Y] = E[X] = b/(1-r) for all t
# 
# So with that mean vector and covariance matrix in hand 
# you can directly calculate mvnorm(Y, mu, Cov(Y))
# 
# Kalman filter approach and direct LL are shown in code below.  
# Sorry  I didn't think of this.  I had this from the early 
# explorations that led to Knape, Besbeas and de Valpine in 
# Ecology last year (or was it this year?).  
# 
# In this setup, you want numSites = 1 (we were considering time-series 
# with data from multiple sites).  The LL functions use the meanData 
# (mean over all the sites, which we were comparing to a more complete 
# site-specific model), which for you will be the same as the fullData.



simData <- function(mu=10, a=0.8, sigPN=0.2, sigOE=0.4, t=20) {
    b <- mu*(1-a)
    X <- mu
    for(i in 1:50) X <- a*X + b + rnorm(1,0,sigPN) ## get to stationary
    states <- numeric(t)
    states[1] <- X
    for(i in 2:t) states[i] <- a*states[i-1]+b + rnorm(1,0,sigPN)
    obsStates <- rnorm(t, states, sigOE)
    list(x = states, y = obsStates, mu=mu, a=a, b=b, sigOE=sigOE, sigPN=sigPN, t=t)
}

mvn_ll <- function(d) {
    require(mvtnorm)
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_vec <- rep( b/(1-a), t)
    cov_mat <- matrix(nrow = t, ncol = t)
    var_x <- sigPN2 / (1-a^2)
    var_y <- var_x + sigOE2
    diag(cov_mat) <- var_y
    for(i in 2:t) {
        for(j in 1:(i-1)) {
            cov_mat[i,j] <- (a^(i-j))*var_x
            cov_mat[j,i] <- cov_mat[i,j]
        }
    }
    dmvnorm(y, mu_vec, cov_mat, log = TRUE)
}

KF_ll <- function(d) {
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_x <- b/(1-a)
    var_x <- sigPN2 / (1-a^2)
    cov_xy <- var_x
    var_y <- var_x + sigOE2
    ll <- dnorm(mu_x, y[1], sqrt(var_y), log = TRUE)
    mu_x <- mu_x + (cov_xy / var_y) * (y[1] - mu_x)
    var_x <- var_x - cov_xy*cov_xy / var_y
    for(i in 2:t) {
        mu_x <- a*mu_x + b
        var_x <- a^2 * var_x + sigPN2
        if(!is.na(y[i])) {
            cov_xy <- var_x
            var_y <- var_x + sigOE2
            ll <- ll + dnorm(mu_x, y[i], sqrt(var_y), log = TRUE)
            mu_x <- mu_x + (cov_xy / var_y) * (y[i] - mu_x)
            var_x <- var_x - cov_xy*cov_xy / var_y
        }
    }
    ll
}

code <- quote({
    x[1] ~ dnorm(mu, sd = sqrt((sigPN^2) / (1 - a^2)))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

d <- simData()

library(nimble)
constants <- list(a=d$a, b=d$b, t=d$t, sigOE=d$sigOE, sigPN=d$sigPN, mu=d$mu)
data <- list(y = d$y)
inits <- list(x = d$x)
Rmodel <- nimbleModel(code = code, constants = constants, data=data, inits=inits)
calculate(Rmodel)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

Cpf$run(100000)

mvn_ll(d)

KF_ll(d)



## exploring how truncation works in our system

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ T(dnorm(a, 1), 1.5, 1.6)
    c ~ T(dnorm(b, 1), 3, Inf)
    d ~ dnorm(c, 1)
    e ~ dnorm(d, 1)
})
constants <- list()
data <- list()
inits <- list()

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)
Rmodel <- md$newModel(data=data, inits=inits)

Rmodel$modelDef$printDI()

node <- 'c'
Rmodel$nodes[[node]]$calculate
Rmodel$nodes[[node]]$simulate
Rmodel$nodes[[node]]$getLogProb

lapply(Rmodel$modelDef$declInfo, function(di) di$truncation)
for(node in c('a', 'b', 'c')) print(Rmodel$isTruncated(node))
for(node in c('a', 'b', 'c')) print(Rmodel$getBounds(node))

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## using rankSample

library(nimble)
##Cmodel <- compileNimble(Rmodel = nimbleModel(quote({a ~ dnorm(0, 1)})))
nfDef <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer(), silent = logical()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp, silent)
        print(samp)
    }
)
nf <- nfDef()
Cnf <- compileNimble(nf)

nfDefOmit <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp)
        print(samp)
    }
)
nfOmit <- nfDefOmit()
CnfOmit <- compileNimble(nfOmit)
x <- 1:10

silent <- FALSE
silent <- TRUE

set.seed(1);   .Call('rankSample', c(.5,.4,.1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,.1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,.1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,.1), 10L,    silent)

set.seed(1);   .Call('rankSample', c( 0, 0, 0), 10L, x, silent)
set.seed(1);            rankSample(c( 0, 0, 0), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c( 0, 0, 0), 10L,    silent)
set.seed(1);               Cnf$run(c( 0, 0, 0), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,-1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,-1), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x);            print(x)
set.seed(1);            nfOmit$run(c(.5,.4,-1), 10L)
set.seed(1);           CnfOmit$run(c(.5,.4,-1), 10L)



## difference in rank sample
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        declare(wts, double(1, 1))
        wts[1] <- 1
        m <- 1L
        declare(samp, integer(1, 1))
        rand <- rnorm(1, 0, 1)
        print('random number before rankSample: ', rand)
        rankSample(wts, m, samp)
        rand <- rnorm(1, 0, 1)
        print('random number after rankSample: ', rand)
    }
)

nf <- nfDef()
Cnf <- compileNimble(nf)
set.seed(0)
nf$run()
set.seed(0)
Cnf$run()



## testing buildPF() particle filter (pf) algorithm
library(nimble)
code <- quote({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
t <- 5
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
Rmodel$mu <- 1/(1-.95)
Rmodel$b <- 1
Rmodel$sigPN <- .2
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)
rm(Rmodel)
Rmodel <- nimbleModel(code, constants, data, inits)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

m <- 100
set.seed(0);    pf$run(m)
set.seed(0);   Cpf$run(m)

1
2
3



## making MCMCsuite work for Ryan using Stan MCMC

allModels <- c('blocker', 'bones', 'dyes', 'line', 'pump', 'rats')
library(nimble)
library(coda)
allModels <- 'blocker'
mcmcs=c('nimble','jags','noConj','stan')
x=list()
for (i in 1:length(allModels)){
    x[[i]]<-readBUGSmodel(model=allModels[i],
                          dir=getBUGSexampleDir(allModels[i]),
                          returnModelComponentsOnly=TRUE)
}
i <- 1
suite_output <- MCMCsuite(
    x[[i]]$model,
    constants = x[[i]]$data,
    inits = x[[i]]$inits,
    MCMCs = mcmcs,
    makePlot=F,
    savePlot=F,
    summaryStats=c('mean','median','sd','CI95_low','CI95_upp','effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) })),
    stan_model = paste('~/temp/', allModels[i], '.stan', sep='')
)


## learning about R's memory management and storage of objects
gcinfo(TRUE)
gc()

library(pryr)
?object_size
object_size(1:100)
object_size(mtcars)

x <- 0:50
s <- sapply(x, function(i) object_size(1:i))
plot(x, s, type='s', ylim=c(0,300))

obj <- complex()
object_size(obj)

mem_used
mem_used()

pryr:::node_size()
pryr:::show_bytes

mem_change(NA)

bytes(1)
bytes(1L)
bytes('ab')

address(2)
x <- 1:10
address(x)
refs(x)
y <- x
refs(x)
z <- x
x[1] <- 16L

tracemem(x)
x[1] <- .1

x <- 15
address(x)

system(paste0('~/temp/t ', address(x)))



## testing of mcmcplots package function mcmcplot()

rm(list=ls())
library(nimble)
model <- 'SSMcorrelated'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)

out <- MCMCsuite(code, constants, data, inits, monitors=c('a','b','p'), MCMCs=c('nimble','nimble_RW', 'nimble_slice', 'jags'), makePlot=FALSE, niter=10000)

library(coda)
samples <- out$samples
dim(samples)
dimnames(samples)
mcmcList <- list()
for(i in 1:dim(samples)[1])     mcmcList[[i]] <- coda::mcmc(t(samples[i,,]))
names(mcmcList) <- dimnames(samples)[[1]]
codaList <- coda::as.mcmc.list(mcmcList)

library(mcmcplots)
mcmcplot(codaList)
traplot(codaList)




## testing conjugacy in newNimbleModel branch

library(nimble)
code <- nimbleCode({
    for(i in 1:10) {
        a[i] ~ dnorm(0, 1)
        b[i] ~ dgamma(1, 1)
    }
    for(i in 1:5)   a_temp[i] <- 3*a[i]
    for(i in 6:10)  a_temp[i] <- a[i]^2
    for(i in 1:10)  x[i] ~ dnorm(a_temp[i], 1)
    for(i in 1:10) {
        y1[i] ~ dpois(2*b[i])
        y2[i] ~ dgamma(1, rate = b[i])
    }
})
Rmodel <- nimbleModel(code)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
for(i in 1:200) spec$addSampler('RW', 'x[2]', print=F)
spec$addSampler('RW', 'x[2]')

spec$setSamplers(c(1:40, 1:40))
spec$setSamplers(c('a', 'y2'))
spec$removeSamplers('a')
spec$setSamplers()


## NEED TO RE-RUN THIS TEST, ONCE PERRY FIXES THE MODELVALUES COPYING ISSUE
## Dirichlet-multinomial conjugacy
## single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- nimbleCode({
    y[1:K] ~ dmulti(p[1:K], n)
    p[1:K] ~ ddirch(alpha[1:K])
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001)
    }
})
constants <- list(n = n, K = K)
data <- list(y = y)
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
Rmodel <- nimbleModel(code, constants=constants, data=data, inits=inits)

mv <- modelValues(Rmodel)
Rmodel$y
mv$y
nimCopy(from=Rmodel, to=mv, nodes = 'y', row = 1)
mv$y


spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

set.seed(0)
Rmcmc$run(10)

conjugate posterior density appears to be wrong, off by Inf
conjugate posterior density appears to be wrong, off by Inf
There were 50 or more warnings (use warnings() to see the first 50)

test_mcmc(model = code, data= data, seed = 0, numItsC = 10000,
          inits = inits,
          results = list(mean = list(p = p)),
          resultsTolerance = list(mean = list(p = rep(.06, K))))




## just a good old fashioned test that it's working
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})

constants <- list()
data <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$checkConjugacy()

spec <- configureMCMC(Rmodel, control = list(scale=0.01), monitors = c('x', 'c'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
set.seed(0)
Rmcmc$run(10)
samples <- as.matrix(Rmcmc$mvSamples)
samples


sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977),
  c = c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))

test_mcmc(model = code, exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01))




## testing if a Ref Class can just return a smaller object (not the whole Ref Class object)
RC <- setRefClass(
    Class = 'RC',
    fields = list(a='ANY'),
    methods = list(
        initialize = function(a) {
            a <<- a
            return(a)
        }
    )
)

rc <- RC(1)
class(rc)
rc$a



## bunch of random stuff

fileToList <- function(file) {
    env <- new.env()
    source(fileName, local = env)
    lst <- list()
    for(name in ls(env))   lst[[name]] <- get(name, env)
    return(lst)
}

library(nimble)

NF <- nimbleFunction(
    setup = function() { },
    run = function(a = double(0)) {
        returnType(double(0))
        return(a)
    }
)

myNF <- NF()
myNF_C <- compileNimble(myNF)

> myNF$run
## function (a) 
##     return(a)

myNF_C$run
## function (a) 
## {
##     ans <- .Call('CALL_nfRefClass60_operator_', a, .basePtr)
##     ans <- ans[[2]]
##     ans
## }


### testing the new autoBlock
rm(list=ls())
library(nimble)
model <- 'litters'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$getNodeNames()

ab <- autoBlock(Rmodel, run = runList)
ab
ab$spec$getSamplers()
Rmcmc <- buildMCMC(ab$spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)


spec <- configureMCMC(Rmodel, autoBlock=TRUE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

Rmcmc <- buildMCMC(Rmodel, autoBlock=TRUE)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)







## March 2015
## timing comparisons for red-state-blue-state
library(nimble)
library(rjags)
rm(list=ls())
load('~/GitHub/autoBlock/data/model_redblue.RData')

## this version of red-state-blue-state code is 'jags friendly'
## reparametrized to use Precision matrix, replaced expit with logit
code <- nimbleCode({
    for (i in 1:2) {
        for (j in 1:2) {
            gamma[i, j] ~ dnorm(0, 1e-04)
        }
    }
    sigmaIntercept ~ dunif(0, 100)
    sigmaSlope ~ dunif(0, 100)
    rho ~ dunif(-1, 1)
    Precision[1, 1] <- 1/(1-rho^2) * 1/(sigmaIntercept^2)
    Precision[2, 2] <- 1/(1-rho^2) * 1/(sigmaSlope^2)
    Precision[1, 2] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    Precision[2, 1] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    for (i in 1:Nstates) {
        stateBetaMeans[i, 1] <- gamma[1, 1] + gamma[1, 2] * stateIncome[i]
        stateBetaMeans[i, 2] <- gamma[2, 1] + gamma[2, 2] * stateIncome[i]
        stateBetas[i, 1:2] ~ dmnorm(stateBetaMeans[i, 1:2], Precision[1:2, 1:2])
    }
    for (i in 1:N) {
        logit(p[i]) <- stateBetas[state[i], 1] + stateBetas[state[i], 2] * income[i]
        y[i] ~ dbern(p[i])
    }
})

niter <- 50000
monitorVars <- c('sigmaIntercept', 'sigmaSlope', 'rho', 'gamma')
constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)


system.time(modelDef <- nimbleModel(code=code, constants=constants, data=data, inits=inits, returnDef=TRUE))[3]/60
26

system.time(Rmodel <- modelDef$newModel(data=data, inits=inits))[3]/60
5

system.time(spec <- configureMCMC(Rmodel))
6 seconds

system.time(Rmcmc <- buildMCMC(spec))[3]
15 seconds

system.time(Cmodel <- compileNimble(Rmodel))[3]/60
11

system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))[3]/60
1

system.time(Cmcmc$run(niter))[3]/60
11 minutes for 50,000 iterations

system.time(jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE))[3]/60
10 seconds

system.time(jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1))[3]/60
2.001833





## displaying calc and sim functions
code <- nimbleCode({ x ~ dnorm(3, sd = 5) })
model <- nimbleModel(code)
model$nodes$x$simulate
## function()
##     model$x <<- rnorm(1, mean = 3, sd = 5)
model$nodes$x$calculate
## function() {
##     model$logProb_x <<- dnorm(model$x, mean = 3, sd = 5, log = 1)
##     return(invisible(model$logProb_x))
## }
model$nodes$x$getLogProb
## function() 
##     return(model$logProb_x)


### March 2015
### trying to figure out new MCMCspec

library(nimble)


#These first few lines are nothing new
myModelCode <- nimbleCode({
	a[1] ~ dnorm(0,1)
	a[2] ~ dexp(a[1]^2)
	a[3] ~ dnorm(0,1)
})
myModel <- nimbleModel(myModelCode)
mySpec1 <- configureMCMC(myModel)
myMCMC1 <- buildMCMC(mySpec1)
cmod <- compileNimble(myModel)
cmcmc1 <- compileNimble(myMCMC1, project = myModel)


# Rebuilding the way we are supposed to: 
# building a new spec with configureMCMC 
# adding a sampler which has already been compiled 
# not altering the monitors
mySpec2 <- configureMCMC(oldSpec = mySpec1)
mySpec2$addSampler('RW', control = list(targetNode = 'a[3]'))
rmcmc2 <- buildMCMC(mySpec2)
cmcmc2 <- compileNimble(rmcmc2, project = myModel)




#### March 2015 working through a full compileNimble call

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 1)
Rmodel <- nimbleModel(code, inits = inits)

debug(compileNimble)
Cmodel <- compileNimble(Rmodel)

debug(project$compileModel)
debug(modelCpp$buildAll)
debug(buildNodes)
debug(nimbleProject$compileNimbleFunctionMulti)
debug(compileNimbleFunction)
debug(buildNimbleFunctionCompilationInfo)


#### March 2015 testing of RW_llFunction or Carl B, submitted as test to testing suite

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 0)
Rmodel <- nimbleModel(code=code, inits=inits)

llFunction <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        ll <- dnorm(1, model$a, 1, log=1)
        returnType(double())
        return(ll)
    }
)
 
myLL <- llFunction(Rmodel)

spec <- configureMCMC(Rmodel, nodes=NULL)
spec$addSampler('RW_llFunction', list(targetNode='a', llFunction=myLL, includesTarget=FALSE))
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
unlist(Rmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

set.seed(0)
Cmcmc$run(10)
unlist(Cmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
## should be within some tolerance of 0.5


###### Feb 2015 example for Cal Poly interview talk

library(nimble)

code <- nimbleCode({
    for( i in 1 : N ) {
        b[i] ~ dnorm(mu, tau)
        r[i] ~ dbin(p[i], n[i])
        logit(p[i]) <- b[i]
    }
    pop.mean <- 1 / (1 + exp(-mu))
    mu ~ dnorm(0, 0.001)
    sigma <- 1 / sqrt(tau)
    tau ~ dgamma(0.001, 0.001)
})

constants <- list(N = 10, n = c(12,20,18,19,10,23,19,15,8,10))
data <- list(r = c(4,7,10,3,6,7,11,9,3,8))
inits <- list(mu = 0, tau = 1)

model <- nimbleModel(code=code, constants=constants, data=data, inits=inits)

spec <- configureMCMC(model)
spec$getSamplers()
spec$addSampler(‘slice’, list(targetNode = ‘mu’))
spec$getSamplers()
mcmc <- buildMCMC(spec)

Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model)

Cmcmc$run(10000)

samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)





## testing new addition to NIMBLE: conf$addSampler('node', 'conjugate')
library(nimble)

code <- nimbleCode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)

conf$printSamplers()
##[1]  conjugate_dgamma_dnorm_dpois_dgamma sampler: x,  dep_dnorm: a, a2,  dep_dpois: b, b2,  dep_dgamma: c
##[2]  conjugate_dnorm_dnorm_dlnorm sampler: a,  dep_dnorm: jNorm[2], jNorm[3],  dep_dlnorm: kLogNorm[2], kLogNorm[3]
##[3]  posterior_predictive sampler: a2
##[4]  posterior_predictive sampler: b
##[5]  posterior_predictive sampler: b2
##[6]  RW sampler: c
##[7]  posterior_predictive sampler: kLogNorm[2]
##[8]  posterior_predictive sampler: kLogNorm[3]
##[9]  posterior_predictive sampler: jNorm[2]
##[10] posterior_predictive sampler: jNorm[3]

pr <- TRUE
pr <- FALSE
nd <- 'x'
nd <- 'a'
nd <- 'c'
nd <- 'a2'
nd <- 'kLogNorm[3]'
conf$addSampler(nd, 'RW',        print=pr)
conf$addSampler(nd, 'slice',     print=pr)
conf$addSampler(nd, 'conjugate', print=pr)

conf$printSamplers()

conf$addSampler(nd, print=TRUE)
conf$addSampler(nd, print=TRUE)
conf$addSampler(nd, print=TRUE)





par(mfrow=c(1,3))
ns <- c(50, 500, 5000)
for(n in ns) {
    x <- rbinom(n=10000, size=n, prob=0.8) / n
    hist(x, xlim=c(0,1))
}



library(nimble)
df <- read.csv('~/Downloads/UsedCars.csv')
code <- nimbleCode({
    b0  ~ dnorm(0, sd=10000)
    bage ~ dnorm(0, sd=10000)
    bhp ~ dnorm(0, sd=10000)
    btype ~ dnorm(0, sd=10000)
    sigma ~ dunif(0, 50000)
    for(i in 1:N) {
        y[i] ~ dnorm(mu[i], sd = sigma)
        mu[i] <- b0 + bage*age[i] + bhp*hp[i] + btype*type[i]
    }
})
constants <- list(N = dim(df)[1], age=df$Age, hp=df$HP, type=df$Type)
data <- list(y = df$Price)
inits <- list(b0=0, bage=0, bhp=0, btype=0, sigma=1)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 20000, nburnin=10000)
dim(samples)

samplesPlot(samples)
samplesPlot(samples, 'btype')

apply(samples, 2, effectiveSize)







1

a <- acf(samples)
a
a$acf
plot(a)
acfplot(as.mcmc(samples), thin=10)

?acf
acfplot)

colnames(samples)
quantile(samples[, 'bage'], c(0.025, 0.975))

dim(samples)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var='btype', burnin=2000)
samplesPlot(samples, var='btype', ind=1:200)
samplesPlot(samples, var='btype', ind=2001:10000)
samplesPlot(samples, var=c('btype', 'bage'), ind=2001:10000)
samplesPlot(samples, var=c('btype', 'bage'), ind=2001:10000, densityplots=FALSE)

samplesPlot(samples, burnin=2000)
samplesPlot(samples, var=1:3, burnin=2000)
samplesPlot(samples, var=4, burnin=5000)
samplesPlot(samples, var=1, burnin=2000)

quantile(samples[-(1:2000), 1], c(0.025, 0.975))








## "seizures" example for STAT365

library(nimble)
load('~/Downloads/seizures.RData')
N <- dim(seizures$Counts)[1]

code <- nimbleCode({
    b0 ~ dnorm(0, sd=1000)
    bbase ~ dnorm(0, sd=1000)
    bage ~ dnorm(0, sd=1000)
    btrt ~ dnorm(0, sd=1000)
    sigma ~ dunif(0, 1000)
    sigma_patient ~ dunif(0, 1000)
    for(i in 1:N) {
        g[i] ~ dnorm(0, sd=sigma_patient)
        for(j in 1:4) {
            eps[i, j] ~ dnorm(0, sd=sigma)
            log(lambda[i,j]) <- b0 + bbase * log(baseline[i]) + bage * age[i] + btrt * treatment[i] + eps[i,j] + g[i]
            y[i,j] ~ dpois(lambda[i,j])
        }
    }
})
constants <- list(N = dim(seizures$Counts)[1], baseline=seizures$Baseline, age=seizures$Age, treatment=seizures$Treatment)
data <- list(y=seizures$Counts)
inits <- list(b0=1, bbase=0, bage=0, btrt=0, sigma=1, sigma_patient=1, g=rep(0,N), eps = array(0,c(N,4)))
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, onlySlice=TRUE)
conf$printSamplers()
##conf$removeSamplers('b0')
##conf$addSampler('b0', 'slice')
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samplesList <- runMCMC(Cmcmc, 500000, nburnin=10000, nchains=3, returnCodaMCMC=TRUE)

samples <- samplesList[[1]]
samplesPlot(samples)
samplesPlot(samples, 'btrt')
samplesPlot(samples, 'b0')
##samplesPlot(samples, 'bbase')
library(coda)
apply(samples, 2, effectiveSize)
##dim(samples)
gelman.diag(samplesList)

apply(samples, 2, effectiveSize)
apply(samples, 2, mean)
sqrt(apply(samples, 2, var)) / sqrt(apply(samples, 2, effectiveSize))


codalist <- runMCMC(Cmcmc, 100000, nchains=3, returnCodaMCMC=TRUE)

gelman.diag(codalist)
geweke.diag(codalist)
geweke.plot(codalist)

apply(codalist[[1]], 2, effectiveSize)
apply(codalist[[1]], 2, var) / apply(codalist[[1]], 2, effectiveSize)


## "seeds" example for STAT365
write.csv(df, '~/Downloads/Seeds.csv')


df <- read.csv('~/Downloads/Seeds.csv')
df

cucumber <- as.numeric(df$plant) - 1
fertB <- as.numeric(df$fertilizer) - 1
y <- df$germinations
n <- df$seeds
N <- dim(df)[1]

library(nimble)

sds <- c(.1, .5, 1, 5, 10)
samps <- array(0, c(10000, length(sds)))

for(i in seq_along(sds)) {
    sdC <- sds[i]
    code <- nimbleCode({
        sigma ~ dunif(0, 100)
        b0 ~ dnorm(0, sd=1000)
        bCuc ~ dnorm(0, sd=sdC)
        bFertB ~ dnorm(0, sd=1000)
        for(i in 1:N) {
            mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
            logit(p[i]) ~ dnorm(mu[i], sd=sigma)
            y[i] ~ dbinom(size=n[i], prob=p[i])
        }
    })
    constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n, sdC=sdC)
    data <- list(y=y)
    inits <- list(sigma=1, b0=1, bCuc=0, bFertB=0)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$mu
    calculate(Rmodel)
    conf <- configureMCMC(Rmodel)
    conf$printSamplers()
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    set.seed(0)
    samples <- runMCMC(Cmcmc, 100000)
    ##Cmcmc$run(10000)
    ##samples <- as.matrix(Cmcmc$mvSamples)
    samps[, i] <- samples[ 50001:60000, 'bCuc']
}

dim(samps)
dimnames(samps)
colnames(samps) <- paste0('sd', as.character(sds))
dimnames(samps)

samplesPlot(samples)

## prior sensitity analysis for different
## parametrizations of the 'sigma' term
niter <- 100000
samps <- array(0, c(niter, 3))
colnames(samps) <- c('sd', 'var', 'tau')


code <- nimbleCode({
    sigma ~ dunif(0, 1000)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], sd=sigma)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(sigma=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 1] <- samples[, 'sigma']
code <- nimbleCode({
    v ~ dunif(0, 1000)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], var=v)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(v=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 2] <- samples[, 'v']
code <- nimbleCode({
    tau ~ dgamma(0.001, 0.001)
    b0 ~ dnorm(0, sd=1000)
    bCuc ~ dnorm(0, sd=1000)
    bFertB ~ dnorm(0, sd=1000)
    for(i in 1:N) {
        mu[i] <- b0 + bCuc * cucumber[i] + bFertB * fertB[i]
        logit(p[i]) ~ dnorm(mu[i], tau=tau)
        y[i] ~ dbinom(size=n[i], prob=p[i])
    }
})
constants <- list(N=N, cucumber=cucumber, fertB=fertB, n=n)
data <- list(y=y)
inits <- list(tau=1, b0=1, bCuc=0, bFertB=0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
samples <- runMCMC(Cmcmc, niter)
samps[, 3] <- samples[, 'tau']

samps[, 'var'] <- sqrt(samps[, 'var'])
samps[, 'tau'] <- 1/sqrt(samps[, 'tau'])


head(samps)

dim(samps)
dimnames(samps)

samplesPlot(samps, burnin=50000)







## "pairs"
## studying estimation of normal variance / covariance,
## from paired observations from a common variance
## (and possibly different mean)


## this shows that for pairs x1,x2 ~ N(mean=[changing], var= constant v)
## using the average of (x2-x1)^2 / 2
## gives an unbiased estimate of the common variance v
estimateV <- function(n, v) {
    pairs <- t(replicate(n, rnorm(2, mean=runif(1,0,100), sd=sqrt(v))))
    vest <- apply(pairs, 1, function(x) ((x[2]-x[1])^2)/2)
    mean(vest)
}
n <- 1000
v <- 2
n.average.over=1000
estimates <- replicate(n.average.over, estimateV(n=n, v=v))
mean(estimates)

## this shows that the above (vest = (x2-x1)^2 / 2)
## is *not* the same as a traditional sd(x), if the x are used to
## generate a vector of n consecutive "pairs"
n <- 10
v <- 2
xs <- rnorm(n, sd=sqrt(v))
sd(xs)^2
pairs <- array(NA, c(n,2))
for(i in 1:(n-1)) pairs[i,] <- xs[i:(i+1)]
pairs[n,] <- xs[c(n,1)]
xs
pairs
vest <- apply(pairs, 1, function(x) ((x[2]-x[1])^2)/2)
sum(vest) / (n-1)
sd(xs)^2


## this shows that for multivariate-normal pairs
## x1,x2 ~ MVN(mean=[changing], Sigma = constant V)
## using the average of dif=x2-x1, dif %*% t(dif) / 2
## gives an unbiased estimate of the common covariance matrix V
estimateV <- function(n, V) {
    pairs <- lapply(1:n, function(x) {
        mus <- runif(2,0,100)
        list(rmvnorm(1, mus, sigma=V)[1,],
             rmvnorm(1, mus, sigma=V)[1,])
    })
    vest <- lapply(pairs, function(x) {
        dif <- x[[2]] - x[[1]]
        (dif %*% t(dif)) / 2
    })
    apply(array(unlist(vest), dim=c(2,2,n)), c(1,2), sum) / n
}
library(mvtnorm)
n <- 200
v1 <- 2
v2 <- 5
rho <- 0.6
V <- array(c(v1^2, rho*v1*v2, rho*v1*v2, v2^2), c(2,2))
n.average.over <- 200
estimates <- replicate(n.average.over, estimateV(n=n, V=V))
V
apply(estimates, c(1,2), mean)









x <- rgamma(10000, 0.001, 0.001)
curve(dgamma(x, 0.001, 0.001), col='blue')
hist(x, prob=TRUE, breaks=50000, add=TRUE)
curve(dgamma(x, 0.001, 0.001), col='blue', add=TRUE)


## testing the new adaptive properties for covariance in RW_block sampler,
## not adapting until acceptance rate >= 0.15 at least once
## uses scaleHistory and propCovHistory

library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

code <- nimbleCode({
    a ~ dnorm(0, sd=100)
    b ~ dnorm(0, sd=100)
    c ~ dnorm(a/2, sd=.1)
    for(i in 1:n1) {
        y1[i] ~ dnorm(a, sd=0.1)
    }
    for(i in 1:n2) {
        y2[i] ~ dnorm(b, sd=0.1)
    }
})
n1 <- 10
n2 <- n1
y1 <- rnorm(n1, 3, 0.1)
y2 <- rnorm(n2, y1+5, 0.01)
constants <- list(n1=n1, n2=n2)
data <- list(y1=y1, y2=y2)
inits <- list(a = 0, b=0, c=0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
##conf$addSampler('a', 'RW')
conf$addSampler(c('a', 'b', 'c'), 'RW_block')
conf$printSamplers()
conf$addMonitors(c('a', 'b', 'c'))
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(4000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

samplesPlot(samples)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])


samplesPlot(samples)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])








## getting model running for Colin Lewis-Beck at Iowa State

## Create Constant W, G and F matrices to pass to NIMBLE
library(nimble)
y1 <- rep(0,2190)
W<-diag(2)
F1<-matrix(c(1,0), nrow = 1, ncol = 2)
G<-matrix(c(cos(2*pi/365),-sin(2*pi/365),sin(2*pi/365),cos(2*pi/365)),nrow=2,ncol=2)
smosCode <- nimbleCode({
    ## Initial Values for States
    c0[1]<-0
    c0[2]<-0
    m0[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]   ## CHANGE 1: add full indexing to m0[...]
    x0[1:2]~dmnorm(c0[1:2],m0[1:2,1:2])        ## CHANGE 2: add full indexing to x0[...]
    ## initial values
    m[1:2,1]<-G[1:2,1:2] %*% x0[1:2]           ## CHANGE 3: add indexing to m[...], inprod() only returns a scalar (not a vector or matrix) so use %*%
    var0[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2] ## CHANGE 4: add full indexing to var0[...]
    x[1:2,1]~dmnorm(m[1:2,1],var0[1:2,1:2])
    y[1]~dnorm(inprod(F1[1,1:2],x[1:2,1]),sigmaSquaredInvv)
    ## Model
    for (t in 2:T){
        mu[1:2, t] <- G[1:2,1:2] %*% x[1:2,t-1]     ## CHANGE 5: something was funny with your mu[] declaration.  In addition to not having
                                                    ## the required indexing, it should should be indexed by t
        sigs[1:2,1:2]<-sigmaSquaredInvw*W[1:2,1:2]  ## CHANGE 6: add full indexing to sigs[...]
        x[1:2,t] ~ dmnorm(mu[1:2,t],sigs[1:2,1:2])  ## CHANGE 7: added appropriate 't' indexing to mu[...]
        y[t] ~ dnorm(inprod(F1[1,1:2],x[1:2,t]),sigmaSquaredInvv)
    }
    ## Priors
    sigmaSquaredInvv~dgamma(5,20)
    sigmaSquaredInvw~dgamma(5,200)
})
smosModel<-nimbleModel(code=smosCode,name='2harm',constants=list(T=2190,pi=pi,W=W,G=G,F1=F1),data=list(y=y1), inits=list(sigmaSquaredInvv=1,sigmaSquaredInvw=1))
Cmodel <- compileNimble(smosModel)




## doing the Surgeries example for STAT 365 homework

library(nimble)
df <- read.csv('~/Downloads/Surgeries.csv')
df

code <- nimbleCode({
    for(i in 1:N) {
        p[i] ~ dbeta(1, 1)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(p = rep(0.5, dim(df)[1]))

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=1)

mean(samples[,1])
median(samples[,1])
sd(samples[,1])
quantile(samples[,1], probs=c(0.025, 0.975))



code <- nimbleCode({
    mu ~ dnorm(0, 0.0001)
    sigma ~ dunif(0, 1000)
    for(i in 1:N) {
        logit(p[i]) ~ dnorm(mu, sd=sigma)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(mu=0, sigma=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=1)

mean(samples[,1])
quantile(samples[,1], probs=c(0.025, 0.975))

expit(mean(samples[,1]))
expit(quantile(samples[,1], probs=c(0.025, 0.975)))

sum(df$Mortalities)/sum(df$Surgeries)


code <- nimbleCode({
    b0 ~ dnorm(0, 0.0001)
    b1 ~ dnorm(0, 0.0001)
    sigma ~ dunif(0, 1000)
    for(i in 1:N) {
        logit(p[i]) ~ dnorm(b0 + b1*x[i], sd=sigma)
        y[i] ~ dbinom(size = n[i], prob = p[i])
    }
    old_procedure <- expit(b0)
    new_procedure <- expit(b0 + b1)
    mortality_difference <- new_procedure - old_procedure
})
constants <- list(N = dim(df)[1], n = df$Surgeries)
data <- list(y = df$Mortalities)
inits <- list(b0=0, b1=0, sigma=1, x=df$Procedure)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$addMonitors(c('old_procedure', 'new_procedure', 'mortality_difference'))
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

apply(samples, 2, mean)
colnames(samples)
samplesPlot(samples)
samplesPlot(samples, var=3:5, burnin=500)

mean(samples[,3])
quantile(samples[,3], probs=c(0.025, 0.975))

expit(mean(samples[,1]))
expit(quantile(samples[,1], probs=c(0.025, 0.975)))

sum(df$Mortalities)/sum(df$Surgeries)


## spatial NIMBLE model example from Abhirup Datta
library(nimble)
library(mvtnorm)
nimbleOptions(showCompilerOutput=TRUE)

## model code ###
gpCode <- nimbleCode({
    for(i in 1:n){ y[i] ~ dnorm(w[i],sd=tau)}
    Vw[,] <- sigma*chol(G[,])
    w[] ~ dmnorm(muw[], cholesky=Vw[,], prec_param=0)
    sigma ~ dunif(0, 100)
    tau ~ dunif(0, 10)
})

## data ###
set.seed(1)
n=100
sigs=4
w=as.vector(rmvnorm(1,rep(0,n),sigs*0.5^as.matrix(dist(1:n))))
Y=w+sqrt(0.5*sigs)*rnorm(n)

### setting up nimble model ####
gpModel <- nimbleModel(gpCode, constants = list(n = n, G=0.5^as.matrix(dist(1:n)), muw=rep(0,n)), 
    dimensions=list(y=n, Vw=c(n, n), w=n, Vy=c(n, n), G = c(n, n)), check=FALSE)

gpModel$setData(list(y = Y))

gpModel$sigma <- 1
gpModel$tau <- 0.1
gpModel$w <- rep(0,n)

CgpModel <- compileNimble(gpModel)

### MCMC ####

gpconf <- configureMCMC(CgpModel, print=TRUE)


gpconf$printSamplers()
gpconf$removeSamplers('w')
gpconf$printSamplers()
for(node in gpModel$expandNodeNames('w', returnScalarComponents=TRUE)) {
    gpconf$addSampler(node, 'RW')
}
gpconf$printSamplers()



gpconf$addMonitors(c('w'))	

gpMCMC <- buildMCMC(gpconf)

CgpMCMC <- compileNimble(gpMCMC, project = CgpModel)

##Nmcmc=100000
Nmcmc=1000
set.seed(1)
CgpMCMC$run(Nmcmc)

gpMCMCsamples <- as.matrix(CgpMCMC$mvSamples)

plot(gpMCMCsamples[,"w[1]"])

dimnames(gpMCMCsamples)

samplesPlot(gpMCMCsamples, var=3:10)


## implementing RJMCMC for a simple 2 models, logistic regression
## yi ~ bernoulli(logit(p) = a + bx)
set.seed(0)
n <- 50
x <- runif(n, -10, 10)
a <- .3
b <- 0.1
p <- 1/(1 + exp(-1*(a + b*x)))
y <- rbinom(n, prob=p, size=1)
##plot(x,y)
m1 <- glm(y~1, family=binomial())
m1
m2 <- glm(y~x, family=binomial())
m2
aic1 <- AIC(m1)  ## for this AIC function, "lower values indicate better fit"
aic2 <- AIC(m2)
aics <- c(aic1, aic2)
aic_min <- min(aics)
daics <- aics - aic_min
w1_aic <- exp(-daics[1]/2) / sum(exp(-daics/2))
w2_aic <- exp(-daics[2]/2) / sum(exp(-daics/2))
bic1 <- BIC(m1)  ## for this AIC function, "lower values indicate better fit"
bic2 <- BIC(m2)
bics <- c(bic1, bic2)
bic_min <- min(bics)
dbics <- bics - bic_min
w1_bic <- exp(-dbics[1]/2) / sum(exp(-dbics/2))
w2_bic <- exp(-dbics[2]/2) / sum(exp(-dbics/2))
w1_aic
w1_bic

## NEED TO KEEP WORKING ON LOGISTIC EXAMPLE FROM HERE
prior_sd <- 10
rj_proposal_sd <- 10
prior_m1 <- 0.5
{
    decide <- function(lMHR) {
        if(is.nan(lMHR)) return(FALSE)
        if(log(runif(1,0,1)) < lMHR) return(TRUE) else return(FALSE)
    }
    prior <- function(x) return(dnorm(x, 0, sd=prior_sd, log=TRUE))
    ll <- function(m, a, b, y, x) {
        if(m==1) return(sum(dnorm(y, a, 1, log=TRUE)))
        if(m==2) return(sum(dnorm(y, a+b*x, 1, log=TRUE)))
    }
    updateM <- function(mab, y, x) {
        ##browser()
        if(mab[1] == 1) {
            bprop <- rnorm(1, 0, rj_proposal_sd)   ## proposal for new slope 'b'
            lMHR <- prior(bprop) + ll(2,mab[2],bprop,y,x)  -  (ll(1,mab[2],0,y,x) + dnorm(bprop,0,rj_proposal_sd,log=TRUE)) + log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(2,mab[2],bprop)) else return(mab)
        }
        if(mab[1] == 2) {
            lMHR <- ll(1,mab[2],0,y,x) + dnorm(mab[3],0,rj_proposal_sd,log=TRUE)  -  (prior(mab[3]) + ll(2,mab[2],mab[3],y,x)) - log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(1,mab[2],0)) else return(mab)
        }
    }
    updateAB <- function(mab, y, x) {
        ## update a
        sd.aprop <- 0.1
        aprop <- rnorm(1,mab[2],sd.aprop)
        lMHR <- prior(aprop) + ll(mab[1],aprop,mab[3],y,x) - prior(mab[2]) - ll(mab[1],mab[2],mab[3],y,x)
        if(decide(lMHR)) mab <- c(mab[1],aprop,mab[3])
        ## update b
        if(mab[1] == 2) {
            sd.bprop <- 0.1
            bprop <- rnorm(1,mab[3],sd.bprop)
            lMHR <- prior(bprop) + ll(mab[1],mab[2],bprop,y,x) - prior(mab[3]) - ll(mab[1],mab[2],mab[3],y,x)
            if(decide(lMHR)) mab <- c(mab[1],mab[2],bprop) else mab
        }
        return(mab)
    }
}

set.seed(0)
iter <- 100000
mab <- c(1, 0, 0)  ## c(1, 0, 0)
##samp <- data.frame(a=rep(NA,iter), b=NA, c=NA)
samp <- array(NA, c(iter,3))
colnames(samp) <- c('m', 'a', 'b')
for(i in 1:iter) {
    mab <- updateM(mab, y, x)
    mab <- updateAB(mab, y, x)
    samp[i,] <- mab
}
df <- as.data.frame(samp)

mean(df$m==1)
##head(df,10)

c(mean(df$a[df$m==1]), m1$coef[1])
c(mean(df$a[df$m==2]), m2$coef[1])
c(mean(df$b[df$m==2]), m2$coef[2])

##plot(df$a, type='l')
##plot(1:iter, df$b)

prod(dnorm(y, 0, sd=sqrt(prior_sd^2 + 1^2)/10))   ## wrong for p(y|m1), don't know why

pym1 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd), 1))))
pym2 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd) + rnorm(1,0,prior_sd)*x, 1))))
pym1
pym2
ppm1 <- prior_m1*pym1 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm2 <- (1-prior_m1)*pym2 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm1
ppm2

w1_aic
w1_bic
mean(df$m==1)
ppm1

library(BMA)
bma <- bic.glm(x=data.frame(x=x), y=y, glm.family='gaussian')
##class(bma)
##str(bma)
bma$postprob[1]
bma$deviance
bma$label
bma$size
bma$which
bma$probne0
bma$postmean
bma$condpostmean
bma$mle





## implementing RJMCMC for a simple 2 models
## yi ~ N(a + bx, 1)
set.seed(0)
n <- 10
x <- runif(n, -1, 1)
a <- .3
b <- 0.1
y <- rnorm(n, a+b*x, 1)
##plot(x,y)
m1 <- lm(y~1)
m1
m2 <- lm(y~x)
m2
aic1 <- AIC(m1)  ## for this AIC function, "lower values indicate better fit"
aic2 <- AIC(m2)
aics <- c(aic1, aic2)
aic_min <- min(aics)
daics <- aics - aic_min
w1_aic <- exp(-daics[1]/2) / sum(exp(-daics/2))
w2_aic <- exp(-daics[2]/2) / sum(exp(-daics/2))
bic1 <- BIC(m1)  ## for this AIC function, "lower values indicate better fit"
bic2 <- BIC(m2)
bics <- c(bic1, bic2)
bic_min <- min(bics)
dbics <- bics - bic_min
w1_bic <- exp(-dbics[1]/2) / sum(exp(-dbics/2))
w2_bic <- exp(-dbics[2]/2) / sum(exp(-dbics/2))
w1_aic
w1_bic

prior_sd <- 10
rj_proposal_sd <- 10
prior_m1 <- 0.5
{
    decide <- function(lMHR) {
        if(is.nan(lMHR)) return(FALSE)
        if(log(runif(1,0,1)) < lMHR) return(TRUE) else return(FALSE)
    }
    prior <- function(x) return(dnorm(x, 0, sd=prior_sd, log=TRUE))
    ll <- function(m, a, b, y, x) {
        if(m==1) return(sum(dnorm(y, a, 1, log=TRUE)))
        if(m==2) return(sum(dnorm(y, a+b*x, 1, log=TRUE)))
    }
    updateM <- function(mab, y, x) {
        ##browser()
        if(mab[1] == 1) {
            bprop <- rnorm(1, 0, rj_proposal_sd)   ## proposal for new slope 'b'
            lMHR <- prior(bprop) + ll(2,mab[2],bprop,y,x)  -  (ll(1,mab[2],0,y,x) + dnorm(bprop,0,rj_proposal_sd,log=TRUE)) + log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(2,mab[2],bprop)) else return(mab)
        }
        if(mab[1] == 2) {
            lMHR <- ll(1,mab[2],0,y,x) + dnorm(mab[3],0,rj_proposal_sd,log=TRUE)  -  (prior(mab[3]) + ll(2,mab[2],mab[3],y,x)) - log((1-prior_m1)/prior_m1)
            if(decide(lMHR)) return(c(1,mab[2],0)) else return(mab)
        }
    }
    updateAB <- function(mab, y, x) {
        ## update a
        sd.aprop <- 0.1
        aprop <- rnorm(1,mab[2],sd.aprop)
        lMHR <- prior(aprop) + ll(mab[1],aprop,mab[3],y,x) - prior(mab[2]) - ll(mab[1],mab[2],mab[3],y,x)
        if(decide(lMHR)) mab <- c(mab[1],aprop,mab[3])
        ## update b
        if(mab[1] == 2) {
            sd.bprop <- 0.1
            bprop <- rnorm(1,mab[3],sd.bprop)
            lMHR <- prior(bprop) + ll(mab[1],mab[2],bprop,y,x) - prior(mab[3]) - ll(mab[1],mab[2],mab[3],y,x)
            if(decide(lMHR)) mab <- c(mab[1],mab[2],bprop) else mab
        }
        return(mab)
    }
}

set.seed(0)
iter <- 100000
mab <- c(1, 0, 0)  ## c(1, 0, 0)
##samp <- data.frame(a=rep(NA,iter), b=NA, c=NA)
samp <- array(NA, c(iter,3))
colnames(samp) <- c('m', 'a', 'b')
for(i in 1:iter) {
    mab <- updateM(mab, y, x)
    mab <- updateAB(mab, y, x)
    samp[i,] <- mab
}
df <- as.data.frame(samp)

mean(df$m==1)
##head(df,10)

c(mean(df$a[df$m==1]), m1$coef[1])
c(mean(df$a[df$m==2]), m2$coef[1])
c(mean(df$b[df$m==2]), m2$coef[2])

##plot(df$a, type='l')
##plot(1:iter, df$b)

prod(dnorm(y, 0, sd=sqrt(prior_sd^2 + 1^2)/10))   ## wrong for p(y|m1), don't know why

pym1 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd), 1))))
pym2 <- mean(replicate(500000, prod(dnorm(y, rnorm(1,0,prior_sd) + rnorm(1,0,prior_sd)*x, 1))))
pym1
pym2
ppm1 <- prior_m1*pym1 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm2 <- (1-prior_m1)*pym2 / (prior_m1*pym1 + (1-prior_m1)*pym2)
ppm1
ppm2

w1_aic
w1_bic
mean(df$m==1)
ppm1

library(BMA)
bma <- bic.glm(x=data.frame(x=x), y=y, glm.family='gaussian')
##class(bma)
##str(bma)
bma$postprob[1]
bma$deviance
bma$label
bma$size
bma$which
bma$probne0
bma$postmean
bma$condpostmean
bma$mle





## testing samplers used by JAGS for multivariate normal (dmnorm) nodes
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    aa <- a*a
    b ~ dnorm(0, aa)
    ab <- a*b
    d ~ dnorm(ab, 1)
    C[1,1] <- 1
    C[1,2] <- 0
    C[2,1] <- 0
    C[2,2] <- 1
    mu[1] <- 0
    mu[2] <- 0
    y[1:2] ~ dmnorm(mu[1:2], C[1:2, 1:2])
    ymu[1] <- exp(y[1])
    ymu[2] <- exp(y[2])
    y2[1:2] ~ dmnorm(ymu[1:2], C[1:2, 1:2])
})
constants <- list()
data <- list()
inits <- list(y = c(0,0), y2=c(0,0), a=1, b=1, d=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()


niter <- 10000
monitorVars <- c('y', 'y2')

constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)

library(rjags)
jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE)

list.samplers(jags_mod)

list.factories(jags_mod)
list.samplers
class(jags_mod)


dimnames(jags_out[[1]])
means <- apply(jags_out[[1]][,], 2, mean)
means
sds <- apply(jags_out[[1]][,], 2, sd)
sds





## doing the midterm question about normal mean hypothesis test
## "weights of sheep raised on a farm"
## now also with predictive distribution

library(nimble)

y <- c(78, 81, 77, 76, 75, 74, 78, 75, 77, 75)
n <- length(y)
np <- 5

code <- nimbleCode({
    mu ~ dnorm(75, sd=10)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, sd=3)
    }
    for(i in 1:np) {
        p[i] ~ dnorm(mu, sd=3)
    }
    pmean <- mean(p[1:np])
})
constants <- list(n=n, np=np)
data <- list(y=y)
inits <- list(mu = 75, p = rep(0,np))

Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$addMonitors(c('mu', 'p', 'pmean'))
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, resetFunctions=TRUE)

set.seed(0)
Cmcmc$run(1000000)
samples <- as.matrix(Cmcmc$mvSamples)
colnames(samples)
dim(samples)
samples[1:20,]
apply(samples, 2, mean)
apply(samples, 2, sd)
apply(samples, 2, var)
1/apply(samples, 2, var)

tau0 <- 1/100
mu0 <- 75
ybar <- mean(y)
tau <- 1/9
tau_n <- tau0 + n*tau
mu_n <- (mu0*tau0 + ybar*n*tau)/(tau0+n*tau)
mu_n
tau_n
1/tau_n
var(samples[,'mu'])

v <- 1/tau
v_n <- 1/tau_n
v_p1 <- v_n + v
v_p1
apply(samples, 2, var)

v_p1/5
var(samples[,'pmean'])
v_n/5 + v
v_n + v/5

var(apply(matrix(sample(as.numeric(samples[,2:6])), ncol=5), 1, mean))




## working on used cars regression example for Bayes course 365
df <- read.csv('~/Downloads/UsedCars.csv')
str(df)
hist(df$Age)
hist(df$HP)
dim(df)
head(df)

m <- lm(Price ~ Age + HP | Type, df=df)
summary(m)

summary(lm(Price ~ Age + HP, df=subset(df, Type==0)))
summary(lm(Price ~ Age + HP, df=subset(df, Type==1)))
summary(lm(Price ~ Age + HP + Type, df=df))


sd(residuals(lm(Price ~ Age + HP, df=subset(df, Type==0))))
sd(residuals(lm(Price ~ Age + HP, df=subset(df, Type==1))))




## recreating and fixing Dao's issue with the correlated SSM,
## where 'a' and 'b' don't mix until after 150,000 iterations
library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
##
code <- nimbleCode({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(1e-04, 1)
    sigOE ~ dunif(1e-04, 1)
    x[1] ~ dnorm(b/(1 - a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for (i in 2:t) {
        x[i] ~ dnorm(x[i - 1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
constants <- list(t = 100)
data <- list(y = c(20.24405,20.57693,20.49357,20.34159,20.45759,20.43326,20.20554,20.12860,20.14756,20.20781,20.23022,20.26766,20.22984,20.37703,20.13641,20.05309,19.95709,20.19303,20.30562,20.54443,20.91010,20.70580,20.42344,20.19795,20.28816,20.31894,20.76939,20.77023,20.83486,20.29335,20.40990,20.19601,20.04083,19.76056,19.80810,19.83129,19.69174,19.90069,19.87623,19.63371,19.62360,19.72630,19.64450,19.86779,20.17104,20.34797,20.32968,20.48027,20.46694,20.47006,20.51676,20.40695,20.18715,19.97552,19.88331,19.67831,19.74702,19.47502,19.24408,19.37179,19.38277,19.15034,19.08723,19.37051,19.14274,19.46433,19.62459,19.77971,19.54194,19.39081,19.61621,19.51307,19.34745,19.17019,19.26829,19.58943,19.77143,19.83582,19.71198,19.67746,19.75053,20.40197,20.49363,20.37079,20.19005,20.55862,20.48523,20.33071,19.97069,19.79758,19.83811,19.79728,19.86277,19.86836,19.92481,19.88095,20.24899,20.55165,20.22707,20.11235))
inits <- list(a = 0.95, b=1, sigPN = 0.2, sigOE=0.05, x = c(20.26036,20.51331,20.57057,20.35633,20.33736,20.47321,20.22002,20.14917,20.19216,20.26969,20.21135,20.22745,20.20466,20.41158,20.13408,20.08023,19.98956,20.13543,20.32709,20.55840,20.88206,20.74740,20.47671,20.14012,20.29953,20.33778,20.80916,20.75773,20.84349,20.35654,20.41045,20.20180,20.02872,19.74226,19.80483,19.81842,19.69770,19.84564,19.88211,19.70559,19.56090,19.73728,19.66545,19.88158,20.13870,20.39163,20.37372,20.47429,20.39414,20.42024,20.55560,20.40462,20.15831,19.89425,19.79939,19.72692,19.74565,19.42233,19.22730,19.36489,19.37289,19.19050,19.00823,19.35738,19.14293,19.48812,19.67329,19.82750,19.58979,19.43634,19.61278,19.56739,19.38584,19.19260,19.32732,19.65500,19.65295,19.84843,19.68285,19.69620,19.77497,20.31795,20.45797,20.32650,20.24045,20.60507,20.51597,20.30076,19.98100,19.86709,19.85965,19.74822,19.86730,19.90523,19.86970,19.87286,20.28417,20.46212,20.22618,20.13689))
##
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## [1] 183.3436
##

conf <- configureMCMC(Rmodel, nodes = NULL)

conf$addSampler(c('a', 'b'), 'RW_block')
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(adaptInterval=100))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(1,-.99,-0.99,1), c(2,2))))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(1,-.99,-0.99,1), c(2,2)), scale=0.01))
##conf$addSampler(c('a', 'b'), 'RW_block', control=list(propCov=array(c(0.001709168, -0.0341986, -0.0341986, 0.6844844), c(2,2))))


conf$printSamplers(c('a','b'))

conf$addSampler('sigOE', 'RW')
conf$addSampler('sigPN', 'RW')
for(node in Rmodel$expandNodeNames('x'))
    conf$addSampler(node, 'RW')
conf$resetMonitors()
conf$addMonitors(c('a', 'b', 'sigOE', 'sigPN'))
##conf$getMonitors()
##
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()   ## [1] 183.3436
nimbleOptions(showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
##
niter <- 500000
system.time(Cmcmc$run(niter))
##
samples <- as.matrix(Cmcmc$mvSamples)
dim(samples)
dimnames(samples)
##

## Plot1
dev.new(width=8, height=5)
par(mfrow=c(1,2))
plot(samples[1:300000,'a'], type='l', ylab='a')
plot(samples[1:300000,'b'], type='l', ylab='b')
par(mfrow=c(1,1))
getwd()
dev.copy2pdf(file='~/Downloads/plot1.pdf')

samplesPlot(samples, col=c('b','a'), ind=1:300000)

##samplesPlot(samples, col=c('sigOE','sigPN'))

##xs <- Rmodel$expandNodeNames('x')
##samplesPlot(samples, col=xs[ 1:10])
##samplesPlot(samples, col=xs[11:20])
##samplesPlot(samples, col=xs[21:30])
##samplesPlot(samples, col=xs[31:40])

##i <-  4:13   ## x[ 1]...x[10]
##i <- 14:23   ## x[11]...x[20]
##i <-  1:8    ## (a,b), sigOE, sigPN, x[1]...x[5]
##i <-  2:3    ## sigOE, sigPN

## all the latent state scales follow sigOE scale
##i <-  c(2,4:8)    ## sigOE, x[1]...x[5]
##nms <- sapply(conf$samplerConfs[i], function(x) x$target)
##nms
##ar <- do.call(cbind, lapply(Cmcmc$samplerFunctions$contentsList[i], function(x) x$scaleHistory))
##colnames(ar) <- nms
##samplesPlot(ar)
##dim(ar)
##samplesPlot(ar, burnin=500)

block_scales <- Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
block_propCovHistory <- Cmcmc$samplerFunctions$contentsList[[1]]$propCovHistory
## create block_propCovScale
block_propCovScale <- block_propCovHistory
for(i in 1:length(block_scales))   block_propCovScale[i,,] <- block_scales[i] * block_propCovHistory[i,,]
##dim(block_propCovScale)
block_scale_a <- apply(block_propCovScale, 1, function(x) sqrt(x[1,1]))
block_scale_b <- apply(block_propCovScale, 1, function(x) sqrt(x[2,2]))
block_cors <- apply(block_propCovHistory, 1, function(x) cov2cor(x)[1,2])
ar <- cbind(block_scales, block_scale_a, block_scale_b, block_cors)
colnames(ar) <- c('scale', 'sig_a', 'sig_b', 'cor')
samplesPlot(ar)

## final constant that scale approaches:
block_scales[length(block_scales)]
##propCov adapts very nicely to true covariance between 'a' and 'b'
cov(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
block_propCovHistory[dim(ar)[1],,]
## final adapted (and scaled) proposal corrleation is very accurate:
cor(samples[(dim(samples)[1]/2):(dim(samples)[1]), c('a','b')])
cov2cor(block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])
## final adapted (and scaled) proposal standard deviations for 'a' and 'b':
sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[1,1])
sqrt((block_scales[length(block_scales)] * block_propCovHistory[dim(ar)[1],,])[2,2])

## expand cor, sig_a, and sig_b by adaptInterval:
length(block_scale_a)
aI <- 200
block_scale_a_ex <- rep(block_scale_a, each=aI)
block_scale_b_ex <- rep(block_scale_b, each=aI)
block_cors_ex    <- rep(block_cors,    each=aI)
block_scales_ex  <- rep(block_scales,  each=aI)
samples_block_info <- cbind(samples[,'a'], samples[,'b'], block_scale_a_ex, block_scale_b_ex, block_cors_ex, block_scales_ex)
dimnames(samples_block_info)[[2]] <- c('a', 'b', 'sig_a', 'sig_b', 'cor', 'scale')

##samplesPlot(samples_block_info, ind=1:300000, col=c('b', 'a'))
##samplesPlot(samples_block_info, ind=1:300000, col=c('sig_a', 'sig_b', 'cor', 'scale'))

## this one is best:
## artificially trim 'b' samples:
samples_block_info_trim <- samples_block_info
samples_block_info_trim[,'b'] <- pmin(samples_block_info_trim[,'b'], 2)
samplesPlot(samples_block_info_trim, ind=1:300000, col=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'))
dev.copy2pdf(file='~/Downloads/plot2.pdf')

## same thing, on the "early" time scale
samplesPlot(samples_block_info_trim, ind=1:2000, col=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
dev.copy2pdf(file='~/Downloads/plot3.pdf')

samplesPlot(samples_block_info_trim, ind=1:20000, col=c('b', 'a', 'sig_a', 'sig_b', 'cor', 'scale'), densityplot=FALSE)
dev.copy2pdf(file='~/Downloads/plot4.pdf')







## testing addition of scaleHistory and propCovHistory
## back into RW and RW_block samplers


library(nimble)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0, b=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel, nodes=NULL)
conf$addSampler('a', 'RW')
conf$addSampler(c('a', 'b'), 'RW_block')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$scaleHistory
a <- Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])

set.seed(0)
Cmcmc$run(1000)

Cmcmc$samplerFunctions$contentsList[[1]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$scaleHistory
Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
a <- Cmcmc$samplerFunctions$contentsList[[2]]$propCovHistory
dim(a)
for(i in 1:dim(a)[1]) print(a[i,,])







## playing with colorRamp and colorRampPalette
## to get gradients of colors in plots

FLcrime <- read.csv("~/github/courses/stat201/data/fl_crime.csv")
head(FLcrime[, 1:4])
FLcrime <- FLcrime[, 2:4]
colnames(FLcrime) <- c("Crime", "Education", "Urbanization")
head(FLcrime)
pairs(FLcrime, panel=panel.smooth)
pairs(FLcrime)
plot(FLcrime$Education, FLcrime$Crime)
hist(FLcrime$Urbanization)
boxplot(FLcrime$Urbanization)
quantile(FLcrime$Urbanization)

plot(FLcrime$Education, FLcrime$Crime, pch=19)
with(subset(FLcrime, Urbanization<20), points(Education, Crime, col='green', pch=19))
with(subset(FLcrime, Urbanization>80), points(Education, Crime, col='red', pch=19))


pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(nrow(FLcrime))[cut(FLcrime$Urbanization, nrow(FLcrime))]
plot(FLcrime$Education, FLcrime$Crime, col=cols, pch=19)
with(subset(FLcrime, Urbanization<20), points(Education, Crime, col='green', pch=19))
with(subset(FLcrime, Urbanization>80), points(Education, Crime, col='red', pch=19))


x <- runif(100)
dat <- data.frame(x = x,y = x^2 + 1)

?colorRampPalette
rbPal <- colorRampPalette(c('red','blue'))

cc <- rbPal(100)[cut(dat$x, 100)]
plot(dat$x, dat$y, col = cc)

cc <- colorRampPalette(1:2)
plot(1:10, 1:10, col=cc(1:10/10))


plot(FLcrime$Education, FLcrime$Crime)

## practice problem 6 for STAT 365

## likelihood: yi ~ Normal(mu, sigma)
n <- 4
y <- c(100, 110, 112, 118)
sigma <- 10

## prior: mu ~ Normal(mu0, sigma0)
mu0 <- 100
sigma0 <- 10

## derive precisions
tau0 <- sigma0^-2
tau <- sigma^-2

## posterior distribution
tau_n <- tau0 + n*tau
mu_n <- (tau0*mu0 + n*tau*mean(y)) / tau_n

sigma_n <- 1/sqrt(tau_n)

mu_n
tau_n
sigma_n

## prior, likelihood, posterior plots

## 95% BCI

## one-sided hypothesis test
## H0: mu|y <= 100
## Ha: mu|y >  100

## two-sided hypothesis test
## H0: mu|y  = 100
## Ha: mu|y != 100

## samples from posterior

## 95% BCI

## one-sided hypothesis test
## H0: mu|y <= 100
## Ha: mu|y >  100


## two-sided hypothesis test
## H0: mu|y  = 100
## Ha: mu|y != 100


.





file <- 'HurricaneDamage.csv'
data <- read.csv(paste0('~/github/courses/stat201/data/', file))
Year <- data$Year
Damage <- data$Damage
Year2 <- Year[-1]
Damage2 <- Damage[-1]
m2 <- lm(Damage2~Year2)

res <- residuals(m2)

par(mfrow=c(1,1))
hist(res, breaks=12)

par(mfrow=c(2,1))
plot(Year2, res)
plot(fitted(m2), res)






## modifying hurricanes.csv data file into Hurricane_Damage.csv,
## for use in the STAT201 miderm
file <- 'hurricanes.csv'
data <- read.csv(paste0('~/Downloads/', file))
str(data)
names(data)
names(data)[1] <- 'Rank'
names(data)[2] <- 'Tropical.Cyclone'
names(data)[3] <- 'Year'
names(data)
dim(data)
data <- data[-(2:5),]
dim(data)
mtemp <- lm(data$Damage ~ data$Year)
coef(mtemp)
mtemp <- lm(data$Damage[-1] ~ data$Year[-1])
coef(mtemp)
fitted(mtemp)
newY <- c(data$Damage[1], 0.4*data$Damage[-1] + 0.6*fitted(mtemp))
data$Damage <- newY
write.csv(data, file='~/Downloads/HurricaneDamage.csv', row.names=FALSE)


file <- 'HurricaneDamage.csv'
data <- read.csv(paste0('~/Downloads/', file))
str(data)
x <- data$Year
y <- data$Damage
m <- lm(y~x)
summary(m)
plot(x, y)
abline(m, col='red')
coef(m)

sort(y)
which(y>40000)
i <- which(y>40000)
x2 <- x[-i]
y2 <- y[-i]
m2 <- lm(y2~x2)
summary(m2)
plot(x2, y2)
abline(m2, col='red')
coef(m2)
cor(x2,y2)
cor(x2,y2)^2





## doing LAX flight departures problem
## from the STAT201 midterm

x <- 3648
y <- 25843
w <- 3407
z <- 26134
N <- x+y+w+z

year <- c(rep(2014, x+y), rep(2015, w+z))
delay <- c(rep('ayes',x), rep('no',y), rep('ayes',w), rep('no',z))
tab <- table(year, delay)
tab

## (a) display contingency table
rbind(cbind(tab, margin.table(tab,1)), c(margin.table(tab,2), N))
##     ayes    no      
##2014 3648 25843 29491
##2015 3183 26284 29467
##     6831 52127 58958

## (b) percentage in 2014?
(x+y)/N * 100
prop.table(margin.table(tab, 1))[1]
## 50.02035 %

## (c) percentage of delayed departures in 2015?
w/(w+x) * 100
## 46.5964 %

## (d) diff. in prop, delays in 2015 relative to 2014?
ptab <- prop.table(tab, 1)
ptab[2,1] - ptab[1,1]
## -0.01567962 = -1.567962 %

## (e) interpret this diff. in prop.
## The fraction of LAX December 2015 departures that were delayed was 1.57 *percentage points* lower than the fraction of LAX December 2016 departures that were delayed.

iter <- 1000
data <- data.frame(year=year, delay=delay)
nr <- nrow(data)
ratio <- numeric(iter)
for(i in 1:iter) {
    ind <- sample(1:nr)
    tab <- table(data$year, data$delay[ind])
    cond.tab <- prop.table(tab, 1)
    ratio[i] <- cond.tab[2,1] - cond.tab[1,1]
}

hist(ratio, breaks=15)
tab <- table(data$year, data$delay)
cond.tab <- prop.table(tab, 1)
observed_ratio <- cond.tab[2,1] - cond.tab[1,1]
abline(v = observed_ratio, col = 'red', lwd = 2)



## demo of coin flips and LLN

n <- 30
n

flips <- rbinom(n, 1, 0.5)
flips

cumsum(flips)

seq_along(flips)

cumsum(flips) / seq_along(flips)

running_avg <- cumsum(flips) / seq_along(flips)

running_avg

par(mfrow=c(1,1))
plot(running_avg, type='l')
abline(h=0.5, col='red', lty=3)

ks <- 2:4
par(mfrow=c(length(ks),1), mar=c(2,2,2,2))

for(i in ks) {
    n <- 10^i
    flips <- rbinom(n, 1, 0.5)
    running_avg <- cumsum(flips) / seq_along(flips)
    plot(running_avg, type='l')
    abline(h=0.5, col='red', lty=3)
}





## testing seq_along in NIMBLE run code
library(nimble)

nfDef <- nimbleFunction(
    setup = function() {
        a <- 1:10
    },
    run = function() {
        for(i in seq_along(a)) {
            print(i)
        }
    }
)

Rnf <- nfDef()

Rnf$run()

Cnf <- compileNimble(Rnf)

Cnf$run()




## trying out MCMC with two Gibbs samplers for Normal

y <- c(4.3, 2.5, 3.2, 3.8, 2.9, 3.1, 4.2, 4.0)
n <- length(y)
n
y
mean(y)
sd(y)

## mu ~ dnorm(mu0=0, tau0=0.001)
## tau ~ dgamma(r0=0.001, v0=0.001)
## yi ~ dnorm(mu, tau)

mu0 <- 0
tau0 <- 0.001
r0 <- 0.001
v0 <- 0.001

## iter sampling iterations to run
iter <- 100000
samp <- cbind(mu = rep(NA,iter), tau = rep(NA,iter))

## inits: mu=0, tau=1
mu <- 0
tau <- 1

set.seed(0)
for(i in 1:iter) {
    mu <- rnorm(1, (tau0*mu0+tau*sum(y))/(tau0+n*tau), (tau0+n*tau)^-0.5)
    tau <- rgamma(1, r0+n/2, v0+0.5*sum((y-mu)^2))
    samp[i, 1] <- mu
    samp[i, 2] <- tau
}


1/2 * sum((y-mu)^2)
n/2 * (mean(y)-mu)^2


head(samp)

samp <- cbind(samp, sd = samp[,'tau']^-0.5)

apply(samp, 2, mean)
apply(samp, 2, median)
mean(y)
sd(y)
sd(y)^-2

samplesPlot(samp)

library(nimble)

code <- nimbleCode({
    mu ~ dnorm(0, 0.001)
    tau ~ dgamma(0.001, 0.001)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, tau)
    }
})
constants <- list(n=n)
data <- list(y=y)
inits <- list(mu=0, tau=1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

head(Rsamples, 10)
head(Csamples, 10)
head(samp, 10)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)




## fixing bug in conjugacy

library(nimble)
set.seed(0)
n <- 10
##y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=1))
y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=2), rnorm(n,20,sd=3))

code <- nimbleCode({
    for(i in 1:3) {
        ##sig[i] ~ dunif(0, 100)
        tau[i] ~ dgamma(0.001, 0.001)   ## using TAU
    }
    mu[1] ~ dnorm(0, sd=1000)
    mu[2] <- mu[1] + delta1
    delta1 ~ T(dnorm(0, sd=1000), 0, 1000)
    mu[3] <- mu[2] + delta2
    delta2 ~ T(dnorm(0, sd=1000), 0, 1000)
    for(i in 1:N) {
        ##means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2]
        means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2] + equals(z[i],3)*mu[3]
        ##sigmas[i] <- equals(z[i],1)*sig[1] + equals(z[i],2)*sig[2] + equals(z[i],3)*sig[3]
        taus[i] <- equals(z[i],1)*tau[1] + equals(z[i],2)*tau[2] + equals(z[i],3)*tau[3]   ## using TAU
        ##y[i] ~ dnorm(means[i], sd=sigmas[i])
        y[i] ~ dnorm(means[i], taus[i])   ## using TAU
        z[i] ~ dcat(pi[1:3])
    }
    for(i in 1:3) {
        pi0[i] ~ dgamma(1, 1)
        pi[i] <- pi0[i] / (pi0[1] + pi0[2] + pi0[3])
    }
})

N <- length(y)
constants <- list(N=N)
data <- list(y=y)
##inits <- list(mu=c(1,2,3), delta1=1, delta2=1, pi0=c(1,1,1), sig=c(1,1,1), z=rep(1:3, each=n))
inits <- list(mu=c(1,2,3), delta1=1, delta2=1, pi0=c(1,1,1), tau=c(1,1,1), z=rep(1:3, each=n))  ## using TAU
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$pi0
Rmodel$pi
Rmodel$z
Rmodel$mu
Rmodel$means
Rmodel$sigmas
Rmodel$taus

##undebug(Rmodel$checkConjugacy)
##Rmodel$checkConjugacy('mu[1]')
## 
##undebug(Rmodel$checkConjugacy2)
##Rmodel$checkConjugacy2('mu[1]')

conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, useConjugacy=FALSE)
conf$resetMonitors()
##conf$addMonitors(c('mu', 'z', 'sig'))
conf$addMonitors(c('mu', 'z', 'tau'))
conf$printSamplers('mu')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
calculate(Rmodel)
calculate(Cmodel)

##iter <- 20
## 
##set.seed(0)
##Rmcmc$run(iter)
## 
##set.seed(0)
##Cmcmc$run(iter)
## 
##Rsamples <- as.matrix(Rmcmc$mvSamples)
##Csamples <- as.matrix(Cmcmc$mvSamples)
## 
##sampNames <- colnames(Rsamples)
## 
##Rsamples[, sampNames]
##Csamples[, sampNames]
##Rsamples[, sampNames] - Csamples[, sampNames]

set.seed(0)
Cmcmc$run(50000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)

mus <- Rmodel$expandNodeNames(c('mu'))
sigs <- Rmodel$expandNodeNames(c('sig'))
taus <- Rmodel$expandNodeNames(c('tau'))
zs <- Rmodel$expandNodeNames(c('z'))

ind <- 1:dim(samples)[1]
ind <- 1:5000
samplesPlot(samples, col = mus, ind = ind)
samplesPlot(samples, col = sigs, ind = ind)
samplesPlot(samples, col = taus, ind = ind)
samplesPlot(samples, col = c('mu[3]','sig[3]'), ind = ind)
samplesPlot(samples, col = c('mu[3]','tau[3]'), ind = ind)
samplesPlot(samples, col = c(        'tau[3]'), ind = ind)

samplesPlot(samples, col = zs, ind = ind)

samplez <- samples[, zs]
samplezsum <- cbind(samplez, z1=apply(samplez, 1, function(x) sum(x==1)), z2=apply(samplez, 1, function(x) sum(x==2)), z3=apply(samplez, 1, function(x) sum(x==3)))
dim(samplez)
dim(samplezsum)

samplesPlot(samplezsum, col = c('z1','z2','z3'), ind=1:5000)
samplesPlot(samplezsum, col = c('z1','z2','z3'), ind=1:15000)

mean(y[1:n])
mean(y[(n+1):(2*n)])
mean(y[(2*n+1):(3*n)])



## doing the midterm question about normal mean hypothesis test
## "weights of sheep raised on a farm"

tau0 <- 1/100
mu0 <- 75

y <- c(78, 81, 77, 76, 75, 74, 78, 75, 77, 75)
n <- length(y)
ybar <- mean(y)
ybar
tau <- 1/9

taup <- tau0 + n*tau
mup <- (mu0*tau0 + ybar*n*tau)/(tau0+n*tau)

mup
taup
pnorm(75, mup, 1/sqrt(taup))


library(nimble)

code <- nimbleCode({
    mu ~ dnorm(75, sd=10)
    for(i in 1:n) {
        y[i] ~ dnorm(mu, sd=3)
    }
})
constants <- list(n=n)
data <- list(y=y)
inits <- list(mu = 75)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(500000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
mup
mean(samples[,1] < 75)
pnorm(75, mup, 1/sqrt(taup))






## inplementing a finite mixture model in NIMBLE

library(nimble)

set.seed(0)
n <- 500
y <- c(rnorm(n,0,sd=1), rnorm(n,10,sd=2), rnorm(n,20,sd=3))
##hist(y, breaks=30)

code <- nimbleCode({
    for(i in 1:3) {
        sig[i] ~ dunif(0, 100)
    }
    mu[1] ~ dnorm(0, sd=1000)
    mu[2] <- mu[1] + delta1
    delta1 ~ T(dnorm(0, sd=1000), 0, 1000)
    mu[3] <- mu[2] + delta2
    delta2 ~ T(dnorm(0, sd=1000), 0, 1000)
    for(i in 1:N) {
        means[i] <- equals(z[i],1)*mu[1] + equals(z[i],2)*mu[2] + equals(z[i],3)*mu[3]
        sigmas[i] <- equals(z[i],1)*sig[1] + equals(z[i],2)*sig[2] + equals(z[i],3)*sig[3]
        y[i] ~ dnorm(means[i], sd=sigmas[i])
        z[i] ~ dcat(pi[1:3])
    }
    ##pi[1:3] ~ ddirch(alpha[1:3])
    for(i in 1:3) {
        pi0[i] ~ dgamma(1, 1)
        pi[i] <- pi0[i] / (pi0[1] + pi0[2] + pi0[3])
    }
})

N <- length(y)
constants <- list(N=N)
data <- list(y=y)
inits <- list(sig=rep(1,3), mu=c(1,2,3), delta1=1, delta2=1, z=rep(1,N), pi0=c(1,1,1))

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$pi0
Rmodel$pi
Rmodel$z
Rmodel$mu
Rmodel$means
Rmodel$sigmas

conf <- configureMCMC(Rmodel)
##conf <- configureMCMC(Rmodel, useConjugacy = FALSE)

conf$printSamplers()
conf$printSamplers('mu')
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

calculate(Rmodel)
calculate(Cmodel)

##set.seed(0)
##Rmcmc$run(2)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


                   

## bootstrapping the slope 

data <- read.csv('~/Downloads/Animals\ (1).csv')
names(data)
head(data)

plot(data$gestation, data$longevity)

m <- lm(data$longevity ~ data$gestation)
m
summary(m)

s <- numeric(10000)
n <- dim(data)[1]
for(i in 1:10000) {
    ind <- sample(1:n, replace = TRUE)
    dat <- data[ind,]
    m <- lm(dat$longevity ~ dat$gestation)
    s[i] <- unname(m$coef[2])
}

hist(s)
m <- lm(data$longevity ~ data$gestation)
abline(v = unname(m$coef[2]), col='red')

              



## doing the bigmac problem

data <- read.delim('~/Downloads/bigmac.txt')
names(data)
head(data)
with(data, plot(EngSal, BigMac))
rownames(data)
colnames(data)
dimnames(data)

data <- read.delim('~/Downloads/bigmac.txt')
with(data, plot(EngSal, BigMac))
with(data, text(EngSal, BigMac, labels=rownames(data), cex=0.5, pos=4))

i <- -c(17,21)

m <- lm(data$EngSal[i] ~ data$BigMac[i])
coef(m)
summary(m)

plot(data$EngSal[i], residuals(m))

dim(data)
length(data$EngSal)
length(unname(residuals(m)))
names(m)
unname(residuals(m))
data$EngSal
data$BigMac

library(MASS)
boxcox
boxcox(m)

plot(data$EngSal[i], (data$BigMac[i])^.01)
plot(data$EngSal[i], log(data$BigMac[i]))

## likelihood: y ~ Normal(mu, sigma)
## what prior for sigma?
## sigma ~ Uniform(0, 1000)

sigma <- runif(100000, 0, 1000)

hist(sigma)



## verifying the change of variable formula

a <- runif(100000, 0, 1000)

par(mfrow=c(2,1))
plot(density(a))
plot(density(sqrt(a)))
curve(2*x/1000, col='red', add=TRUE)



## STAT 365 Monte Carlo Exercise 9.1
## comparing Bayesian and Frequentist estimators of pi

n <- 10
msef <- numeric()
mseb <- numeric()
biasf <- numeric()
biasb <- numeric()
varf <- numeric()
varb <- numeric()
pis <- c(1:5/100, 1:9/10, 95:99/100)
for(i in 1:length(pis)) {
    pi <- pis[i]
    samp <- rbinom(10000, size=n, prob=pi)
    pihatf <- samp / n
    pihatb <- (samp+1)/(n+2)
    biaspif <- mean(pihatf) - pi
    biaspib <- mean(pihatb) - pi
    varpif <- var(pihatf)
    varpib <- var(pihatb)
    msef[i] <- biaspif^2 + varpif
    mseb[i] <- biaspib^2 + varpib
    biasf[i] <- biaspif
    biasb[i] <- biaspib
    varf[i] <- varpif
    varb[i] <- varpib
}
par(mfrow = c(3,1))
plot(pis, biasf, main = 'Bias', col='red', type = 'b', ylim = range(c(biasf,biasb)))
lines(pis, biasb, col='blue', type = 'b')
plot(pis, varf, main = 'Variance', col='red', type = 'b', ylim = range(c(varf,varb)))
lines(pis, varb, col='blue', type = 'b')
plot(pis, msef, main = 'MSE', col='red', type = 'b', ylim = range(c(msef,mseb)))
lines(pis, mseb, col='blue', type = 'b')



## trying HW for STAT201

data <- read.delim('~/Downloads/Saratoga.txt')
names(data)
str(data)
attach(data)
boxplot(Price)
hist(Price)
boxplot(Living.Area)
hist(Living.Area)
plot(Living.Area, Price)

data <- read.delim('~/Downloads/Mauna-Loa-and-DJIA.txt')
attach(data)
ind <- Year>=1982
plot(DJIA, CO2.Avg)
plot(DJIA[ind], CO2.Avg[ind])
cor(DJIA[ind], CO2.Avg[ind])
str(data)

pairs(data, panel = panel.smooth)
?pairs

data <- read.delim('~/Downloads/Kentucky_Derby_2014.txt')
names(data)
str(data)
attach(data)
plot(data$Year, data$Speed..mph.)


## pulling together pieces of not having to recompile models and MCMCs, for Dao

library(nimble)
nimbleOptions(showCompilerOutput = TRUE) ### DELETE THIS later
code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dnorm(0, 1)
    c ~ dnorm(0, 1)
})
Rmodel <- nimbleModel(code, inits = list(a=0, b=0, c=0))
conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

model_orig <- Cmodel

## recover the R (uncompiled) model object
if(inherits(model_orig, 'CmodelBaseClass'))
    model_orig <- model_orig$Rmodel


##md <<- Rmodel_orig$modelDef
Rmodel <- model_orig$newModel(replicate = TRUE, check = FALSE)

conf_initial <- configureMCMC(Rmodel)

monitorsVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
conf_initial$addMonitors(monitorsVector, print=FALSE)

scalarNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
discreteInd <- sapply(scalarNodeVector, function(n) Rmodel$isDiscrete(n), USE.NAMES=FALSE)
scalarNodeVectorContinuous <<- scalarNodeVector[!discreteInd]
scalarNodeVectorContinuous

firstScalarNode <- scalarNodeVectorContinuous[1]
firstScalarNode

conf_initial$printSamplers()

samplersWeMightUse <- c('RW', 'slice', 'RW_block')
for(sampler in samplersWeMightUse)
    conf_initial$addSampler(target = firstScalarNode, type = sampler)

conf_initial$printSamplers()

Rmcmc_initial <- buildMCMC(conf_initial)
Cmodel <- compileNimble(Rmodel)
Cmcmc_initial <- compileNimble(Rmcmc_initial, project = Rmodel)

conf_new <- configureMCMC(oldConf = conf_initial)

conf_new$setSamplers()  ## remove all samplers
conf_new$printSamplers()

nodes <- c('a', 'b', 'c')
for(node in nodes)
    conf_new$addSampler(target = node, type = 'slice')
conf_new$printSamplers()

Rmcmc_new <- buildMCMC(conf_new)
Cmcmc_new <- compileNimble(Rmcmc_new, project = Rmodel)


nimCopy(from = model_orig, to = Cmodel, logProb = TRUE)
calculate(Cmodel)


## doing the random drug abuse survey using NIMBLE MCMC

library(nimble)

N <- 29
y <- 14


code <- nimbleCode({
    theta ~ dunif(0, 1)
    y ~ dbinom(size = N, prob = 0.25 + 0.5*theta)
})

constants <- list(N = N)
data <- list(y = y)
inits <- list(theta = 0.5)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)

samplesPlot(samples, col='theta')
theta <- samples[,'theta']
mean(theta)
median(theta)
quantile(theta, c(0.25, 0.5, 0.75))


## testing random drug abuse survey strategy, the variance of the results

f <- function(n, iter) {
    res <- rbinom(n=iter, size=n, prob=0.25)
    var(res)
}

## expect variance is 3*n/16
iter <- 100000

n <- 100
f(n, iter) - 3*n/16


g <- function(n, iter) {
    res <- numeric(iter)
    for(i in 1:iter) {
        a <- rbinom(n=1, size=n, prob=1/2)
        res[i] <- rbinom(n=1, size=a, prob=1/2)
    }
    var(res)
}

g(n, iter) - 3*n/16



test <- function(iter, n, theta) {
    res <- numeric()
    c <- numeric()
    d <- numeric()
    for(i in 1:iter) {
        a <- rbinom(1, size=n, prob=1/2)
        b <- n - a
        c[i] <- rbinom(1, size=a, prob=1/2)
        d[i] <- rbinom(1, size=b, prob=theta)
        res[i] <- c[i] + d[i]
    }
    var(d)
}

iter <- 10000
n <- 100
theta <- .2
test(iter, n, theta)

3*n/16

theta*n*(2-theta)/4



## solving 5 statisticians problem

one <- function() {
    cur <- 1  ## 1, 2, 3, 4, 5
    while(TRUE) {
        n <- runif(1)  # stay, +1, -1
        if(n < 1/3) {
            return(cur)
        } else if(n < 2/3) {
            cur <- ((cur - 1) - 1) %% 5 + 1
        } else cur <- ((cur + 1) - 1) %% 5 + 1
    }
}

many <- function(n) {
    ret <- numeric(n)
    for(i in 1:n)
        ret[i] <- one()
    return(ret)
}

prop.table(table(many(1000000))) * 11

##        1        2        3        4        5 
## 5.002679 2.000020 0.996589 0.999779 2.000933 
##     5/11,    2/11,    1/11,    1/11,    2/11


## playihg with Rmd

setwd('~/temp/lecTEMP')
getwd()
library(methods)
library(knitr)
library(rmarkdown)

list.files()

render('Lecture1Slides.rmd')

## testing Nick Michaud's question about nimbleFunctionLists
library(nimble)

bigFunction <- nimbleFunction(
    setup = function(N) {
        functions <- nimbleFunctionList(littleFunction_virtual)  ## CHANGE
        for(n in 1:N)
            functions[[n]] <-littleFunction(n)
    },
    run = function() {
        returnType(integer(0))
        sum_N <- 0
        for(n in 1:N) 
            sum_N <- sum_N + functions[[n]]$run()
        return(sum_N)
    }
)

## NEW
littleFunction_virtual <- nimbleFunctionVirtual(
    run = function() {
        returnType(integer(0))
    }
)

littleFunction <- nimbleFunction(
    contains = littleFunction_virtual,
    setup = function(n){},
    run = function(){
        returnType(integer(0))
        return(n)
    }
)

testFunction <- bigFunction(5)
testFunction$run()
CtestFunction <- compileNimble(testFunction)
CtestFunction$run()


## make configureMCMC respect dconstraint()
library(nimble)

code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[j,i] ~ dconstraint(w[j,i] > 0)
            w[j,i] ~ dnorm(theta[j], 1)
        }
        theta[j] ~ dnorm(mu, itau2)
    }
    itau2 ~ dgamma(a, b)
    mu ~ dnorm(0, .00001)
})

J <- 3
n <- rep(2, J)

y <- matrix( sample(c(0,1), sum(n), replace = TRUE), nrow = J)
m <- nimbleModel(code, constants = list(n = n, J = J), data = list(y = y))

conf <- configureMCMC(m)
conf$printSamplers()


## nested sampler function wrapper for Dao
library(nimble)

sampler_record_wrapperNEW <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control){
      numSamples <- 0
      before <- c(0, 0)
      after <- c(0, 0)
      samplerFunctionList <- nimbleFunctionList(sampler_BASE)
    ###### make sure to provide *named* arguments to this function
    ###### shouldn't require anything in control$control, if you don't want
    controlListForNestedSampler <- mcmc_generateControlListArgument(samplerFunction = control$sampler_function, control = control$control)
    samplerFunctionList[[1]] <- eval(call( control$sampler_function, model = model, mvSaved = mvSaved, target = target, control =  controlListForNestedSampler))}, 
    run = function() {
      ## these lines are new:
      numSamples <<- numSamples + 1
      setSize(before, numSamples)
      setSize(after, numSamples)
      before[numSamples] <<- model[[target]]
      ## back to the original sampler function code:
      samplerFunctionList[[1]]$run()
      ## this line new:
      after[numSamples] <<- model[[target]]
    },
    methods = list(
        reset = function() {samplerFunctionList[[1]]$reset()}
    ))

code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)

conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$removeSamplers('sigma')
conf$printSamplers()

## SEE CHANGES HERE
conf$addSampler(target = 'sigma', type = sampler_record_wrapperNEW, control = list(sampler_function = 'sampler_slice', control=list()))
conf$printSamplers()

## SEE CHANGES HERE
conf$addSampler(target = 'sigma', type = 'sampler_record_wrapperNEW', control = list(sampler_function = 'sampler_RW', control = list()))
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



## trying out a recursive nimble function
library(nimble)

Rnf <- nimbleFunction(
    run = function(x = double()) {
        if(x == 0 || x == 1) {
            return(1)
        } else {
            a <- 1
        }
    })



## testing the new NIMBLE function: runMCMC()
library(nimble)
code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)
Rmodel$setData(list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3)))
conf <- configureMCMC(Rmodel)
conf$getMonitors()
conf$setThin(10)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


## testing new functions for addSampler(), 'name', 'libraryTag', etc...
library(nimble)
code <- nimbleCode({
    mu ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 1000)
    for(i in 1:10) {
        x[i] ~ dnorm(mu*mu, sd = sigma)
    }
})
Rmodel <- nimbleModel(code)
conf <- configureMCMC(Rmodel)

debug(conf$addSampler)
undebug(conf$addSampler)
conf$printSamplers()
conf$removeSamplers('sigma')
conf$printSamplers()
conf$addSampler(target = 'sigma', type = sampler_slice, name='slice1')
conf$addSampler(target = 'sigma', type = 'slice', name='slice2')
conf$addSampler(target = 'sigma', type = sampler_slice)
conf$addSampler(target = 'sigma', type = 'slice')
conf$addSampler(target = 'sigma', type = sampler_RW, name='slice1')
conf$addSampler(target = 'sigma', type = 'RW', name='slice2')
conf$addSampler(target = 'sigma', type = sampler_RW)
conf$addSampler(target = 'sigma', type = 'RW')
conf$printSamplers()
conf$controlNamesLibrary



##debug(buildMCMC)
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



            

mc <- Cmcmc
mc <- Rmcmc

pb <- TRUE
pb <- FALSE

si <- TRUE
si <- FALSE

ni <- 10
ni <- 100

nb <- 5
nb <- 9



inits <- function() list(mu = rnorm(1,0,1000), sigma = runif(1,0,10))
##inits <- function() list(mu = 1:2, sigma = runif(1,0,10), x = 3)

initsList <- list(inits(), inits(), inits())
initsList <- list(inits(), inits())
initsList <- list(inits())

debug(runMCMC)
undebug(runMCMC)

debug(mcmc$run)

runMCMC(Rmcmc, niter = 3, nchains = 3, inits = inits())
runMCMC(Cmcmc, niter = 3, nchains = 3, inits = inits())

runMCMC(Rmcmc, niter = 3, nchains = 3, inits = initsList)
runMCMC(Cmcmc, niter = 3, nchains = 3, inits = initsList)

runMCMC(Rmcmc, niter = 300, nchains = 3, inits = inits())
runMCMC(Cmcmc, niter = 300, nchains = 3, inits = ii, nburnin=10)
a <- runMCMC(Rmcmc, niter = 300, nburnin=10, nchains = 4)

runMCMC(Rmcmc, niter = 300, nchains = 3, inits = inits(), nburnin=295)
runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29995)

runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE, returnCodaMCMC = TRUE)

runMCMC(Cmcmc, niter = 30000, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE, returnCodaMCMC = TRUE)

runMCMC(Cmcmc, niter = 30000, inits = inits, nburnin=29995, progressBar=TRUE, silent=TRUE)

runMCMC(Cmcmc, niter = 30000, nchains = 3, inits = inits, nburnin=29999, progressBar=TRUE, silent=TRUE, setSeed=TRUE)

runMCMC(Rmcmc, niter = 30, nchains = 3, inits = inits, nburnin=29, silent=TRUE, setSeed=TRUE)

runMCMC(Cmcmc, niter = 30000, nchains = 3, nburnin=29999, progressBar=TRUE, silent=TRUE, setSeed=TRUE)

runMCMC(Rmcmc, niter = 30, nchains = 3, nburnin=29, silent=TRUE, setSeed=TRUE)



## testing inconsistancy in dgamma() between R and C

library(nimble)

Rnf <- nimbleFunction(
    run = function(x = double(), a = double(), b = double()) {
        lp <- dgamma(x, a, b, log = 1)
        returnType(double())
        return(lp)
    }
)

Cnf <- compileNimble(Rnf)

x <- 6e-100
a <- 0.001
b <- 1.0
Rnf(x, a, b)
Cnf(x, a, b)




## making utility function for MCMC sample traceplots and density histograms

library(nimble)   ## get samples from birats2 model
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(20000)
samples <- as.matrix(Cmcmc$mvSamples)

library(nimble)   ## get samples from Dave Pleydel's multinomial test model
codeTest <- nimbleCode ({
    X[1:nGroups] ~ dmultinom(size=N, prob=pVecX[1:nGroups])
    Y[1:nGroups] ~ dmultinom(size=N, prob=pVecY[1:nGroups])
    for (ii in 1:nGroups)
        Z[ii] ~ dbeta(1 + X[ii], 1 + Y[ii]) })
nGroups   <- 5
N         <- 1E6
pVecX     <- rdirch(1, rep(1, nGroups))
pVecY     <- rdirch(1, rep(1, nGroups))
X         <- rmultinom(1, N, pVecX)[,1]
Y         <- rmultinom(1, N, pVecY)[,1]
Z         <- rbeta(nGroups, 1+X, 1+Y)
Xini      <- rmultinom(1, N, sample(pVecX))[,1]
Yini      <- rmultinom(1, N, sample(pVecY))[,1]
Constants <- list(nGroups=nGroups)
Inits     <- list(X=Xini, Y=Yini, pVecX=pVecX, pVecY=pVecY, N=N)
Data      <- list(Z=Z)
modelTest <- nimbleModel(codeTest, constants=Constants, inits=Inits, data=Data)
mcmcTest  <- buildMCMC(modelTest) 
cModelTest <- compileNimble(modelTest)
cMcmcTest <- compileNimble(mcmcTest, project=modelTest)
cModelTest$N     <- N <- 1E3
cModelTest$pVecX <- sort(rdirch(1, rep(1, nGroups)))
cModelTest$pVecY <- sort(rdirch(1, rep(1, nGroups)))
simulate(cModelTest, c('X','Y','Z'), includeData=TRUE)
niter  <- 1E4
cMcmcTest$run(niter)
samples <- as.matrix(cMcmcTest$mvSamples)

samplesPlot <- function(samples, ind=1:ncol(samples), burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright') {
    ## device window and plotting parameters
    dev.new(height=height, width=width)
    par(mfrow=c(1,2), cex=0.7, cex.main=1.5, lab=c(3,3,7), mgp=c(0,0.6,0), mar=c(2,1,2,1), oma=c(0,0,0,0), tcl=-0.3, yaxt='n', bty='l')
    ## process samples
    samples <- samples[, ind, drop=FALSE]
    if(!is.null(burnin))
        samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
    nparam <- ncol(samples)
    rng <- range(samples)
    ## traceplots
    plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
    for(i in 1:nparam)
        lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
    ## posterior densities
    xMin <- xMax <- yMax <- NULL
    for(i in 1:nparam) {
        d <- density(samples[,i])
        xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
    plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='Posterior Densities', xlab='', ylab='')
    alpha_density <- 0.2
    for(i in 1:nparam)
        polygon(density(samples[,i]), col=rainbow(nparam, alpha=alpha_density)[i], border=rainbow(nparam, alpha=alpha_density)[i])
    if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
        legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
}


dim(samples)
dimnames(samples)
samplesPlot(samples)

apply(samples, 2, mean)

samplesPlot(samples, ind=c(1,2,3,7), burnin=1000, legend.location='topleft')
samplesPlot(samples, ind=c(4,6), burnin=1000)
samplesPlot(samples, ind=c(5), burnin=1000)

## better to just use the plotting functions in coda package!!!
library(coda)
mcmcSamples <- as.mcmc(samples)
acfplot(mcmcSamples)
plot(mcmcSamples)



## playing with plot(density(x))
x <- rnorm(10000)
plot(density(x))
d <- density(x)
class(d)
ls(d)
length(d$x)
length(d$y)
plot(d$x, d$y, type='l')  ## this creates the standard plot(density(x)) plot

      plot(prior ~ param.x, ylim = yLims, type = "l", xlim = range(param.x), 
            xlab = "", ylab = "", main = "", axes = FALSE, ...)
        polygon(param.x, prior, col = "red")
        box()
        r = legend("topleft", legend = "Prior", lty = 1, bty = "n", 
            plot = FALSE)$text
        text(r$x, r$y, "Prior", adj = 0)
        plot(likelihood ~ param.x, type = "l", xlab = "", ylab = "", 
            main = "", axes = FALSE, ...)
        polygon(param.x, likelihood, col = "green")
        box()
        r = legend("topleft", legend = "Prior", lty = 1, bty = "n", 
            plot = FALSE)$text
        text(r$x, r$y, "Likelihood", adj = 0)
        plot(posterior ~ param.x, ylim = yLims, type = "l", xlab = "", 
            ylab = "", main = "", axes = F, ...)
        polygon(param.x, posterior, col = "blue")




## testing funny behavior of DSL round() and nimRound()
library(nimble)

Rnf <- nimbleFunction(
    run = function() {
        for(i in 0:50) {
            x <- i/10
            xRound <- round(x)
            print('x: ', x, ',   round(x): ', xRound)
        }
    }
)

Cnf <- compileNimble(Rnf)

Rnf()
Cnf()



## testing making my own progress bar in NIMBLE DSL nimbleFunction

library(nimble)

rfun <- nimbleFunction(
    run = function(pb = logical(default=TRUE)) {
        ##print('|')
        ##for(i in 1:20)   { a <- i %% 7; print(a) }
        ##print('|')
        ##cat('|')
        ##for(i in 1:3)   cat('-', i)
        ##cat('|')
        a <- 1
        pb <- pb & 0
        print(pb)
    }
)

##rfun()
cfun <- compileNimble(rfun)
##cfun()
cfun()


library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


library(nimble)
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
Rmcmc <- buildMCMC(Rmodel)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(20000)
Cmcmc$run(30007, reset = FALSE)         ## continue previous run
Cmcmc$run(10000, progressBar = FALSE)   ## turn off progress bar
Cmcmc$run(100003, reset = FALSE)        ## sort of slow run
Cmcmc$run(5000)     ## faster
Cmcmc$run(2000)     ## faster still
Cmcmc$run(1000)     ## ...
Cmcmc$run(500)
Cmcmc$run(100)
Cmcmc$run(40)  ## no bar when too few iterations



## playing with R progress bars: txtProgressBar

f <- function(n = 1e6, frac = 0.01) {
    pb <- txtProgressBar(style = 3, char = '-')
    nupdate <- floor(frac * n) 
    for(i in 1:n) {
        a <- rnorm(1)
        if(i %% nupdate == 0) {
            setTxtProgressBar(pb, i/n)
        }
    }
    setTxtProgressBar(pb, 1)
    close(pb)
}

f(n=2e6, frac = .01)



library(Bolstad)

## datasets in library(Bolstad):
## bears (exercise in chapter 3)
## slug  (exercise in chapter 14)

## functions that I might possibly use:
## decomp: makes plots of prior, likelihood, posterior, but only for class=Bolstad.
##         maybe steal a bunch of the code, make one that works for general samples?
##         BUT WAIT -- this only works for sorted (x,y) pairs -- not for samples.


## getting MCMC samples for logProbs of variables

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$getMonitors()
conf$addMonitors('logProb_a')
conf$getMonitors()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##paste0('logProb_', letters)

set.seed(0)
Cmcmc$run(10)
samples <- as.matrix(Cmcmc$mvSamples)
samples
apply(samples, 2, mean)


## tracking down error in conjugacy checking system,
## this model was submitted by: Eduardo Martins <egmartins@gmail.com>
## during the ISEC 2016 NIMBLE workshop
library(nimble)
load('~/Downloads/data_inits.RData')
fec_code <- nimbleCode({
    ## Likelihood
    ## Fecundity
    for (i in 1:ny.fec) {
        y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec) # observation model
        z.fec[i] ~ dnegbin(p.fec[i], size.fec) # sampling model
        p.fec[i] <- size.fec / (size.fec + mu.fec[i])
        log(mu.fec[i]) <- a.fec.raw + b.fec * z.sl[i] + nu.fec.k.raw[idx.k[i]] + nu.fec.l.raw[idx.l[i]]
    }
    a.fec <- a.fec.raw + mean(nu.fec.k.raw[1:nk]) + mean(nu.fec.l.raw[1:nl])
    for (k in 1:nk){
        nu.fec.k.raw[k] ~ dnorm(mu.nu.fec.k, tau.nu.fec.k)
        nu.fec.k[k] <- nu.fec.k.raw[k] - mean(nu.fec.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.fec.l.raw[l] ~ dnorm(mu.nu.fec.l, tau.nu.fec.l)
        nu.fec.l[l] <- nu.fec.l.raw[l] - mean(nu.fec.l.raw[1:nl])
    }
    ## Length
    for(i in 1:ny.fec) {
        z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
        mu.sl[i] <- a.sl.raw + b.sl * z.ocr[i] + nu.sl.k.raw[idx.k[i]] + nu.sl.l.raw[idx.l[i]]
    }
    a.sl <- a.sl.raw + mean(nu.sl.k.raw[1:nk]) + mean(nu.sl.l.raw[1:nl])
    for (k in 1:nk){
        nu.sl.k.raw[k] ~ dnorm(mu.nu.sl.k, tau.nu.sl.k)
        nu.sl.k[k] <- nu.sl.k.raw[k] - mean(nu.sl.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.sl.l.raw[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
        mu.nu.sl.l[l] <- b.nu.sl.l[1] * sst6.std[l] + b.nu.sl.l[2] * pdo6.std[l] + b.nu.sl.l[3] * nsal.std[l]
        nu.sl.l[l] <- nu.sl.l.raw[l] - mean(nu.sl.l.raw[1:nl])
    }
    ## Years of ocean residence
    for(i in 1:ny.fec){
        z.ocr[i] ~ dbern(p.z.ocr[i])
        logit(p.z.ocr[i]) <- a.ocr.raw + nu.ocr.k.raw[idx.k[i]] + nu.ocr.l.raw[idx.l[i]]
    }
    a.ocr <- a.ocr.raw + mean(nu.ocr.k.raw[1:nk]) + mean(nu.ocr.l.raw[1:nl])
    for (k in 1:nk){
        nu.ocr.k.raw[k] ~ dnorm(mu.nu.ocr.k, tau.nu.ocr.k)
        nu.ocr.k[k] <- nu.ocr.k.raw[k] - mean(nu.ocr.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.ocr.l.raw[l] ~ dnorm(mu.nu.ocr.l, tau.nu.ocr.l)
        nu.ocr.l[l] <- nu.ocr.l.raw[l] - mean(nu.ocr.l.raw[1:nl])
    }
    ## Priors
    ## Fecundity
    tau.y.fec <- pow(sig.y.fec, -2)
    sig.y.fec ~ dunif(0, 100)
    num.size.fec ~ dnorm(0, 0.0016)
    den.size.fec ~ dnorm(0, 1)
    size.fec <- abs(num.size.fec / den.size.fec)
    mu.nu.fec.k ~ dnorm(0, 0.0001)
    tau.nu.fec.k <- pow(sig.nu.fec.k, -2)
    sig.nu.fec.k ~ dunif(0, 100)
    mu.nu.fec.l ~ dnorm(0, 0.0001)
    tau.nu.fec.l <- pow(sig.nu.fec.l, -2)
    sig.nu.fec.l ~ dunif(0, 100)
    a.fec.raw ~ dnorm(mean.z.fec, tau.z.fec)
    b.fec ~ dnorm(0, 0.0001)
    ## Length
    tau.sl.i <- pow(sig.sl.i, -2)
    sig.sl.i ~ dunif(0, 100)
    mu.nu.sl.k ~ dnorm(0, 0.0001)
    tau.nu.sl.k <- pow(sig.nu.sl.k, -2)
    sig.nu.sl.k ~ dunif(0, 100)
    tau.nu.sl.l <- pow(sig.nu.sl.l, -2)
    sig.nu.sl.l ~ dunif(0, 100)
    a.sl.raw ~ dnorm(0, 0.0001)
    b.sl ~ dnorm(0, 0.0001)
    for (b in 1:nb.nu.sl.l) {
        b.nu.sl.l[b] ~ dnorm(0, 0.0001)
    }
    ## Years of ocean residence
    mu.nu.ocr.k ~ dnorm(0, 0.0001)     
    tau.nu.ocr.k <- pow(sig.nu.ocr.k, -2)
    sig.nu.ocr.k ~ dunif(0, 100)
    mu.nu.ocr.l ~ dnorm(0, 0.0001)
    tau.nu.ocr.l <- pow(sig.nu.ocr.l, -2)
    sig.nu.ocr.l ~ dunif(0, 100)
    a.ocr.raw ~ dnorm(0, 0.0001)
    ## Residuals and metrics for model assessment
    ## Fecundity
    for(i in 1:ny.fec){
        res.y.fec[i] <- y.fec[i] - z.fec[i]
        new.y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec)
    }
    for(i in 1:ny.fec){
        pres.z.fec[i] <- (z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        new.z.fec[i] ~ dnegbin(p.fec[i], size.fec)
        pres.new.z.fec[i] <- (new.z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        d.z.fec[i] <- pow(pres.z.fec[i], 2)
        d.new.z.fec[i] <- pow(pres.new.z.fec[i], 2)
    }
    fit.z.fec <- sum(d.z.fec[1:ny.fec])
    fit.new.z.fec <- sum(d.new.z.fec[1:ny.fec])
    ## Length
    for (i in 1:ny.fec) {
        res.z.sl[i] <- z.sl[i] - mu.sl[i]
        new.z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
    }
    for (l in 1:nl){
        res.nu.sl.l[l] <- nu.sl.l[l] - mu.nu.sl.l[l]
        new.nu.sl.l[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
    }
})
fec_mod <- nimbleModel(fec_code, constants = dat, inits = inits())
fec_cmod <- compileNimble(fec_mod)

##options(error = recover)

## next line errors out in conjugacy check:
## Error in x[[1]] : subscript out of bounds
mcmcConf <- configureMCMC(fec_mod)

mcmcConf$printSamplers()
mcmcConf$getMonitors()
fecMCMC <- buildMCMC(mcmcConf)
cfecMCMC <- compileNimble(fecMCMC, project = fec_mod)

cfecMCMC$run(10)

## testing of new setData() functionality:
## taking a character vector of variable names

library(nimble)
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')

Rmodel$getNodeNames(includeData = FALSE)

modelLs <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R', returnModelComponents = TRUE)
modelLs$model

Rmodel$beta
Rmodel$tau.c
Rmodel$r

Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('Y'))





## testing new MCMC runtime argument: time
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    c ~ dnorm(b + b^2, 1)
})
constants <- list()
data <- list(b=0, c=0)
inits <- list(a=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000, time = TRUE)
Cmcmc$getTimes()
Cmcmc$run(100000, time = TRUE, reset = FALSE)
Cmcmc$getTimes()



## error from Perry about not finding nimArray
## but I can't reproduce this error

library(nimble)

Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
##spec <- configureMCMC(Rmodel)
##spec$printSamplers()
##spec$getSamplerDefinition(1)
##spec$getSamplerDefinition(4)
##mcmc <- buildMCMC(spec)
mcmc <- buildMCMC(Rmodel)
mcmc$run(1)
##Error in samplerFunctions[[i]]$run() : could not find function "nimArray"

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R', data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)




## Soledad's compilation error

library(nimble)
x <- c(2, 2, 3)

## works
test1 <- nimbleFunction ( run = function (indx = double(1) ) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        aa <- indx[1] 
    } else {
        newindx <- indx
        aa         <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

## doesn't work
test1 <- nimbleFunction(run = function(indx = double(1)) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        newindx <- indx[1]
        aa          <- newindx     ## THIS IS THE MAIN CHANGE IN THIS FUNCTION
    } else {
        newindx <- indx
        aa          <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

Ctest1 <- compileNimble(test1)


test1(x)
Ctest1(x)





## using a custom distribution to have flexible length data
## for floriane plard population models

library(nimble)


dxxx <- nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1), length = integer(), log.p = double()) {
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dnorm(x[i], mu[i], sd=sigma[i], log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)

rxxx <- nimbleFunction(
    run = function(n = integer(), mu = double(1), sigma = double(1), length = integer()) {
        print('this should never run')
        x <- numeric(length)
        return(x)
    }
)

registerDistributions(list(
    dxxx = list(
        BUGSdist = 'dxxx(mu, sigma, length)',
        types    = c('value = double(1)', 'mu = double(1)', 'sigma = double(1)', 'length = integer()')
    )
))

code <- nimbleCode({
    y[1:N] ~ dxxx(mu[1:N], sigma[1:N], length)
})
constants <- list(N = 10)
data <- list()
inits <- list(length = 10)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$length
Rmodel$mu
Rmodel$sigma
Rmodel$y
Rmodel$getNodeNames(dataOnly = TRUE)

length <- 10
length
mu <- ((1:10)/2)[1:length]
mu
sigma <- (1+(1:10)/10)[1:length]
sigma
y <- c(1,0,1,3,2,5,4,6,7,5)[1:length]
y
yaug <- c(y, rep(as.numeric(NA), 10-length))
yaug

Rmodel$length <- length
Rmodel$length
Rmodel$mu[1:length] <- mu
Rmodel$mu
Rmodel$sigma[1:length] <- sigma
Rmodel$sigma
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(list(y=yaug))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$y

sum(dnorm(y, mu, sd=sigma, log=TRUE))
calculate(Rmodel)


debug(exprClasses_setSizes)



## testing MCEM for pump model
library(nimble)

pumpCode <- nimbleCode({ 
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta) 
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i]) 
  }
  alpha ~ dexp(1.0) 
  beta ~ dgamma(0.1,1.0) 
}) 
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                       31.4, 1.05, 1.05, 2.1, 10.5)) 
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)) 
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N)) 
newPump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits) 

pumpMCEM <- buildMCEM(model = newPump,
                      latentNodes = 'theta',
                      burnIn = 100,
                      mcmcControl = list(adaptInterval = 20),
                      boxConstraints = list( list( c('alpha', 'beta'), 
                          limits = c(0, Inf) ) ), 
                      buffer = 1e-6)

pumpMCEM(maxit = 20, m1 = 250, m2 = 500)
pumpMCEM(maxit = 50, m1 = 1000, m2 = 5000)


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
Rfun <- nimbleFunction(run = function() {
    ##ans <- numeric(10, value = 2)
    ##
    #x <- 100 
    #ans <- integer(x) 
    #for(i in 1:x) {
    #    ans[i] <- i 
    #}
    ##
    ##ans <- matrix(1, nrow = 10, ncol = 1)
    ##y <- 3
    ##ans <- numeric(10, value = y) 
    ##ans <- array(y, dim = 10) 
    ##ans <- array(y, dim = c(10))
    ##z <- numeric(10)
    ##z[5] <- 3
    ##x <- 20
    ##y <- 30
    ##ans <- integer(z[5], value = x + y) 
    ##ans <- array(x+y, dim = z[5], type = 'integer')
    ##ans <- array(x+y, dim = c(z[5]), type = 'integer')
    ##
    ##x <- array(0, c(4,5))
    ##ans <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2]) 
    ##ans <- array(0, dim = c(dim(x)[1], dim(x)[2]))
    ##
    x <- 1
    y <- 2
    z <- 3
    ans <- array(0, dim = c(x, y, z))
    ##returnType(double(1))
    ##returnType(integer(1))
    ##returnType(double(2))
    returnType(double(3))
    return(ans)
})
Cfun <- compileNimble(Rfun)

Rfun()
Cfun()
class(Rfun()[1])
class(Cfun()[1])
class(Rfun()[1,1])
class(Cfun()[1,1])
class(Rfun()[1,1,1])
class(Cfun()[1,1,1])

## creates a length-100 integer vector, containing 1, 2, ..., 100 

## creates a 10x1 ones-matrix 

## the following three lines are equivalent 
## each creates a length-10 vector, with all elements equal to y 

## the following two lines are equivalent 
## each creates an integer vector of length z[5], with all elements equal to x+y 

## the following two lines are equivalent 
## each one creates a matrix of 0's of the same size as matrix x 

## the following creates a 3-dimensional array of 0's 



## testing adding numeric(), integer(), array(), matrix()
library(nimble)
library(testthat)
##source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

expected <- numeric(10)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(10)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))


expected <- rep(3, length = 2)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(value = 3, length = 2)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- rep(9, 3)
Rfun <- nimbleFunction(run = function() {
    x <- numeric(10, value = 3)
    ans <- integer(x[2], x[2]+x[3]*2)
    returnType(integer(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1]), 'integer'))

expected <- array(4, c(10,11))
Rfun <- nimbleFunction(run = function() {
    ans <- array(4, c(10,11))
    returnType(double(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('integer', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- matrix(as.integer(0), nrow=4, ncol=5)
Rfun <- nimbleFunction(run = function() {
    x <- 4
    y <- 5
    ans <- matrix(init=FALSE, nrow=x, ncol=y, type='integer')
    returnType(integer(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1,1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1,1]), 'integer'))




expected
Rfun()
Cfun()


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        x <- 2
        bbb <- array(3, x+10, type='integer')
        ##print(bbb)
        bbb2 <- array(value=3, dim=c(3,4), type='double')
        ##print(bbb2)
        ccc <- array(value=3, dim=c(3,4), type='double', init=FALSE)
        ##print(ccc)
        x <- 5
        y <- matrix(2, ncol = 3, nrow = x)
        ##print(y)
        zxc <- array(2, c(x, dim(y)[2], y[2,2]))
        zxc[1,1,1] <- 98
        ##returnType(double(3));  return(zxc)
        zz <- numeric(4)
        zz[1] <- 99
        ##print(zz)
        z <- numeric(value=.5, 4)
        ##print(z)
        zz[2:4] <- z[1:3] + numeric(3, 10)
        ##print(zz)
        mat <- matrix(0, 5, 5)
        mat[1, 1] <- 97
        ##print(mat)
        arr <- array(4, type = 'integer', init = FALSE, dim=c(3,4))
        ##print(arr)
        ar2 <- array(dim = dim(z)[1])
        ##print(ar2)
        ar3 <- array(dim = c(dim(z)[1]))
        ##print(ar3)
        ar4 <- array(dim = c(dim(arr)[1], dim(arr)[2]), type = 'integer')
        ##print(ar4)
    }
)
Rnf <- nfDef()
##Rnf$run
Cnf <- compileNimble(Rnf, dirName)

Rnf$run()

Cnf$run()



tempdir()

##Cnf$run()



## ???

library(nimble)
sampler_dt <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        browser()
        1
        2
        3
        ###  node list generation  ###
        calcNodes  <- model$getDependencies(target)
    },
    run = function() {
        simulate(model, target)
        calculate(model, calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)
spec$getMonitors()
spec$addMonitors(c('theta', 'alpha2'))
spec$getMonitors()
spec$printSamplers()
spec$addSampler(type = 'dt', target = 'alpha')
spec$printSamplers()

Rmcmc <- buildMCMC(spec)

## testing new function in modelBaseClass: model$getNodeFunctions()
library(nimble)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)
Cmodel <- compileNimble(Rmodel)

nodes <- c('theta[3]')
nodes <- c('theta[3]', 'x[10]')
nodes <- c('theta[3]', 'x[9:10]')
Rmodel$getNodeFunctions(nodes)$calculate
Cmodel$getNodeFunctions(nodes)[[1]]$calculate

Rmodel$getNodeNames()
modelDef <- Rmodel$getModelDef()
gids <- mdef$nodeName2GraphIDs(nodes)
gids
mdef$graphIDs2indexedNodeInfo(gids)
Rmodel$nodeFunctions[[2]]$calculate


## making dmnorm conjugate sampler work with dependent
## nodes of different size from target node
library(nimble)
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

code <- nimbleCode({
    x[1:3] ~ dmnorm(mu0[1:3], prec = ident[1:3,1:3])
    mu_y2[1:2] <- asCol(a[1:2]) + B[1:2,1:3] %*% asCol(x[1:3])
    mu_y3[1:3] <- asCol(a[1:3]) + B[1:3,1:3] %*% asCol(x[1:3])
    mu_y5[1:5] <- asCol(a[1:5]) + B[1:5,1:3] %*% asCol(x[1:3])
    y2[1:2] ~ dmnorm(mu_y2[1:2], prec = prec_y[1:2,1:2])
    y3[1:3] ~ dmnorm(mu_y3[1:3], prec = prec_y[1:3,1:3])
    y5[1:5] ~ dmnorm(mu_y5[1:5], prec = prec_y[1:5,1:5])
})

mu0 <- rep(0,3)
ident <- diag(3)
a <- 11:15
B <- matrix(1:15, nrow=5, ncol=3, byrow=TRUE)
prec_y <- diag(1:5)

constants <- list(mu0=mu0, ident=ident, a=a, B=B, prec_y=prec_y)
data <- list(y2=1:2, y3=1:3, y5=1:5)
inits <- list(x=rep(0,3))

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
##spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##         x[1]       x[2]      x[3]
##[1,] 4.966869 -1.1252473 -3.922769
##[2,] 4.972376 -0.5288582 -4.511282


##node <- quote(mu0[1:3])
##node <- quote(prec_y[1:2, 1:2])
##node <- quote(prec_y[1:3, 1:3])
##node <- quote(prec_y[1:5, 1:5])
##nimble:::cc_expandDetermNodesInExpr(Rmodel, node)
##debug(nimble:::cc_expandDetermNodesInExpr)
##undebug(nimble:::cc_expandDetermNodesInExpr)
##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
##a <- spec$getSamplerDefinition(1)
##createNamedObjectsFromList(a, writeToFile = '~/temp/del.R')







## ???
library(nimble)

d <- read.csv("http://personal.bgsu.edu/~albert/data/gateway.csv")
library(dplyr)
d$month <- as.numeric(d$Month)
d$year <- d$Year - 2002

worshipCode <- nimbleCode({ 
  for (i in 1:N){
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- mu + epsilon[i] + alpha[month[i]]
  epsilon[i] ~ dnorm(0, sd=sigma.y)
  }
  for (j in 1:J){
  alpha[j] ~ dnorm(0, sd=sigma.month)
}
mu ~ dnorm(0, sd=1000)
sigma.y ~ dunif(0, 100)
sigma.month ~ dunif(0, 100)
})

worshipConsts <- list(N=484,
                      J=12,
                      K=10,
                      month=d$month)
worshipData <- list(y=d$Count)
worshipInits <- list(mu=5, sigma.y=.4, sigma.month=.4,
                     alpha=rep(0, 12))

worship <- nimbleModel(code = worshipCode, name = 'worship', 
                    constants = worshipConsts,
                    data = worshipData,
                    inits=worshipInits)

Rmodel <- worship

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## error produced in model$check() when using equals()
library(nimble)
code <- nimbleCode({
    for (i in 1:5) {
        x[i] ~ dbern(0.5)
        ##y[i] <- equals(x[i], 0)
    }
    sss <- sum(x[1:5])
    none <- equals(sss, 0)
})
constants <- list()
data <- list()
inits <- list(x = c(0,0,1,1,1))

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$x
Rmodel$sss
Rmodel$none

## testing new log=TRUE option of RW sampler

library(nimble)
code <- nimbleCode({
    sigma ~ dunif(0,10)
    a ~ dnorm(0, sd=sigma)
    b ~ dnorm(sigma, tau=3)
})
constants <- list()
data <- list(a=1, b=4)
inits <- list(sigma=1)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes=NULL)

spec$addSampler('sigma', 'RW_log')
### --- or ---
spec$addSampler('sigma', 'RW', control=list(log=TRUE))

Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)

##         sigma
## [1,] 3.535852
## [2,] 3.535852
## [3,] 3.535852
## [4,] 5.352671
## [5,] 5.352671
## [6,] 3.986347
## [7,] 3.963423
## [8,] 3.963423
## [9,] 3.963423
##[10,] 3.963423




## Github issue #107
## initializeModel issues false warning about RHSonly variable not initialized
## reduced case
library(nimble)

code <- nimbleCode({
    for (n in 1:N){
        a0[n,1:3] ~ dmnorm(mu[1:3], Sigma[1:3,1:3])
        ##lambda[n,1:3] <- exp(a0[n,1:3]) ## alternative, which shifts the problem to lambda
        for (k in 1:K){
            lambda[n,k] <- exp(a0[n,k])
            y[n,k]~dpois(lambda[n,k])
        }
    }
})

constants <- list(N = 3, K = 3)
data <- list(y = matrix(rpois(9, 3), nrow=3))
inits <- list(mu = rep(0,3), Sigma = diag(3))

Rmodel <- nimbleModel(code, constants, data, inits)

Cmodel <- compileNimble(Rmodel)

Rmcmc <- buildMCMC(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$run(1)
Cmcmc$run(1)



## testing new binary sampler
library(nimble)

code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dbern(0.6)
    c ~ dbern(0.05)
    d ~ dbin(prob=0.2, size=1)
    e ~ dbinom(prob=0.9, size=1)
    f ~ dbern(0.5)
    g ~ dbern(0.5)
    h ~ dbern(0.5)
    for(i in 1:10)
        yf[i] ~ dnorm(f, sd = 1)
    for(i in 1:10)
        yg[i] ~ dnorm(g, sd = 1)
    for(i in 1:10)
        yh[i] ~ dnorm(h, sd = 1)
    ##x ~ dnorm(0,1)
    ##y ~ dnorm(x*x, 1)
    ##z ~ dnorm(x*y, 1)
    ##zz ~ dnorm(z + z, 1)
})
constants <- list()
data <- list(yf = c(rep(0,2), rep(1,8)), yg = c(rep(0,8), rep(1,2)), yh = c(rep(0,5), rep(1,5)))
inits <- list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
Rmodel <- nimbleModel(code, constants, data, inits)

##spec <- configureMCMC(Rmodel, autoBlock = TRUE)
##spec$printSamplers()

Rmodel$isBinary('a')
Rmodel$isBinary('b')
Rmodel$isBinary('c')
Rmodel$isBinary('d')
Rmodel$isBinary('e')
Rmodel$isBinary('f')
Rmodel$isBinary('g')
Rmodel$isBinary('h')

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'binary', print=FALSE)
spec$addSampler('b', 'binary', print=FALSE)
spec$addSampler('c', 'binary', print=FALSE)
spec$addSampler('d', 'binary', print=FALSE)
spec$addSampler('e', 'binary', print=FALSE)
spec$addSampler('f', 'binary', print=FALSE)
spec$addSampler('g', 'binary', print=FALSE)
spec$addSampler('h', 'binary', print=FALSE)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)
##Rmcmc$run(3)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)
means


# Slice sampler name conflict with library(bbmle) 'slice' function
library(bbmle)
library(nimble)  # No messages about slice are generated

pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})

pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(2000000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## testing passing a compiled model as a setup argument to a NF

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        simulate(model, node)
    }
)

Rnf <- nfDef(Rmodel, 'a')
Rnf <- nfDef(Cmodel, 'a')

Cnf <- compileNimble(Rnf, project = Rmodel)

set.seed(0)
Rmodel$a
Rnf$run()
Rmodel$a

set.seed(0)
Cmodel$a
Cnf$run()
Cmodel$a



## testing passing a compiled model object to configureMCMC(), or buildMCMC()

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

###
Rmcmc <- buildMCMC(Cmodel)

###
spec <- configureMCMC(Cmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)


## testing problem with covariance matrix in BUGS model causing a chol() error at time of model checking

library(nimble)

code <- nimbleCode({
    C[1,1] <- 1
    C[1,2] <- 0
    C[2,1] <- 0
    C[2,2] <- 1
    mu[1] <- 0
    mu[2] <- 0
    y[1:2] ~ dmnorm(mu[1:2], prec = C[1:2, 1:2])
})
constants <- list()
data <- list(y = c(0,0))
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)





##If you're interested, here's the ecological/scientific background to the project: the potato psyllid spreads a disease called Zebra Chip Disease that is damaging potato and tomato crops in California and western US and Mexico. The disease is new, having only been described about 20 years ago, but the potato psyllid is native to the area. So what's changed to cause the disease outbreaks? I'm using data from where and when potato psyllid museum specimens were collected in California to test if psyllid populations have increased over the last century and if climate change may play a role. If so, this could explain the disease outbreaks. It's complicated though because museum data is messy, specifically, the museum specimens give me presence-only information. So I'm using lists of related species collected at the same time as the psyllids to infer psyllid absences, based on an assumed correlation between number of species collected and collecting effort. This is the "list_length" covariate that's in the model.

##BLOCK SAMPLE PARAMS THAT GO IN THE SAME LINEAR PREDICTOR,
##i.e. all the betapN terms in this one

{
    ## Priors
    ## For the site random effect
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu.alpha, tau.alpha) 
    }
    mu.alpha ~ dnorm(0, 0.001)
    tau.alpha <- 1 / (sigma.alpha * sigma.alpha)
    sigma.alpha ~ dunif(0, 5)
    
    ## Grand mean
    muq ~ dnorm(0, 0.001)
    
    ## For the fixed effect coefficients
    betaq ~ dnorm(0, 0.001)
    betap1 ~ dnorm(0, 0.001)
    betap2 ~ dnorm(0, 0.001)
    betap3 ~ dnorm(0, 0.001)
    betap4 ~ dnorm(0, 0.001)
    betap5 ~ dnorm(0, 0.001)
    betap6 ~ dnorm(0, 0.001)
    betap7 ~ dnorm(0, 0.001)
    betap8 ~ dnorm(0, 0.001)

    ## Likelihood
    for (i in 1:nlist){ # i = events (year-months)
        for(j in 1:nsite) { # j = sites

            logit(p[i,j]) <- betap1*year[i,j] + betap2*pow(year[i,j],2) + # Year quadratic effects
                betap2*month[i,j] + betap3*pow(month[i,j],2) + # month quadratic effects
                    betap4*aet[i,j] + betap5*cwd[i,j] + betap6*tmn[i,j] + betap7*tmx[i,j] + # Climate effects
                        alpha[j] # Random effects

            Y[i,j] ~ dbern(p[i,j]) # Occupancy probability

            logit(q[i,j]) <- Y[i,j] + muq + betaq*list_length[i,j] # Logistic regression for detection

            detectionMatrix[i,j] ~ dbern(q[i,j])  # Distribution for random part; observed presences relating to detection probability

        } #j
    } #i
}



## bug with forwardsolve in nimbleFunctions

library(nimble)
A <- matrix(c(2,1,0,3), nrow = 2)
b <- c(1,5)

nf <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        declare(x, double(1))
        setSize(x, 2)
        x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
        ##x <- forwardsolve(A[1:2,1:2], b[1:2])
        x[1:2,1:2] <- chol(A[1:2,1:2])
        returnType(double(1))
        ##declare(x, double(2))
        ##setSize(x, 2, 2)
        ##returnType(double(2))
        return(x)
    })

cnf <- compileNimble(nf, dirName = '.')

nf(A,b)
cnf(A,b)



## bug with forwardsolve in BUGS models

library(nimble)

code <- nimbleCode({
    A[1,1] <- 2
    A[2,1] <- 1
    A[1,2] <- 0
    A[2,2] <- 3
    b[1] <- 1
    b[2] <- 5
    x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
})

Rmodel <- nimbleModel(code)
Cmodel <- compileNimble(Rmodel)

Cmodel$A
Cmodel$b
Cmodel$x
calculate(Cmodel)
## Process R floating point exception: 8 at Thu Apr  7 14:38:02 2016



## testing using forwardsolve and backsolve in a user-defined function
## in a BUGS model

library(nimble)
A <- array(NA, c(2,2))
A[1,1] <- 2
A[2,1] <- 1
A[1,2] <- 0
A[2,2] <- 3
b <- c(1,5)

myFS <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        returnType(double(1))
        ret <- forwardsolve(A, b)
        return(ret)
    }
)

myChol <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- chol(A)
        return(ret)
    }
)

myInverse <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- inverse(A)
        return(ret)
    }
)

code <- nimbleCode({
    x1[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])  ## forwardsolve() directly here does NOT work
    x2[1:2] <- myFS(A = A[1:2, 1:2], b = b[1:2])
    C1[1:2,1:2] <- chol(A[1:2, 1:2])
    C2[1:2,1:2] <- myChol(A[1:2, 1:2])
    I1[1:2,1:2] <- inverse(A[1:2, 1:2])
    I2[1:2,1:2] <- myInverse(A[1:2, 1:2])
    xx1[1:2,1:2] <- chol(A[1:2, 1:2]) + 2 * inverse(A[1:2, 1:2])
    xx2[1:2,1:2] <- myInverse(A[1:2, 1:2]) * 2 + myChol(A[1:2, 1:2])
})

Rmodel <- nimbleModel(code, constants = list(A=A, b=b))
Cmodel <- compileNimble(Rmodel)

m <- Rmodel
m <- Cmodel

calculate(m)
m$x1; m$x2
m$C1; m$C2
m$I1; m$I2
m$xx1; m$xx2

## another example of MVN conjugate sampler, for test-mcmc.R
## using both cov and prec parametrizaions of MVN,
## and various linear links

library(nimble)

set.seed(0)
mu0 <- rep(0,5)
ident <- diag(5)
a <- array(rnorm(20), c(4,5))
B <- array(NA, c(4, 5, 5))
for(i in c(1,3)) {
    A <- array(rnorm(25,i), c(5,5))
    A <- A + t(A) + 5*i*diag(5)
    B[i,,] <- A
}
B[2,,] <- diag(5)
B[4,,] <- diag(5)
M_y <- array(NA, c(4, 5, 5))
for(i in 1:4)
    M_y[i,,] <- i * diag(5)
x <- array(0, c(4,5))
y <- array(0, c(4,5))

code <- nimbleCode({
    for(i in 1:2)
        x[i,1:5] ~ dmnorm(mu0[1:5], prec = ident[1:5,1:5])
    for(i in 3:4)
        x[i,1:5] ~ dmnorm(mu0[1:5], cov  = ident[1:5,1:5])
    for(i in 1:4)
        mu_y[i,1:5] <- asCol(a[i,1:5]) + B[i,1:5,1:5] %*% asCol(x[i,1:5])
    y[1,1:5] ~ dmnorm(mu_y[1,1:5], prec = M_y[1,1:5,1:5])
    y[2,1:5] ~ dmnorm(mu_y[2,1:5], cov  = M_y[2,1:5,1:5])
    y[3,1:5] ~ dmnorm(mu_y[3,1:5], prec = M_y[3,1:5,1:5])
    y[4,1:5] ~ dmnorm(mu_y[4,1:5], cov  = M_y[4,1:5,1:5])
})
constants <- list(mu0=mu0, ident=ident, a=a, B=B, M_y=M_y)
data <- list(y=y)
inits <- list(x=x)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##        x[1, 1]     x[2, 1]     x[3, 1]    x[4, 1]     x[1, 2]   x[2, 2]
##[1,] -0.7098220 -0.92579244  0.02859492 -0.9271968 -0.03440439 0.1133790
##[2,] -0.6642612  0.09915454 -0.13079585  0.6298495 -0.08432850 0.8806305
##        x[3, 2]    x[4, 2]     x[1, 3]   x[2, 3]     x[3, 3]    x[4, 3]
##[1,] 0.09213596  0.3257091  0.27465408 -1.410726 -0.12937422 -0.2311789
##[2,] 0.01777368 -0.6051182 -0.09322287 -1.723685  0.07759651  0.4119755
##       x[1, 4]     x[2, 4]     x[3, 4]     x[4, 4]   x[1, 5]  x[2, 5]
##[1,] 0.5071505 0.140652780  0.01778143  0.51383003 0.3639871 2.146307
##[2,] 0.4781541 0.002361055 -0.02517465 -0.06995162 0.4202336 1.225317
##         x[3, 5]    x[4, 5]
##[1,] -0.02913630 -0.2563026
##[2,]  0.04303507  0.2716980




## Richard McElreath Statistical Rethinking package,
## might have lots of good data and hierarchical (Bayesian) modeling examples
install.packages('rethinking')
library(rethinking)



## tesing MVN conjugate sampler

library(nimble)

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

Rmodel <- nimbleModel(code, data= data)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
##spec$getSamplerDefinition(1)
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples
##         mu[1]    mu[2]    mu[3]
## [1,] 3.877574 2.201065 1.084173
## [2,] 2.537805 1.916830 2.084021
## [3,] 4.426193 2.326982 1.567054
## [4,] 2.456652 1.921560 1.892777
## [5,] 2.959890 2.133226 1.506502
## [6,] 3.084360 1.560040 1.941621
## [7,] 3.178667 2.072212 2.611881
## [8,] 2.996492 2.210513 2.795391
## [9,] 2.494759 1.539610 2.118232
##[10,] 2.880285 1.826356 1.805385

Rsamples - Csamples
##              mu[1] mu[2] mu[3]
## [1,] -8.881784e-16     0     0
## [2,]  0.000000e+00     0     0
## [3,] -8.881784e-16     0     0
## [4,] -4.440892e-16     0     0
## [5,] -4.440892e-16     0     0
## [6,]  0.000000e+00     0     0
## [7,] -4.440892e-16     0     0
## [8,] -4.440892e-16     0     0
## [9,]  0.000000e+00     0     0
##[10,] -4.440892e-16     0     0


## inconsistent behavior between R and C chol()

n <- 4
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
mat <- A + L

sym <- function(mat, lower) {
    ind <- if(lower) lower.tri(mat) else upper.tri(mat)
    ans <- mat
    ans[!ind] <- 0
    ans <- ans + t(ans) + diag(diag(mat))
    ans
}

## R execution
chol(mat)
chol(sym(mat, FALSE))
chol(mat) - chol(sym(mat, FALSE)) ## R takes the UPPER part of the matrix

## Rnf execution
Rnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Rnf also UPPER part of matrix

## Cnf execution
Cnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Cnf now also uses UPPER part of matrix





## test case for forwardsolve, backsolve in the DSL
## A.triangularView<Eigen::Upper>().solve(b);
## could make a file called nimble_macros.h:
##define ......

library(nimble)
## inputs
n <- 4
tests <- c(
    'chol(A)',
    'forwardsolve(A, b)', 'forwardsolve(A, B)', 'forwardsolve(L, b)', 'forwardsolve(L, B)',
       'backsolve(A, b)',    'backsolve(A, B)',    'backsolve(R, b)',    'backsolve(R, B)',
           'solve(A, b)',        'solve(A, B)',        'solve(R, b)',        'solve(R, B)'
)
## derived objects
testNames <- as.character(sapply(tests, nimble:::Rname2CppName))
testDims <- sapply(grepl('_b_', testNames), function(x) if(x) 1 else 2)
methods <- mapply(function(test, dim) eval(substitute(function() { ans <- TEST; returnType(double(DIM)); return(ans) }, list(TEST = parse(text=test)[[1]], DIM = dim))), tests, testDims)
names(methods) <- testNames
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
nfDef <- nimbleFunction(
    setup = function(A, B, b) {
        R <- chol(A)
        L <- t(R)
    },
    run = function() {},
    methods = methods
)
Rnf <- nfDef(A, B, b)

Cnf <- compileNimble(Rnf)

testOneCase <- function(test, testName, Rnf, Cnf) {
    Rans <- eval(parse(text=test)[[1]])
    Rnfans <- eval(substitute(Rnf$TEST(), list(TEST=as.name(testName))))
    Cnfans <- eval(substitute(Cnf$TEST(), list(TEST=as.name(testName))))
    dif <- max(abs(Rans - Rnfans))
    if(dif > 1E-15) return(1)
    dif <- max(abs(Rans - Cnfans))
    if(dif > 1E-15) return(1)
    return(0)
}
for(i in seq_along(tests)) {
    test <- tests[i]
    testName <- testNames[i]
    testResult <- testOneCase(test, testName, Rnf, Cnf)
    if(testResult != 0) message('failed test: ', test)
}




library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ##ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
        ans4 <- chol(A + A)
    }
)

Cnf <- compileNimble(Rnf)   ## this works




nfDef <- nimbleFunction(
    setup = function(A, B, b) {},
    run = function() {
        ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
    }
)

Rnf <- nfDef(A, B, b)
Cnf <- compileNimble(Rnf)   ## FAIL


library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ans1 <- chol(A)
        print('chol(A)')
        print(ans1)
        ##
        ##ans2 <- forwardsolve(A, b)
        ##print('forwardsolve(A, b)')
        ##print(ans2)
        ####
        ##ans3 <- forwardsolve(A, B)
        ##print('forwardsolve(A, B)')
        ##print(ans3)
        ####
        ##ans4 <- backsolve(A, b)
        ##print('backsolve(A, b)')
        ##print(ans4)
        ####
        ans5 <- backsolve(A, B)
        print('backsolve(A, B)')
        print(ans5)
    }
)

Cnf <- compileNimble(Rnf)


set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
B <- array(1:9, c(3,3))
b <- 1:3

chol(A)
forwardsolve(A, b)
forwardsolve(A, B)
backsolve(A, b)
backsolve(A, B)

Rnf(A, b, B)
Cnf(A, b, B)

##debug(exprClasses_labelForEigenization)
##undebug(exprClasses_labelForEigenization)
##debug(exprClasses_eigenize)
##undebug(exprClasses_eigenize)
##writeCode(nimDeparse(code))
##writeCode(nimDeparse(compileInfo$nimExpr))


## using compareMCMCs, and MCMC comparisons

library(nimble)

code1 <- nimbleCode({ ## really toy model
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    sigma ~ dunif(0.5, 1)
    for(i in 1:3) y[i] ~ dnorm(b, sd = sigma)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5)) ## can also provide constants
results1 <- compareMCMCs(input1, MCMCs = c('nimble')) ## run a single case

results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')) ## run two more cases
results2[[1]]$efficiency ## inspect
results2[[1]]$summary
results2[[1]]$timing
make_MCMC_comparison_pages(results2, 'model1b') ## could generate a comparison of these

## can combine results as follows
## If the same name (“nimble” in this case) was used multiple times one will have to be renamed:
results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]])
results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results')
make_MCMC_comparison_pages(results3, 'model1c')



## testing calc_dmnormConjugacyContributions

nimble:::identityMatrix(3)

calc_dmnormConjugacyContributions(diag(3), diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 7*diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 8*diag(3), 2)
calc_dmnormConjugacyContributions(3*diag(3), 1+diag(3), 1)
calc_dmnormConjugacyContributions(3*diag(3), 2+diag(3), 2)

Cnf <- compileNimble(calc_dmnormConjugacyContributions)

Cnf(diag(3), diag(3), 1)
Cnf(diag(3), 7*diag(3), 1)
Cnf(diag(3), 8*diag(3), 2)
Cnf(3*diag(3), 1+diag(3), 1)
Cnf(3*diag(3), 2+diag(3), 2)



## compiler lesson 101
library(nimble)

nf <- nimbleFunction(
    setup = function() {
        a <- 1
    },
    run = function(arg = double()) {
        b <- a + arg
        returnType(double())
        return(b)
    })

Rnf <- nf()

## nfProc object is created every time we compile a new NF generator
## nfCompileInfo data structure for compilation info RCfunctionCompileClass"
## RCfunProcessing class to manage compilation info for a *single* run/member function
proj <- nimble:::nimbleProjectClass()
nfProc <- nimble:::nfProcessing(Rnf, project = proj)
nfProc$process(control = list(debug = TRUE))

## place to define 1st arg switching at top of genCpp_sizeProcessing
## switch name to nimArr_xxx in genCpp_processSpecificCalls






## testing of new parsed posterior text in conjugacy system
posteriorText = '{a<-1;
b
c
dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                            sd   = (prior_tau + contribution_tau)^(-0.5))
}'

posteriorText = '{ R <- chol(prior_prec + contribution_prec)
                        A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
                        mu <- backsolve(R, forwardsolve(t(R), A))[,1]
                        dmnorm_chol(mean = mu, cholesky = R, prec_param = 1) }'

parsedTotalPosterior <- parse(text = posteriorText)[[1]]
parsedTotalPosterior
if(parsedTotalPosterior[[1]] != '{') parsedTotalPosterior <- substitute({POST}, list(POST = parsedTotalPosterior))
parsedTotalPosterior
prePosteriorCodeBlock <- parsedTotalPosterior[-length(parsedTotalPosterior)]
prePosteriorCodeBlock
posteriorExpr <- parsedTotalPosterior[[length(parsedTotalPosterior)]]
posteriorExpr
rDistribution <- cc_makeRDistributionName(as.character(posteriorExpr[[1]]))
rDistribution
dDistribution <- as.character(posteriorExpr[[1]])
dDistribution
argumentExprs <- as.list(posteriorExpr)[-1]
argumentExprs
argumentNames <- names(argumentExprs)
argumentNames
rCallExpr <- as.call(c(as.name(rDistribution), 1, argumentExprs))
rCallExpr
dCallExpr <- as.call(c(as.name(dDistribution), quote(VALUE), argumentExprs, log = 1))
dCallExpr
posteriorVars <- all.vars(parsedTotalPosterior)
posteriorVars
neededPriorParams <- gsub('^prior_', '', posteriorVars[grepl('^prior_', posteriorVars)])
neededPriorParams
neededContributionNames <- posteriorVars[grepl('^contribution_', posteriorVars)]
neededContributionNames
neededContributionDims <- inferContributionTermDimensions(prior)


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

code <- nimbleCode({ x ~ dnorm(mu, var = v) })
inits <- list(x=0, mu=0, v=1)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
code <- nimbleCode({ x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3]) })
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

md <- nimbleModel(code, inits = inits, returnDef = TRUE, debug=TRUE)

Rmodel <- nimbleModel(code, inits = inits)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node
Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_var
Rmodel$nodes[[nn]]$get_sd
Rmodel$nodes[[nn]]$get_tau
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_cov


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)

## parametrize in terms of prec matrix (Q):
code <- nimbleCode({
    x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3])
})
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

Rmodel <- nimbleModel(code, inits = inits)
Cmodel <- compileNimble(Rmodel)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel[[lifted]]

Q <- A
ch <- chol(Q)
V <- solve(Q)
I <- diag(3)
Q
ch
V
t(ch) %*% ch
Q
backsolve(ch, forwardsolve(t(ch), I))
V

Rmodel$nodes[[nn]]$get_mean
Rmodel$nodes[[nn]]$get_mean()
Rmodel$getParam(nn, 'mean')

Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
Rmodel$getParam(nn, 'cholesky')
ch

Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_prec()
Rmodel$getParam(nn, 'prec')
Q

Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
Rmodel$getParam(nn, 'cov')
V

xx <- forwardsolve(t(ch), I)
t(xx) %*% xx

xx <- backsolve(ch, I)  ## would this work too???
xx %*% t(xx)

t(ch) %*% ch %*% V
t(ch) %*% ch
ch %*% V
Q


nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        print('chol:')
        chh <- model$getParam(node, 'cholesky')
        print(chh)
        print('prec:')
        q <- model$getParam(node, 'prec')
        print(q)
        print('cov:')
        vv <- model$getParam(node, 'cov')
        print(vv)
    }
)

Rnf <- nfDef(Rmodel, nn)
Cnf <- compileNimble(Rnf, project=Rmodel)


ch
Q
V
Rnf$run()
Cnf$run()

## store something into model variable 'Q'
Rmodel$Q <- Q

## calculate deterministic dependents
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node

Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
ch

## request 'prec'
Rmodel$nodes[[nn]]$get_prec  ## direct: model$Q
Rmodel$nodes[[nn]]$get_prec()
Q

## request 'cov'
## here's what we want to get cov: backsolve(ch,  forwardsolve(t(ch), I)  )
Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
V







    


## figuring out R / nimble functions:
## inverse, solve, forwardsolve, backsolve

library(nimble)
set.seed(0)
A <- array(rnorm(9), c(3,3))
A <- t(A) %*% A

nfDef <- nimbleFunction(
    setup <- function(A) {
        D <- diag(3)
        b <- 1:3
    },
    methods = list(
        inv = function() {
            Ainv <- inverse(A)
            print(Ainv)
        },
        ch = function() {
            R <- chol(A)
            print(R)
        },
        fs = function() {
            R <- chol(A)
            L <- t(R)
            Linv <- forwardsolve(L, D)
            print(Linv)
        },
        bs = function() {
            R <- chol(A)
        }
    )
)

Rnf <- nfDef(A)
Cnf <- compileNimble(Rnf)

print('inverse of A')
solve(A)
Rnf$inv()
Cnf$inv()

print('chol of A')
chol(A)
Rnf$ch()
Cnf$ch()

R <- chol(A)
L <- t(R)

print('forwardward solve')
forwardsolve(L, diag(3))

print('backward solve')
backsolve(R, diag(3))



## example of using nf$getDefinition() and nf_getDefinition()
library(nimble)
nfDef <- nimbleFunction(
    setup = function(a) { },
    run = function(b = double()) {
        a <<- a + b
        returnType(double())
        return(a)
    },
    methods = list(
        nothing = function() {
            print(a)
        })
)
Rnf <- nfDef(3)
Cnf <- compileNimble(Rnf)
getDefinition(nfDef)
Rnf$getDefinition()
Cnf$getDefinition()

## example of spec$getSamplerDefinition(..)
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    for(i in 1:140)
        c[i] ~ dnorm(0,1)
})

Rmodel <- nimbleModel(code, check = FALSE)

spec <- configureMCMC(Rmodel, print = TRUE)
spec$printSamplers()
spec$getSamplerDefinition(1)
spec$getSamplerDefinition('b')

spec$getSamplers()
spec$printSamplers()
spec$addSampler(type='RW_block', target=c('a', 'b'))
spec$getSamplers()
spec$printSamplers()


## performance comparison of of new dynamic conjugate samplers, rats example
library(nimble)
rats <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
Rmodel <- nimbleModel(rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits)
spec <- configureMCMC(Rmodel, print = TRUE)
out <- MCMCsuite(
    rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits,
    niter = 100000,
    MCMCs = c('staticConj', 'dynamicConj'),
    MCMCdefs = list(
        staticConj  = quote({ nimbleOptions(useDynamicConjugacy = FALSE)
                              configureMCMC(Rmodel) }),
        dynamicConj = quote({ nimbleOptions(useDynamicConjugacy = TRUE)
                              configureMCMC(Rmodel) })),
    makePlot = FALSE
)
sampdiff <- out$samples['staticConj',,] - out$samples['dynamicConj',,]
all(sampdiff == 0)
## [1] TRUE                      ## samples from each MCMC are identical
out$timing
## staticConj    dynamicConj
##     26.309         17.801     ## ~ 33% decrease in conjugate sampler runtime



## testing proper dmnorm conjugate calculations for Chris

## code = nimbleCode({
##    y[1:3]  ~ dmnorm(mu [1:3], cov = covy [1:3,1:3])
##    mu[1:3] ~ dmnorm(mu0[1:3], cov = covmu[1:3,1:3])
## })
## Concerned about some explicit inverses in conjugate update for mu,
## compared to use of already existing Cholesky factors.  In the R sampler,
## 1) we have: inverse(prior_prec + contribution_prec)
## 2) we use get_prec, which when we specify a dmnorm with a covariance
## has code like the following (in calc_dmnormAltParams()):
##  tmp <- t(cholesky) %*% cholesky
##  return(inverse(tmp))

library(nimble)

nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions()

code = nimbleCode({
   mu[1:3] ~ dmnorm(mu0[1:3], cov=covmu[1:3,1:3])
   y[1:3] ~ dmnorm(mu[1:3], cov = covy[1:3,1:3])
   a ~ dnorm(y[1], y[1])
   b ~ dgamma(a, y[1])
})

cov <- diag(3)
constants <- list(covy = cov, covmu = cov, mu0 = rep(0,3))
data <- list(y = 1:3)
inits <- list(mu = 4:6, a = 1, b = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$printSamplers()
spec$getSamplerDefinition('mu')
Rmcmc <- buildMCMC(spec)


Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$samplerFunctions$contentsList[[2]]$scale
Cmcmc$samplerFunctions$contentsList[[2]]
Cmcmc$samplerFunctions$contentsList[[2]][['scale']]
Cmcmc$samplerFunctions$contentsList[[2]]$scale



set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## testing what works and what doesn't with nimbleFunctionLists, getParam, etc..

library(nimble)

getNimbleOption('checkModel')
nimbleOptions(checkModel = FALSE)
getNimbleOption('checkModel')

code <- quote({
    x ~ dgamma(1, 1)
    y[1] ~ dnorm(1, tau = 1*x)
    y[2] ~ dnorm(2, tau = 2*x)
    y[3] ~ dnorm(3, tau = 3*x)
})
Rmodel <- nimbleModel(code, data=list(), inits=list(x=1, y=11:13))

nfDef <- nimbleFunction(
    setup = function(model) {
        nl <- model$expandNodeNames('y')
        N <- length(nl)
        nfl <- nimbleFunctionList(node_stoch_dnorm)
        for (i in 1:N) nfl[[i]] <- model$nodeFunctions[[nl[i]]]
    },
    run = function() {
        for(i in 1:N) {
            print('i = ', i)
            val <- nfl[[i]]$get_mean()
            print('mean = ', val)
            val <- nfl[[i]]$get_tau()
            print('tau = ', val)
            val <- nfl[[i]]$get_var()
            print('var = ', val)
            val <- nfl[[i]]$get_value()
            print('value = ', val)
        }
        returnType(double(0))
        return(1)
    }
)
Rnf <- nfDef(Rmodel)
Rnf$run()
Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)
Cnf$run()


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
    z ~ dbeta(10, 10)
    zz[1] ~ dbin(z, 10)
    zz[2] ~ dbin(z, 20)
    aaa ~ dbeta(10, 10)
    bbb[1] ~ dbin(aaa, 10)
    bbb[2] ~ dbin(aaa, 20)
    bbb[3] ~ dbin(aaa, 20)
    d ~ dnorm(0, 1)
    dd ~ dnorm(0, 1)
    ddd ~ dnorm(0, 1)
    e ~ dnorm(d + dd + ddd, 1)
    f ~ dlnorm(d + dd + ddd, 1)
    ff ~ dlnorm(d + dd, 1)
    g ~ dnorm(0, 1)
    gg ~ dnorm(5+g, 1)
})
constants <- list()
inits <- list(x = 0.5, z=0.6, aaa=0.3, d=0, dd=0, ddd=0, g=0)
data = list(y = c(3,4), zz=c(5,5), bbb=c(3,3,3), ff=1, e=0, f=1, gg=0)
Rmodel <- nimbleModel(code, constants, data, inits)


spec <- configureMCMC(Rmodel)

spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Rmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(3)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]
##     x z aaa d dd ddd g
##[1,] 0 0   0 0  0   0 0
##[2,] 0 0   0 0  0   0 0
##[3,] 0 0   0 0  0   0 0
Rsamples[, Rnames]
##             x         z       aaa         d         dd        ddd         g
##[1,] 0.3254630 0.3746455 0.3658422 0.6362147 -0.2698403 -1.1333402 -3.156596
##[2,] 0.1975665 0.3995488 0.3210974 0.3695457 -0.2843177 -0.2239394 -2.711577
##[3,] 0.1900000 0.4199454 0.2194685 0.5430496 -0.9140867  0.1178770 -2.233141


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
})
data = list(y = c(3,4))
constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples
as.numeric(Rsamples) - c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)



## testing of new dynamic conjugate samplers
### checkConjugacy_demo3_run.R - various conjugacies
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})
data = list()
constants <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('x', 'c'), control = list(scale = 0.01))
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]

Rsamples[, Rnames] - cbind(c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977), c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))


## testing of new dynamic conjugate samplers
### Dirichlet-multinomial conjugacy
# as of v0.4, exact numerical results here have changed because
# ddirch now sometimes returns NaN rather than -Inf (when an
# alpha is proposed to be negative) -- this changes the RNG
# sequence because NaN values result in no runif() call in decide()
# single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- quote({
    y[1:K] ~ dmulti(p[1:K], n);
    p[1:K] ~ ddirch(alpha[1:K]);
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001);
    }
               })
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
constants <- list(n=n, K=K)
data <- list(y = y)
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)[6:10]
p
apply(samples, 2, mean)[6:10] - p



## testing of new dynamic conjugate samplers
## rats example
library(nimble)
lst <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
code <- lst$model
data <- lst$data[c('Y')]
constants <- lst$data[c('N', 'T', 'x')]
inits <- lst$inits

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Cmcmc$run(5000)
samples <- as.matrix(Cmcmc$mvSamples)

a <- samples
b <- samples
a
b
a-b
any(a-b != 0)




## phase transition regression example for Mark.

## generate data
sigma <- 5
a <- 2
b <- 0.5
t0 <- 30
x <- 1:100
n <- length(x)
y <- rnorm(n, a + (x>t0) * (b*(x-t0)), sd=sigma)
plot(x, y)

## fit NIMBLE model
library(nimble)

code <- nimbleCode({
    a ~ dunif(-100, 100)
    b ~ dunif(-100, 100)
    t0 ~ dunif(0, 100)
    sigma ~ dunif(0, 100)
    for(i in 1:n) {
        phase[i] <- step(x[i] - t0)
        mu[i] <- a + phase[i] * b * (x[i]-t0)
        y[i] ~ dnorm(mu[i], sd = sigma)
    }
})
constants <- list(n=n, x=x)
data <- list(y=y)
inits <- list(sigma=1, a=0, b=1, t0=50)

## MCMC Suite will let us do a few MCMCs quickly, and look at results easily
out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_slice', 'bugs'),
    niter = 20000,
    burnin = 2000,
    makePlot = TRUE,
    savePlot = FALSE,
    calculateEfficiency = TRUE
)

out$summary
out$efficiency
dim(out$samples)
out$timing





## nimble problem with using 'sd' name for 'tau' parameter, e.g.
library(nimble)

code <- nimbleCode({
    sd ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd)
    shape1 ~ dnorm(0, 1)
    f ~ dbeta(mean = y, sd = shape1)
})

m <- nimbleModel(code, inits = list(sd=.1))

m$sd

m$getModelDef()$printDI()

code <- nimbleCode({
    sd2 ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd2)
})
m <- nimbleModel(code, debug=TRUE)

m$sd2






     
## Plot the fake MCMC output
denoverplot(fakemcmc, fakemcmc2)
     denoverplot(fakemcmc, fakemcmc2, style="plain",
                 col=mcmcplotsPalette(3, type="grayscale"),
                 ci=0.95, greek=TRUE)
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data")
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data", greek=TRUE)

setwd('/Users/dturek/GitHub/legacy/autoBlock/')
getwd()

source("autoBlock.R")
load(file.path("data", "model_test.RData"))
saveSamples <- TRUE
niter <- 10000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dftest <- ab$summary
save(dftest, file = file.path("results_samples", "results_test.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_test_samples.RData"))
}

class(burnedSamplesList)
names(burnedSamplesList)

i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])

dimnames(ab$samples[[3]])


ret$summary

dftest <- autoBlock(code, constants, data, inits, 10000, runList)$summary
save(dftest, file = file.path("results_hclust_complete2", "results_test.RData"))


rm(list=ls())
getwd()
load('results_samples/results_test.RData')
ls()
dftest

rm(list=ls())
getwd()
load('results_samples/results_test_samples.RData')
ls()
niter
class(burnedSamplesList)
names(burnedSamplesList)
i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])







## renaming stupid data frame from the autoBlock results
setwd('~/GitHub/legacy/autoBlock')

rm(list=ls())
file <- 'results/results_samplingEfficiency.RData'
load(file)
ls()
dfsamplingEfficiency <- dfSamplingEfficiency
rm(dfSamplingEfficiency)
ls()
save(list='dfsamplingEfficiency', file = file)

rm(list=ls())
file <- 'results/results_computationalRequirement.RData'
load(file)
ls()
dfcomputationalRequirement <- dfComputationalRequirement
rm(dfComputationalRequirement)
ls()
save(list='dfcomputationalRequirement', file = file)





## working on the spatial capture-recapture (SCR) models

tr<-seq(15,85, length=10)

X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. traps

plot(X, xlim=c(0,100), ylim=c(0,100), pch=3, cex=0.75)

set.seed(10)
xlim <- c(0,100); ylim <- c(0,100)      # Area 100*100=1e4
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/10000
mu <- 50                 # Density
N <- rpois(1, mu*A) ;N   # Generate population


s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(s, pch=16, col=2)

sigma <- 5
lambda0 <- 0.4
J <- nrow(X)
K <- 5
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
    dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2]-s[,2])^2)
    lambda <- lambda0*exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        yy[,j,k] <- rpois(N, lambda)
    }
}

n <- apply(yy, c(2,3), sum)

# Plot capture events
tot<-apply(n, 1,sum)
symbols(X, circles=tot, inches=F, bg="#00000022", add=T)
points(X, pch=3, cex=0.75); points(s, pch=16, col=2)


library(nimble)

## define the model
code <- nimbleCode({
    sigma ~ dunif(0,10)
    lam0 ~ dunif(0,5)
    psi ~ dbeta(1,1)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        s[i,1] ~ dunif(xlim[1], xlim[2])
        s[i,2] ~ dunif(ylim[1], ylim[2])
        for(j in 1:J) {# Number of traps
            dist[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
            lam[i,j] <- lam0*exp(-dist[i,j]/(2*sigma^2))*z[i]
        }
    }
    for(j in 1:J){
        bigLambda[j] <- sum(lam[1:M,j])
        for(k in 1:K) {
            n[j,k] ~ dpois(bigLambda[j])
        }
    }
    N <- sum(z[1:M])
})


M<-200

constants <- list(M = M, K=K, J=J)
n1<-apply(n,1,sum)
data<-list(n=n, X=X, xlim=xlim, ylim=ylim)
s<-cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
z<-rep(1,M)
inits <- list (sigma=0.5, lam0=0.1, s=s, z=z)

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE) ## check=FALSE is faster

mcmcspec<-configureMCMC(Rmodel, print=TRUE, monitors = c("N", "lam0", "psi", "sigma"))
scrMCMC <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel) 
CscrMCMC <- compileNimble(scrMCMC, project = Rmodel)

## It's pretty slow to run so I just want an execution time first for a small sample
t1_100 <- system.time(CscrMCMC$run(100))

## And now I want a decent sample
t1_20k <- system.time(CscrMCMC$run(20000))
t1_20k_samples <- as.matrix(CscrMCMC$mvSamples)
save(t1_100, t1_20k, t1_20k_samples, file = "case1results.Rdata")


## checking dcat distribution for Chris

library(nimble)

code <- nimbleCode({
    y ~ dcat((p[1:2,1:2] %*% ones[1:2])[1:2,1])
})
constants <- list(ones = c(1,1))
data <- list()
inits <- list(y = 2, p = diag(c(.5,.5)))

Rmodel <- nimbleModel(code, constants, data, inits)




library(nimble)
Rmodel <- readBUGSmodel('rats', dir = '~/GitHub/nimble/nimble/packages/nimble/inst/classic-bugs/vol1/rats/')

customSpec <- configureMCMC(Rmodel)
customSpec$getSamplers()
customSpec$removeSamplers(1:65)
customSpec$getSamplers()
customSpec$addSampler(target = 'alpha.c', type = 'slice')
customSpec$addSampler(target = 'beta.c', type = 'slice')
customSpec$addSampler(target = 'tau.c', type = 'slice')
customSpec$addSampler(target = 'tau.alpha', type = 'slice')
customSpec$addSampler(target = 'tau.beta', type = 'slice')
customSpec$getSamplers()
alphaNames <- Rmodel$getDependencies('tau.alpha', self = FALSE, stochOnly = TRUE)
alphaNames
length(alphaNames)  ## 30
for(aN in alphaNames) customSpec$addSampler(target = aN, type = 'slice')
customSpec$getSamplers()
betaNames <- Rmodel$getDependencies('tau.beta', self = FALSE, stochOnly = TRUE)
betaNames
length(betaNames)  ## 30
for(bN in betaNames) customSpec$addSampler(target = bN, type = 'slice')
customSpec$getSamplers()
customSpec


## making a simple MCMC traceplot, and density plot

n <- 1000
x <- rnorm(n, 3, 1)
dev.new(width=4, height=3)
plot(1:n, x, type='l')
dev.copy2pdf(file = '~/Downloads/traceplot.pdf')
dev.new(width=4, height=3)
plot(density(x))
dev.copy2pdf(file = '~/Downloads/densityplot.pdf')


## make sure I get the empirical covariance right

x <- matrix(c(1,2,1,2,1,4,6,5,7,4), 5, 2)
timesRan <- 5
statSums <- apply(x, 2, sum)
statSums <- t(statSums)
statProds <- matrix(0, 2, 2)
for(i in 1:timesRan) statProds <- statProds + t(t(x[i,])) %*% t(x[i,])
statSums
statProds
(statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
cov(x)
colMeans <- apply(x, 2, mean)
xCentered <- x - rep(colMeans, each=timesRan)
xCentered
(t(xCentered) %*% xCentered) / (timesRan-1)

set.seed(0)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 4
mu <- 1:d
ch <- array(rnorm(d^2), c(d,d))
Sigma <- ch %*% t(ch)
print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$addSampler('x', 'RW_block_NEW', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(1000)
samples <- as.matrix(Cmcmc$mvSamples)
cov(samples)
           x[1]       x[2]     x[3]       x[4]
x[1]  2.0934846 -0.1465907 1.107342  1.1252451
x[2] -0.1465907  4.1696671 1.644427 -0.5950559
x[3]  1.1073420  1.6444267 2.581730  1.5097074
x[4]  1.1252451 -0.5950559 1.509707  1.8741393
samples[c(1,250,500,750,1000),]
           x[1]      x[2]       x[3]     x[4]
[1,]  1.0000000  1.000000  1.0000000 1.000000
[2,]  0.2162332  2.260594  1.9804354 2.913190
[3,]  1.8118014  3.966497  2.7563534 2.552825
[4,] -1.6015689 -1.004725 -0.9864974 1.065231
[5,]  0.8037007  2.496745  1.2659404 2.208456






## looking into RStudio hanging problem (once again... Nov 2015)

library(nimble)
code <- nimbleCode({
  a ~ dnorm(0, 1)
  sd ~ dunif(0, 100)
  for(i in 1:N) {
    b[i] ~ dnorm(a, sd = sd)
  }
})
N <- 100
constants <- list(N = N)
data <- list()
inits <- list(a = 0, sd = 1, b = rep(0, N))


Rmodel <- nimbleModel(code, constants, data, inits)


## NEED TO FIND THE PROBLEM IN MODEL CHECKING

## EXAMPLE 1:
rm(list = ls())
library(nimble)
classicDyesModel <- readBUGSmodel('dyes', dir = getBUGSexampleDir('dyes'))

## EXAMPLE 2:
set.seed(0)
library(nimble)
n_obs <- 50
### Parameters
mu <- 2
rho_a <- 1
rho_b <- 0.1
monitored_param <- c('mu', 'rho_a', 'rho_b', 'ngis')
SV =  rgamma(n=n_obs, 1, 1)
### Data generation
lambda <-  exp(mu)*SV
ngis  <-  rpois(n=n_obs, lambda=lambda )
Y <- rgamma(n=n_obs, ngis*rho_a, rho_b)
### Remove zero values from ngis
ngis <- ngis[which(ngis>0)]
abs <- which(Y==0)
n_abs <- length(abs)
pres <- which(Y>0)
n_pres <- length(pres)

# Install from source for OS X, Linux, or Windows:
#install.packages("nimble", repos = "http://r-nimble.org", type = "source")
M1_CPG_nimbleCode <- nimbleCode(
    {
        ## Covariates
        for (s in 1:n_obs){
            lambda[s] <- exp(mu)*SV[s]
        }
        ## Observation model
        ## Strictly positive Observation
        for ( s in 1:n_pres){
            ngis[s] ~ T(dpois(lambda[pres[s]]), 1, )  
            Y[pres[s]] ~ dgamma(rho_a * ngis[s], rho_b)
        }
        ## Zero Observation
        for ( s in 1:n_abs){
            proba[s] <- 1 - exp(-lambda[abs[s]])
            Y[abs[s]] ~ dbern(proba[s])
        }
        ## Prior
        rho_a ~ dgamma(0.01, 0.01)
        rho_b ~ dgamma(0.01, 0.01)
        mu ~ dnorm(0, 0.01)
    })

constants_list <- list(n_obs = length(Y), n_abs=n_abs, n_pres=n_pres, abs=abs, pres=pres)
nimble_data <- list(Y = Y, SV=SV)
## Initial values
inits_fc <- function(chain_id = 1){
  mu <- rnorm(1,0,1)
  rho_a <- runif(1,0,10)
  rho_b <- runif(1,0,10)
  ngis <- round(runif(n_pres, 2, 5), digits=0)
  list(mu=mu, rho_a=rho_a, rho_b=rho_b ,ngis=ngis)
}
init_ll <- lapply(1, function(id) inits_fc(chain_id = id))

Rmodel <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll[[1]], check = TRUE)

calculate(Rmodel, 'ngis[1]')

Rmodel$ngis

check <- TRUE
check <- FALSE

md <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll, returnDef = TRUE)

M1_CPG_nimble <- md$newModel(data =nimble_data, inits = init_ll, check = check)

debug(md$newModel)
md$newModel(data =nimble_data, inits = init_ll, check = TRUE)

debug(model$check)





## random stuff for CR paper....

rm(list=ls())
library(nimble)
##source('~/GitHub/legacy/dipper/dipperCode.R')

Rmodel <- nimbleModel(code_dipper, constants, data, inits)
lnodes <- Rmodel$expandNodeNames('x')
length(lnodes) - nind

rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)

con <- orchidDHMM$constants
con$f
con$nind
#### actually 250 individuals in original dataset!!!

con$f
12 - con$f
sum(12 - con$f)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
or <- results$df[results$df$model == 'orchid', ]
or
or[or$mcmc %in% c('jags', 'nimbleDHMM2'), ]
or1 <- or[or$param == 's[1]', ]
data.frame(mcmc = or1$mcmc, runtime_minutes = or1$timing/60)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
go <- results$df[results$df$model == 'goose', ]
go
go1 <- go[go$param == 'p[1]', ]
data.frame(mcmc = go1$mcmc, runtime_minutes = go1$timing/60, runtime_hours = go1$timing/60/60)


rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(gooseDHMM$code, gooseDHMM$constants, gooseDHMM$data, gooseDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)



## testing the MVN conjugacy test from test-mcmc.R
library(nimble)
set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

##### not currently working - see Perry's email of ~ 10/6/14
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3]), Q[1:3, 1:3])
## })

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})

## Simplest version of model w/o 'a' and 'B'
## a = rep(0,3)
## B = diag(rep(1,3))
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(mu[1:3], Q[1:3, 1:3])
## })


mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)

muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
set.seed(0)
Rmcmc$run(5)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0)
Cmcmc$run(5)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples

Rsamples - Csamples

Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
muMeanTrue

apply(samples, 2, mean) - muMeanTrue


source('~/GitHub/nimble/nimble/packages/nimble/inst/tests/test_utils.R')
##test_mcmc
test_mcmc(model = code, name = 'two-level multivariate normal', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))



## working the basketball "streaks" problem emailed out by Aaron Strauss
p <- 0.1
ns <- 1
run <- function(p, ns) {
    out <- numeric()
    n <- 0
    while(n < ns) {
        nmade <- 1
        while(rbinom(1, 1, p) == 1) nmade <- nmade + 1
        out <- c(out, nmade)
        n <- length(out)
    }
    return(out)
}

pOfStreaks <- function(p, streaks, niter=1000) {
    pMissing <- missing(p)
    ns <- length(streaks)
    suc <- 0
    for(i in 1:niter) {
        if(pMissing) p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) suc <- suc+1
    }
    prob <- suc/niter
    return(prob)
}

pOfStreaks(p=0.2, streaks=4, niter=100000)
0.2^3*(.8)

p <- 0.7
n <- 4
m <- 7
pOfStreaks(p=p, streaks=c(n,m), niter=1000000)
p^(n+m-2)*(1-p)^2

n <- 2
m <- 1
pOfStreaks(streaks=c(n,m), niter=1000000)
1/12

pOfPlessthanPROB <- function(streaks, PROB, niter=1000) {
    ns <- length(streaks)
    suc <- 0
    ranSoFar <- 0
    while(ranSoFar < niter) {
        p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) {
            if(p < PROB) suc <- suc+1
            ranSoFar <- ranSoFar+1
        }
    }
    prob <- suc/niter
    return(prob)
}

pOfPlessthanPROB(streaks=c(2,1), PROB=0.5, niter=20000)
11/16

## continusing the "sreaks" basketball problem, now in NIMBLE

library(nimble)

code <- nimbleCode({
    ##p ~ dunif(0, 1)
    p ~ dbeta(1, 1)
    y[1] ~ dnegbin(prob=p, size=1)
    y[2] ~ dnegbin(prob=p, size=1)
})
constants <- list()
data <- list(y = 0:1)  ## the 'streak' history of {2,1} in the original problem
inits <- list(p = 0.5) ## translates into c(1,0) in terms of the negative-binomial distribution

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
#Rmcmc$samplerFunctions$contentsList[[1]]$run
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
debug(Rmcmc$run)
debug(samplerFunctions[[1]]$run)
Rmcmc$run(1)
set.seed(0)
Cmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
Rsamples
Csamples
Rsamples - Csamples

apply(samples, 2, mean)

mean(samples[,1] > 0.5)
5/16

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_noConj'),
    summaryStats = c('mean', 'median', 'sd', 'function(x) mean(x>0.5)'),
    calculateEfficiency = TRUE,
    makePlot = TRUE,
    savePlot = FALSE,
    niter = 100000
)

out$timing
out$summary
11/16


## creating a really small, blank figure, for use in captions in CR-MCMC paper
dev.new(width=3, height=2)
plot(1, 1)
dev.copy2pdf(file='~/GitHub/nimble/nimblePapers/CR-MCMC/blankFigure.pdf')

## testing the problem with naming a variable 'ans' in a nimbleFunction
library(nimble)

rcFunction1 <- nimbleFunction(
    run = function() {
        declare(ans, double(1, 5))
        for(i in 1:5) ans[i] <- i
        print(ans)
    })

rcFunction1()

rcFunction2 <- nimbleFunction(
    run = function() {
        declare(functionAsList, double(1, 5))
        for(i in 1:5) functionAsList[i] <- i
        print(functionAsList)
    })

rcFunction2()

code <- nimbleCode({
    mu ~ dnorm(0, 1)
    mu2 <- timesTwo(mu)
})

Rmodel <- nimbleModel(code, inits = list(mu=1), check=FALSE)

##undebug(calculate)
##debug(nimble:::rCalcNodes)
##Rmodel$nodes$mu2$calculate
##Rmodel$nodes$mu2$calculate()
##Rmodel$nodes$mu2$simulate
##Rmodel$nodes$mu2$simulate()

calculate(Rmodel)
Cmodel <- compileNimble(Rmodel)

Rmodel$mu
Rmodel$mu2
Cmodel$mu
Cmodel$mu2

Rmodel$mu <- 2
calculate(Rmodel)
Rmodel$mu2


## finding the posterior pairwise correlations of
## the blocked parameters for the orchid model
## (multistate CR-MCMC, userDistMCMC)

library(nimble)
rm(list=ls())
source('~/GitHub/userDistMCMC/defs.R')
load('~/GitHub/userDistMCMC/models.RData')

Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$removeSamplers(numeric(0))
spec$removeSamplers(c('b[1:2]','c[1:2]'))
spec$addSampler(c('b[1]','b[2]'), 'RW_block')
spec$addSampler(c('c[1]','c[2]'), 'RW_block')
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
## n <- 10000
n <- 50000
system.time(Cmcmc$run(n)) / 60
samples <- as.matrix(Cmcmc$mvSamples)
samples2 <- samples[(n/2+1):n, ]   ## truncate to second half of samples
apply(samples2, 2, mean)
Cov <- cov(samples2)
Cov
Cor <- cov2cor(Cov)
Cor
sort(as.numeric(abs(Cor)))
Cor['b[1]', 'b[2]']
## [1] 0.9866091
Cor['c[1]', 'c[2]']
## [1] 0.9675666


## printing the numeric performance results for CR-MCMC paper, or userDistMCMC paper
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')
source('~/GitHub/userDistMCMC/defs.R')
options(digits = 3)
df <- results$df
for(mod in unique(df$model)) {
    dfmod <- df[df$model==mod, ]
    message('******************************************')
    message(mod, ' model')
    message('******************************************')
    for(mc in unique(dfmod$mcmc)) {
        dfmodmcmc <- dfmod[dfmod$mcmc==mc, ]
        message(mc, ' MCMC:')
        message('minimum: ',  min(dfmodmcmc$Efficiency))
        message('mean: ',    mean(dfmodmcmc$Efficiency))
    }
    message('******************************************')
}

df[df$model=='orchid',]
df[df$model=='goose',]

results$check()
results$quickplot()


## problem with multivariate conjugate samplers... ?
## fixed.  problem was in vector demoting in the conjugacy definition
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
mu <- rep(0, d)
Sigma <- diag(d)
constants <- list(d=d, mu=mu, Sigma=Sigma)
Rmodel <- nimbleModel(code, constants, data=list(y=rep(0,d)), inits=list(x=rep(0, d)))
Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
solve(cov(samples))

library(nimble)
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dnorm
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dbeta
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dmnorm



## testing different priors for rho in package gpmanagement
library(gpmanagement)
library(ggplot2)
example_data <- function(p, x0, Tobs){
    f <- function (x, h){
        sapply(x, function(x) {
            x <- pmax(0, x - h)
            x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
        })
    }
    sigma_g <- 0.1
    pdfn <- function(x, mu, sigma = sigma_g){
        dlnorm(x, log(mu), sdlog = sigma)
    }
    z_g <- function() rlnorm(1, 0, sigma_g)
    x <- numeric(Tobs)
    x[1] <- x0
    for(t in 1:(Tobs-1))
        x[t+1] = z_g() * f(x[t], h=0)
    obs <- data.frame(x = c(0, 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                      y = c(0, 
                          x[2:Tobs]))
    obs
}
myfit <- function(xObs, yObs, priors = NULL, niter = 100000){    ## default value for priors
  xPred <- seq(0, 1.1 * max(xObs), length = 50)
  fit <- gp_setup(xObs, yObs, xPred, priors)     ## included priors argument
  Cmcmc <- fit$Cmcmc 
  Cpred <- fit$Cpred
  Cmodel <- fit$Cmodel
  system.time(Cmcmc$run(niter))
  samples <- as.matrix(Cmcmc$mvSamples)
  system.time(Cpred$run(samples))
  E <- Cpred$getE()
  V <- sqrt(diag(Cpred$getC()))
  list(samples = tidyr::gather(as.data.frame(samples)),
       pred = data.frame(x = xPred, y = E, ymin = E - V, ymax = E + V))
}

densPlot <- function(fit, paramName, paramPrior) {
        dev.new()
        samp <- fit$samples[fit$samples$key==paramName, 'value']
        plot(density(samp), main=paste0('black posterior, red prior: ', paramName, ' ~ ', deparse(paramPrior)))
        xs <- seq(range(samp)[1], range(samp)[2], length.out = 500)
        yCall <- as.call(c(list(as.name(paramPrior[[1]])), list(quote(xs)), as.list(paramPrior[-1])))
        ys <- eval(yCall, envir = environment())
        lines(xs, ys, col='red')
}

## Carrying capcacity at 10, critical point at 5, x0=6
## 100, 50, x0=60
plotPred <- function(ccap, crit, x0, Tobs=40, plotData=FALSE, rhoPrior, plotRhoDensity=TRUE, niter=100000, plotSigGPDensity=FALSE, plotSigOEDensity=FALSE, plotPredictionIntervals=TRUE) {
    set.seed(0)
    obs <- example_data(p=c(2,ccap,crit), x0=x0, Tobs=Tobs)
    if(plotData) {dev.new(); qplot(seq_along(obs$x[-1]), obs$x[-1]) + geom_line()}
    priors <- expression({
        rho ~ X
        sigGP ~ dunif(0, 1e5)
        sigOE ~ dunif(0, 1e5)
    })
    if(missing(rhoPrior)) stop('must specify a prior for rho')
    rhoPrior[-1] <- lapply(rhoPrior[-1], function(param)
        eval(eval(substitute(substitute(PARAM, list(xObs=obs$x[-1])), list(PARAM=param)))))
    print(rhoPrior)
    priors[[1]][[2]][[3]] <- rhoPrior
    fit <- myfit(obs$x, obs$y, priors = priors, niter=niter)
    if(plotRhoDensity) densPlot(fit, 'rho', rhoPrior)
    if(plotSigGPDensity) densPlot(fit, 'sigGP', priors[[1]][[3]][[3]])
    if(plotSigOEDensity) densPlot(fit, 'sigOE', priors[[1]][[4]][[3]])
    if(plotPredictionIntervals) {
        dev.new()
        print(ggplot2::ggplot(fit$pred) + 
                  geom_ribbon(aes(x = x, y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
                      geom_line(aes(x = x, y = y), size = 1) + 
                          geom_point(data = obs, aes(x, y)) +
                              coord_cartesian(xlim = range(c(obs$x, fit$pred$x)), 
                                              ylim = range(c(obs$y, fit$pred$y))) +
                                                  ggtitle(paste0('rho ~ ', deparse(rhoPrior))))
    }
}

rp <- quote(dunif(0, diff(range(xObs))/sqrt(6)/2))

plotPred(rhoPrior = rp, ccap= 10, crit= 5, x0= 6)
plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60)

plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60, plotSigGPDensity=TRUE, plotSigOEDensity=TRUE)

plotPred(rhoPrior = quote(dunif(0,200)), ccap=10, crit=5, x0=6)

for(uMax in 5:6) {
    rp <- substitute(dunif(0, UMAX), list(UMAX = as.numeric(uMax)))
    plotPred(rhoPrior = rp, ccap=10, crit=5, x0=6, niter=10000)
}



## testing ggplot() inside a loop??
library(ggplot2)
df <- data.frame(a=1:2, b=3:4)
for(i in 1:3) {
    dev.new()
    ## need explict print(...) to make ggplot() plot, when inside a loop!!!
    print(ggplot(df, aes(a, b)) + geom_line() + ggtitle(i))
}



## testing use of substitute...
## which I used to be more comfortable with!

f <- function(a) {
    b <- substitute(a)
    print(b)
    print(class(b))
    vv <- eval(b)
    print(vv)
    print(class(vv))
}

f(dnorm(0,1))
f(4+5)

x <- f(quote(dnorm(0,1)))
x
class(x)


## testing NIMBLE compiler to handle pi
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        a <- pi
        print(a)
        b <- 2 * pi
        print(b)
    }
)
Rnf <- nfDef()
Rnf$run()
Cnf <- compileNimble(Rnf)    ## ERROR
Cnf$run()


## example of trying to use NIMBLE functions in another package
## this uses my minimalist nimtest package; same thing I sent to Duncan
library(devtools)
install_github('danielturek/nimtest')
library(nimtest)
library(nimble)
## First Test:
## a nimbleFunction defined as a package function
nfDefinition
Rnf <- nfDefinition()
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
## Second Test:
## A package function which internally defines and returns a nimbleFunction definition
## This one gives ugly looking warnings
nextFunction
def <- nextFunction()        ## WARNINGS
Rnf <- def()
Cnf <- compileNimble(Rnf)    ## WARNINGS
Rnf$run()
Cnf$run()

## alternate.  from shell:
## R CMD BUILD nimtest
## R CMD install nimtest_0.0.0.9000.tar.gz
library(nimtest)
library(nimble)
nfDefinition    ## finds it -- it created the NIMBLE function at build time
Rnf <- nfDefinition()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
nextFunction
nextFunction()
def <- nextFunction()
Rnf <- def()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()


## testing out the new elliptical slice sampler (ess)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
##mu <- rep(0, d)
mu <- 1:d
##Sigma <- diag(d);  print(Sigma)
ch <- array(c(1, 0.7, 0, 2), c(2,2));  Sigma <- ch %*% t(ch);  print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print=FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0); Cmcmc$run(20000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
Sigma


## conjugate MVN-MVN
set.seed(1)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
})
d <- 3
mu_x <- rnorm(d)
temp <- array(rnorm(d^2), c(d,d))
prec_x <- solve(temp %*% t(temp))
temp <- array(rnorm(d^2), c(d,d))
prec_y <- solve(temp %*% t(temp))
y <- rnorm(d)
constants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
data <- list(y = y)
inits <- list(x = rep(0, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

## block sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)

## ess sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)



## testing elliptical slice sampler on 'spatial' model
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_spatial.RData')

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('block', 'ess'),
    MCMCdefs = list(
        block = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'RW_block', print=FALSE)
            spec
        }),
        ess = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'ess', print=FALSE)
            spec
        })),
    niter = 10000,
    burnin = 2000,
    monitors = c('mu', 'sigma', 'rho', 'g'),
    makePlot = FALSE
)

out$timing
out$summary[,,1:5]    ## WRONG again...
cov(t(out$samples['block', -(1:3), ])) - cov(t(out$samples['ess', -(1:3), ]))






## Taken from Carl's lab notebook:
## Gpdp Via MDPtoolbox Cont; 30 Jul 2015
## http://www.carlboettiger.info/2015/07/30/gpdp-via-mdptoolbox-cont.html
## (Markov Decision Processes (MDP) toolbox)
##knitr::opts_chunk$set(eval=FALSE)
##devtools::install_github("cboettig/gpmanagement@b3b765cbceb51c9b0b8cb2724e395353ec365df9")
library("MDPtoolbox")
library("gpmanagement")
library("tidyr")
library("dplyr")
library("ggplot2")

## True model
p <- c(2, 100, 50)
f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

sigma_g <- 0.1
pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

z_g <- function() rlnorm(1, 0, sigma_g)
set.seed(0)
Tobs <- 40

x <- numeric(Tobs)
x[1] <- 60
for(t in 1:(Tobs-1)) x[t+1] = z_g() * f(x[t], h=0)
obs <- data.frame(x = c(0, pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                  y = c(0, x[2:Tobs]))
xObs <- obs$x
yObs <- obs$y
xPred <- seq(0, 1.1 * max(xObs), length = 50)
qplot(seq_along(x), x) + geom_line()

## Now the GP estimation from NIMBLE. Let’s emphasize shorter length-scales with the prior to compare:

## having this work with output from `nimbleCode` might be more natural than `expression`
priors <- expression({
  rho ~ dgamma(1, 1)
  sigGP ~ dunif(0, 1e5)
  sigOE ~ dunif(0, 1e5)
})

fit <- gp_setup(xObs, yObs, xPred)

Cmcmc <- fit$Cmcmc 
Cpred <- fit$Cpred
Cmodel <- fit$Cmodel
system.time(Cmcmc$run(100000))
samples <- as.matrix(Cmcmc$mvSamples)
## basic sanity check
testthat::expect_identical(Cmodel$getNodeNames(topOnly = TRUE), colnames(samples))

## predict from GP model using posterior MCMC samples
system.time(Cpred$run(samples))

## Posteriors
samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
df <- tidyr::gather(samples)

ggplot(df) + 
  geom_density(aes(value)) + 
  facet_wrap(~key, scale='free')

## extract predictions: E and C
E <- Cpred$getE()
C <- Cpred$getC()

obs <- data.frame(x = xObs, y = yObs)
pred <- data.frame(x = xPred, y = E, ymin = E - sqrt(diag(C)), ymax = E + sqrt(diag(C)))

ggplot2::ggplot(pred) + 
  geom_ribbon(aes(x = x,y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
  geom_line(aes(x = x, y = y), size=1) + 
  geom_point(data = obs, aes(x,y)) +
  coord_cartesian(xlim = range(c(xObs, xPred)), ylim = range(c(yObs,E))) +
  theme_bw()

## Decision theory
states <- xPred # Vector of all possible states
actions <- states # Vector of actions: harvest
## Let’s consider a slight variation of the most trivial utility function:
## one which explicitly adds a cost to completely exhausting the stock
## (or reducing the stock by more than, say 95% in this case.)
## This should be somewhat similar to the impact of no discount rate.

## Utility function
discount = 0.99

#get_utility <- function(x,h) pmin(x,h)
#R <- outer(states, actions, get_utility)

R <- sapply(actions, function(h){
    sapply(states, function(x){
        if(h < x) h else - 1 * max(states)
    })
})

## Implementing policy
z <- function() rlnorm(1, meanlog = 0, sdlog = sigma_g)

simulate_policy <- function(states, actions, policy, f, z, s0, steps = 50, utility = function(s,a) NA, discount = 1){
    s <- numeric(steps)
    a <- numeric(steps)
    u <- numeric(steps)
    s[1] <- s0
    for(t in 1:(steps-1)){
        a[t] <- actions[policy[which.min(abs(states - s[t]))]]
        s[t+1] <- z() * f(s[t], a[t])
        u[t] <- utility(s[t], a[t]) * discount ^ t
    }
    ## Final action determined but not implemented
    a[steps] <- actions[policy[which.min(abs(states - s[t]))]]
    data.frame(time = 1:steps, state = s, action = a, utility = u)
}

## GP model
gp_matrix <- function(states, actions, E, C){
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    sigmas <- sqrt(diag(C))
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            nextpop <- E[k] - actions[i]
            if(nextpop <= 0) {
                transition[k, , i] <- c(1, rep(0, K - 1))
            } else {
                transition[k, , i] <- dnorm(states, nextpop, sigmas[i]) / sum(dnorm(states, nextpop, sigmas[i]))
            }
        }
    }
    transition
}

P_gp <- gp_matrix(states, actions, E, C)
mdp_check(P = P_gp, R = R)

gp <- mdp_value_iteration(P_gp, R, discount = discount, epsilon = 0.00001, max_iter = 5e3, V0 = numeric(length(states)))

plot(states, states - actions[gp$policy],  xlab="Population size", ylab="Escapement")

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")
## With this amount of data, the gp solution is too cautious, and avoids any exploitation.

## Simulate under the true model

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
## (Average utility is approximate here since it does not include penalty;
## since a function and not a matrix is requred by this function at this time.)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")

P <- transition_matrix(states, actions, f, pdfn)
mdp_check(P = P, R = R)

mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Note that the altered award structure has almost no effect on the optimal policy
## given the true model, other than to avoid harvesting directly to zero even when
## the stock cannot persist, due to the explicit penalty for doing so.

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, mdp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")





## Taken from Carl's lab notebook:
## MDPtoolbox Allen Model; 29 Jul 2015
## http://www.carlboettiger.info/2015/07/29/mdptoolbox-allen-model.html
## (Markov Decision Processes (MDP) toolbox)
library("MDPtoolbox", quietly = TRUE)
library("ggplot2", quietly = TRUE)
K <- 10 # state space limit
states <- 0:K # Vector of all possible states
actions <- states # Vector of actions: harvest

sigma_g = 0.1
p <- c(2, 15, 5)

f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

## Utility function
discount = 0.95
get_utility <- function(x,h) {
    pmin(x,h)
}

R <- outer(states, actions, get_utility)

transition_matrix <- function(states, actions, f, pdfn){
    ## Initialize
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            ## Calculate the transition state at the next step, given the 
            ## current state k and action i (harvest H[i])
            nextpop <- f(states[k], actions[i])
            ## Population always extinct if this is negative.
            ## since multiplicitive shock z_t * f(n) < 0 for all f(n) < 0
            if(nextpop <= 0)
                transition[k, , i] <- c(1, rep(0, length(states) - 1))
            ## Implement demographic stochasticity 
            else {
                ## Cts distributions need long-tailed denominator as normalizing factor:
                fine_states <- seq(min(states), 10 * max(states), by = states[2] - states[1])
                N <- sum(pdfn(fine_states, nextpop))  
                transition[k, , i] <-pdfn(states, nextpop) / N
                ## We need to correct this density for the final capping state ("Pile on boundary")
                ## (discrete or cts case)
                ## this can be a tiny but negative value due to floating-point errors.
                ## so we take max(v,0) to avoid
                transition[k, K, i] <- max(1 - sum(transition[k, -K, i]), 0)
            }
        } 
    }
    transition
}

P <- transition_matrix(states, actions, f, pdfn)
apply(P, c(1,3), sum)  ## double-check

## Using toolbox
mdp_check(P = P, R = R)
mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Compare to Reed
## From Reed (1979) we know that the optimal solution is a constant-escapement rule
## when the growth function in convex. Note that this condition is violated by the
## growth function with alternative stable states (Allen/Ricker-Allee model),
## resulting in a very different optimal policy:
## f^prime(s^star) = 1/discount
## For growth-rate function f, where discount is the discount factor and s∗ the stock size
## for the constant escapement. Analytic solutions are clearly possible for certain
## growth functions, but here I’ve just implemented a generic numerical solution.
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)
D <- actions[mdp$policy]
# The difference between Bellman and the analytical solution is small:
plot(states, states - D,  xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)




## Carl's examples using MDPtoolbox
## MDPtoolbox Ex 2; 14 Jul 2015
## (Markov Decision Processes (MDP) toolbox)
## http://www.carlboettiger.info/2015/07/14/mdptoolbox-ex-2.html
## Adapted from Marescot et al. appendix 5, to Reed optimal control problem, including direct comparison against (semi) analytic optimum.

## step 1: define objectives
## This is a conceptual step which does not require coding
## step 2: define states
##K <- 150 # state space limit
K <- 10 # state space limit
states <- 0:K # Vector of all possible states

## step 3: define control actions
## Vector of actions: harvest
H <- states

## step 4: define dynamic model (with demographic parameters)
##p <- c(6,0.05)
p <- c(3, 0.2)
f <- function(x, h){
  A <- p[1] 
  B <- p[2] 
  s <- pmax(x-h, 0)
  A * s/(1 + B * s)   #### Bellman equation
}
sigma_g = 0.1

## step 5: define utility
## Utility function
get_utility <- function(x,h) {
    pmin(x,h)
}

## step 6: solve bellman equation with value iteration
## Initialize transition matrix
transition <- array(0, dim = c(length(states), length(states), length(H)))
## Initialize utility matrix
utility <- array(0, dim = c(length(states), length(H)))
## Fill in the transition and utility matrix
## Loop on all states
for (k in 0:K) {
    ## Loop on all actions
    for (i in 1:length(H)) {
        ## Calculate the transition state at the next step, given the 
        ## current state k and the harvest H[i]
        nextpop <- f(k, H[i])
        if(nextpop <= 0)
            transition[k+1, , i] <- c(1, rep(0, length(states) - 1))
        ## Implement demographic stochasticity by drawing probability from a density function
        else {
            ## We need to correct this density for the final capping state ("Pile on boundary").
            ## For discrete probability distribution,
            ## this is easy if `states` includes all possible
            ## discrete states below the capping state.
            ## (e.g. all non-negative integers less than K).  
            ## For a continuous distribution, this is more problematic
            ## as we have to first normalize the densities.
            ## EDIT: this can be negative, due to floating-point errors.
            ## so we take max(v,0) to avoid.
            ## Get long-tailed denominator as normalizing factor (continuous distributions only):
            fine_states <- seq(min(states), 10 * max(states), by = states[2]-states[1])
            N <- sum(dlnorm(fine_states, log(nextpop), sdlog = sigma_g))
            transition[k+1, , i] <- dlnorm(states, log(nextpop), sdlog = sigma_g) / N
            ## We need to correct this density for the final capping state ("Pile on boundary")
            transition[k+1, K+1, i] <- max(1 - sum(transition[k+1, -(K+1), i]), 0)
        }
        ## Compute utility
        utility[k+1, i] <- get_utility(k, H[i])
    } # end of action loop
} # end of state loop

## Solution calculated explicitly:
## The backward iteration consists in storing action values in the vector
## Vt which is the maximum of utility plus the future action values for all
## possible next states. Knowing the final action values, we can then backwardly
## reset the next action value Vtplus to the new value Vt.  We start The backward
## iteration at time T-1 since we already defined the action value at Tmax.

## Discount factor
discount <- 0.95
## Action value vector at tmax
Vtmax <- numeric(length(states))
## Action value vector at t and t+1
Vt <- numeric(length(states))
Vtplus <- numeric(length(states))
## Optimal policy vector
D <- numeric(length(states))
## Time horizon
##Tmax <- 150
Tmax <- 4

for (t in (Tmax - 1):1) {
    ## We define a matrix Q that stores the updated action values for 
    ## all states (rows)
    ## actions (columns)
    Q <- array(0, dim = c(length(states), length(H)))
    for (i in 1:length(H)) {
        ## For each harvest rate we fill for all states values (row) 
        ## the ith column (Action) of matrix Q
        ## The utility of the ith action recorded for all states is 
        ## added to the product of the transition matrix of the ith 
        ## action by the action value of all states 
        Q[,i] <- utility[, i] + discount * (transition[,,i] %*% Vtplus)
    } # end of the harvest loop
    ## Find the optimal action value at time t is the maximum of Q
    Vt <- apply(Q, 1, max)
    ## After filling vector Vt of the action values at all states, we 
    ## update the vector Vt+1 to Vt and we go to the next step standing 
    ## for previous time t-1, since we iterate backward
    Vtplus <- Vt
} # end of the time loop

## Find optimal action for each state
for (k in 0:K) {
    ## We look for each state which column of Q corresponds to the 
    ## maximum of the last updated value 
    ## of Vt (the one at time t + 1). If the index vector is longer than 1 
    ## (if there is more than one optimal value we chose the minimum 
    ## harvest rate)
    D[k + 1] <- H[(min(which(Q[k + 1, ] == Vt[k + 1])))]
}

## plot solution
plot(states, states - D, xlab="Population size", ylab="Escapement")

## proof of optimality: compare with analytical solution
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)

## The difference between Bellman equation solution and the analytical solution is small:
plot(states, states - D, xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)

## Using MDPtoolbox
library('MDPtoolbox')
mdp_check(P = transition, R = utility)
out <- mdp_value_iteration(transition, utility, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = Vtmax)
plot(states, states - D, xlab="Population size", ylab="Escapement")
lines(states, states - H[out$policy], col="red", lty=2)




## testing using new NIMBLE option: verbose
nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleOptions('verbose')

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)


## testing the new includeEfficiency option for MCMCsuite
## calculates N, ESS, and Eff

library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
load('~/GitHub/legacy/autoBlock/data/model_ice.RData')

out <- MCMCsuite(code, constants, data, inits,
##                 MCMCs = c('nimble', 'nimble_slice', 'autoBlock'),
                 MCMCs = c('nimble', 'nimble_slice'),
                 makePlot = FALSE,
                 calculateEfficiency = TRUE
                 )

out$timing

out$summary

apply(out$summary[, 'efficiency', ], 1, min)
apply(out$summary[, 'efficiency', ], 1, mean)

out$efficiency


## running OpenBUGS installation from OSX (using wine)

##library(R2OpenBUGS)
library(R2WinBUGS)
library(BRugs)

## schools data in the R2OpenBUGS library
data(schools)
schools

## define the model
nummodel <- function() {
    for (j in 1:J) {
        y[j] ~ dnorm(theta[j], tau.y[j])
        theta[j] ~ dnorm(mu.theta, tau.theta)
        tau.y[j] <- pow(sigma.y[j], -2)
    }
    mu.theta ~ dnorm(0, 1E-6)
    tau.theta <- pow(sigma.theta, -2)
    sigma.theta ~ dunif(0, 1000)
}

## write the model code out to a file
write.model(nummodel, 'nummodel.txt')
model.file1 = paste(getwd(), 'nummodel.txt', sep='/')
## and let's take a look:
file.show('nummodel.txt')

## prepare the data for input into OpenBUGS
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list(J = J, y = y, sigma.y = sigma.y)

## initialization of variables
inits <- function() {
    list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100), sigma.theta = runif(1, 0, 100))
}

## set the WINE working directory and the directory to OpenBUGS
## change the OpenBUGS.exe location as necessary
WINE <- '/opt/local/bin/wine'
WINEPATH <- '/opt/local/bin/winepath'
##OpenBUGS.pgm <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe'
bugs.directory <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323'

## these are the parameters to save
parameters = c('theta', 'mu.theta', 'sigma.theta')

## run the model
schools.sim <- bugs(data, inits, model.file = model.file1, parameters=parameters,
                    n.chains = 3, n.iter = 1000,
                    ##OpenBUGS.pgm = OpenBUGS.pgm,   ## syntax for R2OpenBUGS::bugs
                    program = 'OpenBUGS',
                    bugs.directory = bugs.directory,
                    WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE)

## R will pause
## When model is complete a prompt will reappear
print(schools.sim)


## debugging commands in R
## n: next
## c: continue
## s: step into function call
## f: step out of current function call



library(nimble)
modelName <- 'rats'
rats <- readBUGSmodel(model = modelName, dir = getBUGSexampleDir(modelName), returnModelComponentsOnly = TRUE)

code <- rats$model
constants <- rats$data[c('N','T','x')]
data <- rats$data[c('Y')]
inits <- rats$inits

Rmodel <- nimbleModel(code, constants, data, inits, dimensions = rats$dims)

spec <- configureMCMC(Rmodel)
spec$getSamplers()




## creating dendogram figure of autoBlocking for JSM presentation
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
Rmcmc <- buildMCMC(Rmodel, autoBlock = TRUE, makePlots = TRUE)




## working on the conjugacy system, to infer dimensions of 'contributions'
library(nimble)
getDistributionsInfo('dnorm')



## let's manually run all NIMBLE testing
## (to find what's crashing, easier than using travis)
library(nimble)
library(testthat)
setwd('~/GitHub/nimble/nimble/packages/nimble/inst/tests/')
source('test_utils.R')

testFiles <- list.files()
print(testFiles)

##### this listing is out of date!
test_package('nimble', 'copy')      ## pass
test_package('nimble', 'dsl_dists') ## pass
test_package('nimble', 'math')      ## pass
test_package('nimble', 'mcmc')      ## errors/warnings copied below
test_package('nimble', 'meta')      ## pass
test_package('nimble', 'models')    ## 
test_package('nimble', 'trunc')     ## 
test_package('nimble', 'user')      ## 



## fully rebuild NIMBLE user manual

setwd('~/GitHub/nimble/nimble-docs/UserManual/')
library(knitr) 
knit2pdf('NimbleUserManual.Rnw') 
system('open NimbleUserManual.pdf')


## investigating the bug in R M-H sampler

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'RW')
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
niter <- 1000

set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples

d <- Rsamples - Csamples
all(d == 0)



Rmodel$nodes[['a']]$calculateDiff

Rmodel$logProb_a
Cmodel$logProb_a

calculateDiff(Rmodel, 'a')
calculateDiff(Cmodel, 'a')


nimble:::rCalcDiffNodes(Rmodel, 'a')
nimble:::rCalcDiffNodes(Cmodel, 'a')


## testing nimbleOptions() system
## and also the new nimble option: MCMCcontrolDefaultList

library(nimble)

getNimbleOption('MCMCcontrolDefaultList')

environment(sampler_conjugate_dnorm)$nfRefClassDef

getNimbleOption('verifyConjugatePosteriors')
nimble:::setNimbleOption('verifyConjugatePosteriors', TRUE)
getNimbleOption('verifyConjugatePosteriors')
environment(sampler_conjugate_dnorm)$nfRefClassDef
buildConjugateSamplerFunctions()
environment(sampler_conjugate_dnorm)$nfRefClassDef


code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 2)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b = .5, c = 4, d = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

lst <- getNimbleOption('MCMCcontrolDefaultList')
lst
lst$sliceWidth <- 99
lst$scale <- 9999
lst
nimbleOptions(MCMCcontrolDefaultList = lst)
getNimbleOption('MCMCcontrolDefaultList')

spec <- configureMCMC(Rmodel)
spec$getSamplers()


## complete example of using RW_llFunction sampler,
## for the User Manual

library(nimble)

code <- nimbleCode({
    p ~ dunif(0, 1)
    y ~ dbin(p, n)
})

Rmodel <- nimbleModel(code, data = list(y = 3),
                      inits = list(p = 0.5, n = 10))

llFun <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        y <- model$y
        p <- model$p
        n <- model$n
        ll <- lfactorial(n) - lfactorial(y) - lfactorial(n-y)
                 + y*log(p) + (n-y)*log(1-p)
        returnType(double())
        return(ll)
    }
)

RllFun <- llFun(Rmodel)

mcmcspec <- configureMCMC(Rmodel, nodes = NULL)

mcmcspec$addSampler(target = 'p', type = 'RW_llFunction',
    control = list(llFunction = RllFun, includesTarget = FALSE))

Rmcmc <- buildMCMC(spec)





## testing new nimble options system

library(nimble)
nimbleOptions()
ls(nimble:::.nimbleOptions)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b=.5, c=1, d=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## working on dipper models

rm(list=ls())
library(nimble)
trunc <- TRUE
load('~/GitHub/userDistMCMC/dipperData.RData')
## optionally truncate data:
last <- apply(y, 1, function(hist) max(which(hist==1)))
yDHMM <- 2 - y
onesMatrix = matrix(1,nind,k)  ## 'ones' matrix for survival
onesVector = rep(1,nind)       ## 'ones' vector for chi
if(trunc) { ind <- 1:3;   nind<-length(ind);   first<-first[ind];   last<-last[ind];   y<-y[ind,,drop=FALSE];   yDHMM<-yDHMM[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE];   onesMatrix<-onesMatrix[ind,,drop=FALSE];   onesVector<-onesVector[ind] }
code <- quote({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:nind) {
        for(t in (first[i] + 1):last[i]) {
            y[i,t] ~ dbin(p, 1)
            ##onesMatrix[i,t] ~ dbin(phi, 1)
        }
        ##onesVector[i] ~ dbin(chi[last[i]], 1)
    }
    ##chi[k] <- 1
    ##for(t in 1:(k-1)) {
    ##chi[k-t] <- (1-phi) + phi * (1-p) * chi[k-t+1]
    ##}
})
constants <- list(k=k, nind=nind, first=first, last=last)
data      <- list(y=y)##, onesMatrix=onesMatrix, onesVector=onesVector)
inits     <- list(phi=0.6, p=0.9)

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)


Rmodel <- nimbleModel(code, constants, data, inits)



## reproducible example of problem in new igraph package
## emailed this example to Gabor


remove.packages('igraph')
install.packages('igraph')

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
})
constants <- list()
data <- list()
inits <- list(a=0, b=0)

Rmodel <- nimbleModel(code, constants, data, inits)


library(igraph)
graph <- graph.empty()
graph <- add.vertices(graph, 2, name = c('a', 'b'))
graph <- add.edges(graph, c(1, 2))
toposortReturn <- topological.sort(graph, mode = 'out')
print(toposortReturn)
class(toposortReturn)



oldGraphID_2_newGraphID <- sort(newGraphID_2_oldGraphID, index = TRUE)$ix
graph <<- permute.vertices(graph, oldGraphID_2_newGraphID)  # re-label vertices in the graph



## writing stan model file for RSBS
##source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
##install_rstan()

load('~/GitHub/legacy/automated-blocking-examples/data/model_redblue.RData')
constantsAndData <- c(constants, data)

library(rstan)

stan_mod <- stan_model(file = '~/temp/redblue.stan')

stan_out <- sampling(stan_mod, data=constantsAndData, chains=1, iter=10000, thin=1, init=list(inits))

tempArray <- extract(stan_out, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
dimnames(tempArray)[[2]] <- gsub('_', '.', dimnames(tempArray)[[2]])
if(!all(monitorNodesBUGS %in% dimnames(tempArray)[[2]])) {
    missingNames <- setdiff(monitorNodesBUGS, dimnames(tempArray)[[2]])
    warning(paste0('Stan output is missing values for: ', paste0(missingNames,collapse=', ')))
}
samplesArray <- array(0, dim = c(nkeep, length(monitorNodesBUGS)))
dimnames(samplesArray)[[2]] <- monitorNodesBUGS
monitorsWeHave <- intersect(monitorNodesBUGS, dimnames(tempArray)[[2]])
samplesArray[, monitorsWeHave] <- tempArray[(burnin+1):floor(niter/thin), monitorsWeHave, drop=FALSE]
addToOutput('stan', samplesArray, timeResult)

## playing around to build, install, run, and test gpmanagement package

q('no')
library(devtools)
setwd('~/GitHub/forks/gpmanagement')
remove.packages('gpmanagement')
install_github('danielturek/gpmanagement')
##document('gpmanagement')
system('R CMD build gpmanagement')
system('R CMD INSTALL .')

library(nimble)
library(testthat)
library(gpmanagement)

check('.')  ## this gives me the errors 

gp_setup

test_package('gpmanagement')

source('~/GitHub/forks/gpmanagement/tests/testthat/test_gp.R')

##setwd('~/GitHub/forks/gpmanagement/tests/')
##test_check('gpmanagement')

nf <- gpmanagement:::nfd(3)
nf$run()

set.seed(0)
x <- 1:100
y <- sin(x/5) + rnorm(100, 0.1)
ind <- sort(sample(1:100, 40))
xObs <- x[ind]                ## input
yObs <- y[ind]                ## input
xPred <- c(x, 101:120)        ## input

fit <- gp_setup(xObs, yObs, xPred)


## testing of 3-D array arguments, and their automatic demotion to smaller dim arrays

library(nimble)

nf1 <- nimbleFunction(
    run = function(a = double(2), b = double(3)) {
        print('a11: ', a[1,1])
        print('b111: ', b[1,1,1])
        print(dim(a)[1])
        print(dim(a)[2])
        print(dim(b)[1])
        print(dim(b)[2])
        print(dim(b)[3])
    }
)

nf <- nimbleFunction(
    run = function() {
        a <- nimArray(0, 3, 3)
        declare(b, double(3, c(3,3,3)))
        for(i in 1:3)
            for(j in 1:3) {
                a[i,j] <- 10*i + j
                for(k in 1:3)
                    b[i,j,k] <- 100*i+10*j+k
            }
        print(a)
        ##print(b)
        nf1(a[2,1:2], b[,,])
    }
)

nf()
Cnf <- compileNimble(nf)
Cnf()


## testing distribution dDHMM on a multistate example

CdDHMM <- compileNimble(dDHMM)
CrDHMM <- compileNimble(rDHMM)

checkD <- function(ind) {
    y <- yDHMM
    print(y[ind,])
    x <- y[ind, first[ind]:k]
    length <- length(x)
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(dDHMM(x, length, pi, Z, T, condition, 0))
    print(CdDHMM(x, length, pi, Z, T, condition, 0))
    print(exp(dDHMM(x, length, pi, Z, T, condition, 1)))
    print(exp(CdDHMM(x, length, pi, Z, T, condition, 1)))
}
checkD(225)   ## 0.2484
checkD(202)   ## 0.017496
checkD(167)   ## 0.1246882

checkR <- function(length) {
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(rDHMM(1, length, pi, Z, T, condition))
    print(CrDHMM(1, length, pi, Z, T, condition))
}
checkR(1)
checkR(2)
checkR(3)
checkR(10)


Z <- array(c(.9,.05,.05,.1,.8,.1,0,1/2,1/2), c(3,3))
Z
T <- array(c(0,0,1,0,0,1,0,0,1), c(3,3))
T
pi <- c(1, 1, 1)
condition <- c(1, 1, 1)

y <- c(1,1)
length <- length(y)
dDHMM( y, length, pi, Z, T, condition, 0)
CdDHMM(y, length, pi, Z, T, condition, 0)


## testing new DSL functions:
## nimVector()
## nimArray()

library(nimble)

Rnf <- nimbleFunction(
    run = function(val = double(), len = double(), r = double(), c = double()) {
        onesVec <- nimVector(1, 5)
        print('onesVec:');   print(onesVec)
        zerosVec <- nimVector(0, 5)
        print('zerosVec:');   print(zerosVec)
        otherVec <- nimVector(val, len)
        print('otherVec:');   print(otherVec)
        sumVec <- zerosVec + onesVec + otherVec[1:5]
        print('sumVec:');   print(sumVec)
        arr <- nimArray(.4, r, c)
        print('arr:');   print(arr)
        print('sum of otherVec: ', sum(otherVec))
        print('sum of arrL: ', sum(arr))
    }
)

Cnf <- compileNimble(Rnf)

Rnf(10, 7, 3, 4)
Cnf(10, 7, 3, 4)

## trying nimVector() and nimArray() in BUGS code ?!?
## it works!!
code <- nimbleCode({
    a[1:4] <- nimVector(2, 4)
    phi ~ dnorm(0, 1)
    b[2:10] <- nimVector(10, 9)
    b[1] <- phi
})

Rmodel <- nimbleModel(code)
simulate(Rmodel, 'phi')
calculate(Rmodel)

Cmodel <- compileNimble(Rmodel)

model <- Rmodel
model <- Cmodel

model$a
model$phi
model$b

simulate(model)
calculate(model)

## re-creating NA's data example error
## 11 June 2015:
## Perry said he might get to this by July
## running tests/test-copy.R also gives this same error

library(nimble)

code <- nimbleCode({
    y[2:4] ~ dmnorm(mu3[1:3], cov = cov3[1:3, 1:3])
})
constants <- list(mu3=rep(0,3), cov3=diag(3))
data <- list(y = c(NA,0,0,0))

Rmodel <- nimbleModel(code, constants, data)

calculate(Rmodel)
simulate(Rmodel)

Cmodel <- compileNimble(Rmodel)
calculate(Rmodel)
simulate(Rmodel)


## same error from test-copy.R
library(testthat)
library(nimble)
test_package('nimble', 'copy')




## Zhenglei's example, trying MCMC suite on it

library(lme4)
library(ggplot2)
library(scales)
library(nimble)
library(coda)

I <- 25
J <- 5
beta <- c(-0.2, 0.1)
theta <- 1
sigma <-0.5

simfun <- function(substance=letters[1:10],experiment=factor(1:5),log10HQ=seq(-1,3,by=0.2),theta=1,sigma=0.1,beta=c(-0.2, 0.5)){
    expdat <- expand.grid(experiment = experiment,log10HQ=log10HQ,substance=substance)
    expdat$HQ <- 10^expdat$log10HQ
    expdat$obs <- factor(seq(nrow(expdat)))
    n <- nrow(expdat)
    nsub <- length(substance)
    re <- rnorm(nsub,0,theta)
    names(re) <- substance
    fixed <- beta[1]+beta[2]*expdat$log10HQ
    expdat$y <- re[expdat$substance]+rnorm(n,0,sigma)+fixed
    expdat$effect <- gtools::inv.logit(expdat$y)
    return(expdat)
}

expdat <- simfun(substance=factor(paste0('S',1:I)),experiment = factor(1:J),log10HQ=seq(-1,3,by=0.2),theta=0.5,beta=c(-0.2,0.5),sigma=0.3)

code <- nimbleCode({
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm( 0, 0.001 )
    tau.y ~ dgamma(0.001, 0.001)  
    tau.re ~ dgamma(0.001, 0.001)
    for(j in 1:M){
        ranef.v[j]~dnorm(0, tau.re)
    }
    for (i in 1:N){
        y[i] ~ dnorm ( mu.y[i], tau.y )
        mu.y[i] <- alpha+beta*x[i] + ranef.v[unit[i]]
    }
})

constants <- list(N = nrow(expdat),M=nlevels(expdat$substance),x = expdat$log10HQ,unit=as.numeric(expdat$substance))
##data <- data.frame(y = expdat$y)
data <- list(y = expdat$y)
inits <- list(alpha = 0, beta = 0.5, tau.y=1,tau.re=1)

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'jags', 'noConj'),
    niter = 20000,
    makePlot = TRUE,
    savePlot = TRUE,
    summaryStats = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
)

out$summary[, 'effectiveSize', ] / out$timing[c('nimble', 'jags', 'noConj')]


## helping Nick get started with using node$get_mean() function for particle filter

library(nimble)

code <- nimbleCode({
    x ~ dnorm(1, 1)
    y ~ dnorm(x, 1)
})

Rmodel <- nimbleModel(code)
simulate(Rmodel)
Rmodel$x
Rmodel$y
Rmodel$nodeFunctions[['x']]$get_mean()
Rmodel$nodeFunctions[['y']]$get_mean()

## only necessary for the 'higherLevel' example
virtual_NF <- nimbleFunctionVirtual(
    ## don't need to supply prototpye for run(),
    ## since the default is no args, and returnType void
    methods = list(  ## definitely need prototype for return_mean() function
        return_mean = function() {
            returnType(double())
        }
    )
)

nfDef <- nimbleFunction(
    contains = virtual_NF,  ## only necessary for the 'higherLevel' example
    setup = function(model, node) {
        nfList <- nimbleFunctionList(node_stoch_dnorm)
        nfList[[1]] <- model$nodeFunctions[[node]]
    },
    run = function() {
        print('node value is: ', model[[node]])
        print('(assumed to be dnorm) node mean is: ', nfList[[1]]$get_mean())
    },
    methods = list(                        ## only necessary for the 'higherLevel' example
        return_mean = function() {         ##
            returnType(double())           ##
            return(nfList[[1]]$get_mean()) ##
        }                                  ##
    )                                      ##
)

nfX <- nfDef(Rmodel, 'x')
nfY <- nfDef(Rmodel, 'y')
nfX$run()
nfY$run()

compiledList <- compileNimble(list(Rmodel, nfX, nfY))
Cmodel <- compiledList[[1]]
CnfX <- compiledList[[2]]
CnfY <- compiledList[[3]]

CnfX$run()
CnfY$run()

CnfX$return_mean()
CnfY$return_mean()

## this is the 'higherLevel' example
nfDefHigherLevel <- nimbleFunction(
    setup = function(model, nodes) {
        nfList <- nimbleFunctionList(virtual_NF)
        for(i in seq_along(nodes))
            nfList[[i]] <- nfDef(model, nodes[i])
    },
    run = function() {
        declare(meansVector, double(1, length(nfList)))  ## declares a vector of doubles, with length equal to the length of nfList
        for(i in seq_along(nfList)) {
            nfList[[i]]$run()
            meansVector[i] <- nfList[[i]]$return_mean()
        }
        print('vector of the means is: ', meansVector)
    }
)

nfHigherLevel <- nfDefHigherLevel(Rmodel, c('x', 'y'))

nfHigherLevel$run()

CnfHigherLevel <- compileNimble(nfHigherLevel, project = Rmodel)

CnfHigherLevel$run()


## fixing values(model, nodes) <- value reference class '<<-' warning


library(nimble)

nfdef <- nimbleFunction(
    setup = function(model) {},
    run = function(val = double(1)) {
        values(model, 'xxx') <<- val
    }
)

code <- nimbleCode({
    xxx ~ dnorm(0, 1)
})

Rmodel <- nimbleModel(code)

nf <- nfdef(Rmodel)

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(nf, project = Rmodel)

Rmodel$xxx
nf$run(3)
Rmodel$xxx

Cmodel$xxx
Cnf$run(4)
Cmodel$xxx


## testing / trying to find some error in conjugacy system(??)  May 2015

library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')

## run entire MCMC testing system
source('~/GitHub/nimble/packages/nimble/inst/tests/test-mcmc.R')

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)
code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list(mu0 = mu0, Q0 = Q0, Q = Q, B = B, a=a)
data <- list(y=y)
m <- nimbleModel(code, constants = constants, data=data)

m <- nimbleModel(code, constants = data)

spec <- configureMCMC(m)
Rmcmc <- buildMCMC(spec)
spec$getSamplers()

options(error = recover)

test_mcmc(model = code, data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))


## testing distribution of particle filter likelihood estimates

ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))

plot(muvec, varll)
varest <- 1/2*sigx^4/sigy^4 + muvec^2*sigx^2/sigy^4
lines(muvec, varest, col='red')



## library(nimble)
## Rpf <- buildPF(Rmodel, 'x')
## Cmodel <- compileNimble(Rmodel)
## Cpf <- compileNimble(Rpf, project = Rmodel)
## Rpf$run(10)
## Cpf$run(10)
## ll <- numeric()
## for(i in seq_along(muvec)) {
##     mu <- muvec[i]
##     Cmodel$mu <- mu
##     ll[i] <- Cpf$run(10000)
## }
## plot(muvec, ll)


library(nimble)
source('~/GitHub/pfLL/pfLL.R')
m <- 1000
rep <- 1000
n <- m*rep
muvec <- seq(-5, 5, by=0.2)
sigx <- 1
sigy <- 1
taux <- 1/sigx^2; tauy <- 1/sigy^2
y <- 0   ## observation
constants <- list()
inits <- list(mu = y, sigx = sigx, sigy = sigy)
data <- list(y = y)
code <- quote({
    mu ~ dnorm(0, 0.00001)  ## these should not be necessary
    sigx ~ dunif(0, 1000)  ## these should not be necessary
    sigy ~ dunif(0, 1000)  ## these should not be necessary
    x ~ dnorm(mu, sigx)
    y ~ dnorm(x, sigy)
})
Rmodel <- nimbleModel(code, constants, data, inits)

out <- pfLL(Rmodel, 'x', param = data.frame(mu=muvec), m=m, rep=rep, makePlot=FALSE)

ll <- apply(out$ll, 1, mean)
##ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))
plot(muvec, ll)
ll.pred <- -1/2 * muvec^2 / (sigx^2+sigy^2) + log(1/sqrt(2*pi*(sigx^2+sigy^2)))
lines(muvec, ll.pred, col='red')


varll <- apply(out$ll, 1, var)
##varll <- varll - min(varll)  ### TEMPORARY
plot(muvec, varll)
var.pred <- exp(-8.8 + 1.2592*abs(muvec))
lines(muvec, var.pred, col='red')

1



n <- 100000
sigx <- 1
muvec <- seq(-3, 3, by=0.1)
nmu <- length(muvec)
X <- array(NA, c(n, nmu))
for(i in 1:nmu) X[,i] <- rnorm(n, muvec[i], sigx)

Y <- X^2
var <- apply(Y, 2, var)
pred <- muvec^2 + sigx^2
plot(muvec, var, type='p')
lines(muvec, pred, col='red')



## figuring out 3D plotting in R
rm(list=ls())
df <- expand.grid(list(x=1:10, y=1:10))
rep <- 1
ll <- array(NA, c(dim(df)[1], rep))
for(i in 1:dim(df)[1]) {
    for(j in 1:rep) {
        x <- df[i, 'x']
        y <- df[i, 'y']
        dist <- sqrt((x-4)^2 + (y-7)^2)
        ll[i,j] <- rnorm(1, 100-dist^2, sd = dist/10)
    }
}

library(plot3D)
?scatter3D

x <- rep(df$x, rep)
y <- rep(df$y, rep)
z <- as.numeric(ll)

scatter3D(x, y, z)


## minimally reproducible example of getDependencies error for Perry
library(nimble)
code <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
    y[1] ~ dnorm(x[1], 1)
    y[2] ~ dnorm(x[2], 1)
})

Rmodel <- nimbleModel(code)

Rmodel$getDependencies('x[1]')   ## should be x[1], y[1]
## [1] 'x[1]' 'y[1]' 'y[2]'

Rmodel$getDependencies('x[2]')   ## should be x[2], y[2]
## [1] 'x[2]'



## trying to figure out conjugate sampling in 'equiv' BUGS example model
library(nimble)
code <- nimbleCode({
    tau[1] ~ dgamma(0.001, 0.001)
    tau[2] ~ dgamma(0.001, 0.001)
    pi ~ dnorm(0, 1.0E-06)
    phi ~ dnorm(0, 1.0E-06)
    mu ~ dnorm(0, 1.0E-06)
    for(i in 1:N) { 
        d[i] ~ dnorm(0,tau[2])  # Subject random effect
        for (k in 1:2) {
            Treat[i,k] <- group[i]*(k-1.5) + 1.5  # treatment given
            m[i,k] <- mu + pow(-1, Treat[i,k]-1)* phi /2 + pow(-1, k-1)* pi /2 + d[i]
            Y[i,k] ~ dnorm(m[i,k], tau[1])
        }
    }
    ##theta <- exp(phi)
    ##equivalence <- step(theta - 0.8) - step(theta - 1.2)
})
Y <- structure(c(1.4, 1.64, 1.44, 1.36, 1.65, 1.08, 1.09, 1.25, 1.25, 1.3, 1.65, 1.57, 1.58, 1.68, 1.69, 1.31, 1.43, 1.44, 1.39, 1.52), .Dim = c(10, 2))
group <- c(1, 1, -1, -1, -1, 1, 1, 1, -1, -1)
N <- 10
##T <- 2
constants <- list(N=N, group=group)
data <- list(Y=Y)
inits <- list(mu=0, phi=0, pi=0, tau=c(1,1))
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 10

set.seed(0); Rmcmc$run(niter)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0); Cmcmc$run(niter)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples <- Csamples[, dimnames(Rsamples)[[2]]]
Csamples

Rsamples-Csamples


## test equiv model
library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')
test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)


1

## trying Kalman Filter and calculating likelihood for a linear SSM

# Say state process is
# X(t+1) = r * X(t) + b + nu(t)
# nu(t) ~ dnorm(0, sig2nu)
# 
# Obs process is
# Y(t) = X(t) + eps(t)
# eps(t) ~ dnorm(0, sig2eps)
# 
# Time series goes from t = 1:T
# 
# This is the simplest possible model for us.
# 
# Direct mvn approach:
#   1. Stationary variance of states:
#   V[x] = r^2 V[x] + sig2nu
# V[x] = sig2nu / (1-r^2)
# 
# 2. Deriving entries in covariance matrix of Y (from 1:T)
# Cov[X(t), X(t+1)] = r V[x]
# Cov[X(t), X(t+d)] = r^d V[x]
# Therefore Cov[Y(t), Y(t+d)] = r^d V[x]
# 
# And on the diagonal:
#   V[y] = V[x] + sig2eps
# 
# And mu(t) = E[Y] = E[X] = b/(1-r) for all t
# 
# So with that mean vector and covariance matrix in hand 
# you can directly calculate mvnorm(Y, mu, Cov(Y))
# 
# Kalman filter approach and direct LL are shown in code below.  
# Sorry  I didn't think of this.  I had this from the early 
# explorations that led to Knape, Besbeas and de Valpine in 
# Ecology last year (or was it this year?).  
# 
# In this setup, you want numSites = 1 (we were considering time-series 
# with data from multiple sites).  The LL functions use the meanData 
# (mean over all the sites, which we were comparing to a more complete 
# site-specific model), which for you will be the same as the fullData.



simData <- function(mu=10, a=0.8, sigPN=0.2, sigOE=0.4, t=20) {
    b <- mu*(1-a)
    X <- mu
    for(i in 1:50) X <- a*X + b + rnorm(1,0,sigPN) ## get to stationary
    states <- numeric(t)
    states[1] <- X
    for(i in 2:t) states[i] <- a*states[i-1]+b + rnorm(1,0,sigPN)
    obsStates <- rnorm(t, states, sigOE)
    list(x = states, y = obsStates, mu=mu, a=a, b=b, sigOE=sigOE, sigPN=sigPN, t=t)
}

mvn_ll <- function(d) {
    require(mvtnorm)
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_vec <- rep( b/(1-a), t)
    cov_mat <- matrix(nrow = t, ncol = t)
    var_x <- sigPN2 / (1-a^2)
    var_y <- var_x + sigOE2
    diag(cov_mat) <- var_y
    for(i in 2:t) {
        for(j in 1:(i-1)) {
            cov_mat[i,j] <- (a^(i-j))*var_x
            cov_mat[j,i] <- cov_mat[i,j]
        }
    }
    dmvnorm(y, mu_vec, cov_mat, log = TRUE)
}

KF_ll <- function(d) {
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_x <- b/(1-a)
    var_x <- sigPN2 / (1-a^2)
    cov_xy <- var_x
    var_y <- var_x + sigOE2
    ll <- dnorm(mu_x, y[1], sqrt(var_y), log = TRUE)
    mu_x <- mu_x + (cov_xy / var_y) * (y[1] - mu_x)
    var_x <- var_x - cov_xy*cov_xy / var_y
    for(i in 2:t) {
        mu_x <- a*mu_x + b
        var_x <- a^2 * var_x + sigPN2
        if(!is.na(y[i])) {
            cov_xy <- var_x
            var_y <- var_x + sigOE2
            ll <- ll + dnorm(mu_x, y[i], sqrt(var_y), log = TRUE)
            mu_x <- mu_x + (cov_xy / var_y) * (y[i] - mu_x)
            var_x <- var_x - cov_xy*cov_xy / var_y
        }
    }
    ll
}

code <- quote({
    x[1] ~ dnorm(mu, sd = sqrt((sigPN^2) / (1 - a^2)))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

d <- simData()

library(nimble)
constants <- list(a=d$a, b=d$b, t=d$t, sigOE=d$sigOE, sigPN=d$sigPN, mu=d$mu)
data <- list(y = d$y)
inits <- list(x = d$x)
Rmodel <- nimbleModel(code = code, constants = constants, data=data, inits=inits)
calculate(Rmodel)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

Cpf$run(100000)

mvn_ll(d)

KF_ll(d)



## exploring how truncation works in our system

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ T(dnorm(a, 1), 1.5, 1.6)
    c ~ T(dnorm(b, 1), 3, Inf)
    d ~ dnorm(c, 1)
    e ~ dnorm(d, 1)
})
constants <- list()
data <- list()
inits <- list()

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)
Rmodel <- md$newModel(data=data, inits=inits)

Rmodel$modelDef$printDI()

node <- 'c'
Rmodel$nodes[[node]]$calculate
Rmodel$nodes[[node]]$simulate
Rmodel$nodes[[node]]$getLogProb

lapply(Rmodel$modelDef$declInfo, function(di) di$truncation)
for(node in c('a', 'b', 'c')) print(Rmodel$isTruncated(node))
for(node in c('a', 'b', 'c')) print(Rmodel$getBounds(node))

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## using rankSample

library(nimble)

##Cmodel <- compileNimble(Rmodel = nimbleModel(quote({a ~ dnorm(0, 1)})))

nfDef <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer(), silent = logical()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp, silent)
        print(samp)
    }
)
nf <- nfDef()
Cnf <- compileNimble(nf)

nfDefOmit <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp)
        print(samp)
    }
)
nfOmit <- nfDefOmit()
CnfOmit <- compileNimble(nfOmit)
x <- 1:10

silent <- FALSE
silent <- TRUE

set.seed(1);   .Call('rankSample', c(.5,.4,.1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,.1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,.1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,.1), 10L,    silent)

set.seed(1);   .Call('rankSample', c( 0, 0, 0), 10L, x, silent)
set.seed(1);            rankSample(c( 0, 0, 0), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c( 0, 0, 0), 10L,    silent)
set.seed(1);               Cnf$run(c( 0, 0, 0), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,-1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,-1), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x);            print(x)
set.seed(1);            nfOmit$run(c(.5,.4,-1), 10L)
set.seed(1);           CnfOmit$run(c(.5,.4,-1), 10L)



## difference in rank sample
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        declare(wts, double(1, 1))
        wts[1] <- 1
        m <- 1L
        declare(samp, integer(1, 1))
        rand <- rnorm(1, 0, 1)
        print('random number before rankSample: ', rand)
        rankSample(wts, m, samp)
        rand <- rnorm(1, 0, 1)
        print('random number after rankSample: ', rand)
    }
)

nf <- nfDef()
Cnf <- compileNimble(nf)
set.seed(0)
nf$run()
set.seed(0)
Cnf$run()



## testing buildPF() particle filter (pf) algorithm
library(nimble)
code <- quote({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
t <- 5
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
Rmodel$mu <- 1/(1-.95)
Rmodel$b <- 1
Rmodel$sigPN <- .2
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)
rm(Rmodel)
Rmodel <- nimbleModel(code, constants, data, inits)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

m <- 100
set.seed(0);    pf$run(m)
set.seed(0);   Cpf$run(m)

1
2
3



## making MCMCsuite work for Ryan using Stan MCMC

allModels <- c('blocker', 'bones', 'dyes', 'line', 'pump', 'rats')
library(nimble)
library(coda)
allModels <- 'blocker'
mcmcs=c('nimble','jags','noConj','stan')
x=list()
for (i in 1:length(allModels)){
    x[[i]]<-readBUGSmodel(model=allModels[i],
                          dir=getBUGSexampleDir(allModels[i]),
                          returnModelComponentsOnly=TRUE)
}
i <- 1
suite_output <- MCMCsuite(
    x[[i]]$model,
    constants = x[[i]]$data,
    inits = x[[i]]$inits,
    MCMCs = mcmcs,
    makePlot=F,
    savePlot=F,
    summaryStats=c('mean','median','sd','CI95_low','CI95_upp','effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) })),
    stan_model = paste('~/temp/', allModels[i], '.stan', sep='')
)


## learning about R's memory management and storage of objects
gcinfo(TRUE)
gc()

library(pryr)
?object_size
object_size(1:100)
object_size(mtcars)

x <- 0:50
s <- sapply(x, function(i) object_size(1:i))
plot(x, s, type='s', ylim=c(0,300))

obj <- complex()
object_size(obj)

mem_used
mem_used()

pryr:::node_size()
pryr:::show_bytes

mem_change(NA)

bytes(1)
bytes(1L)
bytes('ab')

address(2)
x <- 1:10
address(x)
refs(x)
y <- x
refs(x)
z <- x
x[1] <- 16L

tracemem(x)
x[1] <- .1

x <- 15
address(x)

system(paste0('~/temp/t ', address(x)))



## testing of mcmcplots package function mcmcplot()

rm(list=ls())
library(nimble)
model <- 'SSMcorrelated'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)

out <- MCMCsuite(code, constants, data, inits, monitors=c('a','b','p'), MCMCs=c('nimble','nimble_RW', 'nimble_slice', 'jags'), makePlot=FALSE, niter=10000)

library(coda)
samples <- out$samples
dim(samples)
dimnames(samples)
mcmcList <- list()
for(i in 1:dim(samples)[1])     mcmcList[[i]] <- coda::mcmc(t(samples[i,,]))
names(mcmcList) <- dimnames(samples)[[1]]
codaList <- coda::as.mcmc.list(mcmcList)

library(mcmcplots)
mcmcplot(codaList)
traplot(codaList)




## testing conjugacy in newNimbleModel branch

library(nimble)
code <- nimbleCode({
    for(i in 1:10) {
        a[i] ~ dnorm(0, 1)
        b[i] ~ dgamma(1, 1)
    }
    for(i in 1:5)   a_temp[i] <- 3*a[i]
    for(i in 6:10)  a_temp[i] <- a[i]^2
    for(i in 1:10)  x[i] ~ dnorm(a_temp[i], 1)
    for(i in 1:10) {
        y1[i] ~ dpois(2*b[i])
        y2[i] ~ dgamma(1, rate = b[i])
    }
})
Rmodel <- nimbleModel(code)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
for(i in 1:200) spec$addSampler('RW', 'x[2]', print=F)
spec$addSampler('RW', 'x[2]')

spec$setSamplers(c(1:40, 1:40))
spec$setSamplers(c('a', 'y2'))
spec$removeSamplers('a')
spec$setSamplers()


## NEED TO RE-RUN THIS TEST, ONCE PERRY FIXES THE MODELVALUES COPYING ISSUE
## Dirichlet-multinomial conjugacy
## single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- nimbleCode({
    y[1:K] ~ dmulti(p[1:K], n)
    p[1:K] ~ ddirch(alpha[1:K])
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001)
    }
})
constants <- list(n = n, K = K)
data <- list(y = y)
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
Rmodel <- nimbleModel(code, constants=constants, data=data, inits=inits)

mv <- modelValues(Rmodel)
Rmodel$y
mv$y
nimCopy(from=Rmodel, to=mv, nodes = 'y', row = 1)
mv$y


spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

set.seed(0)
Rmcmc$run(10)

conjugate posterior density appears to be wrong, off by Inf
conjugate posterior density appears to be wrong, off by Inf
There were 50 or more warnings (use warnings() to see the first 50)

test_mcmc(model = code, data= data, seed = 0, numItsC = 10000,
          inits = inits,
          results = list(mean = list(p = p)),
          resultsTolerance = list(mean = list(p = rep(.06, K))))




## just a good old fashioned test that it's working
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})

constants <- list()
data <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$checkConjugacy()

spec <- configureMCMC(Rmodel, control = list(scale=0.01), monitors = c('x', 'c'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
set.seed(0)
Rmcmc$run(10)
samples <- as.matrix(Rmcmc$mvSamples)
samples


sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977),
  c = c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))

test_mcmc(model = code, exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01))




## testing if a Ref Class can just return a smaller object (not the whole Ref Class object)
RC <- setRefClass(
    Class = 'RC',
    fields = list(a='ANY'),
    methods = list(
        initialize = function(a) {
            a <<- a
            return(a)
        }
    )
)

rc <- RC(1)
class(rc)
rc$a



## bunch of random stuff

fileToList <- function(file) {
    env <- new.env()
    source(fileName, local = env)
    lst <- list()
    for(name in ls(env))   lst[[name]] <- get(name, env)
    return(lst)
}

library(nimble)

NF <- nimbleFunction(
    setup = function() { },
    run = function(a = double(0)) {
        returnType(double(0))
        return(a)
    }
)

myNF <- NF()
myNF_C <- compileNimble(myNF)

> myNF$run
## function (a) 
##     return(a)

myNF_C$run
## function (a) 
## {
##     ans <- .Call('CALL_nfRefClass60_operator_', a, .basePtr)
##     ans <- ans[[2]]
##     ans
## }


### testing the new autoBlock
rm(list=ls())
library(nimble)
model <- 'litters'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$getNodeNames()

ab <- autoBlock(Rmodel, run = runList)
ab
ab$spec$getSamplers()
Rmcmc <- buildMCMC(ab$spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)


spec <- configureMCMC(Rmodel, autoBlock=TRUE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

Rmcmc <- buildMCMC(Rmodel, autoBlock=TRUE)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)







## March 2015
## timing comparisons for red-state-blue-state
library(nimble)
library(rjags)
rm(list=ls())
load('~/GitHub/autoBlock/data/model_redblue.RData')

## this version of red-state-blue-state code is 'jags friendly'
## reparametrized to use Precision matrix, replaced expit with logit
code <- nimbleCode({
    for (i in 1:2) {
        for (j in 1:2) {
            gamma[i, j] ~ dnorm(0, 1e-04)
        }
    }
    sigmaIntercept ~ dunif(0, 100)
    sigmaSlope ~ dunif(0, 100)
    rho ~ dunif(-1, 1)
    Precision[1, 1] <- 1/(1-rho^2) * 1/(sigmaIntercept^2)
    Precision[2, 2] <- 1/(1-rho^2) * 1/(sigmaSlope^2)
    Precision[1, 2] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    Precision[2, 1] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    for (i in 1:Nstates) {
        stateBetaMeans[i, 1] <- gamma[1, 1] + gamma[1, 2] * stateIncome[i]
        stateBetaMeans[i, 2] <- gamma[2, 1] + gamma[2, 2] * stateIncome[i]
        stateBetas[i, 1:2] ~ dmnorm(stateBetaMeans[i, 1:2], Precision[1:2, 1:2])
    }
    for (i in 1:N) {
        logit(p[i]) <- stateBetas[state[i], 1] + stateBetas[state[i], 2] * income[i]
        y[i] ~ dbern(p[i])
    }
})

niter <- 50000
monitorVars <- c('sigmaIntercept', 'sigmaSlope', 'rho', 'gamma')
constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)


system.time(modelDef <- nimbleModel(code=code, constants=constants, data=data, inits=inits, returnDef=TRUE))[3]/60
26

system.time(Rmodel <- modelDef$newModel(data=data, inits=inits))[3]/60
5

system.time(spec <- configureMCMC(Rmodel))
6 seconds

system.time(Rmcmc <- buildMCMC(spec))[3]
15 seconds

system.time(Cmodel <- compileNimble(Rmodel))[3]/60
11

system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))[3]/60
1

system.time(Cmcmc$run(niter))[3]/60
11 minutes for 50,000 iterations

system.time(jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE))[3]/60
10 seconds

system.time(jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1))[3]/60
2.001833





## displaying calc and sim functions
code <- nimbleCode({ x ~ dnorm(3, sd = 5) })
model <- nimbleModel(code)
model$nodes$x$simulate
## function()
##     model$x <<- rnorm(1, mean = 3, sd = 5)
model$nodes$x$calculate
## function() {
##     model$logProb_x <<- dnorm(model$x, mean = 3, sd = 5, log = 1)
##     return(invisible(model$logProb_x))
## }
model$nodes$x$getLogProb
## function() 
##     return(model$logProb_x)


### March 2015
### trying to figure out new MCMCspec

library(nimble)


#These first few lines are nothing new
myModelCode <- nimbleCode({
	a[1] ~ dnorm(0,1)
	a[2] ~ dexp(a[1]^2)
	a[3] ~ dnorm(0,1)
})
myModel <- nimbleModel(myModelCode)
mySpec1 <- configureMCMC(myModel)
myMCMC1 <- buildMCMC(mySpec1)
cmod <- compileNimble(myModel)
cmcmc1 <- compileNimble(myMCMC1, project = myModel)


# Rebuilding the way we are supposed to: 
# building a new spec with configureMCMC 
# adding a sampler which has already been compiled 
# not altering the monitors
mySpec2 <- configureMCMC(oldSpec = mySpec1)
mySpec2$addSampler('RW', control = list(targetNode = 'a[3]'))
rmcmc2 <- buildMCMC(mySpec2)
cmcmc2 <- compileNimble(rmcmc2, project = myModel)




#### March 2015 working through a full compileNimble call

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 1)
Rmodel <- nimbleModel(code, inits = inits)

debug(compileNimble)
Cmodel <- compileNimble(Rmodel)

debug(project$compileModel)
debug(modelCpp$buildAll)
debug(buildNodes)
debug(nimbleProject$compileNimbleFunctionMulti)
debug(compileNimbleFunction)
debug(buildNimbleFunctionCompilationInfo)


#### March 2015 testing of RW_llFunction or Carl B, submitted as test to testing suite

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 0)
Rmodel <- nimbleModel(code=code, inits=inits)

llFunction <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        ll <- dnorm(1, model$a, 1, log=1)
        returnType(double())
        return(ll)
    }
)
 
myLL <- llFunction(Rmodel)

spec <- configureMCMC(Rmodel, nodes=NULL)
spec$addSampler('RW_llFunction', list(targetNode='a', llFunction=myLL, includesTarget=FALSE))
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
unlist(Rmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

set.seed(0)
Cmcmc$run(10)
unlist(Cmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
## should be within some tolerance of 0.5


###### Feb 2015 example for Cal Poly interview talk

library(nimble)

code <- nimbleCode({
    for( i in 1 : N ) {
        b[i] ~ dnorm(mu, tau)
        r[i] ~ dbin(p[i], n[i])
        logit(p[i]) <- b[i]
    }
    pop.mean <- 1 / (1 + exp(-mu))
    mu ~ dnorm(0, 0.001)
    sigma <- 1 / sqrt(tau)
    tau ~ dgamma(0.001, 0.001)
})
constants <- list(N = 10, n = c(12,20,18,19,10,23,19,15,8,10))
data <- list(r = c(4,7,10,3,6,7,11,9,3,8))
inits <- list(mu = 0, tau = 1)

model <- nimbleModel(code=code, constants=constants, data=data, inits=inits)

spec <- configureMCMC(model)
spec$getSamplers()
spec$addSampler(‘slice’, list(targetNode = ‘mu’))
spec$getSamplers()
mcmc <- buildMCMC(spec)

Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model)

Cmcmc$run(10000)

samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)

