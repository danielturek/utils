

library(rstan)

suppressWarnings(rm('ppp', 'pppp', 'qqq', 'rrr'))

constants <- list(
    N = 10,
    t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
)

data <- list(
    y = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
)

inits <- list(
    alpha = 1,
    beta = 1,
    theta = rep(0.1, 10)
)

monitors <- c('alpha', 'beta')

stan_code <- '
data {
  int<lower=0> N;
  int y[N];
  vector<lower=0>[N] t;
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector<lower=0>[N] theta;
}

model {
  vector[N] lambda;
  alpha ~ exponential(1);
  beta ~ gamma(0.1, 1);
  for(i in 1:N) {
    theta[i] ~ gamma(alpha, beta);
    lambda[i] = theta[i] * t[i];
  }
  y ~ poisson(lambda);
}

'

stan_mod <- rstan::stan_model(model_code = stan_code)

stan_data <- c(constants, data)

## niter and warmup notes:
## - if you omit 'warmup' argument, then half of niter is used as warmup
## - if you specify 'warmup' and 'niter', then the value for 'niter'
##   also includes the warmup iterations

stan_out <- rstan::sampling(
                       stan_mod,
                       data = stan_data,
                       chains = 1,
                       warmup = 5000,
                       iter = 10000,
                       pars = monitors,
                       init = list(inits)  ## or a list of length 'chains'
                   )

rstan::get_elapsed_time(stan_out)    ## timings

samples_stan <- as.matrix(stan_out)  ## samples

head(samples_stan)

## if using multiple chains, then use:
## rstan::extract(stan_out, permute = FALSE)
## which returns a 3D (niter x nchains x nparams) array



