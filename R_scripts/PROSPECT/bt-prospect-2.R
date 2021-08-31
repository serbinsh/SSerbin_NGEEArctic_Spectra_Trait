library(rrtm)
library(BayesianTools)

data("testspec", package = "PEcAnRTM")

obs <- testspec_ACRU[,1]

likelihood <- function(params) {
  #N, Cab, Car, Cbrown, Cw, Cm
  mod <- prospect5(params[1], params[2], params[3], params[4], params[5], 
                   params[6])
  sum(dnorm(obs, mod[["reflectance"]], log = TRUE))
}

prior <- createPrior(
  density = function(x) {
    # Note: Truncated at zero
    N <- dnorm(x[1], 1, 5, log = TRUE)
    Cab <- dnorm(x[2], 40, 15, log = TRUE)
    Car <- dnorm(x[3], 3, 15, log = TRUE)
    Cbrown <- dnorm(x[4], 0.08, 0.15, log = TRUE)
    Cw <- dlnorm(x[5], -4.456, 1.216, log = TRUE)
    Cm <- dlnorm(x[6], -5.15, 1.328, log = TRUE)
    rsd <- dexp(x[7], 10, log = TRUE)
    N + Cab + Car + Cbrown + Cw + Cm + rsd
  },
  sampler = function(n = 1) {
    N <- rnorm(n, 1, 5)
    Cab <- rnorm(n, 40, 15)
    Car <- rnorm(n, 5, 5)
    Cbrown <- rnorm(n, 0.08, 0.15)
    Cw <- rlnorm(n, -4.456, 1.216)
    Cm <- rlnorm(n, -5.15, 1.328)
    rsd <- rexp(n, 10)
    cbind(N, Cab, Car, Cbrown, Cw, Cm, rsd)
  },
  lower = c(1, 0, 0, 0, 0, 0, 0)
)

setup <- createBayesianSetup(
  likelihood = likelihood,
  prior = prior
)

samps <- runMCMC(setup)

plot(samps$Z)
