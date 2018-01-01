source("chapter8.R")
printf <- function(...) invisible(print(sprintf(...)))

###############################################################################
# Load data
y1 <- y.school1
y2 <- y.school2

n1 <- length(y1)
n2 <- length(y2)

nIter <- 5000

###############################################################################
# Define priors

mu.mu0 <- 50  # prior mean for population mean
mu.s20 <- 625 # Prior variance for population mean

delta.mu0 <- 0   # Prior mean for difference between population means
delta.s20 <- 625 # Prior variance for difference between population means

s2.n0 <- 1   # Prior sample size
s2.p0 <- 100 # Prior s2 estimate

###############################################################################
# Set starting values
mu    <- (mean(y1) + mean(y2))/2
delta <- (mean(y1) - mean(y2))/2

###############################################################################
# Gibbs sampler arrays
MU <- NULL
DELTA <- NULL
S2 <- NULL

###############################################################################
# Run

set.seed(1)

for(i in 1:nIter){

  printf("Iteration %s", i)

  ### Estimate/update s2
  s2.a <- (s2.n0 + n1 + n2) / 2
  s2.b <- (s2.n0 * s2.p0 + ( sum((y1 - (mu + delta))^2) + sum((y1 - (mu - delta))^2) )) / 2
  s2 <- 1/rgamma(1, s2.a, s2.b)
  
  ### Update mu
  mu.var <- 1/ (1/mu.s20 + (n1+n2)/s2 )
  mu.mu  <- mu.var * (mu.mu0/mu.s20 + sum(y1 - delta)/s2 + sum(y2 + delta)/s2 )
  mu <- rnorm(1, mu.mu, sqrt(mu.var))
  
  ### Update delta
  delta.var <- 1/ (1/delta.s20 + (n1+n2)/s2 )
  delta.mu  <- delta.var * (delta.mu0/delta.s20 + sum(y1 - mu)/s2 - sum(y2 - mu)/s2 )
  delta <- rnorm(1, delta.mu, sqrt(delta.var))
  
  MU <- c(MU, mu)
  DELTA <- c(DELTA, delta)
  S2 <- c(S2, s2)
    
}

printf("p(d > 0 | y1, y2) = %f", length(DELTA[DELTA > 0])/nIter)


###############################################################################
# Plot the prior and posterior distribution of mu

mu.priorPlot <- NULL
mu.postPlot  <- NULL
mu.range <- (2000:8000)/100
delta <- 0
s2 <- var(c(y1,y2))
for(mu in mu.range) {
  mu.priorPlot <- c(mu.priorPlot, dnorm(mu, mu.mu0, sqrt(mu.s20)))
  mu.var <- 1/ (1/mu.s20 + (n1+n2)/s2 )
  mu.mu  <- mu.var * (mu.mu0/mu.s20 + sum(y1 - delta)/s2 + sum(y2 + delta)/s2 )
  mu.postPlot  <- c(mu.postPlot, dnorm(mu, mu.mu, sqrt(mu.var)))
}
plot(mu.range, mu.priorPlot, type="l", col="red",  ylim=c(0,0.30))
lines(mu.range, mu.postPlot, col="green")
legend(legend=c("prior", "posterior"))

###############################################################################
# Plot the prior and posterior distribution of d

delta.priorPlot <- NULL
delta.postPlot  <- NULL
delta.range <- (-2500:2500)/100
mu <- 49
s2 <- var(c(y1,y2))
for(delta in delta.range) {
  delta.priorPlot <- c(delta.priorPlot, dnorm(delta, delta.mu0, sqrt(delta.s20)))
  delta.var <- 1/ (1/delta.s20 + (n1+n2)/s2 )
  delta.mu  <- delta.var * (delta.mu0/delta.s20 + sum(y1 - mu)/s2 - sum(y2 - mu)/s2 )
  delta.postPlot  <- c(delta.postPlot, dnorm(delta, delta.mu, sqrt(delta.var)))
}
plot(delta.range, delta.priorPlot, type="l", col="red",  ylim=c(0,0.30))
lines(delta.range, delta.postPlot, col="green")
legend(legend=c("prior", "posterior"))
  