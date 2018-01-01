source("chapter7.r")
printf <- function(...) invisible(print(sprintf(...)))

#Data set with missing values
Y <- Y.pima.miss

###############################################################################
#prior parameters

n <- dim(Y)[1] # Number of samples
p <- dim(Y)[2] # Number of dimensions (features)

mu0 <- c(120, 64, 26, 26) # Prior expected values for our four features
sd0 <- mu0 / 2            # Prior standard deviations of the distributions of our four features
nu0 <- p + 2              # Prior covariance matrix sample size

# Construct an initial covariance matrix for our four features
L0 <- matrix(.1,p,p)
diag(L0) <- 1
L0 <- L0*outer(sd0,sd0)
S0 <- L0

###############################################################################
# Initialize the variables (starting values)
Sigma <- S0         # Our starting Covariance matrix
O <- 1*(!is.na(Y)) # The observed values
Y.full <- Y         # Replace NA values with the mean values of all non-NA values
for (j in 1:p){
  Y.full[is.na(Y.full[,j]), j] <- mean(Y.full[,j], na.rm=TRUE)
}

###############################################################################
# The gibbs sampler

nIter <- 1000

# Create empty lists of all values we are computing
THETA <- NULL
SIGMA <- NULL
Y.MISS <- NULL

set.seed(1) # A consistent seed

for (s in 1:nIter){
  printf("Iteration %s", s)
  # Update Theta
  # (The population means, based on a normal prior and likelihood for the mean)
  ybar <- apply(Y.full, 2, mean) # Calculate the sample mean
  Ln <- solve( solve(L0) + n*solve(Sigma))                           # From equation 7.4
  mun <- Ln %*% ( (solve(L0) %*% mu0) + (n * solve(Sigma) %*% ybar)) # From equation 7.5
  theta <- rmvnorm(1, mun, Ln)
  
  # Update Covariance matrix
  Sn <- S0 + (t(Y.full) - c(theta)) %*% t(t(Y.full) - c(theta)) # From equation 7.9
  sigma <- rwish(1, n+nu0, solve(Sn))
  
  # Update the missing values
  for (i in 1:n){
    i=1
    b <- (O[i,] == 0) # Unobserved values
    a <- (O[i,] == 1) # Observed values
    isA <- solve(sigma[a,a])
    beta.j <- sigma[b,a] %*% isA
    
    sigma.ba <- sigma[b,b] - beta.j %*% sigma[a,b]             # Equation 7.11
    theta.ba <- theta[b] + beta.j %*% t(Y.full[i,a] - theta[a]) # Equation 7.10
    Y.full[i,b] <- rmvnorm(1, theta.ba, sigma.ba)
  }
  
  # Add the sampled values to our list
  THETA <- rbind(THETA, theta)
  SIGMA <- rbind(SIGMA, sigma)
  Y.MISS <- rbind(Y.MISS, Y.full[O == 0])
  
}
