##### R script for group coursework #####
set.seed(123) # Ensure reproducibility
#### Part (a): Normal approximation for ruin probability π1 ####
# Use the normal approximation for the random walk F_n to compute π1 = P(F_N < 0), where N= 40

#=========================#
# 1. Liability Parameters #
#=========================#

## Business Line Parameters (mean & variance of losses)
mu_L1 <- 10; var_L1 <- 4
mu_L2 <- 15; var_L2 <- 9
mu_L3 <- 8; var_L3 <- 2


# Total Quarterly Losses L_n = L_{n,1} + L_{n,2} + L_{n,3}
mu_L <- mu_L1 + mu_L2 + mu_L3 # Total mean loss E[L_n]
var_L <- var_L1 + var_L2 + var_L3 # Total variance Var[L_n]

mu_L
var_L

#=====================#
# 2. Asset Parameters #
#=====================#

# Asset parameters 
lambda <- 0.04 # Premium Loading 
kappa  <- 0.75 # Risk Factor


# Mean and variance of the asset income per quarter
mu_A  <- mu_L * (1 + lambda)
var_A <- var_L * kappa

mu_A
var_A


#=======================#
# 3. Initial Fund Value #
#=======================#

alpha <- 0.5 # Proportion of one quarter's asset income
F0 <- alpha * mu_A # Initial fund at time 0

F0

#==========================#
# 4. Increment Per Quarter #
#==========================#

# Quarterly Increment  Y_n = A_n - L_n
mu_Y <- mu_A - mu_L # Mean of increment
var_Y <- var_A + var_L # Variance of increment

mu_Y
var_Y

#=======================================#
# 5. Fund Distribution After N Quarters #
#=======================================#

N <- 40 # Time horizon (40 quarters = 10 years)

# Random walk with independent increments; F_n is the insurers fund after n quarters
mu_Fn <- F0 + N * mu_Y # Mean of F_n
var_Fn <- N * var_Y # Variance of F_n
sd_Fn <- sqrt(var_Fn) # Standard deviation of F_n

mu_Fn
var_Fn
sd_Fn


#================================================#
# 6. Ruin Probability Using Normal Approximation #
#================================================#

# Approximate F_n using the normal distribution:
pi_1 <- pnorm(q=0, mean = mu_Fn, sd = sd_Fn)

cat("Ruin probability π1 (normal approximation) =", pi_1, "\n")

#### part (b): Monte Carlo estimation of π1 and π2

#==========================================#
# 7. Convert to (meanlog, sdlog) for rlnorm
#==========================================#


# Converts the mean and variance of a lognormal loss into the meanlog sdlog 
log_params <- function(mean,var) {
  sdlog_2 <- log(1 + var/ mean ^ 2)
  sdlog <- sqrt(sdlog_2)
  meanlog <- log(mean) - 0.5 * sdlog_2
  list(meanlog = meanlog, sdlog = sdlog)
}


# Parameters for the 3 independent business lines
L1_par <- log_params(mu_L1, var_L1)
L2_par <- log_params(mu_L2, var_L2)
L3_par <- log_params(mu_L3, var_L3)

L1_par
L2_par
L3_par

# Extract lognormal parameters into simple variables
mean_log1 <- L1_par$meanlog
sd_log1   <- L1_par$sdlog

mean_log2 <- L2_par$meanlog
sd_log2   <- L2_par$sdlog

mean_log3 <- L3_par$meanlog
sd_log3   <- L3_par$sdlog



#===================================#
# 8. Gamma (shape,scale) parameters #   
#===================================#

# A_n ~ Gamma(shape_A, scale_A) with given mean and variance as follows:
shape_A <- mu_A ^2/ var_A
scale_A <- var_A / mu_A

shape_A
scale_A


#===================#
# 9. No.simulations #
#===================#

# Calculate required number of simulations using π1 with 0.1% margin of error at 95% CI
target_margin <- 0.001
z_value <- qnorm(1 - 0.05/2) #1.959964
z_value

# Calculate the required number of simulations needed 
required_n <- ceiling(((z_value / target_margin)^2) * pi_1 * (1-pi_1))
cat(sprintf("Required number of simulations for 0.1%% accuracy: %d\n", required_n))

# Use the exact required number of simulations
nsim <- required_n


#==========================#
# 10. Monte Carlo Function #
#==========================#

# Compute the Monte Carlo estimate, standard error, confidence interval
MonteCarlo <- function(Z, alpha = 0.05)
{
  nsim  <- length(Z)
  theta <- mean(Z)
  SE    <- sd(Z) / sqrt(nsim)
  CI    <- theta + c(-1, 1) * qnorm(1 - alpha / 2) * SE
  return(list(theta = theta, SE = SE, CI = CI))
}


#===========================#
# 11. Simulate Fund Process #
#===========================#

# Simulate the insurer's fund over N quarters for nsim independent paths.
# The update follows the random walk: F_{n} = F_{n-1} + A_n - L_n.
sim_fund <- function(N = 40, nsim,
                     F0,
                     mean_log1, sd_log1,
                     mean_log2, sd_log2,
                     mean_log3, sd_log3,
                     shape_A, scale_A)
{
  F <- matrix(0, nrow = N + 1, ncol = nsim)
  F[1, ] <- F0
  
  for (n in 1:N)
  {
    L1 <- rlnorm(nsim, meanlog = mean_log1, sdlog = sd_log1)
    L2 <- rlnorm(nsim, meanlog = mean_log2, sdlog = sd_log2)
    L3 <- rlnorm(nsim, meanlog = mean_log3, sdlog = sd_log3)
    L_n <- L1 + L2 + L3
    
    A_n <- rgamma(nsim, shape = shape_A, scale = scale_A)
    
    F[n + 1, ] <- F[n, ] + A_n - L_n
  }
  
  return(F)
}

#====================#
# 12. Run Simulation #
#====================#

# Simulate nsim fund paths over N quarters
F_sim <- sim_fund(N = N, nsim = nsim,
                  F0 = F0,
                  mean_log1 = mean_log1, sd_log1 = sd_log1,
                  mean_log2 = mean_log2, sd_log2 = sd_log2,
                  mean_log3 = mean_log3, sd_log3 = sd_log3,
                  shape_A = shape_A, scale_A = scale_A)


#=================#
# 13. Estimate π1 #                  
#=================#

# π1 = P(F_40 < 0) ruin at final horizon
Z1 <- as.numeric(F_sim[N + 1, ] < 0)
# Monte Carlo estimate, standard error and 95% confidence Interval for π1.
MC_pi1  <- MonteCarlo(Z1)

cat(sprintf("Monte Carlo estimate of π1 (ruin at horizon): %.5f\n", MC_pi1$theta))
cat(sprintf("95%% CI for π1: [%.5f, %.5f]\n", MC_pi1$CI[1], MC_pi1$CI[2]))


#=================#
# 14. Estimate π2 #                  
#=================#

# π2 = P(at least one year-end fund < 0)

# Year-ends are at quarters 4, 8, ..., 40 
idx_year_end <- 4 * (1:10) + 1

# Checks whether any year-end fund is negative for each simulated path
Z2 <- apply(F_sim[idx_year_end, ], 2,
            function(x) as.numeric(any(x < 0)))

# Monte Carlo estimate, standard error and 95% confidence Interval for π2.
MC_pi2  <- MonteCarlo(Z2)

cat(sprintf("Monte Carlo estimate of π2 (ruin at a year-end): %.5f\n", MC_pi2$theta))
cat(sprintf("95%% CI for π2: [%.5f, %.5f]\n", MC_pi2$CI[1], MC_pi2$CI[2]))



###################### Graph ######################################
#===================#
# Plot sample paths #                    
#===================#

num_paths <- 300
idx_paths <- sample(1:nsim, size = num_paths)

time <- 0:N

# plot first path and set up axes
plot(time, F_sim[, idx_paths[1]],
     type = "l",
     xlab = "Quarter",
     ylab = "Fund value",
     font.lab = 2,
     font.axis = 1,
     ylim = range(F_sim[, idx_paths]),
     main = "Random walk fund simulations (300 sample paths)",
     col = rgb(0, 0, 1, 0.1))

# add the remaining paths
for (j in idx_paths[-1])
{
  lines(time, F_sim[, j], col = rgb(0, 0, 1, 0.1))
}

# ruin threshold
abline(h = 0, col = "red", lwd = 2, lty = 2)


## Part C: Recalculating pi1 for different values of rho


#===============================#
# 15. Lognormal loss parameters #                  
#===============================#

# Reuse meanlog and sdlog values from Part (b)
m <- c(mean_log1, mean_log2, mean_log3)
s <- c(sd_log1,   sd_log2,   sd_log3)

m
s


#===========================#
# 16. Grid of rho and setup #
#===========================#


# Choose nsim for Part C (using the result from part b)
nsim <- required_n   

rhos <- seq(-1, 1, by = 0.1)

# Data frame to store pi1( rho ) for each correlation level
results <- data.frame(rho = rhos, pi1 = NA)




# mean and variance of total loss L_n    #
# (needed to parameterise gamma A_n)     #

E1 <- exp(m[1] + 0.5 * s[1]^2)
E2 <- exp(m[2] + 0.5 * s[2]^2)
E3 <- exp(m[3] + 0.5 * s[3]^2)

Var1 <- (exp(s[1]^2) - 1) * exp(2 * m[1] + s[1]^2)
Var2 <- (exp(s[2]^2) - 1) * exp(2 * m[2] + s[2]^2)
Var3 <- (exp(s[3]^2) - 1) * exp(2 * m[3] + s[3]^2)


#========================================#
# 17. Loop over rho and estimate pi1(rho)
#========================================#

for (j in seq_along(rhos)) {
  
  rho <- rhos[j]
  cat("Running rho =", rho, "\n")
  
  # Simulate standard normals for each line and quarter
  Z1 <- matrix(rnorm(nsim * N), nsim, N)
  Z2 <- matrix(rnorm(nsim * N), nsim, N)
  Z3 <- matrix(rnorm(nsim * N), nsim, N)
  
  # Correlated driver for L2 with correlation rho with L1
  Z2_corr <- rho * Z1 + sqrt(max(0, 1 - rho^2)) * Z2
  
  # Lognormal losses for each line
  L1 <- exp(m[1] + s[1] * Z1)
  L2 <- exp(m[2] + s[2] * Z2_corr)
  L3 <- exp(m[3] + s[3] * Z3)
  
  Lsum <- L1 + L2 + L3   # total loss per quarter
  
  # Covariance between L1 and L2 due to correlation rho
  Cov12 <- exp(m[1] + m[2] + 0.5 * (s[1]^2 + s[2]^2) + s[1] * s[2] * rho) - E1 * E2
  
  E_L   <- E1 + E2 + E3
  Var_L <- Var1 + Var2 + Var3 + 2 * Cov12
  
  
  #Asset distribution A_n ~ Gamma(shape,scale)
  mu_A_rho <- E_L * (1 + lambda)
  var_A_rho  <- Var_L * kappa
  
  shapeA_rho <- mu_A_rho^2 / var_A_rho
  scaleA_rho <- var_A_rho / mu_A_rho
  
  A <- matrix(rgamma(nsim * N, shape = shapeA_rho, scale = scaleA_rho),
              nsim, N)
  
  # Fund at horizon N
  F0_c <- alpha * mu_A_rho
  FN <- F0_c + rowSums(A - Lsum)
  
  # Ruin probability at time N for this rho
  results$pi1[j] <- mean(FN < 0)
}

#====================#
# 18. Print and plot #
#====================#

print(results)

plot(results$rho, results$pi1, pch = 19, col = "black",
     xlab = "Correlation rho",
     ylab = "Estimated ruin probability π1",
     font.lab = 2, 
     main = "π1 vs correlation")

fit <- smooth.spline(results$rho, results$pi1, spar = 0.7)
lines(fit, col = "red", lwd = 2)

# Part D i)
#==================================#
# 19. Define the transition matrix #
#==================================#

# Define the 4x4 transition matrix P for the Markov chain with states:
# N(neutral), G(good), B(bad), T(terrible)
P <- matrix(
  c(0.75, 0.15, 0.05, 0.05, # row for state N
    0.90, 0.10, 0.00, 0.00, # row for state G
    0.50, 0.00, 0.45, 0.05, # row for state B
    0.50, 0.00, 0.50, 0.00), # row for state T
  nrow = 4, # 4 states so 4 rows
  byrow = TRUE # fill matrix row by row
)


# Give names to rows/columns to make output easier to read
rownames(P) <- c("N", "G", "B", "T")
colnames(P) <- c("N", "G", "B", "T")

P

#======================#
# 20. Ergodicity check #
#======================#

# Helps raise a matrix to a power of n 
powmat <- function(X, n)
{
  eg     <- eigen(X, symmetric = FALSE)
  LAMBDA <- diag(eg$values)
  V      <- eg$vectors
  Vinv   <- solve(V)
  
  V %*% (LAMBDA^n) %*% Vinv
}

# Compute P^4 and check if all entries are strictly positive
P4 <- powmat(P, 4)

# Ergodicity check 
all(P4 > 0)      # TRUE => chain is ergodic


#=============================#
# 21. Stationary Distribution #
#=============================#

# Solve for pi such that pi^T P = pi^T and sum(pi) = 1 using least squares
find_stationary_distr <- function(P)
{
  K <- nrow(P)
  
  
  A <- rbind(t(P - diag(K)), rep(1, K))
  b <- c(rep(0, K), 1)
  
  # least-squares solution
  pi <- t(solve(t(A) %*% A, t(A) %*% b))
  pi
}

# Compute stationary distribution
pi_stat <- find_stationary_distr(P)
pi_stat




# D (ii) Markov fund 

#==================================#
# 22. Basic parameters (as above)  #
#==================================#

# Total mean and variance of quarterly loss L_n
EL <- mu_L
VL <- var_L


# Initial fund 
# F0 = alpha * E[A_n] with base lambda = 0.04

F0_D <- F0


#======================================#
# 23. State–dependent lambda and kappa #
#======================================#

# States in order: N, G, B, T
lambda_vec <- c(0.04, 0.06, 0.02, 0.02)
kappa_vec  <- c(0.75, 0.60, 0.75, 1.00)


#==============================================#
# 24. Simulate Markov chain of economic states #
#==============================================#

# Simulates a sequence of economic states for N quarters from the 4-state Markov chain.
simulate_states <- function(N, P, start_state = 1)
{
  # N          : number of quarters
  # P          : 4 x 4 transition matrix
  # start_state: 1 = N, 2 = G, 3 = B, 4 = T
  
  S <- integer(N)
  S[1] <- start_state
  
  for (n in 2:N)
  {
    current <- S[n - 1]
    S[n] <- sample(1:4, size = 1, prob = P[current, ])
  }
  
  return(S)
}


#===========================================#
# 25. Simulate one total quarterly loss L_n #
#===========================================#

# Simulates the total quarterly loss by generating the 3 independent business line losses from their lognormal distributions.
simulate_Ln <- function()
{
  # use lognormal parameters from Part B
  L1 <- rlnorm(1, meanlog = mean_log1, sdlog = sd_log1)
  L2 <- rlnorm(1, meanlog = mean_log2, sdlog = sd_log2)
  L3 <- rlnorm(1, meanlog = mean_log3, sdlog = sd_log3)
  L  <- L1 + L2 + L3
  return(L)
}


#==================================================#
# 26. Simulate one fund path under the Markov Model #
#==================================================#

# Simulates a single trajectory of the insurer’s fund F_n under the Markov model

simulate_fund_one <- function(N, F0, P, lambda_vec, kappa_vec)
{
  # Simulate state path (1= N, 2= G, 3= B, 4= T)
  states <- simulate_states(N, P, start_state = 1)
  
  # Allocate vector for F_n, n = 0,...,N
  F <- numeric(N + 1)
  F[1] <- F0
  
  # Evolve fund
  for (n in 1:N)
  {
    s <- states[n]        
    
    lambda_n <- lambda_vec[s]
    kappa_n  <- kappa_vec[s]
    
    # Simulate total loss this quarter
    L_n <- simulate_Ln()
    
    # State–dependent asset distribution:
    meanA <- EL * (1 + lambda_n)
    varA  <- VL * kappa_n
    
    shapeA <- meanA^2 / varA
    scaleA <- varA / meanA
    
    A_n <- rgamma(1, shape = shapeA, scale = scaleA)
    
    # update fund - random walk
    F[n + 1] <- F[n] + A_n - L_n
  }
  
  return(F) 
}


#==================================================#
# 27. Monte Carlo estimate of pi1 under Markov Model#
#==================================================#

# Use Monte Carlo from Part b

nsim_D <- required_n    # number of simulated paths for Part (d)

# simulate nsim_D fund paths and store final F_40
F40_D <- replicate(nsim_D, {
  Fpath <- simulate_fund_one(N, F0_D, P, lambda_vec, kappa_vec)
  tail(Fpath, 1)
})

# Indicator of ruin at time 40
ruined_D <- as.numeric(F40_D < 0)


MC_pi1_D <- MonteCarlo(ruined_D)

cat(sprintf("Monte Carlo estimate of pi1 under Markov model: %.5f\n", MC_pi1_D$theta))
cat(sprintf("95%% CI for pi1 (Markov model): [%.5f, %.5f]\n",
            MC_pi1_D$CI[1], MC_pi1_D$CI[2]))
cat(sprintf("Standard error (SE): %.5f\n", MC_pi1_D$SE))


#======================================#
# PIE CHART OF STATIONARY DISTRIBUTION
#======================================#

# Stationary distribution vector 
stationary <- c(0.73, 0.12, 0.11, 0.04)
# labels for states
labels <- c("Neutral (N)", "Good (G)", "Bad (B)", "Terrible (T)")

# colours for states
cols <- c("blue", "green", "orange", "red")

# create pie chart
pie(stationary,
    labels = paste(labels, "\n", round(stationary*100, 1), "%"),
    col = cols,
    main = "Stationary Distribution of Markov States")
