############################################################
#
# R Examples for Survey Block Model and Uncertainty Quantification
# - Focus on Moments (e.g., variance), Shannon Entropy, and the Expected Value of Perfect Information
# - Utilizes symbolic algebra to compute relevant statistics
#
############################################################


library(symengine)
library(ggplot2)
library(magrittr)
library(gridExtra)

# import functions

source("funcs.R")


## Survey Block Model and Relevant Statistics

# Set the model parameters
p <- 0.4         # Probability of pest presence
rho <- -log(1 - p) # Assumed prior mean for X0, given p
delta <- 0.1     # Survey sensitivity level
lambda <- 1.5       # Population growth rate (>1)
tr <- 0.95        # Target probability of pest absence



# Estimate the block size kappa with the target probability of absence tr (T)

# Survey block size calculation
k <- kappa(rho, lambda, delta, tr)
k

# Define model parameters, now including rho
parms <- list(rho = rho, lambda = lambda, delta = delta)

#--------------
# Expected number (and variance) of mop-ups given block size k

# Generate the probability generating function (PGF) for mop-ups
pgf <- make_pgfM(k)
calc_moments(pgf, parms)

# Calculate the probability mass function (PMF) for mop-ups
pmfM <- calc_pmf(pgf, parms, support = 0:15)
# Check if PMF sums to 1 (ensures it's a valid probability distribution)
if (abs(sum(pmfM$p) - 1) > 1e-6) {
  warning("PMF does not sum to 1; consider extending the support.")
}
pmfM

#--------------
# Survey number of first detection given block size k

# Generate the PGF for the first detection
pgf <- make_pgfK(k)
calc_moments(pgf, parms)

# Calculate the PMF for the first detection
pmfK <- calc_pmf(pgf, parms, support = 0:15)
# Check if PMF sums to 1 (validity check)
if (abs(sum(pmfK$p) - 1) > 1e-6) {
  warning("PMF does not sum to 1; consider extending the support.")
}
pmfK

# Generate and calculate moments again for case
pgf <- make_pgfK(k)
calc_moments(pgf, parms)

#--------------
# Expected number of surveys (and variance) given block size k

pgf <- make_pgfN(k)
calc_moments(pgf, parms)


# Calculate the PMF for the number of surveys
pmfN <- calc_pmf(pgf, parms, support = 0:25)
# Check if PMF sums to 1 (validity check)
if (abs(sum(pmfN$p) - 1) > 1e-6) {
  warning("PMF does not sum to 1; consider extending the support.")
}
pmfN

#--------------
# Calculate entropy given prior distributions on p, lambda, and delta
# Here, we use uniform priors with an equidistant of 0.01 over a specified range

# Define the ranges for the priors
p_vals <- seq(0.01, 0.8, by = 0.01)      # Probability values ranging from 0.01 to 0.8
delta_vals <- seq(0.01, 0.1, by = 0.01)  # Survey sensitivity values from 0.01 to 0.1
lambda_vals <- seq(1.0, 1.2, by = 0.01)  # Population growth rates from 1.0 to 1.2
tr <- 0.95                               # Target probability of pest absence
support <- 0:30                          # Support range for calculations

# Calculate the entropy using the defined prior values
val <- calculate_entropy(p_vals, delta_vals, lambda_vals, tr, support)
val



#--------------
# Calculate the Expected Value of Perfect Information (EVPI)
# Here, p and lambda are uncertain, while delta is treated as a control variable.
# The kappa Max strategy is applied to optimize the decision-making process.

# Define the ranges for p, delta, and lambda
p_vals <- seq(0.01, 0.8, by = 0.01)      # Probability values ranging from 0.01 to 0.8
delta_vals <- seq(0.01, 0.1, by = 0.01)  # Survey sensitivity values from 0.01 to 0.1
lambda_vals <- seq(1.0, 1.2, by = 0.01)  # Population growth rates from 1.0 to 1.2

# Define cost parameters
cm <- 200000   # Cost of mop-up
cs <- 112000   # Cost of the survey
tr <- 0.95     # Target probability of pest absence

# Calculate the Expected Value of Perfect Information (EVPI)
evpi <- calculate_evpi_p(p_vals, delta_vals, lambda_vals, tr, cs, cm)
evpi



