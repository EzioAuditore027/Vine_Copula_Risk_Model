# ==============================================================================
# Project: Multivariate Risk Modeling with GARCH-EVT and Vine Copulas
# Author: Himanshu Shelke
# Description: 
#   1. Models marginal volatility using GARCH (GJR or Standard) with Student-t innovations.
#   2. Models dependence using R-Vine Copulas (handling non-linear tail dependence).
#   3. Optimizes a portfolio by minimizing Conditional Value-at-Risk (CVaR).
# ==============================================================================

# --- STEP 0: SETUP & LIBRARIES ---
# Install missing packages if needed:
# install.packages(c("quantmod", "rugarch", "VineCopula", "PerformanceAnalytics", "evd", "network"))

library(quantmod)             # Data fetching
library(rugarch)              # GARCH modeling
library(VineCopula)           # Vine Copula modeling (pobs, RVineStructureSelect)
library(PerformanceAnalytics) # Risk metrics (ES/CVaR)
library(network)              # Dependency for plotting vine trees

# --- STEP 1: DATA EXTRACTION ---
# We select 3 Tech stocks to demonstrate high correlation and potential tail dependence.
tickers <- c("AAPL", "MSFT", "GOOG") 
startDate <- "2020-01-01"
endDate <- "2024-01-01"

# Fetch data and calculate Log Returns
getSymbols(tickers, from = startDate, to = endDate)
prices <- do.call(merge, lapply(tickers, function(x) Ad(get(x))))
returns <- na.omit(Return.calculate(prices, method = "log"))
colnames(returns) <- tickers

# Visual Check: Volatility Clustering
par(mfrow=c(1,1))
chart.RollingPerformance(returns, width = 22, FUN = "sd.annualized", 
                         main = "Rolling Volatility (Evidence for GARCH)")

# --- STEP 2: MARGINAL MODELING (GARCH-EVT) ---
# We fit a GARCH(1,1) with Student-t errors to capture heavy tails.

# Define the GARCH Specification
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std" # Student-t distribution for heavy tails
)

# Initialize storage
u_matrix <- matrix(nrow = nrow(returns), ncol = ncol(returns)) # For Pseudo-Observations
colnames(u_matrix) <- tickers
garch_models <- list()

# Loop through assets: Fit GARCH -> Extract Residuals -> Transform to Uniform
print("Fitting GARCH models...")
for(i in 1:ncol(returns)){
  # 1. Fit GARCH
  fit <- ugarchfit(spec = spec, data = returns[,i])
  garch_models[[i]] <- fit
  
  # 2. Extract Standardized Residuals (The "Shocks")
  z <- residuals(fit, standardize = TRUE)
  
  # 3. Transform to Uniform(0,1) using Empirical CDF (pobs)
  #    This serves as the input for the Copula.
  u_matrix[,i] <- pobs(z)
}

# --- STEP 3: DEPENDENCE MODELING (VINE COPULA) ---
# We use AIC to select the best tree structure and families (e.g., Clayton, t, Gumbel).

print("Fitting Vine Copula...")
# type = 0 selects "R-Vine" (General structure)
vine_fit <- RVineStructureSelect(data = u_matrix, 
                                 type = 0, 
                                 selectioncrit = "AIC")

# Output the Summary (Look for 't', 'Clayton', or 'BB8' families)
print("--- Vine Copula Model Summary ---")
summary(vine_fit)

# Visualize the Tree Structure (Tree 1 = Strongest dependencies)
par(mfrow=c(1,1))
RVineTreePlot(vine_fit, tree = 1, edge.labels = "family")
title("Tree 1: Strongest Direct Dependencies")

# --- STEP 4: SIMULATION & OPTIMIZATION (CVaR) ---
# Simulate future market scenarios to estimate downside risk.

print("Simulating Scenarios...")
N_sim <- 10000
sim_u <- RVineSim(N_sim, vine_fit) # Simulate Uniforms from Vine

# Inverse Transform: Uniforms -> Residuals -> Returns
sim_returns <- matrix(nrow = N_sim, ncol = ncol(returns))
colnames(sim_returns) <- tickers

for(i in 1:ncol(returns)){
  # Forecast volatility for T+1
  fc <- ugarchforecast(garch_models[[i]], n.ahead = 1)
  mu_next <- fitted(fc)[1]
  sigma_next <- sigma(fc)[1]
  
  # Convert Uniforms back to residuals using the fitted Student-t shape
  shape_param <- coef(garch_models[[i]])["shape"]
  sim_residuals <- qdist("std", p = sim_u[,i], shape = shape_param)
  
  # Calculate simulated return
  sim_returns[,i] <- mu_next + sigma_next * sim_residuals
}

# Portfolio Optimization: Calculate 95% CVaR for an Equal-Weighted Portfolio
weights <- rep(1/ncol(returns), ncol(returns)) 
port_returns <- sim_returns %*% weights

# Fix: Convert to vector to avoid 'xts' date error in ES function
port_vec <- as.vector(port_returns)

# Calculate Conditional Value at Risk (Expected Shortfall)
cvar_95 <- ES(port_vec, p = 0.95, method = "historical")

# --- FINAL OUTPUT ---
print(paste("Predicted 1-Day Portfolio CVaR (95%):", round(cvar_95, 5)))

# Visualization of Risk
hist(port_vec, breaks=50, main="Simulated Portfolio Returns (GARCH-Vine)", 
     col="lightblue", border="white", xlab="Returns")
abline(v = cvar_95, col="red", lwd=2, lty=2)
legend("topleft", legend=paste("CVaR 95%:", round(cvar_95, 4)), 
       col="red", lty=2, lwd=2)
