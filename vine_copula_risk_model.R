# ==============================================================================
# Project: Multivariate Risk Modeling with GARCH-EVT and Vine Copulas
# Methodology: Based on Karmakar & Paul (2023) "Downside Risk of Energy Stocks"
#
# Description:
#   1. Marginal Modeling: GJR-GARCH (Leverage Effect) + Skewed-t (Fat Tails).
#   2. Tail Modeling: Extreme Value Theory (EVT) using Generalized Pareto (GPD).
#   3. Dependence: R-Vine Copula with restricted financial families.
#   4. Optimization: Min-CVaR Portfolio Optimization using Linear Programming.
# ==============================================================================

# --- STEP 0: SETUP & LIBRARIES ---
# install.packages(c("quantmod", "rugarch", "VineCopula", "PerformanceAnalytics", "evir", "Rglpk"))

library(quantmod)             # Data fetching
library(rugarch)              # GARCH modeling
library(VineCopula)           # Vine Copula modeling
library(evir)                 # [UPGRADE] Extreme Value Theory (GPD fitting)
library(Rglpk)                # [UPGRADE] Linear Programming Solver (Min-CVaR)
library(PerformanceAnalytics) # Risk metrics

# --- STEP 1: DATA EXTRACTION ---
tickers <- c("AAPL", "MSFT", "GOOG") 
startDate <- "2020-01-01"
endDate <- "2024-01-01"

getSymbols(tickers, from = startDate, to = endDate)
prices <- do.call(merge, lapply(tickers, function(x) Ad(get(x))))
returns <- na.omit(Return.calculate(prices, method = "log"))
colnames(returns) <- tickers

# --- STEP 2: MARGINAL MODELING (GJR-GARCH + SKEWED-T) ---
# [UPGRADE 1]: Switched to 'gjrGARCH' to capture Leverage Effects (asymmetry).
# [UPGRADE 2]: Switched to 'sstd' (Skewed Student-t) for non-symmetric fat tails.

spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), 
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd" 
)

z_matrix <- matrix(nrow = nrow(returns), ncol = ncol(returns)) # Standardized Residuals
garch_models <- list()

print("STEP 2: Fitting GJR-GARCH-sstd models...")
for(i in 1:ncol(returns)){
  fit <- ugarchfit(spec = spec, data = returns[,i])
  garch_models[[i]] <- fit
  z_matrix[,i] <- residuals(fit, standardize = TRUE)
}

# --- STEP 3: EVT TAIL MODELING (THE "MISSING LINK") ---
# [UPGRADE 3]: Inserted EVT step. Instead of assuming the whole distribution fits one curve,
# we fit a Generalized Pareto Distribution (GPD) specifically to the tails (Peak Over Threshold).

# Helper Function: Semi-Parametric CDF (Empirical Center + GPD Tails)
p_semi_parametric <- function(x, u_lower, u_upper, gpd_lower, gpd_upper) {
  n <- length(x)
  p <- numeric(n)
  
  # Lower Tail (x < u_lower): Use GPD
  idx_lower <- x < u_lower
  if(any(idx_lower)) {
    excess <- u_lower - x[idx_lower]
    # GPD CDF formula adaptation for lower tail
    p[idx_lower] <- (sum(idx_lower)/n) * (1 + gpd_lower$par.ests['xi'] * excess / gpd_lower$par.ests['beta'])^(-1/gpd_lower$par.ests['xi'])
  }
  
  # Upper Tail (x > u_upper): Use GPD
  idx_upper <- x > u_upper
  if(any(idx_upper)) {
    excess <- x[idx_upper] - u_upper
    # GPD Survival Function for upper tail
    p[idx_upper] <- 1 - (sum(idx_upper)/n) * (1 + gpd_upper$par.ests['xi'] * excess / gpd_upper$par.ests['beta'])^(-1/gpd_upper$par.ests['xi'])
  }
  
  # Interior: Use Empirical Rank (pobs)
  idx_mid <- !idx_lower & !idx_upper
  p[idx_mid] <- rank(x)[idx_mid] / (n + 1)
  
  return(p)
}

# Storage for EVT models and Uniforms
gpd_models <- list()
u_matrix <- matrix(nrow = nrow(returns), ncol = ncol(returns))
threshold_prob <- 0.10 # Top/Bottom 10% are considered "Extreme"

print("STEP 3: Fitting EVT (GPD) to tails...")
for(i in 1:ncol(returns)){
  z <- z_matrix[,i]
  u_lower <- quantile(z, threshold_prob)
  u_upper <- quantile(z, 1 - threshold_prob)
  
  # Fit GPD (Note: 'evir' fits positive excesses, so flip signs for lower tail)
  fit_lower <- gpd(-z, threshold = -u_lower)
  fit_upper <- gpd(z, threshold = u_upper)
  
  gpd_models[[i]] <- list(u_l = u_lower, u_u = u_upper, fit_l = fit_lower, fit_u = fit_upper)
  u_matrix[,i] <- p_semi_parametric(z, u_lower, u_upper, fit_lower, fit_upper)
}

# --- STEP 4: DEPENDENCE MODELING (RESTRICTED VINE) ---
# [UPGRADE 4]: Restricted families to standard financial types (t, Clayton, Gumbel)
# to avoid overfitting with obscure copulas and improve computation speed.

family_set <- c(1, 2, 3, 4, 5) # Gaussian, t, Clayton, Gumbel, Frank
print("STEP 4: Fitting Vine Copula...")
vine_fit <- RVineStructureSelect(data = u_matrix, type = 0, selectioncrit = "AIC", familyset = family_set)
print(summary(vine_fit))

# --- STEP 5: SIMULATION (INVERSE EVT) ---
print("STEP 5: Simulating 10,000 scenarios...")
N_sim <- 10000
sim_u <- RVineSim(N_sim, vine_fit)
sim_returns <- matrix(0, nrow=N_sim, ncol=ncol(returns))
colnames(sim_returns) <- tickers

# Helper Function: Inverse Semi-Parametric CDF (Uniform -> Residuals)
q_semi_parametric <- function(p, u_lower, u_upper, gpd_lower, gpd_upper, original_z) {
  n <- length(original_z)
  k_l <- sum(original_z < u_lower)
  k_u <- sum(original_z > u_upper)
  res <- numeric(length(p))
  
  # Inverse Lower Tail
  idx_l <- p < (k_l/n)
  if(any(idx_l)) {
    xi <- gpd_lower$par.ests['xi']; beta <- gpd_lower$par.ests['beta']
    res[idx_l] <- u_lower - (beta/xi) * (( (p[idx_l] / (k_l/n))^(-xi) ) - 1)
  }
  
  # Inverse Upper Tail
  idx_u <- p > (1 - k_u/n)
  if(any(idx_u)) {
    xi <- gpd_upper$par.ests['xi']; beta <- gpd_upper$par.ests['beta']
    res[idx_u] <- u_upper + (beta/xi) * (( ((1 - p[idx_u]) / (k_u/n))^(-xi) ) - 1)
  }
  
  # Inverse Interior (Quantile of original data)
  idx_m <- !idx_l & !idx_u
  res[idx_m] <- quantile(original_z, probs = p[idx_m], type = 1)
  
  return(res)
}

for(i in 1:ncol(returns)){
  # 1. Forecast GARCH Volatility
  fc <- ugarchforecast(garch_models[[i]], n.ahead = 1)
  mu <- fitted(fc)[1]; sig <- sigma(fc)[1]
  
  # 2. Inverse EVT (Uniform -> Z)
  gpd <- gpd_models[[i]]
  sim_z <- q_semi_parametric(sim_u[,i], gpd$u_l, gpd$u_u, gpd$fit_l, gpd$fit_u, z_matrix[,i])
  
  # 3. Simulate Returns
  sim_returns[,i] <- mu + sig * sim_z
}

# --- STEP 6: MIN-CVAR OPTIMIZATION (LINEAR PROGRAMMING) ---
# [UPGRADE 5]: Replaced equal-weights with a Solver to find optimal Min-CVaR weights.

print("STEP 6: Optimizing Portfolio (Min-CVaR)...")
alpha <- 0.95
S <- nrow(sim_returns)
J <- ncol(sim_returns)

# LP Formulation for Min-CVaR (Rockafellar & Uryasev)
# Obj: Min [gamma + (1 / ((1-alpha)*S)) * sum(z_s)]
obj_coeffs <- c(rep(0, J), 1, rep(1/((1-alpha)*S), S))

# Constraints:
# 1. Scenarios: Returns*w + gamma + z_s >= 0
mat_scenarios <- cbind(sim_returns, 1, diag(S))
# 2. Budget: sum(w) = 1
mat_budget <- c(rep(1, J), 0, rep(0, S))

mat <- rbind(mat_scenarios, mat_budget)
dir <- c(rep(">=", S), "==")
rhs <- c(rep(0, S), 1)

# Solve
res <- Rglpk_solve_LP(obj = obj_coeffs, mat = mat, dir = dir, rhs = rhs, max = FALSE)
opt_weights <- res$solution[1:J]
names(opt_weights) <- tickers

print("--- Optimal Min-CVaR Weights ---")
print(round(opt_weights, 4))

# Validation
ew_cvar <- ES(as.vector(sim_returns %*% rep(1/J, J)), p=0.95, method="historical")
opt_cvar <- ES(as.vector(sim_returns %*% opt_weights), p=0.95, method="historical")

print(paste("Risk Reduction vs Equal Weight:", round((1 - opt_cvar/ew_cvar)*100, 2), "%"))

# --- VISUALIZATION ---
# Compare Equal Weight vs Optimized
par(mfrow=c(1,1))
d_ew <- density(sim_returns %*% rep(1/J, J))
d_opt <- density(sim_returns %*% opt_weights)

plot(d_ew, main="Portfolio Return Density: Equal vs Optimized", col="red", lwd=2, ylim=c(0, max(d_opt$y)))
lines(d_opt, col="green", lwd=2)
abline(v = -opt_cvar, col="green", lty=2)
abline(v = -ew_cvar, col="red", lty=2)
legend("topleft", legend=c("Equal Weight", "Min-CVaR Optimized"), col=c("red", "green"), lwd=2)