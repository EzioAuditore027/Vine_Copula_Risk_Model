# Multivariate Risk Modeling with GARCH-EVT and Vine Copulas

This repository contains an R implementation of a robust portfolio risk model. It combines **GARCH-EVT** (Generalized Autoregressive Conditional Heteroskedasticity with Extreme Value Theory) for marginal filtering and **Vine Copulas** for modeling complex, non-linear dependence structures between assets.

The project demonstrates a robust pipeline for estimating **Conditional Value-at-Risk (CVaR)** in portfolios where standard Gaussian assumptions fail to capture tail dependence.

## 📌 Project Overview

Traditional Mean-Variance optimization often underestimates risk during market stress because it assumes linear correlations (Gaussian Copula). This project addresses that limitation by:

1.  **Filtering Volatility:** Using GARCH(1,1) models with Student-t innovations to handle volatility clustering and heavy tails in marginal distributions.
2.  **Modeling Dependence:** Constructing an **R-Vine Copula** (Regular Vine) to capture asymmetric tail dependence (e.g., assets crashing together) without forcing a specific parametric family.
3.  **Risk Optimization:** Using Monte Carlo simulations (N=10,000) from the fitted Vine structure to estimate and minimize **95% CVaR** (Expected Shortfall).

## 🛠️ Methodology

The statistical pipeline implemented in `Vine_Copula_Risk_Model.R` follows these steps:

### 1. Marginal Estimation (GARCH-EVT)
* **Model:** GARCH(1,1) with `std` (Student-t) distribution.
* **Goal:** Remove serial correlation and volatility clustering.
* **Transformation:** Standardized residuals are transformed into Pseudo-Observations (Uniform $[0,1]$) using the Empirical CDF (Probability Integral Transform).

### 2. Dependence Modeling (Vine Copulas)
* **Structure:** R-Vine (Regular Vine) decomposition.
* **Selection:** Structure and pair-copula families (e.g., Clayton, Gumbel, t-Copula) are selected automatically via **AIC** (Akaike Information Criterion).
* **Advantage:** Captures complex dependency patterns like **Tail Dependence** that simple correlation matrices miss.

### 3. Simulation & Risk Assessment
* **Simulation:** 10,000 scenarios are simulated from the fitted Vine Copula.
* **Inverse Transform:** Uniform simulations are mapped back to returns using the forecasted GARCH volatility and estimated shape parameters.
* **Metric:** 95% Conditional Value-at-Risk (CVaR/ES) is calculated for the portfolio.

## 📦 Requirements

The project is written in **R** and requires the following packages for statistical modeling and optimization:

```r
install.packages(c("quantmod", "rugarch", "VineCopula", "PerformanceAnalytics", "evd", "network"))