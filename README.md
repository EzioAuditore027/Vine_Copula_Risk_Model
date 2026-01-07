# Multivariate Risk Modeling with GARCH-EVT and Vine Copulas

This repository contains an R implementation of a robust portfolio risk model. It combines **GARCH-EVT** (Generalized Autoregressive Conditional Heteroskedasticity with Extreme Value Theory) for marginal filtering and **Vine Copulas** for modeling complex, non-linear dependence structures between assets.

The project demonstrates a rigorous pipeline for **Min-CVaR Portfolio Optimization** in energy/tech markets where standard Gaussian assumptions fail to capture tail dependence and asymmetry.

## 📌 Project Overview

Traditional Mean-Variance optimization often underestimates risk during market stress because it assumes linear correlations (Gaussian Copula) and symmetric volatility. This project addresses those limitations by:

1.  [cite_start]**Capturing Asymmetry:** Using **GJR-GARCH** to model "Leverage Effects" (where volatility spikes more after crashes than booms)[cite: 1596].
2.  [cite_start]**Modeling Extremes:** Applying **Extreme Value Theory (EVT)** with the Generalized Pareto Distribution (GPD) to strictly model tail risks ("Black Swans") beyond historical data[cite: 1345, 1348].
3.  [cite_start]**Modeling Dependence:** Constructing an **R-Vine Copula** to capture asymmetric tail dependence (e.g., assets crashing together) using a restricted set of financial copula families[cite: 1287].
4.  [cite_start]**Optimizing Risk:** Implementing **Min-CVaR Optimization** via Linear Programming to mathematically find the portfolio weights that minimize 95% Expected Shortfall[cite: 1864, 1865].

## 🛠️ Methodology

The statistical pipeline implemented in `Vine_Copula_Risk_Model.R` follows the framework of Karmakar & Paul (2023):

### 1. Marginal Estimation (GJR-GARCH + Skewed-t)
* [cite_start]**Model:** **GJR-GARCH(1,1)** with **Skewed Student-t (`sstd`)** innovations[cite: 1373, 1606].
* [cite_start]**Why:** Captures **Volatility Clustering** (periods of high/low variance), **Leverage Effects** (asymmetric reaction to bad news), and **Skewness** (negative return bias) inherent in energy and tech stocks[cite: 1259, 1596].

### 2. Tail Modeling (EVT - Peak Over Threshold)
* [cite_start]**Technique:** **Peak Over Threshold (POT)** method using the **Generalized Pareto Distribution (GPD)**[cite: 1346].
* [cite_start]**Threshold:** Applied to the upper and lower 10% of standardized residuals[cite: 1621].
* **Why:** Standard empirical distributions cannot simulate events worse than history. [cite_start]EVT fits a parametric curve to the tails, allowing the model to extrapolate extreme "Black Swan" scenarios[cite: 1274, 1348].

### 3. Dependence Modeling (Vine Copulas)
* **Structure:** Regular Vine (R-Vine) decomposition.
* [cite_start]**Selection:** Pair-copula families (Gaussian, Student-t, Clayton, Gumbel, Frank) selected via **AIC**[cite: 1547, 1548].
* [cite_start]**Advantage:** Captures non-linear dependence structures, such as **Tail Dependence** (the likelihood of simultaneous extreme losses), which simple correlation matrices miss[cite: 1281].

### 4. Optimization (Min-CVaR)
* [cite_start]**Objective:** Minimize 95% CVaR (Conditional Value-at-Risk) subject to full investment constraints[cite: 1865].
* [cite_start]**Method:** **Linear Programming (LP)** solver based on the Rockafellar & Uryasev (2000) framework.
* [cite_start]**Simulation:** 10,000 scenarios simulated from the fitted Vine Copula and re-transformed via Inverse-EVT and GARCH forecasts[cite: 1728, 1744].

## 📦 Requirements

The project is written in **R** and requires the following packages for statistical modeling, EVT, and optimization:

```r
install.packages(c("quantmod", "rugarch", "VineCopula", "PerformanceAnalytics", "evir", "Rglpk"))