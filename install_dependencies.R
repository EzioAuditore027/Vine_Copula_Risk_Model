# install_dependencies.R
packages <- c(
  "quantmod",
  "rugarch",
  "VineCopula",
  "copula",
  "PerformanceAnalytics",
  "ggplot2",
  "reshape2",
  "evir",
  "Rglpk"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

print("All dependencies installed successfully.")