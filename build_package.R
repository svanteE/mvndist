#!/usr/bin/env Rscript

# Set working directory to package root
setwd("c:/Users/svante/OneDrive - Aalborg Universitet/geodesic/mvndist")

# Load required packages quietly
suppressMessages({
  library(roxygen2)
})

# Generate documentation using roxygen2 directly
cat("Generating documentation...\n")
roxygen2::roxygenize(".", roclets = c("rd", "namespace"))

cat("\nDocumentation generated successfully!\n")
cat("Check the man/ folder for .Rd files.\n")
cat("\nTo check and install the package, run:\n")
cat("  rcmdcheck::rcmdcheck()\n")
cat("  install.packages('.', repos = NULL, type = 'source')\n")
