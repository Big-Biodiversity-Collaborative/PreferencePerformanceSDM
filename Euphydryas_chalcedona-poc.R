# Test case on E. chalcedona
# Jeff Oliver
# jcoliver@arizona.edu
# 2022-06-08

# Download data from GBIF for E. chalcedona
# Download data from GBIF for two host plants
# Do geographic filtering and maybe date filtering...

# Download climate data from...Swallowtail project?

# Run SDM GLM on two host plants
# Run SDM GLM on E. chalcedona ~ climate data only
# Maybe do some lassoing here?
# Run SDM GLM on E. chalcedona ~ climate + hosts

# See if model with hosts is better? Why?
# Compare model loadings for two hosts
# try Wald test
# https://stats.stackexchange.com/questions/478408/compare-regression-coefficients-within-the-same-model
# A cool alternative would be a LRT:
# Complex model (different slopes): y ~ b0 + b1*x1 + b2*x2
# Simple model (identical slopes):  y ~ b0 + b3*(x1 + x2)
# https://andrewpwheeler.com/2016/10/19/testing-the-equality-of-two-regression-coefficients/
# https://stats.stackexchange.com/questions/211584/testing-linear-restriction-in-r/211597#211597
# or car::linear.hypothesis()
# https://stats.stackexchange.com/questions/228351/how-to-compare-coefficients-within-the-same-multiple-regression-model

# Run SDM GLM on E. chalcedona ~ hosts?
