# Test case on E. chalcedona
# Jeff Oliver
# jcoliver@arizona.edu
# 2022-06-08

require(spocc)   # Downloading from gbif; loaded by function...
require(dplyr)   # Data wrangling
require(ggplot2) # Checking data, exploratory
require(terra)   # Raster things

source(file = "src/download_gbif.R")
num_background <- 10000
set.seed(20220613)

# TODO: Need some (manual?) QC for observation points. There are three 
# questionable observations of D. aurantiacus (northern WA, UT or CO, and one 
# off the coast of Baja California)

# Download data from GBIF for Euphydryas chalcedona
# Download data from GBIF for two host plants
# Diplacus aurantiacus, Scrophularia californica
species <- c("Euphydryas chalcedona", "Diplacus aurantiacus", 
             "Scrophularia californica")
obs_list <- list()
for (one_sp in species) {
  nice_name <- tolower(x = gsub(pattern = " ",
                                replacement = "_",
                                x = one_sp))
  data_file <- paste0("data/gbif/", nice_name, "-gbif-raw.csv")
  if (!file.exists(data_file)) {
  download_gbif(species_name = one_sp, gbif_name = one_sp,
                verbose = TRUE, restrict_n_amer = TRUE)
  } else {
    message(paste0("Data for ", one_sp, " already on disk."))
  }
  # host or insect?
  sp_type <- "host"
  if (one_sp == "Euphydryas chalcedona") {
    sp_type <- "insect"
  }
  obs_list[[nice_name]][["sp_type"]] <- sp_type
  obs_list[[nice_name]][["obs"]] <- read.csv(file = data_file)
}

# Do geographic filtering and maybe date filtering...
# Restrict by D. aurantiacus, which has smallest range?
nice_names <- tolower(x = gsub(pattern = " ",
                               replacement = "_",
                               x = species))

# Drop anything east of -80 lon
for (nice_name in nice_names) {
  obs_list[[nice_name]][["obs"]] <- obs_list[[nice_name]][["obs"]] %>%
    filter(longitude < -80)
}

# While iterating, can keep track of lat/lon bounds
min_lon <- NA
max_lon <- NA
min_lat <- NA
max_lat <- NA
# Only retain records from 2000 to present
for (nice_name in nice_names) {
  total_obs <- nrow(obs_list[[nice_name]][["obs"]])
  obs_list[[nice_name]][["obs"]] <- obs_list[[nice_name]][["obs"]] %>%
    filter(year >= 2000)
  filtered_obs <- nrow(obs_list[[nice_name]][["obs"]])
  message(filtered_obs, " ", nice_name, " records after filtering (",
          total_obs, " before filtering)")
  # Get min/max lat/lon while iterating here
  # So ugly, JCO
  if (is.na(min_lon)) {
    min_lon <- min(obs_list[[nice_name]][["obs"]]$longitude)
  } else {
    min_lon <- min(min_lon, obs_list[[nice_name]][["obs"]]$longitude)
  }
  if (is.na(max_lon)) {
    max_lon <- max(obs_list[[nice_name]][["obs"]]$longitude)
  } else {
    max_lon <- max(max_lon, obs_list[[nice_name]][["obs"]]$longitude)
  }
  if (is.na(min_lat)) {
    min_lat <- min(obs_list[[nice_name]][["obs"]]$latitude)
  } else {
    min_lat <- min(min_lat, obs_list[[nice_name]][["obs"]]$latitude)
  }
  if (is.na(max_lat)) {
    max_lat <- max(obs_list[[nice_name]][["obs"]]$latitude)
  } else {
    max_lat <- max(max_lat, obs_list[[nice_name]][["obs"]]$latitude)
  }
}

# Download climate data from Swallowtail project if necessary
# "https://github.com/Big-Biodiversity-Collaborative/SwallowtailClimateChange/blob/main/data/wc2-1/bio1.tif"
base_url <- "https://github.com/Big-Biodiversity-Collaborative/SwallowtailClimateChange/raw/main/data/wc2-1/"
bio_vars <- paste0("bio", 1:19)
for (bio in bio_vars) {
  dest_file <- paste0("data/climate/", bio, ".tif")
  if (!dir.exists("data/climate")) {
    dir.create("data/climate")
  }
  if (!file.exists(dest_file)) {
    url = paste0(base_url, bio, ".tif")
    download.file(url = url,
                  destfile = dest_file)
  } else {
    message(paste0(bio, " climate data already on disk"))
  }
}

# Create geographic extent of *ALL* observations
# bio1 <- terra::rast(x = "data/climate/bio1.tif")
# geo_ext <- terra::ext(bio1)
geo_ext <- terra::ext(c(min_lon, max_lon, min_lat, max_lat))

# Create background points that will be used for SDM evaluation; use one of the 
# climate files for resolution
mask <- terra::rast(x = "data/climate/bio1.tif")

# Use random sampling to generate pseudo-absence points
# Extend background points a little beyond the 
# geographic extent of observations
absence_points <- terra::spatSample(x = mask,
                                    size = num_background,
                                    method = "random",
                                    ext = geo_ext * 1.25,
                                    as.points = TRUE,
                                    na.rm = TRUE)
# plot(absence_points)
absence_extent <- terra::ext(absence_points)

# Read in climate data into single stack
bio_files <- paste0("data/climate/", bio_vars, ".tif")
predictors <- terra::rast(x = bio_files)

# Predictors won't be needed beyond the extent of the (pseudo)absence points
model_predictors <- terra::crop(x = predictors, y = absence_extent)

# We can use same absence object for both hosts
predictors_absence <- terra::extract(x = model_predictors, y = absence_points)

# Will want to create one model for each of the host plants
host_models <- list()
for (host_nice_name in nice_names[-1]) {
  presence_points <- obs_list[[host_nice_name]][["obs"]][, c("longitude", "latitude")]
  predictors_presence <- terra::extract(x = model_predictors, y = presence_points)

  # TODO: Anything with absence should be taken out of this for loop
  
  # Make a vector of appropriate length with 0/1 values for 
  # (pseudo)absence/presence
  pa_data <- c(rep(x = 1, times = nrow(presence_points)), 
               rep(x = 0, times = nrow(absence_points)))

  # Create a vector of folds for easier splitting into testing/training
  num_folds <- 5 # for 20/80 split
  fold <- c(rep(x = 1:num_folds, length.out = nrow(presence_points)),
            rep(x = 1:num_folds, length.out = nrow(absence_points)))

  # Combine our presence / absence and fold vectors with environmental data we 
  # extracted
  full_data <- data.frame(cbind(pa = pa_data,
                                fold = fold,
                                rbind(predictors_presence, predictors_absence)))

  # Before doing test/train split, drop any presence rows with missing (NA) 
  # climate data (will also drop any absence points missing climate data, 
  # although there should not be any of those)
  full_data <- na.omit(full_data)
  
  # Create separate data frames for testing and training presence data
  presence_train <- full_data %>%
    filter(pa == 1) %>%
    filter(fold != 1)
  presence_test <- full_data %>%
    filter(pa == 1) %>%
    filter(fold == 1)
  # Create separate data frames for testing and training (pseudo)absence data
  absence_train <- full_data %>%
    filter(pa == 0) %>%
    filter(fold != 1)
  absence_test <- full_data %>%
    filter(pa == 0) %>%
    filter(fold == 1)

  # Add presence and pseudoabsence training data into single data frame
  sdmtrain <- rbind(presence_train, absence_train)
  sdmtest <- rbind(presence_test, absence_test)

  message(paste0("Running generalized linear model: ", host_nice_name))

  # Run an GLM model, specifying model with standard formula syntax
  # Exclude bio3 (a function of bio2 & bio7) and bio7 (a function of bio5 and 
  # bio6)
  # May throw a warning, but not sure if this is so bad  
  glm_model <- stats::glm(pa ~ bio1 + bio2 + bio4 + bio5 + bio6 +
                            bio8 + bio9 + bio10 + bio11 + bio12 +
                            bio13 + bio14 + bio15 + bio16 + bio17 + bio18 +
                            bio19,
                          data = sdmtrain,
                          family = binomial(link = "logit"))
  
  message(paste0("Model complete. Calculating predicted presence probabilities: ",
                 host_nice_name))

  # We want now to use that model to create a raster of presence probabilities
  presence_probs <- predict(model_predictors,
                            glm_model,
                            type = "response")
  
  host_models[[host_nice_name]] <- presence_probs
}

# Now want to take those predicted presence probabilities and do an SDM model
# for the insect, where the only predictors are those host probabilities

# Start by creating a raster stack of the host predicted probabilities
host_predictors <- terra::rast(host_models)

presence_points <- obs_list[[1]][["obs"]][, c("longitude", "latitude")]
predictors_presence <- terra::extract(x = host_predictors, y = presence_points)

# Need to grab predicted probabilities for our absence points, too
predictors_absence <- terra::extract(x = host_predictors, y = absence_points)

# Make a vector of appropriate length with 0/1 values for 
# (pseudo)absence/presence
pa_data <- c(rep(x = 1, times = nrow(presence_points)), 
             rep(x = 0, times = nrow(absence_points)))

# Create a vector of folds for easier splitting into testing/training
num_folds <- 5 # for 20/80 split
fold <- c(rep(x = 1:num_folds, length.out = nrow(presence_points)),
          rep(x = 1:num_folds, length.out = nrow(absence_points)))

# Combine our presence / absence and fold vectors with environmental data we 
# extracted
full_data <- data.frame(cbind(pa = pa_data,
                              fold = fold,
                              rbind(predictors_presence, predictors_absence)))

# Before doing test/train split, drop any presence rows with missing (NA) 
# climate data (will also drop any absence points missing climate data, 
# although there should not be any of those)
full_data <- na.omit(full_data)

# Create separate data frames for testing and training presence data
presence_train <- full_data %>%
  filter(pa == 1) %>%
  filter(fold != 1)
presence_test <- full_data %>%
  filter(pa == 1) %>%
  filter(fold == 1)
# Create separate data frames for testing and training (pseudo)absence data
absence_train <- full_data %>%
  filter(pa == 0) %>%
  filter(fold != 1)
absence_test <- full_data %>%
  filter(pa == 0) %>%
  filter(fold == 1)

# Add presence and pseudoabsence training data into single data frame
sdmtrain <- rbind(presence_train, absence_train)
sdmtest <- rbind(presence_test, absence_test)

message(paste0("Running generalized linear model: ", names(obs_list)[1]))

# Run an GLM model, specifying model with standard formula syntax
# Including host predicted probabilities as the predictors
pred_string <- paste(names(host_models), collapse = " + ")
formula_string <- paste0("pa ~ ", pred_string)

glm_model <- stats::glm(formula = eval(expr = formula_string),
                        data = sdmtrain,
                        family = binomial(link = "logit"))
summary(glm_model)

# A quick and dirty LRT
# Make a new column that is the sum of the probabilities for all the hosts
sdmtrain$host_sum <- rowSums(sdmtrain[, names(obs_list)[-1]])

# Run the glm with only that host_sum as predictor
glm_simple_model <- stats::glm(formula = pa ~ host_sum,
                               data = sdmtrain,
                               family = binomial(link = "logit"))

lnL_simple <- as.numeric(stats::logLik(glm_simple_model))
lnL_full <- as.numeric(stats::logLik(glm_model))
delta_lnL <- lnL_full - lnL_simple
# Calculate p-value assuming chi-square distribution and 1 df
p_value <- pchisq(q = delta_lnL, df = 1, lower.tail = FALSE)

########################################
# OLD BELOW

# Run SDM GLM on two host plants
# Run SDM GLM on E. chalcedona ~ climate data only
# Maybe do some lassoing here?
# Run SDM GLM on E. chalcedona ~ climate + hosts
# See if model with hosts is better? Why?
# Run SDM GLM on E. chalcedona ~ hosts

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
