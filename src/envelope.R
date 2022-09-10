#' Return two-dimensional data outside the contour of kernel density estimate
#' 
#' @param data data frame or matrix with longitude and latitude columns
#' @param climate_data raster used to thin data before defining envelope
#' @param coord_cols vector indicating names of x and y columns
#' @param cutoff scalar value between 0 and 1; e.g. 0.95 discards observations 
#' outside the 95% density envelope
#' 
#' @return integer vector indexed by original data frame where 0 indicates data 
#' outside the contour and 1 indicates data inside the contour
envelope <- function(data, climate_data, cutoff = 0.95) {
  if (!require(dismo)) {
    stop("envelope requires dismo package, but it could not be loaded")
  }
  if (!require(dplyr)) {
    stop("envelope requires dplyr package, but it could not be loaded")
  }
  # Thin observations, to avoid lots of observations within individual cells 
  # driving the envelope
  # terra::rast(x = data[, c("longitude", "latitude")], type = "xy")
  thinned <- data %>%
    dplyr::select(longitude, latitude) %>%
    dismo::gridSample(r = climate_data, n = 1)

  # Calculate density envelope; using a 0.5 degree resolution
  n_points <- c(abs(max(thinned$longitude) - min(thinned$longitude)) * 2,
                abs(max(thinned$latitude) - min(thinned$latitude)) * 2)
  obs_kde <- MASS::kde2d(x = thinned$longitude,
                         y = thinned$latitude, 
                         n = n_points)
  
  # The projection string for raster conversion
  wgs_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # Transform this to a raster; not quite sure why, but need to do 90 degree 
  # counter-clockwise rotation of the z matrix...
  kde_raster <-raster::raster(x = apply(X = t(obs_kde$z),
                                        MARGIN = 2,
                                        FUN = rev),
                              xmn = min(obs_kde$x), 
                              xmx = max(obs_kde$x),
                              ymn = min(obs_kde$y), 
                              ymx = max(obs_kde$y), 
                              crs = wgs_crs)
  # From https://mhallwor.github.io/_pages/activities_GenerateTerritories
  # Set zeros to NA
  kde_raster[kde_raster == 0] <- NA
  # Get the values as a vector
  kde_values <- raster::getValues(kde_raster)
  # Sort all the not missing values
  sorted_values <- sort(kde_values[!is.na(kde_values)], 
                        decreasing = TRUE)
  # Create cumulative sum of those sorted values
  summed_values <- cumsum(x = sorted_values)
  # Find index of those sorted values for the cutoff
  cutoff_index <- sum(summed_values <= cutoff * summed_values[length(summed_values)])
  # Set the values of the raster to 0 or 1 based on that cutoff
  kde_envelope <- raster::setValues(kde_raster, 
                                    kde_values >= sorted_values[cutoff_index])
  
  # Pull out values from the kde_envelope raster
  # 0 = outside envelope, 1 = inside envelope
  in_out <- raster::extract(x = kde_envelope, 
                                   y = data[,c("longitude", "latitude")])
  
  # Sometimes values for envelope are NA if point is on edge of the defined 
  # envelope; count those as outliers, too
  in_out[is.na(in_out)] <- 0

  # Reality check
  # plot(kde_envelope, col = c("white", "grey75"))
  # points(x = data$longitude, y = data$latitude, pch = 16,
  #        cex = 0.75, col = "blue")
  # points(x = data$longitude[in_out == 0],
  #        y = data$latitude[in_out == 0],
  #        cex = 1, col = "red", pch = 16)
  return(in_out)
}