#' Create presence / absence data from buffered minimum convex polygon
#' 
#' @param observations data.frame or matrix of observations in x, y coordinates
#' @param tif_file path to climate data file that provides sampling resolution 
#' for pseudo-absence data
#' @param num_absence integer number of pseudo-absence points to sample; the 
#' actual number of pseudo-absence points may be lower than this number
#' @param buffer_mult factor by which to multiply nearest-neighbor distance for
#' value of buffer
#' @param num_folds integer number of folds for splitting data into training / 
#' testing sets
#' @param coord_cols character vector indicating column names of x and y 
#' coordinates, respectively; if data in \code{observations} has longitude and 
#' latitude columns, pass \code{coord_cols = c("longitude", "latitude")}
#' 
#' @details It is assumed the data passed through \code{observations} has 
#' already been de-duplicated
#' 
#' @return list with two elements:
#' \describe{
#'   \item{data}{data.frame of presence/absence data with columns for: presence 
#'   / absence, fold, longitude, latitude}
#'   \item{buffered_mcp}{buffered minimum convex polygon used to generate 
#'   pseudo-absence points}
#' }
pa_mcp <- function(observations, num_absence = 5000, buffer_mult = 1000, 
                   num_folds = 5, coord_cols = c("x", "y")) {
  # TODO: Most of this is taken from E. Zylstra's approach at 
  # https://github.com/Big-Biodiversity-Collaborative/SwallowtailClimateChange/blob/main/src/data/gbif-3-presence-absence.R
  # Load in all default or passed values
  num_absence <- match.args(num_absence)
  buffer_mult <- match.args(buffer_mult)
  num_folds <- match.args(num_folds)
  coord_cols <- match.args(coord_cols)
  
  # Extract function name messages
  function_name <- as.character(match.call())[1]
  
  # Check incoming data to make sure we have the correct columns
  if (!all(coord_cols %in% colnames(observations))) {
    stop(function_name, " requires x, y columns.")
  }
  
  # Load dependencies
  dependencies <- c("terra", "sp")
  if (!all(unlist(lapply(X = dependencies, FUN = require, character.only = TRUE)))) {
    stop("At least one package required by ", function_name, 
         " could not be loaded: ", paste(dependencies, collapse = ", "),
         " are required.")
  }

  # Load tif file used for sampling resolution
  if (!file.exists(tif_file)) {
    stop(function_name, " could not locate tif file ", tif_file)
  }
  
  # Extract just those x, y columns for our use
  presence <- observations[, coord_cols]

  # Load the tif file for sampling and MCP
  predictor <- terra::rast(tif_file)
  
  # Sample one observation from each of the climate data grid cells
  thinned <- presence %>%
    dismo::gridSample(r = predictor, n = 1)
  
  # Convert to a SpatialPoints object using the WGS84 CRS
  thinned_sp <- SpatialPoints(coords = thinned,
                                proj4string = CRS("+init=epsg:4326"))
  
  # Calculate the GreatCircle distance (in km) between points (can take
  # minutes for larger data sets)
  gc_dist <- sp::spDists(thinned_sp, longlat = TRUE) 
  
  # Calculate the maximum of nearest neighbor distances (in km)
  buffer <- gc_dist %>%
    apply(., 1, function(x) min(x[x > 0])) %>%
    max %>%
    round
  rm(gc_dist)
  
  # Garbage can pile up at this point. Clean it up.
  gc(verbose = FALSE)
  
  # Create a minimum convex polygon (MCP) for observations
  ch <- chull(presence)
  ch_coords <- presence[c(ch, ch[1]), ]
  ch_polygon <- SpatialPolygons(list(Polygons(list(Polygon(ch_coords)), ID = 1)),
                                proj4string = CRS("+init=epsg:4326"))
  
  
  
  
  
  
  
  
  



  # Create a minimum convex polygon (MCP) for observations
  ch <- chull(presence)
  ch_coords <- presence[c(ch, ch[1]), ]
  ch_polygon <- SpatialPolygons(list(Polygons(list(Polygon(ch_coords)), ID = 1)),
                                proj4string = CRS("+init=epsg:4326"))
  
  # Convert MCP to sf object and project to NA Albers Equal Area Conic 
  ch_poly_proj <- st_as_sf(ch_polygon) %>%
    st_transform(crs = "ESRI:102008")
  
  # Create a polygon = MCP + buffer
  ch_buffer <- st_buffer(ch_poly_proj,
                         dist = buffer * 1000)
  
  # Transform the buffered MCP back to lat/long 
  ch_buffer_latlong <- st_transform(ch_buffer, 4326)
  
  # Convert the buffered MCP to a SpatVector
  ch_buffer_sv <- terra::vect(ch_buffer_latlong)
  
  # Crop and mask climate data to the buffered MCP polygon
  pred_mask <- predictor %>%
    terra::crop(ch_buffer_sv) %>%
    terra::mask(ch_buffer_sv)
  rm(predictor)

  # TODO: Find number of cells, in case it is less than num_absence
  num_cells <- c()
  
  # Generate pseudo-absence points
  absence <- terra::spatSample(x = pred_mask,  
                               size = min(num_absence, num_cells),
                               method = "random",
                               na.rm = TRUE,
                               values = FALSE,
                               xy = TRUE)
  # Note: if buffered area is small, function may not be able to generate the 
  # indicated number of points (size) because the default is to select cells
  # without replacement. If this happens, R will return a warning:
  # [spatSample] fewer cells returned than requested
  
  # Reality check:
  # plot(ch_buffer_sv)
  # plot(ch_polygon, add = TRUE)
  # points(y ~ x, data = absence, cex = 0.5, col = "gray")
  # points(y ~ x, data = presence, cex = 0.5, col = "blue")
  
  # Make a vector of appropriate length with 0/1 values for 
  # (pseudo)absence/presence
  pa_data <- c(rep(x = 1, times = nrow(presence)), 
               rep(x = 0, times = nrow(absence)))  
  
  # Create a vector of folds for easier splitting into testing/training
  fold <- c(rep(x = 1:num_folds, length.out = nrow(presence)),
            rep(x = 1:num_folds, length.out = nrow(absence)))
  
  # Combine our presence / absence data
  full_data <- data.frame(cbind(pa = pa_data,
                                fold = fold,
                                rbind(presence, absence)))
  
  return(list(data = full_data,
              buffered_mcp = ch_buffer_latlong))
}