
#' Prep inputs: Extract values from spatial data source
#'
#' Reproject and align input to a `masterRaster` template
#' while retaining categorical values by using "mode" resampling.
#' Return a vector of the values for each cell of `masterRaster`.
#'
#' @param input terra SpatRaster, one or more raster files, or sf polygons.
#' @param masterRaster terra SpatRaster. Raster template for alignment
#' @param outPath character. Path for writing aligned raster with overwrite = TRUE.
#'
#' @return vector
prepInputsExtractSpatial <- function(input, masterRaster, outPath = NULL, verbose = TRUE){

  # If saving output: make sure path can be written to
  if (!is.null(outPath) && file.exists(outPath)){
    unlink(outPath)
    if (file.exists(outPath)) stop("Failed to remove file: ", outPath)
  }

  # Align input with masterRaster
  if (verbose) message("Aligning input to masterRaster")
  alignRast <- prepInputsToMasterRaster(input, masterRaster)

  # Extract raster values
  if (verbose) message("Extracting values from aligned raster")
  alignVals <- terra::values(alignRast, mat = FALSE)

  # Set values as text
  alignCats <- terra::cats(alignRast)[[1]]
  if (!is.null(alignCats)) alignVals <- alignCats[[2]][alignVals]

  # Write output to file
  if (!is.null(outPath)){

    if (verbose) message("Writing aligned raster to path: ", outPath)

    dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(alignRast, outPath, overwrite = TRUE)
  }

  return(alignVals)
}


#' Prep inputs to masterRaster
#'
#' Reproject and align input to a `masterRaster` template
#' while retaining categorical values by using "mode" resampling.
#'
#' @param input terra SpatRaster, one or more raster files, or sf polygons.
#' @param masterRaster terra SpatRaster. Raster template for alignment
#'
#' @return SpatRaster
prepInputsToMasterRaster <- function(input, masterRaster, useCache = FALSE){

  if (inherits(input, "sf")){

    prepInputsToMasterRaster_vect(input, masterRaster) |>
      reproducible::Cache(
        useCache = useCache,
        verbose  = ifelse(useCache, getOption("reproducible.verbose", 1), -1))

  }else{

    prepInputsToMasterRaster_rast(input, masterRaster) |>
      reproducible::Cache(
        useCache = useCache,
        verbose  = ifelse(useCache, getOption("reproducible.verbose", 1), -1))
  }
}

# Prep inputs to masterRaster: raster
prepInputsToMasterRaster_rast <- function(input, masterRaster){

  # Read as SpatRaster; mosaic tiles if need be
  if (!inherits(input, "SpatRaster")){

    if (length(input) > 1 && is.character(input) &&
        all(tryCatch(file.exists(input), error = function(e) FALSE))){
      input <- do.call(terra::mosaic, lapply(input, terra::rast))

    }else{
      input <- terra::rast(input)
    }
  }

  # Keep only the first band
  if (terra::nlyr(input) > 1) input <- input[[1]]

  # Crop input
  input <- terra::crop(
    input,
    terra::project(terra::as.polygons(masterRaster, extent = TRUE), terra::crs(input)),
    snap = "out")

  # Assign NAs to -1 so they are considered values in mode calculation
  ## TODO: find a better method
  input <- terra::classify(input, cbind(NA, -1))

  reproj <- !terra::compareGeom(
    input, masterRaster, lyrs = FALSE,
    crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE,
    warncrs = FALSE, stopOnError = FALSE, messages = FALSE)

  if (reproj){

    input <- terra::project(input, masterRaster, method = "mode")

  }else{

    terra::crs(masterRaster) <- terra::crs(input)
    input <- exactextractr::exact_resample(input, masterRaster, fun = "mode")
  }

  terra::classify(input, cbind(-1, NA))

}

# Prep inputs to masterRaster: vector
prepInputsToMasterRaster_vect <- function(input, masterRaster){

  # Crop and reproject
  cropBBOX <- sf::st_buffer(
    sf::st_as_sfc(sf::st_transform(sf::st_bbox(masterRaster), sf::st_crs(input))),
    terra::res(masterRaster)[[1]])

  input <- withCallingHandlers(
    sf::st_crop(input, sf::st_bbox(cropBBOX)),
    warning = function(w){
      if (w$message == "attribute variables are assumed to be spatially constant throughout all geometries"){
        invokeRestart("muffleWarning")
      }
    })
  input <- sf::st_transform(input, sf::st_crs(masterRaster))

  ## Observed to be slow in some cases
  # input <- postProcess(
  #     input,
  #     cropTo    = masterRaster,
  #     projectTo = masterRaster
  #   )

  # Rasterize
  cellIdxRast <- exactextractr::rasterize_polygons(
    input, masterRaster, min_coverage = 0.5)

  rclTable <- data.frame(from = 1:nrow(input), to = input[[1]])
  if (is.character(rclTable$to)){

    levels(cellIdxRast) <- rclTable
    cellIdxRast

  }else{
    terra::classify(cellIdxRast, rcl = rclTable)
  }
}


