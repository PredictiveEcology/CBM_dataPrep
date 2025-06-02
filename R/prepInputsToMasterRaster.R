
#' Prep inputs to masterRaster
#'
#' Reproject and align input to a `masterRaster` template
#' while retaining categorical values by using "mode" resampling.
#'
#' @param input terra SpatRaster, one or more raster files, or sf polygons.
#' @param masterRaster terra SpatRaster. Raster template for alignment
#'
#' @return SpatRaster
prepInputsToMasterRaster <- function(input, masterRaster){

  # If input is file(s): read as SpatRaster; mosaic tiles if need be
  if (is.character(input) &&
      all(tryCatch(file.exists(input), error = function(e) FALSE))){
    input <- do.call(terra::mosaic, lapply(input, terra::rast))
  }

  if (inherits(input, "sf")){

    # Crop and reproject
    input <- postProcess(
      input,
      cropTo    = masterRaster,
      projectTo = masterRaster
    )

    # Rasterize
    cellIdxRast <- exactextractr::rasterize_polygons(
      input, masterRaster, min_coverage = 0.5)

    terra::classify(
      cellIdxRast,
      rcl = data.frame(from = 1:nrow(input), to = input[[1]]))

  }else{

    # Read as SpatRaster
    if (!inherits(input, "SpatRaster")){
      input <- tryCatch(
        terra::rast(input),
        error = function(e) stop(
          "input could not be converted to SpatRaster: ", e$message,
          call. = FALSE))
    }

    reproj <- !terra::compareGeom(
      input, masterRaster, lyrs = FALSE,
      crs = TRUE, ext = FALSE, rowcol = FALSE, res = FALSE,
      warncrs = FALSE, stopOnError = FALSE, messages = FALSE)

    # Assign NAs to -1 so they are considered values in mode calculation
    input <- terra::classify(input, cbind(NA, -1))

    if (reproj){

      input <- postProcess(
        input,
        to     = masterRaster,
        method = "mode"
      )

    }else{

      terra::crs(masterRaster) <- terra::crs(input)
      input <- exactextractr::exact_resample(input, masterRaster, fun = "mode")
    }

    terra::classify(input, cbind(-1, NA))
  }
}


