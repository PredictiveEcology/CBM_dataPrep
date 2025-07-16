
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

  if (inherits(input, "sf")){

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
    # input <- withCallingHandlers(
    #   postProcess(
    #     input,
    #     cropTo    = masterRaster,
    #     projectTo = masterRaster
    #   ),
    #   warning = function(w){
    #     if (w$message == "attribute variables are assumed to be spatially constant throughout all geometries"){
    #       invokeRestart("muffleWarning")
    #     }
    #   })

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

  }else{

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
}


