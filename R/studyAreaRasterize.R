
studyAreaRasterize <- function(studyArea, template, res = NULL){

  if (inherits(studyArea, "bbox")) studyArea <- sf::st_as_sfc(studyArea)

  studyArea <- sf::st_transform(sf::st_geometry(studyArea), sf::st_crs(template))
  studyArea <- terra::vect(studyArea)
  terra::crs(studyArea) <- terra::crs(template)

  masterRaster <- terra::crop(terra::extend(terra::rast(template), studyArea), studyArea)

  if (!is.null(res)){
    resTempl <- terra::res(template)[[1]]
    if (res > resTempl){
      fact <- res / resTempl
      if (fact %% 1 != 0) stop(
        "masterRaster resolution of ", res, " cannot produce a raster aligned with the ", resTempl, " resolution template")
      masterRaster <- terra::aggregate(masterRaster, fact)
    }else{
      fact <- resTempl / res
      if (fact %% 1 != 0) stop(
        "masterRaster resolution of ", res, " cannot produce a raster aligned with the ", resTempl, " resolution template")
      masterRaster <- terra::disagg(masterRaster, fact)
    }
  }

  masterRaster <- terra::rasterize(studyArea, masterRaster, touches = FALSE)

  return(masterRaster)
}

studyAreaRasterTemplate <- function(dataSource){

  mrTemplates <- list(
    SCANFI = terra::rast(
      xmin = -2341500, ymin = 5863500,
      xmax =  3010500, ymax = 9436500,
      res = 30,
      crs = paste(
        "PROJCRS[\"NAD_1983_Canada_Lambert\",",
        "    BASEGEOGCRS[\"NAD83\",",
        "        DATUM[\"North American Datum 1983\",ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,LENGTHUNIT[\"metre\",1]]],",
        "        PRIMEM[\"Greenwich\",0,ANGLEUNIT[\"degree\",0.0174532925199433]],",
        "        ID[\"EPSG\",4269]],",
        "    CONVERSION[\"Lambert Conic Conformal (2SP)\",",
        "        METHOD[\"Lambert Conic Conformal (2SP)\",ID[\"EPSG\",9802]],",
        "        PARAMETER[\"Latitude of false origin\",0,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8821]],",
        "        PARAMETER[\"Longitude of false origin\",-95,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8822]],",
        "        PARAMETER[\"Latitude of 1st standard parallel\",49,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8823]],",
        "        PARAMETER[\"Latitude of 2nd standard parallel\",77,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8824]],",
        "        PARAMETER[\"Easting at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8826]],",
        "        PARAMETER[\"Northing at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8827]]],",
        "    CS[Cartesian,2],",
        "        AXIS[\"easting\",east,ORDER[1],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]],",
        "        AXIS[\"northing\",north,ORDER[2],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]]]"
      ))
  )

  if (!dataSource %in% names(mrTemplates)) stop(
    "No master raster template found for ", sQuote(dataSource),
    ". Choose from: ", paste(sQuote(names(mrTemplates)), collapse = ", "))

  return(mrTemplates[[dataSource]])
}

