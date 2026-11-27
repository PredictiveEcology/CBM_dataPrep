
studyAreas <- function(dataPath, studyAreaAdmin = NULL, studyAreaEco = NULL, crs = NULL, buffer = NULL){

  adminLocator <- reproducible::prepInputs(
    destinationPath = dataPath,
    archive         = file.path(dataPath, "adminLocator.zip"),
    targetFile      = "adminLocator.shp",
    fun             = sf::st_read(targetFile, quiet = TRUE)
  )
  names(adminLocator)[1] <- "admin_name"
  adminLocator$admin_abbrev <- as.character(adminAbbrev(adminLocator$admin_name))

  if (!is.null(studyAreaAdmin)){

    if (!is.vector(studyAreaAdmin)) stop("studyAreaAdmin must be a vector")

    adminLocator <- adminLocator[
      tolower(adminLocator$admin_name) %in% tolower(as.character(studyAreaAdmin)) |
        tolower(adminLocator$admin_abbrev) %in% tolower(as.character(studyAreaAdmin)),]

    if (nrow(adminLocator) == 0) stop(
      "studyAreaAdmin invalid. Choose from: ", paste(sQuote(unique(adminLocator$admin_name)), collapse = ", "))
  }

  ecoLocator <- reproducible::prepInputs(
    destinationPath = dataPath,
    archive         = file.path(dataPath, "ecoLocator.zip"),
    targetFile      = "ecoLocator.shp",
    fun             = sf::st_read(targetFile, quiet = TRUE)
  )
  names(ecoLocator)[1:2] <- c("eco_id", "eco_name")

  if (!is.null(studyAreaEco)){

    if (!is.vector(studyAreaEco)) stop("studyAreaEco must be a vector")

    ecoLocator <- ecoLocator[
      ecoLocator$eco_id %in% studyAreaEco |
        tolower(ecoLocator$eco_name) %in% tolower(as.character(studyAreaEco)),]

    if (nrow(ecoLocator) == 0) stop(
      "studyAreaEco invalid. Choose from: ", paste(sQuote(unique(ecoLocator$eco_name)), collapse = ", "))
  }

  if (is.null(crs)) crs <- sf::st_crs(adminLocator)

  studyAreas <- suppressWarnings(sf::st_intersection(
    sf::st_transform(adminLocator, crs), sf::st_transform(ecoLocator, crs)
  ))

  if (nrow(studyAreas) == 0) stop("Study area studyAreaAdmin and ecozone(s) boundaries do not overlap")

  if (!is.null(buffer)){

    studyAreasBuffer <- sf::st_buffer(studyAreas, buffer)
    studyAreasBuffer <- suppressWarnings(sf::st_collection_extract(sf::st_difference(studyAreasBuffer), "POLYGON"))
    studyAreasBuffer <- studyAreasBuffer[sf::st_geometry_type(studyAreasBuffer) %in% c("POLYGON", "MULTIPOLYGON"),]
    studyAreasBuffer <- suppressWarnings(sf::st_difference(studyAreasBuffer, sf::st_union(sf::st_geometry(studyAreas))))

    studyAreas <- rbind(studyAreas, studyAreasBuffer)
  }

  studyAreas <- aggregate(
    studyAreas[attr(studyAreas, "sf_column")],
    by = as.list(sf::st_drop_geometry(studyAreas)),
    FUN = unique
  )
  studyAreas <- studyAreas[order(studyAreas$admin_name, studyAreas$eco_id), ]
  row.names(studyAreas) <- NULL

  studyAreas

}

