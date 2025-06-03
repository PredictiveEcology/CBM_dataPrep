
## DEFINE MODULE ----

defineModule(sim, list(
  name = "CBM_dataPrep",
  description = "A data preparation module to format and prepare user-provided input to the SpaDES forest-carbon modelling family.",
  keywords = NA,
  authors = c(
    person("Céline", "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Susan",  "Murray",    email = "murray.e.susan@gmail.com",           role = c("ctb"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.2", CBM_dataPrep = "1.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("CBM_dataPrep.Rmd"),
  reqdPkgs = list(
    "data.table", "sf", "terra", "exactextractr",
    "reproducible (>=2.1.2)" ,
    "PredictiveEcology/CBMutils@development (>=2.0.3.0005)"
  ),
  parameters = rbind(
    defineParameter("saveRasters", "logical", FALSE, NA, NA, "Save rasters of inputs aligned to the `masterRaster`"),
    defineParameter(".useCache", "character", c(".inputObjects", "Init"), NA, NA, "Cache module events")
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "masterRaster", objectClass = "SpatRaster",
      desc = "Raster template defining the study area. NA cells will be excluded from analysis."),
    expectsInput(
      objectName = "masterRasterURL", objectClass = "character", desc = "URL for `masterRaster`"),
    expectsInput(
      objectName = "ecoLocator", objectClass = "sf|SpatRaster",
      sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
      desc = paste(
        "Spatial data source of Canada ecozone IDs.",
        "Default is Canada's ecozones as polygon features.")),
    expectsInput(
      objectName = "ecoLocatorURL", objectClass = "character", desc = "URL for `ecoLocator`"),
    expectsInput(
      objectName = "spuLocator", objectClass = "sf|SpatRaster",
      sourceURL = "https://drive.google.com/file/d/1D3O0Uj-s_QEgMW7_X-NhVsEZdJ29FBed",
      desc = paste(
        "Spatial data source of CBM-CFS3 spatial unit IDs.",
        "Default is Canada's spatial units as polygon features.")),
    expectsInput(
      objectName = "spuLocatorURL", objectClass = "character", desc = "URL for `spuLocator`"),
    expectsInput(
      objectName = "ageLocator", objectClass = "sf|SpatRaster",
      desc = "Spatial data source of stand ages."),
    expectsInput(
      objectName = "ageLocatorURL", objectClass = "character", desc = "URL for `ageLocator`"),
    expectsInput(
      objectName = "ageDataYear", objectClass = "numeric",
      desc = paste(
        "Year that the ages in `ageLocator` represent.",
        "If omitted, ages are assumed to represent the simulation start year.")),
    expectsInput(
      objectName = "ageSpinupMin", objectClass = "numeric",
      desc = paste(
        "Minimum age for cohorts during spinup.",
        "Temporary fix to CBM_core issue: https://github.com/PredictiveEcology/CBM_core/issues/1")),
    expectsInput(
      objectName = "cohortLocators", objectClass = "list",
      desc = paste(
        "List of spatial data sources for additional columns in `cohortDT`.",
        "Each item may be a terra SpatRaster object,",
        "a character vector of multiple raster tiles,",
        "or an sf polygons object.")),
    expectsInput(
      objectName = "cohortLocatorURLs", objectClass = "list", desc = "URLs for `cohortLocators`"),
    expectsInput(
      objectName = "gcKey", objectClass = "character",
      desc = paste(
        "Column names defining unique growth curves in `gcMeta` nad `userGcM3`.",
        "Each column must have a corresponding named spatial data source in `cohortLocators`")),
    expectsInput(
      objectName = "gcMeta", objectClass = "data.frame",
      desc = "Growth curve metadata with all columns in `gcKey`",
      columns = list(
        species       = paste(
          "Optional. Species name.",
          "If no other columns are provided, names will be matched with known species to set them."),
        species_id    = "CBM-CFS3 species ID",
        sw_hw         = "'sw' or 'hw'",
        canfi_species = "CanFI species codes",
        genus         = "NFI species genus"
      )),
    expectsInput(
      objectName = "gcMetaURL", objectClass = "character", desc = "URL for `gcMeta`"),
    expectsInput(
      objectName = "userGcM3", objectClass = "data.frame",
      desc = "Growth curve volumes with all columns in `gcKey` and `Age` and `MerchVolume`."),
    expectsInput(
      objectName = "userGcM3URL", objectClass = "character", desc = "URL for `userGcM3`"),
    expectsInput(
      objectName = "disturbanceRasters", objectClass = "list",
      desc = paste(
        "Set of spatial data sources containing locations of disturbance events for each year.",
        "List items must be named by disturbance event IDs found in `disturbanceMeta`.",
        "Within each event's list, items must be named by the 4 digit year the disturbances occured in.",
        "For example, event type 1 disturbance locations for 2025 can be accessed with `disturbanceRasters[[\"1\"]][[\"2025\"]]`.",
        "Each disturbance item can be one of the following:",
        "a terra SpatRaster layer, one or more raster file paths, or sf polygons.",
        "All non-NA areas will be considered events."
      )),
    expectsInput(
      objectName = "disturbanceMeta", objectClass = "data.table",
      desc = "Table defining the disturbance event types",
      columns = c(
        eventID             = "Event type ID",
        disturbance_type_id = "Optional. CBM-CFS3 disturbance type ID. If not provided, the user will be prompted to choose IDs.",
        name                = "Optional. Disturbance name (e.g. 'Wildfire'). Required if 'disturbance_type_id' absent.",
        objectName          = "Optional. Name of the object in the `simList` to retrieve the disturbances from annually.",
        delay               = "Optional. Delay (in years) of when the disturbance will take effect"
      )),
    expectsInput(
      objectName = "disturbanceMetaURL", objectClass = "character", desc = "URL for `disturbanceMeta`"),
    expectsInput(
      objectName = "dbPath", objectClass = "character",
      sourceURL = "https://raw.githubusercontent.com/cat-cfs/libcbm_py/main/libcbm/resources/cbm_defaults_db/cbm_defaults_v1.2.8340.362.db",
      desc = paste(
        "Path to the CBM-CBM3 defaults database. ",
        "Required if `disturbanceMeta` is missing the 'disturbance_type_id' column"))
  ),
  outputObjects = bindrows(
    createsOutput(
      objectName = "standDT", objectClass = "data.table",
      desc = "Table of stand attributes.",
      columns = c(
        pixelIndex      = "`masterRaster` cell index",
        area            = "`masterRaster` cell area in meters",
        ecozone         = "Canada ecozone ID extracted from input `ecoLocator`",
        spatial_unit_id = "CBM-CFS3 spatial unit ID extracted from input `spuLocator`"
      )),
    createsOutput(
      objectName = "cohortDT", objectClass = "data.table",
      desc = "Table of cohort attributes.",
      columns = c(
        cohortID   = "`masterRaster` cell index",
        pixelIndex = "`masterRaster` cell index",
        ages       = "Cohort ages extracted from input `ageLocator`",
        ageSpinup  = "Cohort ages raised to minimum of `ageSpinupMin` to use in the spinup"
      )),
    createsOutput(
      objectName = "curveID", objectClass = "character",
      desc = paste(
        "Name of the column uniquely defining each growth curve in `gcMeta` and `userGcM3`.",
        "This column will also be present in `cohortDT`.")),
    createsOutput(
      objectName = "gcMeta", objectClass = "data.frame", # TODO
      desc = "Growth curve metadata.",
      columns = list(
        curveID         = "`curveID`",
        ecozone         = "Canada ecozone ID",
        spatial_unit_id = "CBM-CFS3 spatial unit ID",
        species_id      = "CBM-CFS3 species ID",
        sw_hw           = "'sw' or 'hw'",
        canfi_species   = "CanFI species codes",
        genus           = "NFI species genus"
      )),
    createsOutput(
      objectName = "userGcM3", objectClass = "data.frame", # TODO
      desc = "Growth curve volumes.",
      columns = list(
        curveID         = "`curveID`",
        ecozone         = "Canada ecozone ID",
        spatial_unit_id = "CBM-CFS3 spatial unit ID",
        species_id      = "CBM-CFS3 species ID",
        sw_hw           = "'sw' or 'hw'",
        canfi_species   = "CanFI species codes",
        genus           = "NFI species genus"
      )),
    createsOutput(
      objectName = "disturbanceEvents", objectClass = "data.table",
      desc = paste(
        "Table with disturbance events for each simulation year.",
        "Input `disturbanceRasters` are aligned with the `masterRaster`",
        "and the events are summarized into this table.")),
    createsOutput(
      objectName = "disturbanceMeta", objectClass = "data.frame",
      desc = "Table defining the disturbance event types.")
  )
))


## MODULE EVENTS ----

doEvent.CBM_dataPrep <- function(sim, eventTime, eventType, debug = FALSE) {

  switch(
    eventType,

    init = {

      # Prepare data
      sim <- Init(sim)

      # Read disturbances
      sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "readDisturbances", eventPriority = 8)
    },

    readDisturbances = {

      # Get disturbances for the year
      if (length(sim$disturbanceRasters) > 0){

        distRasts <- lapply(sim$disturbanceRasters, function(d){
          if (as.character(time(sim)) %in% names(d)) d[[as.character(time(sim))]]
        })

      }else distRasts <- list()

      # Retrieve disturbances from simList
      for (i in which(!is.na(sim$disturbanceMeta$objectName))){

        distRasts[[as.character(sim$disturbanceMeta[i,]$eventID)]] <- get(
          sim$disturbanceMeta[i,]$objectName, envir = sim)
      }

      # Summarize year events into a table
      distRasts <- distRasts[!sapply(distRasts, is.null)]
      if (length(distRasts) > 0){

        if (is.null(names(distRasts)))    stop("disturbanceRasters list names must be disturbance event IDs")
        if (any(is.na(names(distRasts)))) stop("disturbanceRasters list names contains NAs")
        if (any(names(distRasts) == ""))  stop("disturbanceRasters list names contains empty strings")

        eventIDs <- suppressWarnings(tryCatch(
          as.integer(names(distRasts)),
          error = function(e) stop("disturbanceRasters list names must be coercible to integer")))

        for (i in 1:length(distRasts)){

          distAlign <- prepInputsToMasterRaster(
            distRasts[[i]],
            sim$masterRaster
          ) |> Cache()

          eventIndex <- which(!is.na(terra::values(distAlign)[,1]))
          if (length(eventIndex) > 0){

            distDelay <- sim$disturbanceMeta$delay[sim$disturbanceMeta$eventID == eventIDs[[i]]]
            if (length(na.omit(distDelay)) != 1) distDelay <- 0

            sim$disturbanceEvents <- rbind(sim$disturbanceEvents, data.table::data.table(
              pixelIndex = eventIndex,
              year       = as.integer(time(sim) + distDelay),
              eventID    = eventIDs[[i]]
            ))
          }

          if (P(sim)$saveRasters){
            outPath <- file.path(
              outputPath(sim), "CBM_dataPrep",
              sprintf("distRast_%s-%s_%s.tif", i, eventIDs[[i]], time(sim)))
            dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
            terra::writeRaster(distAlign, outPath, overwrite = TRUE)
          }
        }
      }

      # Schedule for next year
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_dataPrep", "readDisturbances", eventPriority = 8)

    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

Init <- function(sim) {

  ## Read spatial inputs ----

  if (length(sim$cohortLocators) > 0){
    if (!is.list(sim$cohortLocators) || is.null(names(sim$cohortLocators))) stop(
      "'cohortLocators' must be a named list")
    if (any(is.na(names(sim$cohortLocators)))) stop("'cohortLocators' names contains NAs")
  }

  # Set which columns come from which input object
  colInputs <- c(
    list(
      ecozone         = sim$ecoLocator,
      spatial_unit_id = sim$spuLocator,
      ages            = sim$ageLocator
    ),
    sim$cohortLocators
  )
  colInputs <- colInputs[!sapply(colInputs, is.null)]

  # Read master raster
  if (is.null(sim$masterRaster)) stop("masterRaster not found")
  if (!inherits(sim$masterRaster, "SpatRaster")){
    sim$masterRaster <- tryCatch(
      terra::rast(sim$masterRaster),
      error = function(e) stop(
        "masterRaster could not be converted to SpatRaster: ", e$message,
        call. = FALSE))
  }

  # Initiate table
  allPixDT <- data.table::data.table(
    pixelIndex = 1:terra::ncell(sim$masterRaster),
    area       = terra::values(terra::cellSize(sim$masterRaster, unit = "m", mask = TRUE, transform = FALSE))[,1]
  )
  data.table::setkey(allPixDT, pixelIndex)

  # Convert to SpatRaster, align with masterRaster, and add columns to table
  for (colName in names(colInputs)){

    inAlign <- prepInputsToMasterRaster(
      colInputs[[colName]],
      sim$masterRaster
    ) |> Cache()

    if (P(sim)$saveRasters){
      outPath <- file.path(outputPath(sim), "CBM_dataPrep", paste0(colName, ".tif"))
      dir.create(dirname(outPath), recursive = TRUE, showWarnings = FALSE)
      terra::writeRaster(inAlign, outPath, overwrite = TRUE)
    }

    allPixDT[[colName]] <- terra::values(inAlign)[, 1]
  }

  # Subset table to cells where masterRaster is not NA
  allPixDT <- allPixDT[!is.na(terra::values(sim$masterRaster)[,1]),]
  if (nrow(allPixDT) == 0) stop("all masterRaster values are NA")

  # Remove pixels that are missing key attributes
  if (length(colInputs) > 0){
    allPixDT_isNA <- is.na(allPixDT)
    if (any(allPixDT_isNA)){

      rmRow <- apply(allPixDT_isNA, 1, any)
      rmCol <- apply(allPixDT_isNA, 2, any)

      if (all(rmRow)) stop(
        "All pixels invalid due to NAs in: ",
        paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

      message(
        sum(rmRow), "/", nrow(allPixDT),
        " pixel(s) removed due to NAs in: ",
        paste(shQuote(names(rmCol)[rmCol]), collapse = ", "))

      allPixDT <- allPixDT[!rmRow,]

      rm(rmRow)
      rm(rmCol)
    }
    rm(allPixDT_isNA)
  }

  # Create unique growth curve ID for every spatial_unit_id
  if (all(c("gcIndex", "spatial_unit_id") %in% names(allPixDT))){
    allPixDT$gcids <- CBMutils::gcidsCreate(allPixDT[, .(spatial_unit_id, gcIndex)])
  }

  # Adjust stand ages
  if ("ages" %in% names(allPixDT) && !is.null(sim$ageDataYear) && sim$ageDataYear != start(sim)){

    # TODO: add step to adjust stand age to the simulation start year
    warning("Cohort age data is from ", sim$ageDataYear,
            " instead of the simulation start year")
  }

  # Set spinup age
  if ("ages" %in% names(allPixDT) && !is.null(sim$ageSpinupMin)){
    allPixDT[, ageSpinup := ages]
    allPixDT[ageSpinup < sim$ageSpinupMin, ageSpinup := sim$ageSpinupMin]
  }


  ## Prepare growth curves ----

  if (!is.null(sim$gcMeta) && !inherits(sim$gcMeta, "data.table")){
    sim$gcMeta <- tryCatch(
      data.table::as.data.table(sim$gcMeta),
      error = function(e) stop(
        "gcMeta could not be converted to data.table: ", e$message, call. = FALSE))
  }
  if (!is.null(sim$userGcM3) && !inherits(sim$userGcM3, "data.table")){
    sim$userGcM3 <- tryCatch(
      data.table::as.data.table(sim$userGcM3),
      error = function(e) stop(
        "userGcM3 could not be converted to data.table: ", e$message, call. = FALSE))
  }

  if (!is.null(sim$gcKey)){

    if (is.null(sim$gcMeta))   stop("'gcMeta' not found")
    if (is.null(sim$userGcM3)) stop("'userGcM3' not found")

    if (!all(sim$gcKey %in% names(sim$gcMeta)))   stop("'gcMeta' must contain all columns in `gcKey`")
    if (!all(sim$gcKey %in% names(sim$userGcM3))) stop("'userGcM3' must contain all columns in `gcKey`")
    if (!all(sim$gcKey %in% names(allPixDT)))     stop("'cohortLocators' must contain all columns in `gcKey`")

    if ("gcids" %in% sim$gcKey)           stop("'gcKey' cannot contain \"gcids\"")
    if ("gcids" %in% names(sim$gcMeta))   stop("'gcMeta' cannot contain \"gcids\"")
    if ("gcids" %in% names(sim$userGcM3)) stop("'userGcM3' cannot contain \"gcids\"")
    if ("gcids" %in% names(allPixDT))     stop("'cohortLocators' cannot contain \"gcids\"")

    # Create a new unique key defining each growth curve
    sim$curveID <- "gcids"

    allPixDT$gcids <- CBMutils::gcidsCreate(allPixDT[, .SD, .SDcols = sim$gcKey])

    sim$gcMeta <- cbind(
      gcids = CBMutils::gcidsCreate(sim$gcMeta[, .SD, .SDcols = sim$gcKey]),
      sim$gcMeta)
    data.table::setkey(sim$gcMeta, gcids)

    sim$userGcM3 <- cbind(
      gcids = CBMutils::gcidsCreate(sim$userGcM3[, .SD, .SDcols = sim$gcKey]),
      sim$userGcM3)
    data.table::setkeyv(sim$userGcM3, c("gcids", "Age"))

    # Add ecozone and spatial_unit_id to sim$gcMeta
    if (any(!c("ecozones", "spatial_unit_id") %in% names(sim$gcMeta))){

      ## NOTE: this could produce different results for GCIDs in more than one ecozone or SPU.
      gcSPU <- unique(allPixDT[, .(gcids, ecozone, ecozones = ecozone, spatial_unit_id)])
      # if (any(duplicated(gcSPU$gcids))) warning(
      #   "Growth curves are not unique to every CBM-CFS spatial_unit_id in the simulation. ",
      #   "Add \"ecozone\" and/or \"spatial_unit_id\" to sim$gcKey to make them unique.")

      sim$gcMeta <- merge(sim$gcMeta, gcSPU, by = intersect(names(gcSPU), names(sim$gcMeta)), all.x = TRUE)

      # Only include curves with this information available
      sim$gcMeta <- subset(sim$gcMeta, gcids %in% gcSPU$gcids)
    }
  }

  # Get species attributes
  if (!is.null(sim$gcMeta) && any(!c("species_id", "sw_hw", "canfi_species", "genus") %in% names(sim$gcMeta))){

    if (!"species" %in% names(sim$gcMeta)) stop(
      "gcMeta requires the 'species' column with species names to retrieve species data with CBMutils::sppMatch")

    sppMatchTable <- CBMutils::sppMatch(
      sim$gcMeta$species,
      return     = c("CBM_speciesID", "Broadleaf", "CanfiCode", "NFI"),
      otherNames = list(
        "White birch" = "Paper birch"
      ))[, .(
        species_id    = CBM_speciesID,
        sw_hw         = data.table::fifelse(Broadleaf, "hw", "sw"),
        canfi_species = CanfiCode,
        genus         = sapply(strsplit(NFI, "_"), `[[`, 1)
      )]

    sim$gcMeta <- cbind(
      sim$gcMeta[, .SD, .SDcols = setdiff(names(sim$gcMeta), names(sppMatchTable))],
      sppMatchTable)
    rm(sppMatchTable)
  }


  ## Prepare sim$disturbanceMeta ----

  if (!is.null(sim$disturbanceMeta) && !"disturbance_type_id" %in% names(sim$disturbanceMeta)){

    if (is.null(sim$dbPath)) stop("'dbPath' input required to set disturbanceMeta 'disturbance_type_id'")

    if (!inherits(sim$disturbanceMeta, "data.table")){
      sim$disturbanceMeta <- tryCatch(
        data.table::as.data.table(sim$disturbanceMeta),
        error = function(e) stop(
          "'disturbanceMeta' could not be converted to data.table: ", e$message, call. = FALSE))
    }

    # Match user disturbances with CBM-CFS3 disturbance type IDs
    askUser <- interactive() & !identical(Sys.getenv("TESTTHAT"), "true")
    if (askUser) message(
      "Prompting user to match input disturbances with CBM-CFS3 disturbances:")

    data.table::setnames(
      sim$disturbanceMeta, c("name", "description"), c("nameUser", "descUser"),
      skip_absent = TRUE)

    sim$disturbanceMeta <- cbind(
      sim$disturbanceMeta, CBMutils::distMatch(
        sim$disturbanceMeta$nameUser,
        dbPath = sim$dbPath,
        ask    = askUser
      ) |> Cache()
    )
  }


  ## Prepare sim$standDT and sim$cohortDT ----

  sim$standDT <- allPixDT[, .SD, .SDcols = intersect(
    c("pixelIndex", "area", "ecozone", "spatial_unit_id"), names(allPixDT))]
  data.table::setkey(sim$standDT, pixelIndex)

  if (is.null(sim$cohortDT)){
    allPixDT[, cohortID := pixelIndex]
    sim$cohortDT <- allPixDT[, .SD, .SDcols = intersect(
      c("cohortID", "pixelIndex", "ages", "ageSpinup", "gcids", names(sim$cohortLocators)), names(allPixDT))]
    data.table::setkey(sim$cohortDT, cohortID)
  }


  ## Return simList ----

  return(invisible(sim))

}

.inputObjects <- function(sim){

  ## Define study area ----

  # Master raster
  if (!suppliedElsewhere("masterRaster", sim) & suppliedElsewhere("masterRasterURL", sim)){

    sim$masterRaster <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$masterRasterURL
    )
  }

  # Canada ecozones
  if (!suppliedElsewhere("ecoLocator", sim)){

    if (suppliedElsewhere("ecoLocatorURL", sim) &
        !identical(sim$ecoLocatorURL, extractURL("ecoLocator"))){

      sim$ecoLocator <- prepInputs(
        destinationPath = inputPath(sim),
        url = sim$ecoLocatorURL
      )

    }else{

      ## 2024-12-04 NOTE:
      ## Multiple users had issues downloading and extracting this file via prepInputs.
      ## Downloading the ZIP directly and saving it in the inputs directory works OK.
      sim$ecoLocator <- tryCatch(

        prepInputs(
          destinationPath = inputPath(sim),
          url         = extractURL("ecoLocator"),
          filename1   = "ecozone_shp.zip",
          targetFile  = "ecozones.shp",
          alsoExtract = "similar",
          fun         = sf::st_read(targetFile, agr = "constant", quiet = TRUE)
        ),

        error = function(e) stop(
          "Canada ecozones Shapefile failed be downloaded and extracted:\n", e$message, "\n\n",
          "If this error persists, download the ZIP file directly and save it to the inputs directory.",
          "\nDownload URL: ", extractURL("ecoLocator"),
          "\nInputs directory: ", normalizePath(inputPath(sim), winslash = "/"),
          call. = FALSE))

      # Drop other fields
      sim$ecoLocator <- sim$ecoLocator[, "ECOZONE"]
    }
  }

  # Canada CBM-CFS3 spatial units
  if (!suppliedElsewhere("spuLocator", sim)){

    if (suppliedElsewhere("spuLocatorURL", sim) &
        !identical(sim$spuLocatorURL, extractURL("spuLocator"))){

      sim$spuLocator <- prepInputs(
        destinationPath = inputPath(sim),
        url = sim$spuLocatorURL
      )

    }else{

      sim$spuLocator <- prepInputs(
        destinationPath = inputPath(sim),
        url         = extractURL("spuLocator"),
        filename1   = "spUnit_Locator.zip",
        targetFile  = "spUnit_Locator.shp",
        alsoExtract = "similar",
        fun         = sf::st_read(targetFile, agr = "constant", quiet = TRUE)
      )

      # Drop other fields
      sim$spuLocator <- sim$spuLocator[, "spu_id"]
    }
  }

  # Stand age
  if (!suppliedElsewhere("ageLocator", sim) & suppliedElsewhere("ageLocatorURL", sim)){

    sim$ageLocator <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$ageLocatorURL
    )
  }

  if (!suppliedElsewhere("ageDataYear", sim) & !is.null(sim$ageLocator)){

    warning("'ageDataYear' not provided by user; `ageLocator` ages assumed to represent cohort age at simulation start")

    sim$ageDataYear <- as.numeric(start(sim))
  }


  ## Growth curves ---

  # Growth curve IDs
  if (!suppliedElsewhere("cohortLocators", sim) & suppliedElsewhere("cohortLocatorURLs", sim)){

    sim$cohortLocators <- lapply(sim$cohortLocatorURLs, function(url){
      prepInputs(
        destinationPath = inputPath(sim),
        url = url
      )
    })
  }

  # Growth curve metadata
  if (!suppliedElsewhere("gcMeta", sim) & suppliedElsewhere("gcMetaURL", sim)){

    sim$gcMeta <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$gcMetaURL,
      fun = data.table::fread
    )
  }

  # Growth curve increments
  if (!suppliedElsewhere("userGcM3", sim) & suppliedElsewhere("userGcM3URL", sim)){

    sim$userGcM3 <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$userGcM3URL,
      fun = data.table::fread
    )
  }


  ## Disturbances ----

  if (!suppliedElsewhere("disturbanceMeta") & suppliedElsewhere("disturbanceMetaURL", sim)){

    sim$disturbanceMeta <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$disturbanceMetaURL,
      fun = data.table::fread
    )
  }

  # CBM-CFS3 defaults database
  if (!suppliedElsewhere("dbPath", sim)){
    if (suppliedElsewhere("dbPathURL", sim)){

      sim$dbPath <- prepInputs(
        destinationPath = inputPath(sim),
        url = sim$dbPathURL
      )

    }else{

      sim$dbPath <- file.path(inputPath(sim), "cbm_defaults_v1.2.8340.362.db")

      if (!file.exists(sim$dbPath)) prepInputs(
        destinationPath = inputPath(sim),
        url         = extractURL("dbPath"),
        targetFile  = basename(sim$dbPath),
        dlFun       = download.file(extractURL("dbPath"), sim$dbPath, mode = "wb", quiet = TRUE),
        fun         = NA
      )
    }
  }


  ## Return simList ----

  return(invisible(sim))
}

