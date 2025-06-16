
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
    "data.table", "RSQLite", "sf", "terra", "exactextractr",
    "reproducible (>=2.1.2)" ,
    "PredictiveEcology/CBMutils@development (>=2.0.3.0005)",
    "PredictiveEcology/LandR@development"
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
      objectName = "adminLocator", objectClass = "sf|SpatRaster",
      sourceURL = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
      desc = paste(
        "Spatial data source of Canada's administrative boundary names or CBM-CFS3 'admin_boundary_id'.",
        "Default is Canada's provinces and territories as polygon features.")),
    expectsInput(
      objectName = "adminLocatorURL", objectClass = "character", desc = "URL for `adminLocator`"),
    expectsInput(
      objectName = "ageLocator", objectClass = "sf|SpatRaster",
      desc = "Spatial data source of cohort ages."),
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
      objectName = "gcIndexLocator", objectClass = "sf|SpatRaster",
      desc = paste(
        "Spatial data source of growth curve locations.",
        "If provided, IDs will be added to the 'curveID' column of `cohortDT` and `curveID` will be set to 'curveID'")),
    expectsInput(
      objectName = "gcIndexLocatorURL", objectClass = "character", desc = "URL for `gcIndexLocator`"),
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
      objectName = "curveID", objectClass = "character",
      desc = paste(
        "Column(s) uniquely defining each growth curve in `cohortDT`, `userGcMeta`, and `userGcM3`.",
        "Each column must have a corresponding named spatial data source in `cohortLocators`")),
    expectsInput(
      objectName = "userGcMeta", objectClass = "data.table",
      desc = paste(
        "Growth curve metadata for CBM_vol2biomass.",
        "If provided, species names will be matched with known species get additional attributes."),
      columns = list(
        species = "Species name"
      )),
    expectsInput(
      objectName = "gcMeta", objectClass = "data.table",
      desc = paste(
        "Growth curve metadata.",
        "If provided, species names will be matched with known species get additional attributes."),
      columns = list(
        species = "Species name"
      )),
    expectsInput(
      objectName = "disturbanceRasters", objectClass = "list",
      desc = paste(
        "Set of spatial data sources containing locations of disturbance events for each year.",
        "List items must be named by disturbance event IDs found in `disturbanceMeta`.",
        "Within each event's list, items must be named by the 4 digit year the disturbances occured in.",
        "For example, event type 1 disturbance locations for 2025 can be accessed with `disturbanceRasters[[\"1\"]][[\"2025\"]]`.",
        "Each disturbance item can be one of the following:",
        "a terra SpatRaster layer, one or more raster file paths, or sf polygons.",
        "All non-NA areas will be considered events unless the 'sourceValue' column is set."
      )),
    expectsInput(
      objectName = "disturbanceMeta", objectClass = "data.table",
      desc = "Table defining the disturbance event types",
      columns = c(
        eventID             = "Event type ID",
        disturbance_type_id = "CBM-CFS3 disturbance type ID. If not provided, the user will be prompted to choose IDs.",
        name                = "Disturbance name (e.g. 'Wildfire'). Required only if 'disturbance_type_id' absent.",
        sourceValue         = "Optional. Value in `disturbanceRasters` to include as events",
        sourceDelay         = "Optional. Delay (in years) of when the `disturbanceRasters` will take effect",
        sourceObjectName    = "Optional. Name of the object in the `simList` to retrieve the `disturbanceRasters` from annually."
      )),
    expectsInput(
      objectName = "disturbanceMetaURL", objectClass = "character", desc = "URL for `disturbanceMeta`"),
    expectsInput(
      objectName = "dbPath", objectClass = "character",
      sourceURL = "https://raw.githubusercontent.com/cat-cfs/libcbm_py/main/libcbm/resources/cbm_defaults_db/cbm_defaults_v1.2.8340.362.db",
      desc = "Path to the CBM-CBM3 defaults database")
  ),
  outputObjects = bindrows(
    createsOutput(
      objectName = "standDT", objectClass = "data.table",
      desc = "Table of stand attributes.",
      columns = c(
        pixelIndex         = "`masterRaster` cell index",
        area               = "`masterRaster` cell area in meters",
        admin_abbrev       = "Canada administrative abbreviation extracted from `adminLocator`",
        admin_boundary_id  = "CBM-CFS3 administrative boundary ID",
        ecozone            = "Canada ecozone ID extracted from `ecoLocator`",
        spatial_unit_id    = "CBM-CFS3 spatial unit ID"
      )),
    createsOutput(
      objectName = "cohortDT", objectClass = "data.table",
      desc = "Table of cohort attributes.",
      columns = c(
        cohortID   = "`masterRaster` cell index",
        pixelIndex = "`masterRaster` cell index",
        ages       = "Cohort ages extracted from input `ageLocator`",
        ageSpinup  = "Cohort ages raised to minimum of `ageSpinupMin` to use in the spinup",
        gcids      = "Growth curve ID unique to every spatial unit and `curveID`"
      )),
    createsOutput(
      objectName = "userGcSPU", objectClass = "data.table",
      desc = "Table of growth curves and spatial unit combinations in the cohorts.",
      columns = list(
        curveID         = "Growth curve ID",
        spatial_unit_id = "CBM-CFS3 spatial unit ID"
      )),
    createsOutput(
      objectName = "userGcMeta", objectClass = "data.table",
      desc = "Growth curve metadata with additional species attributes.",
      columns = list(
        species_id      = "CBM-CFS3 species ID",
        sw_hw           = "'sw' or 'hw'",
        canfi_species   = "CanFI species codes",
        genus           = "NFI species genus"
      )),
    createsOutput(
      objectName = "gcMeta", objectClass = "data.table",
      desc = "Growth curve metadata with additional species attributes.",
      columns = list(
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
      objectName = "disturbanceMeta", objectClass = "data.table",
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
      for (i in which(!is.na(sim$disturbanceMeta$sourceObjectName))){

        distRasts[[as.character(sim$disturbanceMeta[i,]$eventID)]] <- get(
          sim$disturbanceMeta[i,]$sourceObjectName, envir = sim)
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

          distMeta <- if (!is.null(sim$disturbanceMeta)) subset(sim$disturbanceMeta, eventID == eventIDs[[i]])

          distAlign <- prepInputsToMasterRaster(
            distRasts[[i]],
            sim$masterRaster
          ) |> Cache()

          distValues <- terra::values(distAlign)[,1]
          if (length(na.omit(distMeta$sourceValue)) == 1){
            distValues[!distValues %in% distMeta$sourceValue] <- NA
          }

          eventIndex <- which(!is.na(distValues))
          if (length(eventIndex) > 0){

            sim$disturbanceEvents <- rbind(sim$disturbanceEvents, data.table::data.table(
              pixelIndex = eventIndex,
              year       = as.integer(time(sim) + c(na.omit(distMeta$sourceDelay), 0)[[1]]),
              eventID    = eventIDs[[i]]
            ))
          }

          rm(distValues)
          rm(eventIndex)

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

  ## Prepare sim$standDT and sim$cohortDT ----

  if (length(sim$cohortLocators) > 0){
    if (!is.list(sim$cohortLocators) || is.null(names(sim$cohortLocators))) stop(
      "'cohortLocators' must be a named list")
    if (any(is.na(names(sim$cohortLocators)))) stop("'cohortLocators' names contains NAs")
  }

  # Set which columns come from which input object
  colInputs <- c(
    list(
      admin_name = sim$adminLocator,
      ecozone    = sim$ecoLocator,
      ages       = sim$ageLocator,
      curveID    = sim$gcIndexLocator
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
    if (!is.null(terra::cats(inAlign)[[1]])){
      allPixDT[[colName]] <- merge(allPixDT[[colName]], terra::cats(inAlign)[[1]], by = 1, all.x = TRUE)[[2]]
    }
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

  # Set CBM-CFS3 spatial_unit_id
  if ("admin_name" %in% names(allPixDT)){

    cbmDBcon <- DBI::dbConnect(RSQLite::dbDriver("SQLite"), sim$dbPath)
    cbmDB <- list(
      admin_boundary_tr = data.table::data.table(DBI::dbReadTable(cbmDBcon, "admin_boundary_tr")) |> subset(locale_id == 1),
      spatial_unit      = data.table::data.table(DBI::dbReadTable(cbmDBcon, "spatial_unit"))
    )
    RSQLite::dbDisconnect(cbmDBcon)

    ## Add a row for "Yukon"
    cbmDB$admin_boundary_tr <- rbind(
      cbmDB$admin_boundary_tr,
      cbind(cbmDB$admin_boundary_tr[name == "Yukon Territory", !"name"], name = "Yukon")
    )

    if (is.character(allPixDT$admin_name)){

      allPixDT <- merge(
        allPixDT,
        cbmDB$admin_boundary_tr[, .(admin_boundary_id, admin_name = name)],
        by = "admin_name", all.x = TRUE)

      if (any(is.na(allPixDT$admin_boundary_id))) stop(
        "adminLocator name(s) not found in admin_boundary_tr:",
        paste(shQuote(unique(subset(allPixDT, is.na(admin_boundary_id))$admin_name)),
              collapse = ", "))

    }else{

      # Input is admin_boundary_id
      data.table::setnames(allPixDT, "admin_name", "admin_boundary_id")
      allPixDT <- merge(
        allPixDT,
        cbmDB$admin_boundary_tr[, .(admin_boundary_id, admin_name = name)],
        by = "admin_boundary_id", all.x = TRUE)
    }

    # Set admin abbreviation
    allPixDT$admin_abbrev <- c(
      "Newfoundland"              = "NL",
      "Labrador"                  = "NL",
      "Newfoundland and Labrador" = "NL",
      "Prince Edward Island"      = "PE",
      "Nova Scotia"               = "NS",
      "New Brunswick"             = "NB",
      "Quebec"                    = "QC",
      "Ontario"                   = "ON",
      "Manitoba"                  = "MB",
      "Alberta"                   = "AB",
      "Saskatchewan"              = "SK",
      "British Columbia"          = "BC",
      "Yukon"                     = "YT",
      "Yukon Territory"           = "YT",
      "Northwest Territories"     = "NT",
      "Nunavut"                   = "NU")[allPixDT$admin_name] |> unname()

    # Set CBM-CFS3 spatial_unit_id
    allPixDT <- merge(
      allPixDT,
      cbmDB$spatial_unit[, .(admin_boundary_id, ecozone = eco_boundary_id, spatial_unit_id = id)],
      by = c("admin_boundary_id", "ecozone"), all.x = TRUE)
    data.table::setkey(allPixDT, pixelIndex)

    if (any(is.na(allPixDT$spatial_unit_id))){
      noMatch <- unique(allPixDT[is.na(allPixDT$spatial_unit_id), .(ecozone, admin_name, spatial_unit_id)])
      data.table::setkey(noMatch, admin_name, ecozone)
      stop("spatial_unit_id not found for: ",
           paste(paste(noMatch$admin_name, "ecozone", noMatch$ecozone), collapse = "; "))
    }
  }

  # Adjust cohort ages
  if ("ages" %in% names(allPixDT) && !is.null(sim$ageDataYear) && sim$ageDataYear != start(sim)){

    # TODO: add step to adjust cohort ages to the simulation start year
    warning("Cohort age data is from ", sim$ageDataYear,
            " instead of the simulation start year")
  }

  # Set spinup age
  if ("ages" %in% names(allPixDT) && !is.null(sim$ageSpinupMin)){
    allPixDT[, ageSpinup := ages]
    allPixDT[ageSpinup < sim$ageSpinupMin, ageSpinup := sim$ageSpinupMin]
  }

  # Set growth curve ID
  if (!is.null(sim$gcIndexLocator)) sim$curveID <- "curveID"
  if (!is.null(sim$curveID)){

    if (!all(sim$curveID %in% names(allPixDT))) stop("'cohortLocators' must contain all columns in `curveID`")

    # Define unique growth curves with spatial_unit_id
    allPixDT$gcids <- factor(
      CBMutils::gcidsCreate(allPixDT[, .SD, .SDcols = c("spatial_unit_id", sim$curveID)])
    )
    sim$userGcSPU <- unique(allPixDT[, .SD, .SDcols = c("spatial_unit_id", sim$curveID)])
  }

  # Create sim$standDT and sim$cohortDT
  sim$standDT <- allPixDT[, .SD, .SDcols = intersect(
    c("pixelIndex", "area", "admin_abbrev", "admin_boundary_id", "ecozone", "spatial_unit_id"),
    names(allPixDT))]
  data.table::setkey(sim$standDT, pixelIndex)

  if (is.null(sim$cohortDT)){
    allPixDT[, cohortID := pixelIndex]
    sim$cohortDT <- allPixDT[, .SD, .SDcols = intersect(
      c("cohortID", "pixelIndex", "ages", "ageSpinup", "gcids", names(sim$cohortLocators)),
      names(allPixDT))]
    data.table::setkey(sim$cohortDT, cohortID)
  }


  ## Prepare sim$userGcMeta and sim$gcMeta ----

  # Get species attributes
  for (gcMetaTable in intersect(c("gcMeta", "userGcMeta"), objects(sim))){
    if (any(!c("species_id", "sw_hw", "canfi_species", "genus") %in% names(sim[[gcMetaTable]]))){

      if (!"species" %in% names(sim[[gcMetaTable]])) stop(
        gcMetaTable, " requires the 'species' names column to retrieve species data with CBMutils::sppMatch")

      ## TEMPORARY: Add species missing from LandR::sppEquivalencies_CA
      sppEquiv <- LandR::sppEquivalencies_CA
      if (!177 %in% sppEquiv$CBM_speciesID){
        sppEquiv <- data.table::rbindlist(list(
          sppEquiv, data.frame(
            EN_generic_full = "Balsam poplar, largetooth aspen and eastern cottonwood",
            CBM_speciesID = 177,
            Broadleaf     = TRUE,
            CanfiCode     = 1211,
            NFI           = "POPU_"
          )), fill = TRUE)
      }

      sppMatchTable <- CBMutils::sppMatch(
        sim[[gcMetaTable]]$species,
        sppEquivalencies = sppEquiv,
        return     = c("CBM_speciesID", "Broadleaf", "CanfiCode", "NFI"),
        otherNames = list(
          "White birch" = "Paper birch"
        ))[, .(
          species_id    = CBM_speciesID,
          sw_hw         = data.table::fifelse(Broadleaf, "hw", "sw"),
          canfi_species = CanfiCode,
          genus         = sapply(strsplit(NFI, "_"), `[[`, 1)
        )]

      sim[[gcMetaTable]] <- cbind(
        sim[[gcMetaTable]][, .SD, .SDcols = setdiff(names(sim[[gcMetaTable]]), names(sppMatchTable))],
        sppMatchTable)
      rm(sppMatchTable)
    }
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


  ## Return simList ----

  return(invisible(sim))

}

.inputObjects <- function(sim){

  ## Databases ----

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


  ## Define stands and cohorts ----

  # Master raster
  if (!suppliedElsewhere("masterRaster", sim) & suppliedElsewhere("masterRasterURL", sim)){

    sim$masterRaster <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$masterRasterURL
    )
  }

  # Canada admin boundaries
  if (!suppliedElsewhere("adminLocator", sim)){

    if (suppliedElsewhere("adminLocatorURL", sim) &
        !identical(sim$adminLocatorURL, extractURL("adminLocator"))){

      sim$adminLocator <- prepInputs(
        destinationPath = inputPath(sim),
        url = sim$adminLocatorURL
      )

    }else{

      sim$adminLocator <- prepInputs(
        destinationPath = inputPath(sim),
        url         = extractURL("adminLocator"),
        filename1   = "lpr_000a21a_e.zip",
        targetFile  = "lpr_000a21a_e.shp",
        alsoExtract = "similar",
        fun         = sf::st_read(targetFile, agr = "constant", quiet = TRUE)
      )

      # Split Newfoundland and Labrador
      sim$adminLocator <- cbind(name = sim$adminLocator$PRENAME, sim$adminLocator)

      adminSplit <- terra::split(
        terra::vect(sim$adminLocator),
        terra::vect(sf::st_sfc(sf::st_linestring(rbind(c(8476500, 2297500), c(8565300, 2451300))),
                               crs = sf::st_crs(sim$adminLocator))))
      adminSplit <- sf::st_as_sf(adminSplit, agr = "constant")
      sf::st_agr(adminSplit) <- "constant"

      nl_cd <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(
        adminSplit[adminSplit$name == "Newfoundland and Labrador",])))
      adminSplit[adminSplit$name == "Newfoundland and Labrador", "name"] <- sapply(
        nl_cd[, "X"] == min(nl_cd[, "X"]), ifelse, "Labrador", "Newfoundland")

      sim$adminLocator <- adminSplit
    }
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

  # Cohort ages
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

  # Growth curve locations
  if (!suppliedElsewhere("gcIndexLocator", sim) & suppliedElsewhere("gcIndexLocatorURL", sim)){

    sim$gcIndexLocator <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$gcIndexLocatorURL
    )
  }

  # Other cohort data
  if (!suppliedElsewhere("cohortLocators", sim) & suppliedElsewhere("cohortLocatorURLs", sim)){

    sim$cohortLocators <- lapply(sim$cohortLocatorURLs, function(url){
      prepInputs(
        destinationPath = inputPath(sim),
        url = url
      )
    })
  }


  ## Disturbances ----

  if (!suppliedElsewhere("disturbanceMeta") & suppliedElsewhere("disturbanceMetaURL", sim)){

    sim$disturbanceMeta <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$disturbanceMetaURL,
      fun = data.table::fread
    )
  }


  ## Return simList ----

  return(invisible(sim))
}

