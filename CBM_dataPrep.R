
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
    "data.table", "RSQLite", "sf", "terra", "exactextractr", "gstat",
    "reproducible (>=2.1.2)", "digest", "googledrive",
    "PredictiveEcology/CBMutils@development (>=2.4.1.9001)",
    "PredictiveEcology/LandR@development"
  ),
  parameters = rbind(
    defineParameter("saveRasters", "logical", FALSE, NA, NA, "Save rasters of inputs aligned to the `masterRaster`"),
    defineParameter("ageBacktrack", "list", NA, NA, NA, "Age backtracking parameters"),
    defineParameter(".useCache", "character", "init", NA, NA, "Cache module events")
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "masterRaster", objectClass = "SpatRaster|character",
      desc = paste(
        "Raster template defining the study area. NA cells will be excluded from analysis.",
        "This can be provided as a SpatRaster or URL.")),
    expectsInput(
      objectName  = "adminLocator",
      objectClass = "sf|SpatRaster|sourceID|URL|character",
      sourceID    = "StatCan-admin",
      desc = paste(
        "Canada administrative boundary name(s).",
        "This can be provided as a spatial object, a `CBMutils::CBMsources` sourceID, a URL, or a single value for all cohorts.")),
    expectsInput(
      objectName  = "ecoLocator",
      objectClass = "sf|SpatRaster|sourceID|URL|numeric",
      sourceID    = "CanSIS-ecozone",
      desc = paste(
        "Canada ecozone ID(s).",
        "This can be provided as a spatial object, a `CBMutils::CBMsources` sourceID, a URL, or a single value for all cohorts.")),
    expectsInput(
      objectName  = "ageLocator",
      objectClass = "sf|SpatRaster|sourceID|URL|numeric",
      desc = paste(
        "Cohort ages at the simulation start year.",
        "This can be provided as a spatial object, a `CBMutils::CBMsources` sourceID, a URL, or a single value for all cohorts.")),
    expectsInput(
      objectName = "ageDataYear", objectClass = "numeric",
      desc = "Year that the ages in `ageLocator` represent."),
    expectsInput(
      objectName = "ageBacktrackSplit", objectClass = "character",
      desc = "Optional. If backtracking ages, split the age layer by these `cohortDT` or `gcMeta` columns when interpolating ages."),
    expectsInput(
      objectName = "ageSpinupMin", objectClass = "numeric",
      desc = "Minimum age for cohorts during spinup. Temporary fix to CBM_core issue #1: https://github.com/PredictiveEcology/CBM_core/issues/1"),
    expectsInput(
      objectName = "gcIndexLocator", objectClass = "sf|SpatRaster|character",
      desc = paste(
        "Growth curve ID(s).",
        "This can be provided as a spatial object, a URL, or a single value for all cohorts.",
        "If provided, IDs will be added to the 'curveID' column of `cohortDT` and `curveID` will be set to 'curveID'")),
    expectsInput(
      objectName  = "cohortLocators",
      objectClass = "list",
      desc = paste(
        "Named list of data sources defining cohorts.",
        "Each item may be a spatial object, a `CBMutils::CBMsources` sourceID, a URL, or a single value for all cohorts.")),
    expectsInput(
      objectName  = "CBMsourceIDs",
      objectClass = "character",
      desc = "`CBMutils::CBMsources` sourceID(s) to use as cohort locators."),
    expectsInput(
      objectName = "curveID", objectClass = "character",
      desc = "Column(s) uniquely defining each growth curve in `cohortDT` and `userGcMeta`."),
    expectsInput(
      objectName = "userGcMeta", objectClass = "data.table",
      desc = paste(
        "Growth curve metadata. An input to CBM_vol2biomass.",
        "If provided, species names or LandR codes will be matched with known species get additional attributes.")),
    expectsInput(
      objectName = "gcMeta", objectClass = "data.table",
      desc = paste(
        "Growth curve metadata. An input to CBM_core.",
        "If provided, species names or LandR codes will be matched with known species get additional attributes.")),
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
      objectName = "disturbanceSource", objectClass = "character",
      desc = paste(
        "Names of known disturbance sources to use. Can be one or more of: 'NTEMS'.",
        "NTEMS source: CA Forest Fires 1985-2020 and CA Forest Harvest 1985-2020 GeoTIFF layers.",
          "Hermosilla, T., M.A. Wulder, J.C. White, N.C. Coops, G.W. Hobart, L.B. Campbell, 2016.",
          "Mass data processing of time series Landsat imagery: pixels to data products for forest monitoring.",
          "International Journal of Digital Earth 9(11), 1035-1054 (Hermosilla et al. 2016)."
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
        age        = "Cohort ages extracted from `ageLocator`",
        ageSpinup  = "Cohort ages raised to >= `ageSpinupMin`",
        gcids      = "Growth curve ID unique to every spatial unit and `curveID`"
      )),
    createsOutput(
      objectName = "ageDataYear", objectClass = "numeric",
      desc = paste(
        "Year that the ages in `ageLocator` represent.",
        "If `ageLocator` is a `CBMutils::CBMsources` sourceID this will be automatically set.",
        "Otherwise, if omitted, ages are assumed to represent the simulation start year.")),
    expectsInput(
      objectName = "curveID", objectClass = "character",
      desc = paste(
        "Column(s) uniquely defining each growth curve in `cohortDT` and `userGcMeta`.",
        "Defaults to the nmes of the columns created by `cohortLocators` and `CBMsourceIDs`.")),
    createsOutput(
      objectName = "userGcSPU", objectClass = "data.table",
      desc = "Table of growth curve and spatial unit combinations `cohortDT`."),
    createsOutput(
      objectName = "userGcMeta", objectClass = "data.table",
      desc = "Growth curve metadata with additional species attributes.",
      columns = list(
        species_id    = "CBM-CFS3 species ID",
        sw_hw         = "'sw' or 'hw'",
        LandR         = "LandR species code",
        canfi_species = "CanFI species codes",
        genus         = "NFI species genus"
      )),
    createsOutput(
      objectName = "gcMeta", objectClass = "data.table",
      desc = "Growth curve metadata with additional species attributes.",
      columns = list(
        species_id    = "CBM-CFS3 species ID",
        sw_hw         = "'sw' or 'hw'",
        LandR         = "LandR species code",
        canfi_species = "CanFI species codes",
        genus         = "NFI species genus"
      )),
    createsOutput(
      objectName = "disturbanceMeta", objectClass = "data.table",
      desc = "Table defining `disturbanceEvents` event types."),
    createsOutput(
      objectName = "disturbanceEvents", objectClass = "data.table",
      desc = "Table of disturbance events.")
  )
))


## MODULE EVENTS ----

doEvent.CBM_dataPrep <- function(sim, eventTime, eventType, debug = FALSE) {

  switch(
    eventType,

    init = {

      # Prepare master raster
      sim <- PrepMasterRaster(sim)

      # Prepare cohorts
      sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "prepCohorts", eventPriority = 2)

      # Prepare species data
      sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "matchSpecies", eventPriority = 1)

      # Prepare disturbances
      sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "matchDisturbances", eventPriority = 8)
      sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "readDisturbances",  eventPriority = 8)

      if ("NTEMS" %in% sim$disturbanceSource){
        sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "readDisturbancesNTEMS", eventPriority = 1)
      }

      # CBM_vol2biomass prep
      if (!is.null(sim$curveID)){
        sim <- scheduleEvent(sim, start(sim), "CBM_dataPrep", "prepVol2Biomass", eventPriority = 4)
      }
    },

    prepCohorts = {

      # Read cohort data
      sim <- PrepCohorts(sim)

      # Adjust cohort ages
      if ("age" %in% names(sim$cohortDT) && !is.null(sim$ageDataYear) && start(sim) != sim$ageDataYear){

        # Read disturbances
        distYears <- sort(c(sim$ageDataYear, start(sim)))
        for (year in distYears[[1]]:(distYears[[2]]-1)){
          sim <- ReadDisturbances(sim, year = year)
        }

        # Step ages forward or backwards
        if (start(sim) > sim$ageDataYear) sim <- AgeStepForward(sim)
        if (start(sim) < sim$ageDataYear) sim <- AgeStepBackward(sim)
      }

      # Convert ages to integer; set spinup age
      if ("age" %in% names(sim$cohortDT)){

        if (!is.integer(sim$cohortDT$age)){
          sim$cohortDT[, age := as.integer(round(age))]
        }
        if (!is.null(sim$ageSpinupMin)){
          sim$cohortDT[, ageSpinup := age]
          sim$cohortDT[ageSpinup < sim$ageSpinupMin, ageSpinup := sim$ageSpinupMin]
        }
      }

      # Subset cohorts
      sim <- SubsetCohorts(sim)
    },

    prepVol2Biomass = {
      sim <- PrepVol2Biomass(sim)
    },

    matchSpecies = {
      sim <- MatchSpecies(sim)
    },

    matchDisturbances = {
      sim <- MatchDisturbances(sim)
    },
    readDisturbances = {
      sim <- ReadDisturbances(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "CBM_dataPrep", "readDisturbances", eventPriority = 8)
    },
    readDisturbancesNTEMS = {
      sim <- ReadDisturbancesNTEMS(sim)
    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

PrepMasterRaster <- function(sim){

  if (is.null(sim$masterRaster)) stop("masterRaster not found")

  if (!inherits(sim$masterRaster, "SpatRaster")){

    if (isURL(sim$masterRaster)){
      sim$masterRaster <- reproducible::prepInputs(
        destinationPath = inputPath(sim),
        url = sim$masterRaster,
        fun = terra::rast
      )

    }else{
      sim$masterRaster <- tryCatch(
        terra::rast(sim$masterRaster),
        error = function(e) stop(
          "masterRaster could not be converted to SpatRaster: ", e$message,
          call. = FALSE))
    }
  }

  # Mask cells outside of admin boundary
  if (is.character(sim$adminLocator) && length(sim$adminLocator) == 1 &&
      !terra::global(sim$masterRaster, "anyNA")[1, 1]){

    adminBoundaries <- CBMutils::CBMsourcePrepInputs("StatCan-admin")$source
    if (sim$adminLocator %in% adminBoundaries$admin){

      adminMask <- subset(adminBoundaries, admin == sim$adminLocator) |>
        sf::st_segmentize(10000) |>
        sf::st_transform(sf::st_crs(sim$masterRaster))
      sim$masterRaster <- terra::mask(sim$masterRaster, adminMask, touches = FALSE)
    }
  }

  return(invisible(sim))
}

PrepCohorts <- function(sim){

  # Initiate pixel table
  allPixDT <- data.table::data.table(
    pixelIndex = 1:terra::ncell(sim$masterRaster),
    key = "pixelIndex")

  # Set cell area
  if (!terra::is.lonlat(sim$masterRaster)){
    allPixDT$area <- prod(terra::res(sim$masterRaster) * terra::linearUnits(sim$masterRaster))

  }else{
    masterRasterCellSize <- terra::cellSize(
      sim$masterRaster, unit = "m", mask = FALSE, transform = FALSE) |> Cache()
    allPixDT$area <- terra::values(masterRasterCellSize, mat = FALSE) |> Cache()
  }

  # Set cohort attributes from input sources
  colInputs <- list(
    admin_name = sim$adminLocator,
    ecozone    = sim$ecoLocator,
    age        = sim$ageLocator,
    curveID    = sim$gcIndexLocator
  )

  if (length(sim$cohortLocators) > 0){

    if (!is.list(sim$cohortLocators) || is.null(names(sim$cohortLocators))) stop(
      "'cohortLocators' must be a named list")
    if (any(is.na(names(sim$cohortLocators)))) stop("'cohortLocators' names contains NAs")

    colInputs <- c(colInputs, sim$cohortLocators)
  }

  if (length(sim$CBMsourceIDs) > 0){

    if (!all(sim$CBMsourceIDs %in% CBMutils::CBMsources$sourceID)) stop(
      "sourceID(s) not found in `CBMutils::CBMsources$sourceID`: ",
      paste(shQuote(setdiff(sim$CBMsourceIDs, CBMutils::CBMsources$sourceID)), collapse = ", "))

    colInputs <- c(
      colInputs,
      with(subset(CBMutils::CBMsources, sourceID %in% sim$CBMsourceIDs), setNames(sourceID, attr)))
  }

  colInputs <- colInputs[!sapply(colInputs, is.null)]
  for (colName in names(colInputs)){

    if (isValue(colInputs[[colName]])){

      # Set column as a single value
      allPixDT[[colName]] <- colInputs[[colName]]

    }else{

      if (isCBMsource(colInputs[[colName]])){

        message("Extracting CBM source '", colInputs[[colName]], "' into column '", colName, "'")

        sourceCBM <- CBMutils::CBMsourceExtractToRast(colInputs[[colName]], sim$masterRaster) |> Cache()

        allPixDT[[colName]] <- sourceCBM$extractToRast

        if (colName == "age") sim$ageDataYear <- sourceCBM$year

        rm(sourceCBM)

      }else{

        message("Extracting spatial input data into column '", colName, "'")

        if (isURL(allPixDT[[colName]])){
          allPixDT[[colName]] <- reproducible::prepInputs(
            destinationPath = inputPath(sim),
            url             = allPixDT[[colName]])
        }

        allPixDT[[colName]] <- CBMutils::extractToRast(colInputs[[colName]], sim$masterRaster) |> Cache()
      }

      if (P(sim)$saveRasters){
        outPath <- file.path(outputPath(sim), "CBM_dataPrep", paste0("input_", colName, ".tif"))
        message("Writing aligned raster to path: ", outPath)
        tryCatch(
          CBMutils::writeRasterWithValues(sim$masterRaster, allPixDT[[colName]], outPath, overwrite = TRUE),
          error = function(e) warning(e$message, call. = FALSE))
      }

      if (is.factor(allPixDT[[colName]])) allPixDT[[colName]] <- as.character(allPixDT[[colName]])
    }
  }

  # Set cohort age data year if not set
  if ("age" %in% names(allPixDT) & is.null(sim$ageDataYear)){
    warning("'ageDataYear' not provided by user; `ageLocator` ages assumed to represent cohort age at simulation start")
    sim$ageDataYear <- as.numeric(start(sim))
  }

  # Set CBM-CFS3 spatial_unit_id
  if ("admin_name" %in% names(allPixDT)){

    cbmDBcon <- DBI::dbConnect(RSQLite::dbDriver("SQLite"), sim$dbPath)
    cbmDB <- list(
      admin_boundary_tr = DBI::dbReadTable(cbmDBcon, "admin_boundary_tr"),
      spatial_unit      = DBI::dbReadTable(cbmDBcon, "spatial_unit")
    ) |> lapply(data.table::data.table)
    RSQLite::dbDisconnect(cbmDBcon)

    # Set admin_abbrev & admin_boundary_id
    cbmDB$admin_boundary_tr <- cbmDB$admin_boundary_tr[locale_id == 1, .(admin_name = name, admin_boundary_id)]

    ## Add a row for "Yukon"
    cbmDB$admin_boundary_tr <- rbind(cbmDB$admin_boundary_tr, data.frame(
      admin_name = "Yukon",
      admin_boundary_id = cbmDB$admin_boundary_tr[admin_name == "Yukon Territory", admin_boundary_id]))

    cbmDB$admin_boundary_tr[, admin_abbrev := sapply(
      admin_name, switch,
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
      "Yukon Territory"           = "YT",
      "Yukon"                     = "YT",
      "Northwest Territories"     = "NT",
      "Nunavut"                   = "NU",
      NA_character_
    )]

    if (is.character(allPixDT$admin_name)){
      allPixDT <- merge(allPixDT, cbmDB$admin_boundary_tr, by = "admin_name", all.x = TRUE)

    }else{

      # Input is admin_boundary_id
      data.table::setnames(allPixDT, "admin_name", "admin_boundary_id")
      allPixDT <- merge(allPixDT, cbmDB$admin_boundary_tr, by = "admin_boundary_id", all.x = TRUE)
    }

    # Set CBM-CFS3 spatial_unit_id
    allPixDT <- merge(
      allPixDT,
      cbmDB$spatial_unit[, .(admin_boundary_id, ecozone = eco_boundary_id, spatial_unit_id = id)],
      by = c("admin_boundary_id", "ecozone"), all.x = TRUE)
  }

  # Create sim$standDT and sim$cohortDT
  data.table::setkey(allPixDT, pixelIndex)
  allPixDT[, admin_name := NULL]

  sim$standDT <- allPixDT[, .SD, .SDcols = intersect(
    c("pixelIndex", "area", "admin_abbrev", "admin_boundary_id", "ecozone", "spatial_unit_id"),
    names(allPixDT))]

  if (is.null(sim$cohortDT)){
    allPixDT[, cohortID := pixelIndex]
    sim$cohortDT <- allPixDT[, .SD, .SDcols = intersect(
      c("cohortID", "pixelIndex", "age", "ageSpinup", setdiff(names(colInputs), names(sim$standDT))),
      names(allPixDT))]
    data.table::setkey(sim$cohortDT, cohortID)
  }

  return(invisible(sim))
}

SubsetCohorts <- function(sim){

  # Subset stands and cohorts to cells where masterRaster is not NA
  sim$standDT  <- sim$standDT[terra::cells(sim$masterRaster),]
  sim$cohortDT <- sim$cohortDT[sim$standDT$pixelIndex,]
  if (nrow(sim$cohortDT) == 0) stop("all masterRaster values are NA")

  # Check spatial unit IDs
  for (col in c("admin_abbrev", "ecozone")){
    if (any(is.na(sim$standDT[[col]]))) stop(col, " NA in ", sum(is.na(sim$standDT[[col]])), " pixels")
  }
  if (any(is.na(sim$standDT$spatial_unit_id))){
    noMatch <- unique(sim$standDT[is.na(spatial_unit_id), .(admin_abbrev, ecozone, spatial_unit_id)])
    data.table::setkey(noMatch, admin_abbrev, ecozone)
    if (nrow(noMatch) > 0) stop(
      "spatial_unit_id not found for: ",
      paste(paste(noMatch$admin_abbrev, "ecozone", noMatch$ecozone), collapse = "; "))
  }

  # Remove cohorts that are missing key attributes
  cohortCols <- setdiff(names(sim$cohortDT), c("cohortID", "pixelIndex"))
  if (length(cohortCols) > 0){

    isNA  <- is.na(sim$cohortDT[, .SD, .SDcols = cohortCols])
    hasNA <- colSums(isNA) > 0

    if (any(hasNA)){

      sim$cohortDT <- sim$cohortDT[rowSums(isNA[, hasNA, drop = FALSE]) == 0,]

      rmMsg <- paste0(
        round((1 - nrow(sim$cohortDT) / nrow(isNA)) * 100, 2),
        "% of pixels excluded due to NAs in one or more of: ",
        paste(shQuote(names(hasNA)[hasNA]), collapse = ", "))
      if (nrow(sim$cohortDT) == 0) stop(rmMsg)
      message(rmMsg)
    }
    rm(isNA)
    rm(hasNA)
  }

  return(invisible(sim))
}

AgeStepForward <- function(sim){

  # WORK IN PROGRESS
  warning("Cohort age data is from ", sim$ageDataYear, " instead of the simulation start year",
          call. = FALSE)

  return(invisible(sim))
}

AgeStepBackward <- function(sim){

  # Skip process if there are multiple cohorts per pixel
  if (any(duplicated(sim$cohortDT$pixelIndex))){
    warning("Cohort age data is from ", sim$ageDataYear, " instead of the simulation start year",
            call. = FALSE)
    return(invisible(sim))
  }

  # Set function to backtrack ages
  cacheExtra <- list(
    masterRaster = digest::digest(sim$masterRaster),
    distEvents   = digest::digest(sim$disturbanceEvents)
  )
  ageStepBack <- function(pixelIndex, pixelAges, yearIn, yearOut, params, msgPrefix = NULL){

    ageRast <- terra::rast(sim$masterRaster, vals = NA)
    terra::set.values(ageRast, pixelIndex, pixelAges)

    ageRast <- withCallingHandlers(
      do.call(
        CBMutils::ageStepBackward,
        c(list(
          ageRast    = ageRast,
          yearIn     = yearIn,
          yearOut    = yearOut,
          distEvents = sim$disturbanceEvents
        ),
        if (!identical(P(sim)$ageBacktrack, NA)) params)
      ),
      message = function(m){
        message(msgPrefix, gsub("\\n", "", conditionMessage(m)))
        invokeRestart("muffleMessage")
      }
    )

    terra::extract(ageRast, pixelIndex)[,1]
  }

  # Backtrack ages
  if (is.null(sim$ageBacktrackSplit)){

    sim$cohortDT[, ageBacktrack := ageStepBack(
      pixelIndex = sim$cohortDT$pixelIndex,
      pixelAges  = sim$cohortDT$age,
      yearIn     = sim$ageDataYear,
      yearOut    = start(sim),
      params     = P(sim)$ageBacktrack
    ) |>
      Cache(omitArgs = "distEvents", .cacheExtra = cacheExtra)
    ]

  }else{

    # Split age data by splitting columns
    ageTable <- sim$cohortDT[!is.na(age),]

    colMissing <- setdiff(sim$ageBacktrackSplit, names(sim$cohortDT))
    if (length(colMissing) > 0){

      joinCol <- intersect(c("species", "LandR"), sim$curveID)[1]

      spsJoin <- lapply(c("userGcMeta", "gcMeta"), function(tbl){
        if (colMissing %in% names(sim[[tbl]])) unique(
          sim[[tbl]][, .SD, .SDcols = c(joinCol, colMissing)])
      })
      spsJoin <- spsJoin[!sapply(spsJoin, is.null)]

      if (length(spsJoin) == 0) stop(
        "ageBacktrackSplit column(s) not found: ",
        paste(shQuote(colMissing), collapse = ", "))

      ageTable <- merge(ageTable, spsJoin[[1]], by = joinCol, all.x = TRUE)
    }

    ageTable[, setdiff(names(ageTable), c("pixelIndex", "age", sim$ageBacktrackSplit)) := NULL]
    ageTable <- ageTable[rowSums(is.na(ageTable[, .SD, .SDcols = sim$ageBacktrackSplit])) == 0,]
    ageTable[, split := .GRP, by = eval(sim$ageBacktrackSplit)]
    ageTable <- split(ageTable, ageTable$split)

    # Backtrack ages
    for (spl in names(ageTable)){

      ageTable[[spl]][, age := ageStepBack(
        pixelIndex = ageTable[[spl]]$pixelIndex,
        pixelAges  = ageTable[[spl]]$age,
        yearIn     = sim$ageDataYear,
        yearOut    = start(sim),
        params     = P(sim)$ageBacktrack,
        msgPrefix  = paste0("Group ", spl, "/", length(ageTable), ": ")
      ) |>
        Cache(omitArgs = "distEvents", .cacheExtra = cacheExtra)
      ]
      ageTable[[spl]] <- ageTable[[spl]][, .(pixelIndex, ageBacktrack = age)]
    }

    ageTable <- data.table::rbindlist(ageTable)
    sim$cohortDT <- merge(sim$cohortDT, ageTable, "pixelIndex", all.x = TRUE)
    rm(ageTable)

    data.table::setkey(sim$cohortDT, cohortID)
  }

  # Rename ages columns
  data.table::setnames(sim$cohortDT, c("age", "ageBacktrack"), c(paste0("age", sim$ageDataYear), "age"))

  if (P(sim)$saveRasters){
    outPath <- file.path(outputPath(sim), "CBM_dataPrep", paste0("input_age_", start(sim), ".tif"))
    message("Writing backtracked age raster to path: ", outPath)
    tryCatch(
      CBMutils::writeRasterWithValues(sim$masterRaster, sim$cohortDT$age, outPath, overwrite = TRUE),
      error = function(e) warning(e$message, call. = FALSE))
  }

  return(invisible(sim))
}

PrepVol2Biomass <- function(sim){

  if (!all(sim$curveID %in% names(sim$cohortDT))) stop("cohortDT does not contain all columns in `curveID`")

  # Define unique growth curves with spatial_unit_id
  userGcSPU <- cbind(sim$standDT[sim$cohortDT$pixelIndex, .(spatial_unit_id)], sim$cohortDT[, .SD, .SDcols = sim$curveID])

  sim$userGcSPU <- unique(userGcSPU)

  sim$cohortDT[, gcids := factor(CBMutils::gcidsCreate(userGcSPU))]

  return(invisible(sim))
}

MatchSpecies <- function(sim){

  ## TEMPORARY: Add species missing from LandR::sppEquivalencies_CA
  sppEquiv <- LandR::sppEquivalencies_CA
  if (!177 %in% sppEquiv$CBM_speciesID){
    sppEquiv <- data.table::rbindlist(list(
      sppEquiv, data.frame(
        EN_generic_full = "Balsam poplar, largetooth aspen and eastern cottonwood",
        CBM_speciesID = 177,
        Broadleaf     = TRUE,
        CanfiCode     = 1211,
        NFI           = "POPU_",
        LandR         = "POPU_BAL"
      )), fill = TRUE)
  }

  # Get species attributes
  for (gcMetaTable in intersect(c("gcMeta", "userGcMeta"), objects(sim))){
    if (any(!c("species_id", "sw_hw", "canfi_species", "genus") %in% names(sim[[gcMetaTable]]))){

      if (!data.table::is.data.table(sim[[gcMetaTable]])){
        sim[[gcMetaTable]] <- data.table::as.data.table(sim[[gcMetaTable]])
      }

      matchCol <- intersect(c("LandR", "species"), names(sim[[gcMetaTable]]))[1]
      if (length(matchCol) == 0) stop(
        gcMetaTable, " requires column(s) 'species' and/or 'LandR' to retrieve species metadata")

      sppMatchTable <- CBMutils::sppMatch(
        sim[[gcMetaTable]][[matchCol]],
        sppEquivalencies = sppEquiv,
        return     = c("EN_generic_full", "CBM_speciesID", "Broadleaf", "CanfiCode", "NFI", "LandR"),
        otherNames = list(
          "White birch" = "Paper birch"
        ))[, .(
          species       = EN_generic_full,
          species_id    = CBM_speciesID,
          sw_hw         = data.table::fifelse(Broadleaf, "hw", "sw"),
          canfi_species = CanfiCode,
          genus         = sapply(strsplit(NFI, "_"), `[`, 1),
          LandR
        )]

      sim[[gcMetaTable]] <- cbind(
        sim[[gcMetaTable]][, .SD, .SDcols = setdiff(names(sim[[gcMetaTable]]), names(sppMatchTable))],
        sppMatchTable)
      rm(sppMatchTable)
    }
  }

  return(invisible(sim))
}

MatchDisturbances <- function(sim){

  if (isURL(sim$disturbanceMeta)){
    sim$disturbanceMeta <- prepInputs(
      destinationPath = inputPath(sim),
      url = sim$disturbanceMeta,
      fun = data.table::fread
    )
  }

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

  return(invisible(sim))
}

ReadDisturbances <- function(sim, year = time(sim)){

  # Get disturbances for the year
  distRasts <- lapply(sim$disturbanceRasters, function(d){
    if (as.character(time(sim)) %in% names(d)) d[[as.character(year)]]
  })

  # Retrieve disturbances from simList
  for (i in which(!is.na(sim$disturbanceMeta$sourceObjectName))){

    distRasts[[as.character(sim$disturbanceMeta[i,]$eventID)]] <- get(
      sim$disturbanceMeta[i,]$sourceObjectName, envir = sim)
  }

  # Summarize year events into a table
  distRasts <- distRasts[!sapply(distRasts, is.null)]
  if (length(distRasts) == 0) return(invisible(sim))

  if (is.null(names(distRasts)))    stop("disturbanceRasters list names must be disturbance event IDs")
  if (any(is.na(names(distRasts)))) stop("disturbanceRasters list names contains NAs")
  if (any(names(distRasts) == ""))  stop("disturbanceRasters list names contains empty strings")

  eventIDs <- suppressWarnings(tryCatch(
    as.integer(names(distRasts)),
    error = function(e) stop("disturbanceRasters list names must be coercible to integer")))

  newEvents <- lapply(1:length(distRasts), function(i){

    distMeta <- if (!is.null(sim$disturbanceMeta)){
      x <- as.list(subset(sim$disturbanceMeta, eventID == eventIDs[[i]]))
      x[!sapply(x, is.na)]
    }else list(eventID = eventIDs[[1]])

    with(distMeta, message(
      year, ": ",
      "Reading disturbances for eventID = ", eventID,
      if (exists("disturbance_type_id")) paste("; CBM type ID =", disturbance_type_id),
      if (exists("name"))                paste("; name =", shQuote(name))))

    distValues <- CBMutils::extractToRast(distRasts[[i]], sim$masterRaster) |> Cache()

    if (P(sim)$saveRasters){
      outPath <- file.path(outputPath(sim), "CBM_dataPrep", sprintf(
        "distEvents-%s_%s-%s.tif", eventIDs[[i]], year, i))
      message("Writing aligned raster to path: ", outPath)
      tryCatch(
        CBMutils::writeRasterWithValues(sim$masterRaster, distValues, outPath, overwrite = TRUE),
        error = function(e) warning(e$message, call. = FALSE))
    }

    if (length(na.omit(distMeta$sourceValue)) == 1){
      eventIndex <- which(distValues %in% distMeta$sourceValue)
    }else{
      eventIndex <- which(!is.na(distValues))
    }

    data.table::data.table(
      pixelIndex = eventIndex,
      year       = as.integer(year + c(na.omit(distMeta$sourceDelay), 0)[[1]]),
      eventID    = eventIDs[[i]]
    )
  })

  sim$disturbanceEvents <- data.table::rbindlist(c(list(sim$disturbanceEvents), newEvents))

  return(invisible(sim))
}

ReadDisturbancesNTEMS <- function(sim){

  if (any(c(1001, 1002) %in% sim$disturbanceMeta$eventID)) stop(
    "NTEMS disturbances reserve eventIDs 1001 and 1002")

  newDist <- rbind(
    data.table(
      eventID             = 1001L,
      disturbance_type_id = 1, # Wildfire
      name                = "NTEMS CA Forest Fires 1985-2020",
      url                 = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip"
    ),
    data.table(
      eventID             = 1002L,
      disturbance_type_id = 204, # Clearcut harvesting without salvage
      name                = "NTEMS CA Forest Harvest 1985-2020",
      url                 = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Harvest_1985-2020.zip"
    )
  )

  sim$disturbanceMeta <- data.table::rbindlist(list(
    sim$disturbanceMeta, newDist[, 1:3]), fill = TRUE)

  newEvents <- lapply(1:nrow(newDist), function(i){

    url <- newDist[i,]$url
    sourceTIF <- reproducible::prepInputs(
      url,
      destinationPath = inputPath(sim),
      filename1   = basename(url),
      targetFile  = paste0(tools::file_path_sans_ext(basename(url)), ".tif"),
      alsoExtract = "similar",
      fun         = NA
    ) |> Cache()

    with(newDist[i,], message(
      "Reading NTEMS disturbances for eventID = ", eventID,
      "; CBM type ID = ", disturbance_type_id,
      "; name = ", shQuote(name)))

    distValues <- CBMutils::extractToRast(sourceTIF, sim$masterRaster) |> Cache()

    if (P(sim)$saveRasters){
      outPath <- file.path(outputPath(sim), "CBM_dataPrep", paste0(newDist[i,]$name, '.tif'))
      message("Writing aligned raster to path: ", outPath)
      tryCatch(
        CBMutils::writeRasterWithValues(sim$masterRaster, distValues, outPath, overwrite = TRUE),
        error = function(e) warning(e$message, call. = FALSE))
    }

    data.table::data.table(
      pixelIndex = 1:length(distValues),
      year       = as.integer(distValues),
      eventID    = newDist[i,]$eventID
    ) |> subset(!is.na(year)) |> subset(year != 0)
  })

  sim$disturbanceEvents <- data.table::rbindlist(c(list(sim$disturbanceEvents), newEvents))

  return(invisible(sim))
}

.inputObjects <- function(sim){

  # CBM-CFS3 defaults database
  if (!suppliedElsewhere("dbPath", sim)){

    sim$dbPath <- file.path(inputPath(sim), "cbm_defaults_v1.2.8340.362.db")

    if (!file.exists(sim$dbPath)) prepInputs(
      destinationPath = inputPath(sim),
      url         = extractURL("dbPath"),
      targetFile  = basename(sim$dbPath),
      dlFun       = download.file(extractURL("dbPath"), sim$dbPath, mode = "wb", quiet = TRUE),
      fun         = NA
    )
  }

  # Canada admin boundaries & ecozones
  defaultSourceIDs <- with(
    list(x = inputObjects(sim, "CBM_dataPrep")),
    sapply(split(x$sourceID, x$objectName), unlist))

  if (!suppliedElsewhere("adminLocator", sim)) sim$adminLocator <- defaultSourceIDs[["adminLocator"]]
  if (!suppliedElsewhere("ecoLocator",   sim)) sim$ecoLocator   <- defaultSourceIDs[["ecoLocator"]]

  # Growth curve ID
  if (!suppliedElsewhere("curveID", sim)){

    if (suppliedElsewhere("gcIndexLocator")){
      sim$curveID <- "curveID"

    }else{
      curveID <- c(
        names(sim$cohortLocators)[!sapply(sim$cohortLocators, is.null)],
        setdiff(subset(CBMutils::CBMsources, sourceID %in% sim$CBMsourceIDs)$attr, "age")
      )
      if (length(curveID) > 0) sim$curveID <- unique(curveID)
    }
  }

  # Return simList
  return(invisible(sim))
}
