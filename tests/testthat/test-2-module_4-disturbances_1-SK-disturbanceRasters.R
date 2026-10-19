
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: SK with disturbanceRasters", {

  ## Run simInit and spades ----

  # Set up project
  projectName <- "SK-disturbanceRasters"
  times       <- list(start = 1998, end = 1999)

  simInitInput <- SpaDES.project::setupProject(

    modules = "CBM_dataPrep",
    times   = times,
    paths   = list(
      projectPath = spadesTestPaths$projectPath,
      modulePath  = spadesTestPaths$modulePath,
      packagePath = spadesTestPaths$packagePath,
      inputPath   = spadesTestPaths$inputPath,
      cachePath   = spadesTestPaths$cachePath,
      outputPath  = file.path(spadesTestPaths$temp$outputs, projectName)
    ),

    # Set required packages for project set up
    require = "terra",

    # Set study area
    masterRaster = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
      resolution = 30,
      vals       = 1
    ),

    # Set disturbances
    ## Test matching user disturbances with CBM-CFS3 disturbances
    disturbanceMeta = rbind(

      # Events are only in pixels with a value of 1
      # Events are delayed by 1 year
      data.frame(eventID = 1, name = "Wildfire",
                 sourceValue = 1, sourceDelay = 1, sourceObjectName = NA_character_),

      # Events come from the 'clearcut' object
      data.frame(eventID = 2, name = "Clearcut harvesting without salvage",
                 sourceValue = NA_integer_, sourceDelay = NA_integer_, sourceObjectName = "clearcut")
    ),
    disturbanceRasters = list(
      `1` = list(
        `1998` = terra::rast(
          crs        = "EPSG:3979",
          extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
          resolution = 30,
          vals       = c(1, 2, rep(NA, 31300))
        ))),
    clearcut = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
      resolution = 30,
      vals       = c(1, rep(NA, 31301))
    )
  )

  # Run simInit
  simTestInit <- SpaDES.core::simInit2(simInitInput)
  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDES.core::spades(simTestInit)
  expect_s4_class(simTest, "simList")


  ## Check outputs

  # disturbanceMeta
  expect_true(!is.null(simTest$disturbanceMeta))
  expect_true(inherits(simTest$disturbanceMeta, "data.table"))

  expect_equal(nrow(simTest$disturbanceMeta), 2)

  ## Check that disturbances have been matched with CBM-CFS3 disturbances correctly
  expect_equal(simTest$disturbanceMeta$disturbance_type_id, c(1, 204))

  # disturbanceEvents
  expect_true(!is.null(simTest$disturbanceEvents))
  expect_true(inherits(simTest$disturbanceEvents, "data.table"))

  for (colName in c("pixelIndex", "year", "eventID")){
    expect_true(colName %in% names(simTest$disturbanceEvents))
    expect_true(is.integer(simTest$disturbanceEvents[[colName]]))
    expect_true(all(!is.na(simTest$disturbanceEvents[[colName]])))
  }

  expect_equal(
    simTest$disturbanceEvents[order(eventID, year, pixelIndex)],
    rbind(
      data.table(pixelIndex = 1, year = 1999, eventID = 1),
      data.table(pixelIndex = 1, year = 1998, eventID = 2),
      data.table(pixelIndex = 1, year = 1999, eventID = 2)
    )
  )


  ## Update disturbances and continue simulation ----

  # Other modules may add disturbance rasters during the simulation
  simTest$disturbanceRasters$`1`$`2000` <- terra::rast(
    crs        = "EPSG:3979",
    extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
    resolution = 30,
    vals       = c(3, 1, rep(NA, 31300))
  )

  # Other modules may update the disturbance raster during the simulation
  simTest$clearcut <- terra::rast(
    crs        = "EPSG:3979",
    extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
    resolution = 30,
    vals       = c(NA, 1, rep(NA, 31300))
  )

  # Run spades
  SpaDES.core::end(simTest) <- 2000
  simTest <- SpaDES.core::spades(simTest)
  expect_s4_class(simTest, "simList")

  # Check disturbanceEvents
  expect_equal(
    simTest$disturbanceEvents[order(eventID, year, pixelIndex)],
    rbind(
      data.table(pixelIndex = 1, year = 1999, eventID = 1),
      data.table(pixelIndex = 2, year = 2001, eventID = 1),
      data.table(pixelIndex = 1, year = 1998, eventID = 2),
      data.table(pixelIndex = 1, year = 1999, eventID = 2),
      data.table(pixelIndex = 2, year = 2000, eventID = 2)
    )
  )
})


