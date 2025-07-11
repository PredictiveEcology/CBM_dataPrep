
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: RIA-small", {

  ## Run simInit and spades ----

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "2-module_2-RIA-small")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

      times = list(start = 2020, end = 2021),

      modules = "CBM_dataPrep",
      paths   = list(
        projectPath = projectPath,
        modulePath  = spadesTestPaths$modulePath,
        packagePath = spadesTestPaths$packagePath,
        inputPath   = spadesTestPaths$inputPath,
        cachePath   = spadesTestPaths$cachePath,
        outputPath  = file.path(projectPath, "outputs")
      ),

      # Set required packages for project set up
      require = "terra",

      # Set study area
      masterRaster = terra::rast(
        crs        = "EPSG:3979",
        extent     = c(xmin = -1653000, xmax = -1553000, ymin = 1180000, ymax = 1280000),
        resolution = 250,
        vals       = 1
      ),

      # Set disturbances
      disturbanceSource = "NTEMS"
    )
  )

  # Run simInit
  simTestInit <- SpaDEStestMuffleOutput(
    SpaDES.core::simInit2(simInitInput)
  )

  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDEStestMuffleOutput(
    SpaDES.core::spades(simTestInit)
  )

  expect_s4_class(simTest, "simList")


  ## Check output 'standDT' ----

  expect_true(!is.null(simTest$standDT))
  expect_true(inherits(simTest$standDT, "data.table"))

  for (colName in c("pixelIndex", "area", "admin_abbrev", "admin_boundary_id", "ecozone", "spatial_unit_id")){
    expect_true(colName %in% names(simTest$standDT))
    expect_true(all(!is.na(simTest$standDT[[colName]])))
  }

  expect_identical(data.table::key(simTest$standDT), "pixelIndex")

  expect_equal(nrow(simTest$standDT), 160000)
  expect_equal(simTest$standDT$pixelIndex, 1:160000)
  expect_in(simTest$standDT$area,              250*250)
  #expect_in(simTest$standDT$admin_name,        "British Columbia") # Column excluded from result
  expect_in(simTest$standDT$admin_abbrev,      "BC")
  expect_in(simTest$standDT$admin_boundary_id, 11)
  expect_in(simTest$standDT$ecozone,           c(4, 9, 12, 14))
  expect_in(simTest$standDT$spatial_unit_id,   c(38, 39, 40, 42))


  ## Check output 'cohortDT' ----

  expect_true(!is.null(simTest$cohortDT))
  expect_true(inherits(simTest$cohortDT, "data.table"))

  for (colName in c("cohortID", "pixelIndex")){
    expect_true(colName %in% names(simTest$cohortDT))
    expect_true(all(!is.na(simTest$cohortDT[[colName]])))
  }

  expect_identical(data.table::key(simTest$cohortDT), "cohortID")

  expect_equal(nrow(simTest$cohortDT), 160000)
  expect_equal(simTest$cohortDT$pixelIndex, 1:160000)
  expect_equal(simTest$cohortDT$cohortID, simTest$cohortDT$pixelIndex)


  ## Check output 'disturbanceMeta' ----

  expect_true(!is.null(simTest$disturbanceMeta))
  expect_true(inherits(simTest$disturbanceMeta, "data.table"))

  expect_equal(nrow(simTest$disturbanceMeta), 2)

  # Check that disturbances have been matched correctly
  expect_equal(simTest$disturbanceMeta$disturbance_type_id, c(1, 204))


  ## Check output 'disturbanceEvents' ----

  expect_true(!is.null(simTest$disturbanceEvents))
  expect_true(inherits(simTest$disturbanceEvents, "data.table"))

  for (colName in c("pixelIndex", "year", "eventID")){
    expect_true(colName %in% names(simTest$disturbanceEvents))
    expect_true(is.integer(simTest$disturbanceEvents[[colName]]))
    expect_true(all(!is.na(simTest$disturbanceEvents[[colName]])))
  }

  distEventCount <- simTest$disturbanceEvents[, .(N = .N), by = c("eventID")]
  expect_equal(distEventCount, rbind(
    data.table(eventID = 1001, N = 2138),
    data.table(eventID = 1002, N = 4681)
  ))
})


