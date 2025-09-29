
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: vector inputs", {

  ## Run simInit and spades ----

  # Set up project
  projectName <- "vectorInputs"
  times       <- list(start = 2025, end = 2025)

  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

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
        crs        = "EPSG:32610",
        extent     = c(xmin = 475595 - 50, xmax = 475595 + 50, ymin = 5361627 - 50, ymax = 5361627 + 50),
        resolution = 1,
        vals       = 1
      ),

      adminLocator   = "Nova Scotia",
      ecoLocator     = 7,
      ageLocator     = 10,
      ageDataYear    = 2025,
      cohortLocators = list(
        curveID = "GC-1"
      )
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

  expect_equal(nrow(simTest$standDT), 10000)
  expect_equal(simTest$standDT$pixelIndex, 1:10000)
  expect_in(simTest$standDT$area,              1)
  #expect_in(simTest$standDT$admin_name,        "Nova Scotia") # Column excluded from result
  expect_in(simTest$standDT$admin_abbrev,      "NS")
  expect_in(simTest$standDT$admin_boundary_id, 3)
  expect_in(simTest$standDT$ecozone,           7)
  expect_in(simTest$standDT$spatial_unit_id,   5)


  ## Check output 'cohortDT' ----

  expect_true(!is.null(simTest$cohortDT))
  expect_true(inherits(simTest$cohortDT, "data.table"))

  for (colName in c("cohortID", "pixelIndex", "age", "curveID")){
    expect_true(colName %in% names(simTest$cohortDT))
    expect_true(all(!is.na(simTest$cohortDT[[colName]])))
  }
  expect_identical(data.table::key(simTest$cohortDT), "cohortID")

  expect_equal(nrow(simTest$cohortDT), 10000)
  expect_equal(simTest$cohortDT$cohortID,   1:10000)
  expect_equal(simTest$cohortDT$pixelIndex, 1:10000)
  expect_in(simTest$cohortDT$age,           10)
  expect_in(simTest$cohortDT$curveID,       "GC-1")

})


