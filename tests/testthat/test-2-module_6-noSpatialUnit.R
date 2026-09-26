
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: pixels with no CBM spatial unit are excluded", {

  ## Run simInit and spades ----

  # Set up project
  projectName <- "noSpatialUnit"
  times       <- list(start = 2025, end = 2025)

  # Northwest Territories x eco_id 3 (Southern Arctic) has no CBM spatial unit.
  # Northwest Territories x eco_id 4 does have a CBM spatial unit.
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

    # Set study area: 1 masked (NA) pixel + 2 pixels in an invalid pair + 2 in a valid pair
    masterRaster = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = 0, xmax = 5, ymin = 0, ymax = 1),
      resolution = 1,
      vals       = c(NA, 1, 1, 1, 1)
    ),

    adminLocator = "Northwest Territories",
    ecoLocator   = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = 0, xmax = 5, ymin = 0, ymax = 1),
      resolution = 1,
      vals       = c(NA, 3, 3, 4, 4)
    ),
    ageLocator   = 10,
    ageDataYear  = 2025
  )

  # Run simInit
  simTestInit <- SpaDES.core::simInit2(simInitInput)
  expect_s4_class(simTestInit, "simList")

  # Run spades: expect the invalid (admin_name, eco_id) pair pixels to be
  # excluded, with a message explaining why
  expect_message(
    simTest <- SpaDES.core::spades(simTestInit),
    "no CBM spatial unit"
  )
  expect_s4_class(simTest, "simList")


  ## Check output 'standDT' ----

  expect_true(!is.null(simTest$standDT))
  expect_equal(nrow(simTest$standDT), 2)
  expect_equal(simTest$standDT$pixelIndex, c(4L, 5L))
  expect_in(simTest$standDT$eco_id, 4)


  ## Check output 'cohortDT' ----

  expect_true(!is.null(simTest$cohortDT))
  expect_equal(nrow(simTest$cohortDT), 2)
  expect_equal(simTest$cohortDT$pixelIndex, c(4L, 5L))

})
