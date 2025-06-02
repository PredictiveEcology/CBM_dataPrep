
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: SK-small", {

  ## Run simInit and spades ----

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "2-module_1-SK-small")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(

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
        extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
        resolution = 30,
        vals       = 1
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

  for (colName in c("pixelIndex", "area", "spatial_unit_id")){
    expect_true(colName %in% names(simTest$standDT))
    expect_true(all(!is.na(simTest$standDT[[colName]])))
  }

  expect_identical(data.table::key(simTest$standDT), "pixelIndex")

  expect_equal(nrow(simTest$standDT), 31302)
  expect_equal(simTest$standDT$pixelIndex,      1:31302)
  expect_equal(simTest$standDT$area,            rep(900, 31302))
  expect_equal(simTest$standDT$spatial_unit_id, rep(28, 31302))


  ## Check output 'cohortDT' ----

  expect_true(!is.null(simTest$cohortDT))
  expect_true(inherits(simTest$cohortDT, "data.table"))

  for (colName in c("cohortID", "pixelIndex")){
    expect_true(colName %in% names(simTest$cohortDT))
    expect_true(all(!is.na(simTest$cohortDT[[colName]])))
  }

  expect_identical(data.table::key(simTest$cohortDT), "cohortID")

  expect_equal(nrow(simTest$cohortDT), 31302)
  expect_equal(simTest$cohortDT$pixelIndex, 1:31302)
  expect_equal(simTest$cohortDT$cohortID, simTest$cohortDT$pixelIndex)

})

