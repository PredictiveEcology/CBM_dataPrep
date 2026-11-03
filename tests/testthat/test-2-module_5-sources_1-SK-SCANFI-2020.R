
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: SK with SCANFI 2020 data", {

  ## Run simInit and spades ----

  ## Test only to be run manually
  ## Source data is too large to download for a test
  testthat::skip()

  # Set up project
  simInitInput <- setupProjectTest(
    projectName = "SK-SCANFI-2020",
    times       = list(start = 2019, end = 2020), # Check age backtracking

    # Set study area
    masterRaster = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716183),
      resolution = 30,
      vals       = 1
    ),

    # Set input data sources
    CBMsourceIDs = c("SCANFI-2020-age", "SCANFI-2020-LandR")
  )

  # Run simInit
  simTestInit <- SpaDES.core::simInit2(simInitInput)
  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDES.core::spades(simTestInit)
  expect_s4_class(simTest, "simList")


  ## Check outputs ----

  # curveID
  expect_equal(simTest$curveID, c("LandR"))

  # ageDataYear
  expect_equal(simTest$ageDataYear, 2020)

  # cohortDT
  expect_true(!is.null(simTest$cohortDT))
  expect_true(inherits(simTest$cohortDT, "data.table"))

  for (colName in c("cohortID", "pixelIndex", "age", "LandR")){
    expect_true(colName %in% names(simTest$cohortDT))
    expect_true(all(!is.na(simTest$cohortDT[[colName]])))
  }
  expect_identical(data.table::key(simTest$cohortDT), "cohortID")

  expect_true(is.integer(simTest$cohortDT$age))
  expect_true("Pice_mar" %in% simTest$cohortDT$LandR)

})


