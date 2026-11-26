
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: vector inputs", {

  ## Run simInit and spades ----

  # Set up project
  simInitInput <- setupProjectTest(
    projectName = "vectorInputs",
    times       = list(start = 2025, end = 2025),

    # Set study area
    masterRaster = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = 2437830 - 50, xmax = 2437830 + 50, ymin = 126710 - 50, ymax = 126710 + 50),
      resolution = 10,
      vals       = 10
    ),

    # Set vector inputs
    ageLocator     = 10,
    ageDataYear    = 2025,
    cohortLocators = list(
      curveID = "GC-1"
    )
  )

  # Run simInit
  simTestInit <- SpaDES.core::simInit2(simInitInput)
  expect_s4_class(simTestInit, "simList")

  # Run spades
  simTest <- SpaDES.core::spades(simTestInit)
  expect_s4_class(simTest, "simList")


  ## Check output 'cohortDT' ----

  expect_true(!is.null(simTest$cohortDT))
  expect_true(inherits(simTest$cohortDT, "data.table"))

  for (colName in c("cohortID", "pixelIndex", "age", "curveID")){
    expect_true(colName %in% names(simTest$cohortDT))
    expect_true(all(!is.na(simTest$cohortDT[[colName]])))
  }
  expect_identical(data.table::key(simTest$cohortDT), "cohortID")

  expect_equal(nrow(simTest$cohortDT), 100)
  expect_equal(simTest$cohortDT$cohortID,   1:100)
  expect_equal(simTest$cohortDT$pixelIndex, 1:100)
  expect_in(simTest$cohortDT$age,           10)
  expect_in(simTest$cohortDT$curveID,       "GC-1")

})

