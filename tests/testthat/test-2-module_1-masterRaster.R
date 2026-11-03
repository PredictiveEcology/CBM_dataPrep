
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: masterRaster missing", {

  # Expect error due to master raster missing
  expect_error(
    setupProjectTest(
      projectName = "masterRaster-missing",
      times       = list(start = 2026, end = 2026)
    ) |> SpaDES.core::simInitAndSpades2(),
    "masterRaster"
  )
})


