
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

tempDir <- file.path(spadesTestPaths$temp$outputs, "prepInputsExtractSpatial")
dir.create(tempDir, showWarnings = FALSE)

test_that("Function: prepInputsExtractSpatial: raster upsampling", {

  outPath <- tempfile("rast-upsample_", fileext = ".tif", tmpdir = tempDir)

  input <- terra::rast(file.path(spadesTestPaths$testdata, "SaskDist_1987_crop.tif"))

  masterRaster <- terra::rast(
    res = 5, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))

  inAlign <- prepInputsExtractSpatial(input, masterRaster)

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(27859, 3091, 12723, 316327)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: raster downsampling", {

  outPath <- tempfile("rast-downsample_", fileext = ".tif", tmpdir = tempDir)

  input <- terra::rast(file.path(spadesTestPaths$testdata, "SaskDist_1987_crop.tif"))

  masterRaster <- terra::rast(
    res = 50, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(280, 32, 121, 3167)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: raster reprojecting", {

  outPath <- tempfile("rast-reproject_", fileext = ".tif", tmpdir = tempDir)

  input <- terra::rast(file.path(spadesTestPaths$testdata, "tile1.tif"))

  masterRaster <- terra::rast(
    ncols = 213, nrows = 215,
    vals = 1, crs = "EPSG:4326",
    ext = c(xmin = -105.6567825386380974, xmax = -105.6294401406081107,
            ymin =   55.1008597705739831, ymax =   55.1284589047357017))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(16369, 29426)
    ), tolerance = 50, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: TIF file", {

  outPath <- tempfile("rast-TIF_", fileext = ".tif", tmpdir = tempDir)

  input <- file.path(spadesTestPaths$testdata, "tile1.tif")

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28, NaN),
      N   = c(48548, 76452, 125000)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: TIF tiles", {

  outPath <- tempfile("rast-TIF-tiles_", fileext = ".tif", tmpdir = tempDir)

  input <- file.path(spadesTestPaths$testdata, c("tile1.tif", "tile2.tif"))

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(108952, 141048)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: sf polygons with numeric field", {

  outPath <- tempfile("sf-numeric_", fileext = ".tif", tmpdir = tempDir)

  input <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)[, "id"]

  ## Create an NA area
  input <- subset(input, id != 8)

  masterRaster <- terra::rast(
    res = 100, vals = 1, crs = terra::crs(input),
    ext = round(terra::ext(input)))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_is(alignVals, "numeric")
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(1, 2, 3, 4, 5, 6, 7, NaN),
      N   = c(58876, 17794, 1939, 161340, 26220, 26934, 49538, 22779)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: sf polygons with numeric field: reproject", {

  outPath <- tempfile("sf-numeric-reproject_", fileext = ".tif", tmpdir = tempDir)

  input <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)[, "id"]

  ## Create an NA area
  input <- subset(input, id != 8)

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))
  expect_is(alignVals, "numeric")
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(1, 4, 5, NaN),
      N   = c(69058, 116667, 19457, 44818)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsExtractSpatial: sf polygons with text field: reproject", {

  outPath <- tempfile("sf-text-reproject_", fileext = ".tif", tmpdir = tempDir)

  input <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)[, "id"]
  input$id <- paste("Id", input$id)

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  alignVals <- prepInputsExtractSpatial(input, masterRaster, outPath = outPath, verbose = FALSE)
  alignRas  <- terra::rast(outPath)

  expect_true(terra::compareGeom(alignRas, masterRaster, stopOnError = FALSE))

  ## Make sure spaces and letter case are honored
  expect_is(alignVals, "character")
  expect_equal(
    data.table::data.table(val = alignVals)[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c("Id 1", "Id 4", "Id 5", "Id 8"),
      N   = c(69058, 116667, 19457, 44818)
    ), tolerance = 10, scale = 1)
})


