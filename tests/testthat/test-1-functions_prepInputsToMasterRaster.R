
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))
if (interactive()){
  tempDir <- file.path(spadesTestPaths$temp$outputs, "prepInputsToMasterRaster")
  dir.create(tempDir)
}

test_that("Function: prepInputsToMasterRaster: raster upsampling", {

  input <- terra::rast(file.path(spadesTestPaths$testdata, "SaskDist_1987_crop.tif"))

  masterRaster <- terra::rast(
    res = 5, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("rast-upsample_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(27859, 3091, 12723, 316327)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: raster downsampling", {

  input <- terra::rast(file.path(spadesTestPaths$testdata, "SaskDist_1987_crop.tif"))

  masterRaster <- terra::rast(
    res = 50, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("rast-downsample_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(280, 32, 121, 3167)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: raster reprojecting", {

  input <- terra::rast(file.path(spadesTestPaths$testdata, "tile1.tif"))

  masterRaster <- terra::rast(
    ncols = 213, nrows = 215,
    vals = 1, crs = "EPSG:4326",
    ext = c(xmin = -105.6567825386380974, xmax = -105.6294401406081107,
            ymin =   55.1008597705739831, ymax =   55.1284589047357017))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("rast-reproject_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(16369, 29426)
    ), tolerance = 50, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: TIF file", {

  input <- file.path(spadesTestPaths$testdata, "tile1.tif")

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("rast-TIF_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28, NaN),
      N   = c(48548, 76452, 125000)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: TIF tiles", {

  input <- file.path(spadesTestPaths$testdata, c("tile1.tif", "tile2.tif"))

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("rast-TIF-tiles_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(108952, 141048)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: sf polygons with numeric field", {

  input <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)[, "id"]

  ## Create an NA area
  input <- subset(input, id != 8)

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("sf-numeric_", fileext = ".tif", tmpdir = tempDir))

  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(1, 4, 5, NaN),
      N   = c(69058, 116667, 19457, 44818)
    ), tolerance = 10, scale = 1)
})

test_that("Function: prepInputsToMasterRaster: sf polygons with text field", {

  input <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)[, "id"]
  input$id <- paste("Id", input$id)

  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  inAlign <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(inAlign, tempfile("sf-text_", fileext = ".tif", tmpdir = tempDir))

  ## Make sure spaces and letter case are honored
  expect_true(terra::compareGeom(inAlign, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)]$N,
    data.table::data.table(val = terra::values(inAlign)[, 1])[, .N, by = "val"][order(val)]$N
  )
  expect_setequal(terra::cats(inAlign)[[1]][[2]], c("Id 1", "Id 4", "Id 5", "Id 8"))

})




