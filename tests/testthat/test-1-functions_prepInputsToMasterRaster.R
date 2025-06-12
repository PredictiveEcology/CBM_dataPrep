
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Function: prepInputsToMasterRaster", {

  withr::local_libpaths(spadesTestPaths$packagePath, action = "prefix")
  source(list.files(file.path(spadesTestPaths$RProj, "R"), pattern = "\\.R$", full = TRUE))

  # Prep SpatRaster: test upsampling
  input <- terra::rast(file.path(spadesTestPaths$testdata, "SaskDist_1987_crop.tif"))

  masterRaster <- terra::rast(
    res = 5, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))
  prepRast <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(
    prepRast, tempfile("prepRast-upsample_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepRast, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepRast)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(27859, 3091, 12723, 316327)
    ), tolerance = 10, scale = 1)

  # Prep SpatRaster: test downsampling
  masterRaster <- terra::rast(
    res = 50, vals = 1, crs = "EPSG:3979",
    ext = c(xmin = -674500, xmax = -671500, ymin =  702000, ymax =  705000))
  prepRast <- prepInputsToMasterRaster(input, masterRaster)

  if (interactive()) terra::writeRaster(
    prepRast, tempfile("prepRast-downsample_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepRast, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepRast)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(2, 3, 5, NA),
      N   = c(280, 32, 121, 3167)
    ), tolerance = 10, scale = 1)

  # Prep SpatRaster: test reprojecting
  masterRaster <- terra::rast(
    ncols = 213, nrows = 215,
    vals = 1, crs = "EPSG:4326",
    ext = c(xmin = -105.6567825386380974, xmax = -105.6294401406081107,
            ymin =   55.1008597705739831, ymax =   55.1284589047357017))

  prepRast <- prepInputsToMasterRaster(
    input = terra::rast(file.path(spadesTestPaths$testdata, "tile1.tif")),
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepRast, tempfile("prepRast-reproject_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepRast, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepRast)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(16404, 29391)
    ), tolerance = 50, scale = 1)

  # Prep raster file
  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))
  prepFile <- prepInputsToMasterRaster(
    input = file.path(spadesTestPaths$testdata, "tile1.tif"),
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepFile, tempfile("prepFile_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepFile, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepFile)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28, NaN),
      N   = c(48548, 76452, 125000)
    ), tolerance = 10, scale = 1)

  # Prep raster tiles
  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  458500, xmax =  463500, ymin = 6105000, ymax = 6110000))
  prepTiles <- prepInputsToMasterRaster(
    input = file.path(spadesTestPaths$testdata, c("tile1.tif", "tile2.tif")),
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepTiles, tempfile("prepTiles_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepTiles, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepTiles)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(27, 28),
      N   = c(108952, 141048)
    ), tolerance = 10, scale = 1)

  # Prep sf polygons
  inSF <- sf::st_read(
    file.path(spadesTestPaths$testdata, "spuLocator.shp"), agr = "constant",
    quiet = TRUE)
  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  prepSF <- prepInputsToMasterRaster(
    input = inSF[, "id"],
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepSF, tempfile("prepSF_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepSF, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepSF)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(1, 4, 5, 8),
      N   = c(69052, 116674, 19458, 44816)
    ), tolerance = 10, scale = 1)

  # Prep sf polygons with text field
  prepSF_text <- prepInputsToMasterRaster(
    input = cbind(text = as.character(inSF[["id"]]), inSF),
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepSF_text, tempfile("prepSF-text_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepSF_text, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepSF_text)[, 1])[, .N, by = "val"][order(val)]$N,
    data.table::data.table(val = terra::values(prepSF)[, 1])[, .N, by = "val"][order(val)]$N
  )
  expect_setequal(terra::cats(prepSF_text)[[1]][[2]], c("1", "4", "5", "8"))

  ## Allow for NA areas
  masterRaster <- terra::rast(
    res = 10, vals = 1, crs = "EPSG:32613",
    ext = c(xmin =  456000,  xmax = 461000, ymin = 6105000, ymax = 6110000))

  prepSF_NAs <- prepInputsToMasterRaster(
    input = subset(inSF[, "id"], id == 1),
    masterRaster = masterRaster)

  if (interactive()) terra::writeRaster(
    prepSF_NAs, tempfile("prepSF-NAs_", fileext = ".tif", tmpdir = spadesTestPaths$temp$outputs))

  expect_true(terra::compareGeom(prepSF_NAs, masterRaster, stopOnError = FALSE))
  expect_equal(
    data.table::data.table(val = terra::values(prepSF_NAs)[, 1])[, .N, by = "val"][order(val)],
    data.table::data.table(
      val = c(1, NaN),
      N   = c(69052, 250000 - 69052)
    ), tolerance = 10, scale = 1)
})


