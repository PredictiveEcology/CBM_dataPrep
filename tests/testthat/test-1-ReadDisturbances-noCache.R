
# ReadDisturbances() must not Cache() the extraction of rasters taken from the
# simList (a new raster every year, so the Cache never hits but is searched and written).

test_that("ReadDisturbances: simList-sourced rasters are not Cached; events unchanged", {

  # Load module functions without running simInit
  modEnv <- new.env(parent = globalenv())
  modFile <- file.path(dirname(dirname(getwd())), "CBM_dataPrep.R")
  if (!file.exists(modFile)) modFile <- file.path(getwd(), "CBM_dataPrep.R")
  exprs <- parse(modFile)
  isFn <- vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                   is.name(e[[2]]) && as.character(e[[2]]) %in% c("ReadDisturbances", "masterRasterDigest"), logical(1))
  for (e in exprs[isFn]) eval(e, modEnv)
  eval(parse(file.path(dirname(modFile), "R", "utils-cache.R")), modEnv)
  modEnv$P <- function(sim) list(saveRasters = FALSE)
  modEnv$time <- function(sim) 2000

  tmpl <- terra::rast(crs = "EPSG:3979", extent = c(xmin = 0, xmax = 300, ymin = 0, ymax = 300),
                      resolution = 30, vals = 1)
  clearcut <- terra::rast(tmpl, vals = c(1, 2, rep(NA, 98)))
  fileRast <- terra::rast(tmpl, vals = c(rep(NA, 50), 1, 1, rep(NA, 48)))

  makeSim <- function(){
    sim <- new.env()
    sim$masterRaster <- tmpl
    sim$disturbanceMeta <- data.table::data.table(
      eventID = 1:2, sourceValue = c(1L, NA), sourceDelay = c(1L, NA),
      sourceObjectName = c(NA, "clearcut"))
    sim$disturbanceRasters <- list(`1` = list(`2000` = fileRast))
    sim$disturbanceEvents <- NULL
    sim$clearcut <- clearcut
    sim
  }

  # Count Cache() calls by extraction input
  nCache <- 0L
  testthat::local_mocked_bindings(
    Cache = function(...) { nCache <<- nCache + 1L; ..1 },
    .package = "reproducible")

  sim <- modEnv$ReadDisturbances(makeSim(), year = 2000)

  # Only the file-based raster (event 1) goes through Cache; the simList raster (event 2) does not
  expect_equal(nCache, 1L)

  # Events equal those computed directly from extractToRast()
  vals2 <- CBMutils::extractToRast(clearcut, templateRast = tmpl)
  vals1 <- CBMutils::extractToRast(fileRast, templateRast = tmpl)
  expected <- data.table::rbindlist(list(
    data.table::data.table(pixelIndex = which(vals1 %in% 1), year = 2001L, eventID = 1L),
    data.table::data.table(pixelIndex = which(!is.na(vals2)), year = 2000L, eventID = 2L)))
  expect_equal(as.data.frame(sim$disturbanceEvents)[order(sim$disturbanceEvents$eventID, sim$disturbanceEvents$pixelIndex),],
               as.data.frame(expected)[order(expected$eventID, expected$pixelIndex),],
               ignore_attr = TRUE)
})
