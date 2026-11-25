
ecozones <- reproducible::prepInputs(
  destinationPath = tempdir(),
  url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
  targetFile = "ecozones.shp",
  alsoExtract = "similar",
  fun = sf::st_read(targetFile, quiet = TRUE)
)
ecozones <- ecozones[, c("ECOZONE", "ZONE_NAME")]
names(ecozones)[1:2] <- c("id", "name")

zipDir <- file.path(tempdir(), "ecoLocator")
zipFile <- file.path(getwd(), "data/ecoLocator.zip")

unlink(zipDir, recursive = TRUE)
unlink(zipFile)
dir.create(zipDir)
sf::st_write(ecozones, file.path(zipDir, "ecoLocator.shp"))
withr::with_dir(zipDir, zip(zipFile, files = list.files(zipDir)))



