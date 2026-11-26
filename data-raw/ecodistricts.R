
ecodistricts <- reproducible::prepInputs(
  destinationPath = tempdir(),
  url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecodistrict_shp.zip",
  targetFile = "ecodistricts.shp",
  alsoExtract = "similar",
  fun = sf::st_read(targetFile, quiet = TRUE)
)
ecodistricts <- ecodistricts[, c("ECODISTRIC", "ECOREGION")]
names(ecodistricts)[1:2] <- c("ecodist_id", "ecoreg_id")

zipDir <- file.path(tempdir(), "ecodistricts")
zipFile <- file.path(getwd(), "data/ecodistricts.zip")

unlink(zipDir, recursive = TRUE)
unlink(zipFile)
dir.create(zipDir)
sf::st_write(ecodistricts, file.path(zipDir, "ecodistricts.shp"))
withr::with_dir(zipDir, zip(zipFile, files = list.files(zipDir)))



