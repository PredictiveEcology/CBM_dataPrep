
admins <- reproducible::prepInputs(
  destinationPath = "C:/Users/sumurray/Documents/Github/spadesCBM/inputs/StatCan",
  url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
  targetFile = "lpr_000a21a_e.shp",
  alsoExtract = "similar",
  fun = sf::st_read(targetFile, quiet = TRUE)
)
admins <- admins[, c("PRENAME")]
names(admins)[1] <- c("name")

# Split Newfoundland and Labrador
admins <- terra::split(
  terra::vect(admins),
  terra::vect(sf::st_sfc(sf::st_linestring(rbind(c(8486903, 2323176), c(8537414, 2449589))),
                         crs = sf::st_crs(admins)))) |>
  sf::st_as_sf(agr = "constant")

nl_cd <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(
  admins[admins$name == "Newfoundland and Labrador",])))
admins[admins$name == "Newfoundland and Labrador", "name"] <- sapply(
  nl_cd[, "X"] == min(nl_cd[, "X"]), ifelse, "Labrador", "Newfoundland")

zipDir <- file.path(tempdir(), "adminLocator")
zipFile <- file.path(getwd(), "data/adminLocator.zip")
unlink(zipDir, recursive = TRUE)
unlink(zipFile)
dir.create(zipDir)
sf::st_write(admins, file.path(zipDir, "adminLocator.shp"))
withr::with_dir(zipDir, zip(zipFile, files = list.files(zipDir)))

