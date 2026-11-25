
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: set masterRaster", {

  ## Expect error: missing inputs ----

  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-missing",
      times                = list(start = 2026, end = 2026)
    ) |> SpaDES.core::simInitAndSpades2()
  )

  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-missing",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "SCANFI",
    ) |> SpaDES.core::simInitAndSpades2()
  )

  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-missing",
      times                = list(start = 2026, end = 2026),
      studyAreaReporting   = sf::st_sf(sf::st_buffer(sf::st_sfc(sf::st_point(c(-63.57688, 44.64269)), crs = 4326), 1000)),
    ) |> SpaDES.core::simInitAndSpades2()
  )


  ## masterRasterTemplate as SpatRaster ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-template",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = terra::rast(
      crs        = "EPSG:3979",
      extent     = c(xmin = -687696, xmax = -681036, ymin = 711955, ymax = 716185),
      resolution = 30
    ),
    studyAreaReporting   = sf::st_sf(sf::st_buffer(sf::st_sfc(sf::st_point(c(-63.57688, 44.64269)), crs = 4326), 1000))
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs("EPSG:3979"))
  expect_equal(terra::res(simTest$masterRaster), rep(30, 2))
  expect_equal(terra::ncell(simTest$masterRaster), 4830)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 3688)

  expect_equal(nrow(simTest$standDT), 3688)
  expect_setequal(as.character(simTest$standDT$admin_name),   "Nova Scotia")
  expect_setequal(as.character(simTest$standDT$admin_abbrev), "NS")
  expect_setequal(simTest$standDT$eco_id,                     7)


  ## masterRasterTemplate as data source ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-dataSource",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    studyAreaReporting   = sf::st_sf(sf::st_buffer(sf::st_sfc(sf::st_point(c(-63.57688, 44.64269)), crs = 4326), 1000))
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), terra::res(masterRasterTemplate("SCANFI")))
  expect_equal(terra::ncell(simTest$masterRaster), 4761)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 3690)

  expect_equal(nrow(simTest$standDT), 3690)
  expect_setequal(as.character(simTest$standDT$admin_name),   "Nova Scotia")
  expect_setequal(as.character(simTest$standDT$admin_abbrev), "NS")
  expect_setequal(simTest$standDT$eco_id,                     7)

  # Expect error: invalid data source
  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-dataSource-invalid",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "invalidTemplate",
      studyAreaAdmin       = "SK"
    ) |> SpaDES.core::simInitAndSpades2(),
    "invalidTemplate"
  )


  ## masterRasterTemplate with data source and masterRasterRes ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-res",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 240,
    studyAreaReporting   = sf::st_sf(sf::st_buffer(sf::st_sfc(sf::st_point(c(-63.57688, 44.64269)), crs = 4326), 1000))
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), rep(240, 2))
  expect_equal(terra::ncell(simTest$masterRaster), 81)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 59)

  expect_equal(nrow(simTest$standDT), 59)
  expect_setequal(as.character(simTest$standDT$admin_name),   "Nova Scotia")
  expect_setequal(as.character(simTest$standDT$admin_abbrev), "NS")
  expect_setequal(simTest$standDT$eco_id,                     7)

  # Expect error: invalid res
  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-res-invalid",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "SCANFI",
      masterRasterRes      = 250,
      studyAreaReporting   = sf::st_sf(sf::st_buffer(sf::st_sfc(sf::st_point(c(-63.57688, 44.64269)), crs = 4326), 1000))
    ) |> SpaDES.core::simInitAndSpades2(),
    "resolution"
  )


  ## masterRasterTemplate with studyAreaAdmin ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-admin-name",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    studyAreaAdmin       = "Prince Edward Island"
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), terra::res(masterRasterTemplate("SCANFI")))
  expect_equal(terra::ncell(simTest$masterRaster), 20049120)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 9942278)

  expect_equal(nrow(simTest$standDT), 9942278)
  expect_setequal(as.character(simTest$standDT$admin_name),   "Prince Edward Island")
  expect_setequal(as.character(simTest$standDT$admin_abbrev), "PE")
  expect_setequal(simTest$standDT$eco_id,                     7)

  expect_is(simTest$studyAreaReporting, "sf")
  expect_setequal(simTest$studyAreaReporting$admin_name,   as.character(simTest$standDT$admin_name))
  expect_setequal(simTest$studyAreaReporting$admin_abbrev, as.character(simTest$standDT$admin_abbrev))
  expect_setequal(simTest$studyAreaReporting$eco_id,       simTest$standDT$eco_id)
  expect_equal(as.numeric(sf::st_area(simTest$studyAreaReporting)), 8948051787, tolerance = 100, scale = 1)

  # Try using abbreviation
  simTestEquiv <- setupProjectTest(
    projectName          = "masterRaster-admin-abbrev",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    studyAreaAdmin       = "PE"
  ) |> SpaDES.core::simInitAndSpades2()

  expect_true(terra::compareGeom(simTestEquiv$masterRaster, simTest$masterRaster))
  expect_equal(simTestEquiv$standDT, simTest$standDT)

  # Expect error: invalid admin
  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-admin-name-invalid",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "SCANFI",
      studyAreaAdmin       = "USA"
    ) |> SpaDES.core::simInitAndSpades2(),
    "studyAreaAdmin"
  )


  ## masterRasterTemplate with multiple studyAreaAdmin ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-admin-multi",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 3000,
    studyAreaAdmin       = c("Nova Scotia", "Prince Edward Island")
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), rep(3000, 2))
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 9403)

  expect_equal(nrow(simTest$standDT), 9403)
  expect_setequal(as.character(simTest$standDT$admin_name),   c("Nova Scotia", "Prince Edward Island"))
  expect_setequal(as.character(simTest$standDT$admin_abbrev), c("NS", "PE"))
  expect_setequal(unique(simTest$standDT$eco_id),              7)

  expect_is(simTest$studyAreaReporting, "sf")
  expect_setequal(simTest$studyAreaReporting$admin_name,   as.character(simTest$standDT$admin_name))
  expect_setequal(simTest$studyAreaReporting$admin_abbrev, as.character(simTest$standDT$admin_abbrev))
  expect_setequal(simTest$studyAreaReporting$eco_id,       simTest$standDT$eco_id)
  expect_equal(sum(as.numeric(sf::st_area(simTest$studyAreaReporting))), 84629730786, tolerance = 100, scale = 1)


  ## masterRasterTemplate with studyAreaEco ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-eco-id",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 3000,
    studyAreaEco         = 7,
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), rep(3000, 2))
  expect_equal(terra::ncell(simTest$masterRaster), 76953)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 26733)

  expect_equal(nrow(simTest$standDT), 26733)
  expect_setequal(as.character(simTest$standDT$admin_name),   c("Nova Scotia", "Prince Edward Island", "New Brunswick", "Quebec"))
  expect_setequal(as.character(simTest$standDT$admin_abbrev), c("NS", "PE", "NB", "QC"))
  expect_setequal(simTest$standDT$eco_id,                     7)

  expect_is(simTest$studyAreaReporting, "sf")
  expect_setequal(simTest$studyAreaReporting$admin_name,   as.character(simTest$standDT$admin_name))
  expect_setequal(simTest$studyAreaReporting$admin_abbrev, as.character(simTest$standDT$admin_abbrev))
  expect_setequal(simTest$studyAreaReporting$eco_id,       simTest$standDT$eco_id)

  # Try using ecozone name
  simTestEquiv <- setupProjectTest(
    projectName          = "masterRaster-eco-name",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 3000,
    studyAreaEco         = "Atlantic Maritime",
  ) |> SpaDES.core::simInitAndSpades2()

  expect_true(terra::compareGeom(simTestEquiv$masterRaster, simTest$masterRaster))
  expect_equal(simTestEquiv$standDT, simTest$standDT, check.attributes = FALSE)

  # Expect error: invalid ecozone
  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-eco-id",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "SCANFI",
      studyAreaEco         = "USA"
    ) |> SpaDES.core::simInitAndSpades2(),
    "studyAreaEco"
  )


  ## masterRasterTemplate with multiple studyAreaEco ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-eco-multi",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 3000,
    studyAreaEco         = c(13, 14),
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), rep(3000, 2))
  expect_equal(terra::ncell(simTest$masterRaster), 219672)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 80330)

  expect_equal(nrow(simTest$standDT), 80330)
  expect_setequal(as.character(simTest$standDT$admin_name),   c("British Columbia", "Alberta", "Yukon"))
  expect_setequal(as.character(simTest$standDT$admin_abbrev), c("BC", "AB", "YT"))
  expect_setequal(simTest$standDT$eco_id,                     c(13, 14))

  expect_is(simTest$studyAreaReporting, "sf")
  expect_setequal(simTest$studyAreaReporting$admin_name,   as.character(simTest$standDT$admin_name))
  expect_setequal(simTest$studyAreaReporting$admin_abbrev, as.character(simTest$standDT$admin_abbrev))
  expect_setequal(simTest$studyAreaReporting$eco_id,       simTest$standDT$eco_id)


  ## masterRasterTemplate with studyAreaAdmin and studyAreaEco ----

  simTest <- setupProjectTest(
    projectName          = "masterRaster-admin-eco",
    times                = list(start = 2026, end = 2026),
    masterRasterTemplate = "SCANFI",
    masterRasterRes      = 3000,
    studyAreaAdmin       = "Quebec",
    studyAreaEco         = 7,
  ) |> SpaDES.core::simInitAndSpades2()

  expect_s4_class(simTest, "simList")

  expect_is(simTest$masterRaster, "SpatRaster")
  expect_equal(terra::crs(simTest$masterRaster), terra::crs(masterRasterTemplate("SCANFI")))
  expect_equal(terra::res(simTest$masterRaster), rep(3000, 2))
  expect_equal(terra::ncell(simTest$masterRaster), 52210)
  expect_equal(terra::global(simTest$masterRaster, "notNA")[1,1], 8435)

  expect_equal(nrow(simTest$standDT), 8435)
  expect_setequal(as.character(simTest$standDT$admin_name),   c("Quebec"))
  expect_setequal(as.character(simTest$standDT$admin_abbrev), c("QC"))
  expect_setequal(simTest$standDT$eco_id,                     7)

  expect_is(simTest$studyAreaReporting, "sf")
  expect_setequal(simTest$studyAreaReporting$admin_name,   as.character(simTest$standDT$admin_name))
  expect_setequal(simTest$studyAreaReporting$admin_abbrev, as.character(simTest$standDT$admin_abbrev))
  expect_setequal(simTest$studyAreaReporting$eco_id,       simTest$standDT$eco_id)

  # Expect error: no admin and ecozone overlap
  expect_error(
    setupProjectTest(
      projectName          = "masterRaster-admin-eco-invalid",
      times                = list(start = 2026, end = 2026),
      masterRasterTemplate = "SCANFI",
      studyAreaAdmin       = "SK",
      studyAreaEco         = 7
    ) |> SpaDES.core::simInitAndSpades2(),
    "overlap"
  )
})

