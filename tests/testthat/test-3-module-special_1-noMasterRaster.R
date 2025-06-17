
if (!testthat::is_testing()) source(testthat::test_path("setup.R"))

test_that("Module: masterRaster missing", {

  # Set project path
  projectPath <- file.path(spadesTestPaths$temp$projects, "3-module-special_1-noMasterRaster")
  dir.create(projectPath)
  withr::local_dir(projectPath)

  # Set up project
  simInitInput <- SpaDEStestMuffleOutput(

    SpaDES.project::setupProject(
      modules = "CBM_dataPrep",
      paths   = list(
        projectPath = projectPath,
        modulePath  = spadesTestPaths$modulePath,
        packagePath = spadesTestPaths$packagePath,
        inputPath   = spadesTestPaths$inputPath,
        cachePath   = spadesTestPaths$cachePath,
        outputPath  = file.path(projectPath, "outputs")
      )
    )
  )

  # Run simInit
  simTestInit <- SpaDEStestMuffleOutput(
    SpaDES.core::simInit2(simInitInput)
  )

  expect_s4_class(simTestInit, "simList")

  # Run spades: expect error due to master raster missing
  expect_error(
    SpaDEStestMuffleOutput(
      SpaDES.core::spades(simTestInit)
    ))
})


