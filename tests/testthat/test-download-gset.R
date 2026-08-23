test_that("download_gset validates data_id", {
  valid_dir <- tempdir()

  expect_error(download_gset(NULL, valid_dir), "data_id")
  expect_error(download_gset(character(), valid_dir), "data_id")
  expect_error(download_gset(c("GSE1009", "GPL96"), valid_dir), "data_id")
  expect_error(download_gset(NA_character_, valid_dir), "data_id")
  expect_error(download_gset(1009, valid_dir), "data_id")
  expect_error(download_gset("GSE", valid_dir), "data_id")
  expect_error(download_gset("gse1009", valid_dir), "data_id")
  expect_error(download_gset("ABC1009", valid_dir), "data_id")
})

test_that("download_gset validates destdir", {
  missing_dir <- file.path(tempdir(), "bioBulkR-directory-does-not-exist")

  expect_error(download_gset("GSE1009", NULL), "destdir")
  expect_error(download_gset("GSE1009", character()), "destdir")
  expect_error(download_gset("GSE1009", c(tempdir(), tempdir())), "destdir")
  expect_error(download_gset("GSE1009", NA_character_), "destdir")
  expect_error(download_gset("GSE1009", ""), "destdir")
  expect_error(download_gset("GSE1009", missing_dir), "existing directory")
})

test_that("download_gset delegates to GEOquery and returns its result", {
  calls <- list()
  expected <- structure(list("mock expression set"), class = "mock_gset")

  local_mocked_bindings(
    getGEO = function(GEO, destdir, AnnotGPL, getGPL) {
      calls <<- list(
        GEO = GEO,
        destdir = destdir,
        AnnotGPL = AnnotGPL,
        getGPL = getGPL,
        timeout = getOption("timeout")
      )
      expected
    },
    .package = "GEOquery"
  )

  expect_message(
    result <- download_gset("GSE1009", tempdir()),
    "Download successful"
  )

  expect_identical(result, expected)
  expect_identical(calls$GEO, "GSE1009")
  expect_identical(calls$destdir, tempdir())
  expect_false(calls$AnnotGPL)
  expect_false(calls$getGPL)
  expect_identical(calls$timeout, 120)
})

test_that("download_gset restores the user's timeout option", {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 37)

  local_mocked_bindings(
    getGEO = function(...) list(mock = TRUE),
    .package = "GEOquery"
  )

  expect_message(download_gset("GPL96", tempdir()), "Download successful")
  expect_identical(getOption("timeout"), 37)
})
