SPECTRUM_UTILS_TOLERANCE <- 1e-7
df <- DataFrame(msLevel = c(1L, 2L, 1L), rtime = c(1.2, 3.2, 5.2))
df$mz <- list(
    c(34.3, 45.6, 199.1),
    c(45.6, 89.3, 100.1, 200.5),
    c(56.1, 56.2)
)
df$intensity <- list(
    1:3,
    1:4,
    1:2
)
sps <- Spectra(df)

test_that("spectraVariableMapping, spectraVariableMapping<- works", {
    res <- spectraVariableMapping()
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_true(length(names(res)) > 0)
    expect_equal(res, .SPECTRA_2_MATCHMS)

    setSpectraVariableMapping(c(a = "b", d = "e"))
    expect_equal(spectraVariableMapping(), c(a = "b", d = "e"))
    setSpectraVariableMapping(.SPECTRA_2_MATCHMS)
    expect_equal(spectraVariableMapping(),
                 .SPECTRA_2_MATCHMS)
    expect_equal(defaultSpectraVariableMapping(), .SPECTRA_2_MATCHMS)

    res <- spectraVariableMapping("spectrum_utils")
    expect_equal(res, .SPECTRA_2_SPECTRUM_UTILS)

    expect_error(spectraVariableMapping("other"), "Supported values")

    res <- spectraVariableMapping("matchms")
    expect_equal(res, .SPECTRA_2_MATCHMS)

    expect_error(setSpectraVariableMapping(3), "named character vector")
    expect_error(setSpectraVariableMapping(c("a", "b")), "named character")

    res <- spectraVariableMapping("matchms", c("some", "other"))
    expect_equal(res, c(spectraVariableMapping("matchms"),
                        c(some = "some", other = "other")))

    expect_warning(
        res <- spectraVariableMapping("spectrum_utils", c("some", "other")),
        "Provided spectra variables")
    expect_equal(res, spectraVariableMapping("spectrum_utils"))
})

#############
##    R to Py

test_that(".single_rspec_to_pyspec works", {
    res <- .single_rspec_to_pyspec(sps[1L])
    expect_true(is(res, "matchms.Spectrum.Spectrum"))
    expect_true(all(.SPECTRA_2_MATCHMS %in% names(res$metadata)))
    expect_equal(sps$mz[[1L]], as.vector(py_to_r(res$mz)))
    expect_equal(sps$intensity[[1L]], as.vector(py_to_r(res$intensities)))
    expect_equal(sps$rtime[1L], py_to_r(res$metadata["retention_time"]))
    expect_equal(sps$msLevel[1L], py_to_r(res$metadata["ms_level"]))
    expect_equal(sps$precursorMz[1L], py_to_r(res$metadata["precursor_mz"]))

    res <- .single_rspec_to_pyspec(sps[2L], character())
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res$intensities)))
    expect_equal(names(res$metadata), character())
})

test_that(".rspec_to_matchms_pyspec works", {
    res <- .rspec_to_matchms_pyspec(sps)
    expect_true(is(res, "python.builtin.list"))
    expect_true(length(res) == length(sps))

    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_true(all(.SPECTRA_2_MATCHMS %in% names(res[[1]]$metadata)))

    res <- .rspec_to_matchms_pyspec(sps, character())
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_equal(names(res[[1]]$metadata), character())

    ## passing a `DataFrame` or `data.frame` instead
    x <- spectraData(sps@backend)
    res <- .rspec_to_matchms_pyspec(x)
    expect_true(is(res, "python.builtin.list"))
    expect_true(length(res) == length(sps))
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_true(all(.SPECTRA_2_MATCHMS %in% names(res[[1]]$metadata)))

    res <- .rspec_to_matchms_pyspec(x, character())
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_equal(names(res[[1]]$metadata), character())
})

test_that(".rspec_to_spectrum_utils_pyspec works", {
    res <- expect_warning(.rspec_to_spectrum_utils_pyspec(sps),
                          "Ignoring variables")
    expect_equal(length(res), length(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_true(is(res[[1L]], "spectrum_utils.spectrum.MsmsSpectrum"))
    expect_true(is.na(py_to_r(res[[1L]]$precursor_mz)))
    expect_equal(py_to_r(res[[1L]]$retention_time), sps$rtime[2L])

    res <- .rspec_to_spectrum_utils_pyspec(sps, c())
    expect_equal(length(res), length(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_true(is(res[[1L]], "spectrum_utils.spectrum.MsmsSpectrum"))
    expect_true(is.na(py_to_r(res[[1L]]$precursor_mz)))
    expect_true(is.na(py_to_r(res[[1L]]$retention_time)))

    sps$my_id <- c("a", "b", "c")
    res <- .rspec_to_spectrum_utils_pyspec(sps, c(my_id = "identifier"))
    expect_equal(length(res), length(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_true(is(res[[1L]], "spectrum_utils.spectrum.MsmsSpectrum"))
    expect_true(is.na(py_to_r(res[[1L]]$precursor_mz)))
    expect_true(is.na(py_to_r(res[[1L]]$retention_time)))
    expect_equal(py_to_r(res[[1L]]$identifier), "b")

    res <- .rspec_to_spectrum_utils_pyspec(sps, .SPECTRA_2_SPECTRUM_UTILS)
    expect_equal(length(res), length(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_true(is(res[[1L]], "spectrum_utils.spectrum.MsmsSpectrum"))
    expect_true(is.na(py_to_r(res[[1L]]$precursor_mz)))

    ## passing a `DataFrame` or `data.frame` instead
    x <- spectraData(sps@backend)
    res <- .rspec_to_spectrum_utils_pyspec(x)
    expect_true(is(res, "python.builtin.list"))
    expect_true(length(res) == length(sps))
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensity)))
    expect_equal(sps$rtime[2], py_to_r(res[[1]]$retention_time))
    expect_equal(sps$precursorMz[2], py_to_r(res[[1]]$precursor_mz))

    res <- .rspec_to_spectrum_utils_pyspec(x, character())
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensity)))
    expect_true(is.na(py_to_r(res[[1]]$retention_time)))
    expect_true(is.na(py_to_r(res[[1]]$precursor_mz)))
})

test_that("rspec_to_pyspec works", {
    expect_error(rspec_to_pyspec(4), "Spectra object")

    res <- r_to_py(sps)
    expect_equal(res, rspec_to_pyspec(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_equal(length(res), length(sps))
    expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(sps)[[3]])
    expect_equal(rtime(sps)[1], py_to_r(res[0]$metadata$retention_time))
    expect_equal(msLevel(sps)[1], py_to_r(res[0]$metadata$ms_level))
    expect_equal(msLevel(sps)[2], py_to_r(res[1]$metadata$ms_level))
    expect_equal(msLevel(sps)[3], py_to_r(res[2]$metadata$ms_level))
    expect_equal(precursorMz(sps)[1], py_to_r(res[0]$metadata$precursor_mz))

    ## mapping only specific spectra variables.
    sps$new_col <- c(9, 3, 1)
    setSpectraVariableMapping(c(collisionEnergy = "collision_energy",
                                new_col = "new_col"))
    res <- r_to_py(sps)
    expect_true(is(res, "python.builtin.list"))
    expect_equal(length(res), length(sps))
    expect_equal(sort(names(res[0]$metadata)),
                 sort(c("collision_energy", "new_col")))
    expect_equal(py_to_r(res[0]$metadata$new_col), 9)

    s <- sps
    s$rtime <- NULL
    setSpectraVariableMapping(c(precursorMz = "precursor_mz",
                                rtime = "retention_time"))
    res <- r_to_py(s)
    expect_equal(sort(names(res[0]$metadata)),
                 c("precursor_mz", "retention_time"))
    expect_equal(py_to_r(res[0]$metadata$precursor_mz), NA_real_)
    ## It's a bit odd - a retention time of NA is converted to a NULL
    expect_equal(py_to_r(res[0]$metadata$retention_time), NULL)

    s$new_col <- NULL
    setSpectraVariableMapping(c(precursorMz = "precursor_mz",
                                new_col = "new_col"))
    expect_error(res <- r_to_py(s), "requested spectra variables")

    ## No mapping at all.
    setSpectraVariableMapping(character())
    res <- rspec_to_pyspec(s)
    expect_equal(names(res[0]$metadata), character())

    setSpectraVariableMapping(.SPECTRA_2_MATCHMS)

    ## spectrum_utils
    expect_error(rspec_to_pyspec(sps, pythonLibrary = "other"))
    res <- rspec_to_pyspec(sps, spectraVariableMapping("spectrum_utils"),
                           "spectrum_utils")
    expect_equal(length(res), length(sps))
    expect_true(is(res, "python.builtin.list"))
    expect_true(is(res[[1L]], "spectrum_utils.spectrum.MsmsSpectrum"))
})

#############
##    Py to R

test_that(".py_matchms_spectrum_spectra_data works", {
    p <- r_to_py(sps)
    res <- .py_matchms_spectrum_spectra_data(p[[1]])
    expect_true(is.data.frame(res))
    expect_true(all(names(.SPECTRA_2_MATCHMS) %in% colnames(res)))
    expect_equal(res$msLevel, sps$msLevel[2])
    expect_equal(res$rtime, sps$rtime[2])
    expect_equal(res$precursorCharge, sps$precursorCharge[2])
    expect_equal(res$precursorMz, sps$precursorMz[2],
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_equal(res$precursorIntensity, sps$precursorIntensity[2])
    expect_equal(res$collisionEnergy, sps$collisionEnergy[2])

    res <- .py_matchms_spectrum_spectra_data(
        p[[1]], mapping = c(msLevel = "ms_level", other_col = "other_col"))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, sps$msLevel[2])

    res <- .py_matchms_spectrum_spectra_data(p[[1]], mapping = character())
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, NA_integer_)
})

test_that(".py_matchms_spectrum_peaks_data works", {
    p <- r_to_py(sps)
    res <- .py_matchms_spectrum_peaks_data(p[[1]])
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("mz", "intensity"))
    expect_equal(res, peaksData(sps)[[2L]])
})

test_that(".py_matchms_spectrum_peaks_data_columns works", {
    p <- r_to_py(sps)
    res <- .py_matchms_spectrum_peaks_data_columns(p[[1]])
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("mz", "intensity"))
    expect_equal(res, peaksData(sps)[[2L]])
    res <- .py_matchms_spectrum_peaks_data_columns(
                           p[[1]], columns = c("intensity", "mz"))
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("intensity", "mz"))
    res <- .py_matchms_spectrum_peaks_data_columns(
                           p[[1]], columns = c("intensity"))
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("intensity"))

    res <- .py_matchms_spectrum_peaks_data_columns(
                           p[[1]], columns = c("intensity"), drop = TRUE)
    expect_true(is.numeric(res))
    expect_false(is.matrix(res))
})

test_that("pyspec_to_rspec and .single_pyspec_to_rspec work", {
    p <- r_to_py(sps)
    res <- pyspec_to_rspec(p[0])
    expect_s4_class(res, "Spectra")
    expect_equal(peaksData(res), peaksData(sps[1L]))
    expect_equal(res$msLevel, sps[1L]$msLevel)
    expect_equal(res$rtime, sps[1L]$rtime)

    res <- pyspec_to_rspec(p)
    expect_s4_class(res, "Spectra")
    expect_equal(peaksData(res), peaksData(sps))
    expect_equal(res$msLevel, sps$msLevel)
    expect_equal(res$rtime, sps$rtime)

    res2 <- pyspec_to_rspec(py_to_r(p))
    expect_s4_class(res2, "Spectra")
    expect_equal(peaksData(res2), peaksData(sps))
    expect_equal(res2$msLevel, sps$msLevel)
    expect_equal(res2$rtime, sps$rtime)

    res <- .single_pyspec_to_rspec(
        p[1L], mapping = c(rtime = "retention_time"))
    expect_equal(res$rtime, sps[2L]$rtime)
    expect_true(is.na(res$msLevel))

    res <- .single_pyspec_to_rspec(p[1L], mapping = character())
    expect_true(is.na(res$rtime))
    expect_true(is.na(res$msLevel))

    ## spectrum_utils
    map <- spectraVariableMapping("spectrum_utils")
    p <- rspec_to_pyspec(sps, map, "spectrum_utils")
    res <- expect_warning(pyspec_to_rspec(p, map, "spectrum_utils"), "NAs")
    expect_s4_class(res, "Spectra")
    expect_s4_class(res@backend, "MsBackendMemory")
    expect_equal(rtime(res), rtime(sps), tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_true(validObject(res))
    expect_equal(peaksData(res, return.type = "list"),
                 peaksData(sps, return.type = "list"),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
})

test_that("pyspec_to_rspec works", {
    p <- r_to_py(sps)
    res <- pyspec_to_rspec(p)
    expect_s4_class(res, "Spectra")
    expect_equal(peaksData(res), peaksData(sps))
    expect_equal(spectraData(res, names(spectraVariableMapping())),
                 spectraData(sps, names(spectraVariableMapping())))

    res_1 <- pyspec_to_rspec(p[0])
    expect_equal(res_1, res[1])

    ## custom mapping.
    sps$new_col <- seq_along(sps)
    maps <- c(rtime = "retention_time", new_col = "r_new_col")
    p <- rspec_to_pyspec(sps, mapping = maps)
    res <- pyspec_to_rspec(p, mapping = maps)
    expect_equal(peaksData(res), peaksData(sps))
    expect_equal(res$rtime, sps$rtime)
    expect_equal(res$new_col, sps$new_col)
    expect_equal(spectraVariableMapping(), .SPECTRA_2_MATCHMS)

    p <- rspec_to_pyspec(sps, character())
    res <- pyspec_to_rspec(p)
    expect_true(all(is.na(res$msLevel)))
    expect_true(all(is.na(res$rtime)))
    expect_equal(spectraVariableMapping(), .SPECTRA_2_MATCHMS)
})

test_that(".py_matchms_peaks_data works", {
    ## Python variable in R
    p <- r_to_py(sps)
    res <- .py_matchms_peaks_data("r.p", (1:3 - 1L))
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(is.matrix(res[[1L]]))
    res <- lapply(res, function(z) {
        colnames(z) <- c("mz", "intensity")
        z
    })
    expect_equal(res, peaksData(sps@backend))

    res <- .py_matchms_peaks_data("r.p", (c(3, 1, 3) - 1L))
    res <- lapply(res, function(z) {
        colnames(z) <- c("mz", "intensity")
        z
    })
    expect_equal(res, peaksData(sps@backend)[c(3, 1, 3)])

    ## Python variable in Python
    rm("p")
    py_set_attr(py, "p", r_to_py(sps))
    res <- .py_matchms_peaks_data("p", (1:3 - 1L))
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(is.matrix(res[[1L]]))
    res <- lapply(res, function(z) {
        colnames(z) <- c("mz", "intensity")
        z
    })
    expect_equal(res, peaksData(sps@backend))
})

test_that(".py_matchms_peaks_data_cmd works", {
    res <- .py_matchms_peaks_data_cmd("AAA")
    expect_match(res, "AAA")
})

## spectrum_utils

test_that(".py_spectrum_utils_peaks_data_cmd works", {
    res <- .py_spectrum_utils_peaks_data_cmd("BB")
    expect_match(res, "BB")
})

test_that(".py_spectrum_utils_peaks_data works", {
    py_set_attr(py, "p", rspec_to_pyspec(sps, pythonLibrary = "spectrum_utils"))

    res <- .py_spectrum_utils_peaks_data("p", 0:2)
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_true(ncol(res[[1L]]) == 2)
    expect_equal(res[[1L]][, 2L], 1:3)
    expect_equal(res[[1L]], unname(peaksData(sps)[[1L]]),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)

    a <- .py_spectrum_utils_peaks_data("p", 0)
    expect_equal(a[[1L]], res[[1L]])
    a <- .py_spectrum_utils_peaks_data("p", 1)
    expect_equal(a[[1L]], res[[2L]])
    a <- .py_spectrum_utils_peaks_data("p", 2)
    expect_equal(a[[1L]], res[[3L]])
})

test_that(".py_spectrum_utils_spectrum_peaks_data works", {
    tmp <- rspec_to_pyspec(sps, spectraVariableMapping("spectrum_utils"),
                           "spectrum_utils")
    res <- .py_spectrum_utils_spectrum_peaks_data(tmp[0])
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2)
    expect_equal(colnames(res), c("mz", "intensity"))
    expect_equal(res, peaksData(sps)[[1L]],
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    res <- .py_spectrum_utils_spectrum_peaks_data(tmp[1])
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2)
    expect_equal(colnames(res), c("mz", "intensity"))
    expect_equal(res, peaksData(sps)[[2L]],
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    rm(tmp)
})

test_that(".py_spectrum_utils_spectrum_spectra_data works", {
    map <- spectraVariableMapping("spectrum_utils")
    tmp <- rspec_to_pyspec(sps, map, "spectrum_utils")
    res <- .py_spectrum_utils_spectrum_spectra_data(tmp[0])
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 1)
    expect_true(ncol(res) == 2)
    expect_equal(colnames(res), c("precursorMz", "rtime"))
    expect_equal(res$rtime, sps$rtime[1L], tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_equal(res$precursorMz, sps$precursorMz[1L])

    res <- .py_spectrum_utils_spectrum_spectra_data(tmp[1], map)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 1)
    expect_true(ncol(res) == 4)
    expect_equal(colnames(res), c("precursorMz", "precursorCharge", "rtime",
                                  "scanIndex"))
    ## expect_identical(res$precursorCharge, NA_integer_)
    expect_equal(res$rtime, sps$rtime[2L], tolerance = SPECTRUM_UTILS_TOLERANCE)

    res <- .py_spectrum_utils_spectrum_spectra_data(tmp[0], mapping = c())
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, NA_integer_)
})
