library(Spectra)
library(MsBackendMgf)
fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
s <- Spectra(fl, source = MsBackendMgf())
s$spectrum_index <- seq_along(s)

## Convert the data to Python MS data structures; seems we need to store
## the Python object in Python, otherwise it can not be seen by the unit
## tests below.
py_set_attr(py, "s_p", rspec_to_pyspec(s, c(defaultSpectraVariableMapping(),
                                            spectrum_index = "spctrm_idx")))
expect_true(any(names(py$s_p[[1]]$metadata) == "spctrm_idx"))

## Convert the data also to spectrum_utils
py_set_attr(py, "su_p", rspec_to_pyspec(
                            s, spectraVariableMapping("spectrum_utils"),
                            "spectrum_utils"))
SPECTRUM_UTILS_TOLERANCE <- 1e-7

test_that("MsBackendPy constructor works", {
    res <- MsBackendPy()
    expect_s4_class(res, "MsBackendPy")
    expect_true(res@is_in_py)
    expect_equal(res@spectraVariableMapping, defaultSpectraVariableMapping())
})

test_that(".check_spectra_variable_mapping works", {
    expect_true(.check_spectra_variable_mapping(character()))
    expect_true(.check_spectra_variable_mapping(
        defaultSpectraVariableMapping()))
    expect_error(.check_spectra_variable_mapping("hello"),
                 "named character vector")
})

test_that(".get_py works", {
    res <- .get_py()
    expect_true(is(res, "python.builtin.module"))

    ## expect_error(
    ##     with_mocked_bindings(
    ##         "base::get" = function(x, pos = -1L, envir = as.environment(pos),
    ##                          mode = "any", inherits = TRUE) stop("aaaa"),
    ##         code = .get_py()
    ##     ),
    ##     "Failed to get variable"
    ## )
})

test_that(".exists_py_var works", {
    expect_false(.exists_py_var("what_"))
    expect_true(.exists_py_var("r"))
    expect_true(.exists_py_var("s_p"))
})

test_that(".check_py_var_exists works", {
    expect_true(.check_py_var_exists("r"))
    expect_error(.check_py_var_exists("r", FALSE), "No variable")
    expect_true(.check_py_var_exists("py", FALSE))
    expect_true(.check_py_var_exists("s_p", TRUE))
    ## expect_true(.check_py_var_exists("s", FALSE))
})

test_that(".check_py_var works", {
    expect_true(.check_py_var("s_p", TRUE))
    expect_error(.check_py_var("py", FALSE), "supposed to be")
})

test_that(".get_py_var works", {
    ## For whatever reason I can't get something from the global environment.
    ## res <- .get_py_var("s", FALSE)
    ## expect_s4_class(res, "Spectra")
    py$tmp <- 34
    res <- .get_py_var("tmp", TRUE)
    expect_equal(py_to_r(res), 34)
})

test_that("backendInitialize,MsBackendPy works", {
    expect_error(backendInitialize(MsBackendPy()), "has to be provided")
    expect_error(backendInitialize(MsBackendPy(), 3), "name of the")
    expect_error(backendInitialize(MsBackendPy(), "not exists"),
                 "No variable of name")
    expect_error(backendInitialize(MsBackendPy(), "r"), "Python list")
    res <- backendInitialize(MsBackendPy(), "s_p")
    expect_s4_class(res, "MsBackendPy")
})

test_that("backendInitialize,MsBackendPy works with providing data", {
    expect_error(
        backendInitialize(MsBackendPy(), pythonVariableName = "ttt", data = 4),
        "'DataFrame'")
    d <- spectraData(s)
    expect_error(
        backendInitialize(MsBackendPy(), pythonVariableName = "ttt", data = d),
        "are required")
    d <- spectraData(s@backend)
    res <- backendInitialize(
        MsBackendPy(), pythonVariableName = "ttt", data = d)
    expect_s4_class(res, "MsBackendPy")
    expect_equal(res@py_var, "ttt")
    expect_true(res@is_in_py)
    expect_equal(nrow(d), length(res))
    expect_equal(rtime(s), rtime(res))
    expect_equal(intensity(s), intensity(res))
    expect_equal(mz(s), mz(res))
    ## Repeat with spectrum_utils
    py_del_attr(py, "ttt")
    res <- backendInitialize(MsBackendPy(), pythonVariableName = "ttt",
                             data = d, pythonLibrary = "spectrum_utils")
    expect_s4_class(res, "MsBackendPy")
    expect_equal(res@py_var, "ttt")
    expect_true(res@is_in_py)
    expect_equal(nrow(d), length(res))
    expect_equal(rtime(s), rtime(res))
    expect_equal(intensity(s), intensity(res),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    expect_equal(mz(s), mz(res), tolerance = SPECTRUM_UTILS_TOLERANCE)
})

test_that("show,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_output(show(be), "MsBackendPy")
})

test_that(".py_var_length, length,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(.py_var_length(be), 100L)
    expect_equal(.py_var_length(MsBackendPy()), 0L)
    expect_equal(length(be), 100L)
    expect_equal(length(MsBackendPy()), 0L)

    ## dddd <- 1:10
    ## be2 <- MsBackendPy()
    ## be2@is_in_py <- FALSE
    ## be2@py_var <- "dddd"
    ## expect_equal(SpectriPy:::.py_var_length(be2), 10L)
})

test_that("spectraVariables and .py_get_metadata_names works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    vars <- SpectriPy:::.py_get_metadata_names(be)
    expect_true(is.character(vars))
    expect_true(
        all(c("charge", "collision_energy", "ms_level",
              "precursor_intensity", "precursor_mz", "retention_time") %in%
            vars))
    expect_true(any(vars == "spctrm_idx"))
    expect_true(any(names(vars) == "spctrm_idx"))
    res <- spectraVariables(be)
    expect_true(all(names(coreSpectraVariables()) %in% res))
    expect_true("spctrm_idx" %in% res)

    m <- c(defaultSpectraVariableMapping(), spectrum_index = "spctrm_idx")
    be@spectraVariableMapping <- m
    res <- spectraVariables(be)
    expect_true(all(names(coreSpectraVariables()) %in% res))
    expect_true("spectrum_index" %in% res)

    ## spectrum_utils
    be <- backendInitialize(MsBackendPy(), "su_p",
                            pythonLibrary = "spectrum_utils")
    m <- .py_get_metadata_names(be)
    expect_equal(unname(m), c("precursor_mz", "retention_time"))
    be <- backendInitialize(
        MsBackendPy(), "su_p",
        spectraVariableMapping = spectraVariableMapping("spectrum_utils"),
        pythonLibrary = "spectrum_utils")
    m <- .py_get_metadata_names(be)
    expect_equal(unname(m), c("precursor_mz", "precursor_charge",
                      "retention_time", "identifier"))
    res <- spectraVariables(be)
    expect_true(all(names(m) %in% res))

    be2 <- MsBackendPy()
    be2@py_var <- "hey"
    be2@is_in_py <- FALSE
    hey <- integer()
    expect_equal(.py_get_metadata_names(be2), character())
})

test_that("peaksData,MsBackendPy works", {
    ## Empty data
    be <- MsBackendPy()
    res <- peaksData(be)
    expect_equal(res, list())

    ## Read data
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(be@i, 1:100)
    res <- peaksData(be)
    expect_equal(res, peaksData(s@backend))

    res <- peaksData(be, c("intensity", "mz"))
    expect_equal(res, peaksData(s@backend, c("intensity", "mz")))
    expect_equal(colnames(res[[1L]]), c("intensity", "mz"))

    res <- peaksData(be, "mz")
    expect_equal(res, peaksData(s@backend, "mz"))
    expect_equal(colnames(res[[1L]]), c("mz"))

    res <- peaksData(be, "mz", drop = TRUE)
    expect_true(is.list(res))
    expect_equal(length(res), 100L)
    expect_true(is.numeric(res[[1L]]))
    expect_false(is.matrix(res[[1L]]))

    ## spectrum_utils
    be <- backendInitialize(MsBackendPy(), "su_p",
                            pythonLibrary = "spectrum_utils")
    res <- peaksData(be)
    expect_true(is.list(res))
    expect_equal(length(res), length(be))
    expect_equal(res, peaksData(s, return.type = "list"),
                 tolerance = SPECTRUM_UTILS_TOLERANCE)
    a <- peaksData(be[4])
    expect_equal(a, res[4])
})

test_that("spectraData,MsBackendPy works", {
    be <- MsBackendPy()
    res <- spectraData(be)
    expect_true(nrow(res) == 0)
    expect_true(is(res, "DataFrame"))

    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- spectraData(be)
    expect_true(nrow(res) == 100)
    expect_s4_class(res, "DataFrame")
    expect_true(all(names(coreSpectraVariables()) %in% colnames(res)))
    expect_true(all(names(be@spectraVariableMapping) %in% colnames(res)))
    expect_true(any(colnames(res) == "spctrm_idx"))

    ## arbitrary order
    cols <- c("rtime", "msLevel", "spctrm_idx", "precursorMz", "polarity")
    res <- spectraData(be, cols)
    expect_equal(colnames(res), cols)
    expect_equal(res$precursorMz, s$precursorMz)
    expect_equal(res$msLevel, s$msLevel)

    ## individual columns
    res <- spectraData(be, "msLevel")
    expect_s4_class(res, "DataFrame")
    expect_equal(res$msLevel, s$msLevel)

    res <- spectraData(be, "spctrm_idx")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "spctrm_idx")
    expect_equal(res$spctrm_idx, s$spectrum_index)

    res <- spectraData(be, "smoothed")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "smoothed")
    expect_equal(res$smoothed, s$smoothed)

    ## only mz
    res <- spectraData(be, "mz")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "mz")
    expect_equal(res$mz, s$mz)

    ## only intensity
    res <- spectraData(be, "intensity")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "intensity")
    expect_equal(res$intensity, s$intensity)

    ## only core variables not in backend
    cols <- c("polarity", "smoothed", "centroided")
    res <- spectraData(be, cols)
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), cols)
    expect_true(all(is.na(res)))

    ## only variables that are mapped
    cols <- c("precursorMz", "precursorCharge", "msLevel")
    res <- spectraData(be, cols)
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), cols)
    expect_equal(res$precursorMz, s$precursorMz)
    expect_equal(res$precursorCharge, s$precursorCharge)
    expect_equal(res$msLevel, s$msLevel)

    ## only one with drop TRUE
    res <- spectraData(be, "msLevel", drop = TRUE)
    expect_equal(res, s$msLevel)

    ## non-existant column
    expect_error(spectraData(be, c("msLevel", "other_col")), "not available")

    ## custom spectravariablemapping
    m <- c(defaultSpectraVariableMapping(),
           spectrum_index = "spctrm_idx")
    be@spectraVariableMapping <- m
    res <- spectraData(be)
    expect_true(all(names(m) %in% colnames(res)))
    expect_equal(res$spectrum_index, seq_along(be))

    ## spectrum_utils
    be <- backendInitialize(
        MsBackendPy(), "su_p", pythonLibrary = "spectrum_utils",
        spectraVariableMapping("spectrum_utils"))
    res <- spectraData(be)
    expect_s4_class(res, "DataFrame")
    expect_equal(res$dataStorage, rep("su_p", length(be)))
    expect_equal(res$precursorMz, s$precursorMz)
})

test_that(".check_i works", {
    expect_true(.check_i(MsBackendPy()))
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_true(.check_i(be))
    be@i <- 1233:3431
    expect_error(.check_i(be), "out of bound")
})

test_that("reindex works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(be@i, 1:100)
    be@i <- c(4L, 3L)
    be <- reindex(be)
    expect_equal(be@i, 1:100)
})

test_that("extractByIndex and [,MsBackendPy work", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- be[c(43, 2)]
    expect_equal(res@i, c(43L, 2L))

    expect_equal(peaksData(res), peaksData(s@backend)[c(43, 2)])
    cols <- c("msLevel", "precursorMz", "rtime")
    expect_equal(spectraData(res, cols),
                 spectraData(s@backend, cols)[c(43, 2), ])

    ## Single element
    res <- be[3L]
    expect_equal(peaksData(res), peaksData(s@backend)[3L])
})

test_that("$,MsBackendPy works", {
    be <- MsBackendPy()
    expect_equal(be$msLevel, integer())
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(be$msLevel, s$msLevel)
    expect_equal(be$precursorMz, s$precursorMz)
    expect_error(be$not_exists, "not available")

    be <- backendInitialize(
        MsBackendPy(), "su_p",
        spectraVariableMapping = spectraVariableMapping("spectrum_utils"),
        pythonLibrary = "spectrum_utils")
    expect_equal(be$precursorMz, s$precursorMz)
})

test_that("lengths,MsBackendPy works", {
    be <- MsBackendPy()
    expect_equal(lengths(be), integer(0))

    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(lengths(be), lengths(s))
})

test_that("isEmpty,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_false(any(isEmpty(be)))
})

test_that("acquisitionNum,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- acquisitionNum(be)
    expect_true(is.integer(res))
    expect_true(all(is.na(res)))
})

test_that("centroided,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- centroided(be)
    expect_true(is.logical(res))
    expect_true(all(is.na(res)))
})

test_that("collisionEnergy,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- collisionEnergy(be)
    expect_true(is.numeric(res))
    expect_true(all(is.na(res)))
})

test_that("dataOrigin,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- dataOrigin(be)
    expect_true(is.character(res))
    expect_true(all(is.na(res)))
})

test_that("intensity,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- intensity(be)
    expect_equal(res, intensity(s@backend))
    expect_s4_class(res, "NumericList")
})

test_that("isolationWindowLowerMz,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- isolationWindowLowerMz(be)
    expect_true(is.numeric(res))
    expect_equal(res, isolationWindowLowerMz(s@backend))
    expect_true(all(is.na(res)))
})

test_that("isolationWindowUpperMz,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- isolationWindowUpperMz(be)
    expect_true(is.numeric(res))
    expect_equal(res, isolationWindowUpperMz(s@backend))
    expect_true(all(is.na(res)))
})

test_that("isolationWindowTargetMz,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- isolationWindowTargetMz(be)
    expect_true(is.numeric(res))
    expect_equal(res, isolationWindowTargetMz(s@backend))
    expect_true(all(is.na(res)))
})

test_that("msLevel,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- msLevel(be)
    expect_true(is.integer(res))
    expect_equal(res, msLevel(s@backend))
})

test_that("mz,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- mz(be)
    expect_equal(res, mz(s@backend))
    expect_s4_class(res, "NumericList")
})

test_that("polarity,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- polarity(be)
    expect_true(is.integer(res))
    expect_equal(res, polarity(s@backend))
    expect_true(all(is.na(res)))
})

test_that("precScanNum,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- precScanNum(be)
    expect_true(is.integer(res))
    expect_equal(res, precScanNum(s@backend))
    expect_true(all(is.na(res)))
})

test_that("precursorCharge,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- precursorCharge(be)
    expect_true(is.integer(res))
    expect_equal(res, precursorCharge(s@backend))
})

test_that("precursorIntensity,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- precursorIntensity(be)
    expect_true(is.numeric(res))
    expect_equal(res, precursorIntensity(s@backend))
})

test_that("precursorMz,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- precursorMz(be)
    expect_true(is.numeric(res))
    expect_equal(res, precursorMz(s@backend))
})

test_that("rtime,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- rtime(be)
    expect_true(is.numeric(res))
    expect_equal(res, rtime(s@backend))

    s2 <- s[1:4]
    s2$rtime <- c(1.2, 3.4, 5.3, 5.9)
    py_set_attr(py, "s_p2",
                rspec_to_pyspec(s2, defaultSpectraVariableMapping()))
    be <- backendInitialize(MsBackendPy(), "s_p2")
    expect_equal(rtime(be), rtime(s2))
})

test_that("scanIndex,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- scanIndex(be)
    expect_true(is.integer(res))
    expect_equal(res, scanIndex(s@backend))
})

test_that("scanIndex,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- smoothed(be)
    expect_true(is.logical(res))
    expect_equal(res, smoothed(s@backend))
})

test_that("spectraNames,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- spectraNames(be)
    expect_equal(res, spectraNames(s@backend))
})

test_that("tic,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    res <- tic(be, initial = TRUE)
    expect_true(is.numeric(res))
    expect_equal(res, tic(s@backend, initial = TRUE))
    res <- tic(be, initial = FALSE)
    expect_true(is.numeric(res))
    expect_equal(res, tic(s@backend, initial = FALSE))

    s2 <- s[1:4]
    s2$totIonCurrent <- 123.2
    py_set_attr(py, "s_p2",
                rspec_to_pyspec(s2, c(defaultSpectraVariableMapping(),
                                      totIonCurrent = "totIonCurrent")))
    be <- backendInitialize(
        MsBackendPy(), "s_p2",
        spectraVariableMapping = c(defaultSpectraVariableMapping(),
                                   totIonCurrent = "totioncurrent"))
    expect_equal(tic(be, initial = TRUE), c(123.2, 123.2, 123.2, 123.2))
})

test_that("spectraVariableMapping,spectraVariableMapping<- works", {
    s2 <- s[1:4]
    s2$totIonCurrent <- 123.2
    py_set_attr(py, "s_p2",
                rspec_to_pyspec(s2, c(defaultSpectraVariableMapping(),
                                      totIonCurrent = "tot_ion_current")))
    be <- backendInitialize(MsBackendPy(), "s_p2")
    expect_equal(spectraVariableMapping(be), defaultSpectraVariableMapping())

    expect_false("totIonCurrent" %in% spectraVariables(be))
    spectraVariableMapping(be) <- c(spectraVariableMapping(be),
                                    totIonCurrent = "tot_ion_current")
    expect_true("totIonCurrent" %in% spectraVariables(be))

    s3 <- Spectra(be)
    expect_error(spectraVariableMapping(s3) <- "rtime", "needs to be a named")
    spectraVariableMapping(s3) <- c(rtime = "retention_time")
    expect_equal(spectraVariableMapping(s3@backend),
                 c(rtime = "retention_time"))
})
