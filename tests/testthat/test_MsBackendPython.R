library(Spectra)
library(MsBackendMgf)
fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
s <- Spectra(fl, source = MsBackendMgf())
s$spectrum_index <- seq_along(s)

## Convert the data to Python MS data structures; seems we need to store
## the Python object in Python, otherwise it can not be seen by the unit
## tests below.
py_set_attr(py, "s_p", rspec_to_pyspec(s, c(defaultSpectraVariableMapping(),
                                            spectrum_index = "spectrum_index")))

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

    expect_error(backendInitialize(
        MsBackendPy(), "some", data = data.frame(msLevel = 1L)),
        "not yet implemented")
})

test_that("show,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_output(show(be), "MsBackendPy")
})

test_that("length,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    expect_equal(length(be), 100L)
    expect_equal(length(MsBackendPy()), 0L)
})

test_that("spectraVariables and .py_get_metadata_names works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
    vars <- SpectriPy:::.py_get_metadata_names(be)
    expect_true(is.character(vars))
    expect_true(
        all(c("charge", "collision_energy", "ms_level",
              "precursor_intensity", "precursor_mz", "retention_time") %in%
            vars))
    expect_true(any(vars == "spectrum_index"))
    expect_true(any(names(vars) == "spectrum_index"))
    res <- spectraVariables(be)
    expect_true(all(names(coreSpectraVariables()) %in% res))
})

test_that("peaksData,MsBackendPy works", {
    be <- backendInitialize(MsBackendPy(), "s_p")
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
    expect_true(any(colnames(res) == "spectrum_index"))

    ## arbitrary order
    cols <- c("rtime", "msLevel", "spectrum_index", "precursorMz", "polarity")
    res <- spectraData(be, cols)
    expect_equal(colnames(res), cols)
    expect_equal(res$precursorMz, s$precursorMz)
    expect_equal(res$msLevel, s$msLevel)

    ## individual columns
    res <- spectraData(be, "msLevel")
    expect_s4_class(res, "DataFrame")
    expect_equal(res$msLevel, s$msLevel)

    res <- spectraData(be, "spectrum_index")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "spectrum_index")
    expect_equal(res$spectrum_index, s$spectrum_index)

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

    ## non-existant column
    expect_error(spectraData(be, c("msLevel", "other_col")), "not available")
})
