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

test_that("spectraVariableMapping works", {
    res <- spectraVariableMapping()
    expect_true(is.character(res))
    expect_true(length(res) > 0)
})

test_that("rspec_to_pyspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    basiliskRun(cl, function(x) {
        res <- rspec_to_pyspec(x)
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(x))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(x)[[3]])
        expect_equal(rtime(x)[1], py_to_r(res[0]$metadata$retention_time))
    }, x = sps)

    basiliskRun(cl, function(x) {
        x$new_col <- c(9, 3, 1)
        res <- rspec_to_pyspec(x, mapping = c(rtime = "retention_time",
                                              new_col = "new_col"))
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(x))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(x)[[3]])
        expect_equal(rtime(x)[1], py_to_r(res[0]$metadata$retention_time))
        expect_equal(x$new_col[1], py_to_r(res[0]$metadata$new_col))
    }, x = sps)

    basiliskRun(cl, function(x) {
        x$rtime <- NULL
        res <- rspec_to_pyspec(x)
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(x))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(x)[[3]])
        ## Seems NA for rtime is changed to NULL in Spectrum
        expect_equal(py_to_r(res[0]$metadata$retention_time), NULL)
    }, x = sps)

    basiliskStop(cl)
})

test_that("pyspec_to_rspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(x)
        res <- pyspec_to_rspec(p, mapping = spectraVariableMapping())
        expect_true(is(res, "Spectra"))
        expect_equal(spectraData(x), spectraData(res))
    }, x = sps)

    ## Map only selected values back.
    basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(x)
        res <- pyspec_to_rspec(p, mapping = c(rtime = "retention_time",
                                              what = "not_exists"))
        expect_true(is(res, "Spectra"))
        expect_equal(mz(x), mz(res))
        expect_equal(rtime(x), rtime(res))
        expect_equal(intensity(x), intensity(res))
        expect_true(all(is.na(msLevel(res))))
    }, x = sps)

    basiliskRun(cl, function(x) {
        x$rtime <- NULL
        p <- rspec_to_pyspec(x)
        res <- pyspec_to_rspec(p)
        expect_equal(mz(x), mz(res))
        expect_equal(rtime(x), rtime(res))
        expect_equal(intensity(x), intensity(res))
    }, x = sps)
    ## errors
    expect_error(pyspec_to_rspec(5), "Python list")

    basiliskStop(cl)
})

test_that(".single_rspec_to_pyspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping()
    basiliskRun(cl, function(x) {
        res <- .single_rspec_to_pyspec(sps[1L], vars)
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), sort(unname(vars)))
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## No metadata
    basiliskRun(cl, function(x) {
        res <- .single_rspec_to_pyspec(sps[1L], character())
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), character())
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## Only msLevel
    basiliskRun(cl, function(x) {
        res <- .single_rspec_to_pyspec(sps[1L], c(msLevel = "msLevel"))
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), "mslevel")
    }, x = sps)
    basiliskStop(cl)
})

with_parameters_test_that(".single_pyspec_to_rspec works with parameters", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping()

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[index])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(mz(res), mz(sps[index]))
    expect_equal(intensity(res), intensity(sps[index]))
    expect_equal(rtime(res), rtime(sps[index]))
    expect_equal(msLevel(res), msLevel(sps[index]))

    basiliskStop(cl)
}, cases(
       first = list(index = 1L),
       second = list(index = 2L),
       third = list(index = 3L)
   ))

test_that(".single_pyspec_to_rspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping()

    ## Request single spectra variable
    vars <- c(rtime = "retention_time")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))

    ## Request spectra variables that don't exist.
    vars <- c(rtime = "retention_time", other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))

    ## Request spectra variables that don't exist.
    vars <- c(other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_true(is.na(rtime(res)))
    expect_true(is.na(msLevel(res)))
    basiliskStop(cl)
})
