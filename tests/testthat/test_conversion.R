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
    expect_equal(res, SpectriPy:::.SPECTRA_2_MATCHMS)

    setSpectraVariableMapping(c(a = "b", d = "e"))
    expect_equal(spectraVariableMapping(), c(a = "b", d = "e"))
    setSpectraVariableMapping(SpectriPy:::.SPECTRA_2_MATCHMS)
    expect_equal(spectraVariableMapping(),
                 SpectriPy:::.SPECTRA_2_MATCHMS)
})

test_that(".single_rspec_to_pyspec works", {
    res <- .single_rspec_to_pyspec(sps[1L])
    expect_true(is(res, "matchms.Spectrum.Spectrum"))
    expect_true(all(SpectriPy:::.SPECTRA_2_MATCHMS %in% names(res$metadata)))
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

test_that(".rspec_to_pyspec works", {
    res <- .rspec_to_pyspec(sps)
    expect_true(is(res, "python.builtin.list"))
    expect_true(length(res) == length(sps))

    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_true(all(.SPECTRA_2_MATCHMS %in% names(res[[1]]$metadata)))

    res <- .rspec_to_pyspec(sps, character())
    expect_equal(sps$mz[[2L]], as.vector(py_to_r(res[[1]]$mz)))
    expect_equal(sps$intensity[[2L]], as.vector(py_to_r(res[[1]]$intensities)))
    expect_equal(names(res[[1]]$metadata), character())
})

test_that("rspec_to_pyspec works", {
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

    setSpectraVariableMapping(SpectriPy:::.SPECTRA_2_MATCHMS)
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
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], vars)
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), sort(unname(vars)))
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## No metadata
    basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], character())
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), character())
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## Only msLevel
    basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], c(msLevel = "msLevel"))
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), "ms_level")
    }, x = sps)
    basiliskStop(cl)
})

test_that(".single_pyspec_to_rspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping()

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(mz(res), mz(sps[1L]))
    expect_equal(intensity(res), intensity(sps[1L]))
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_equal(msLevel(res), msLevel(sps[1L]))

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[2L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(mz(res), mz(sps[2L]))
    expect_equal(intensity(res), intensity(sps[2L]))
    expect_equal(rtime(res), rtime(sps[2L]))
    expect_equal(msLevel(res), msLevel(sps[2L]))

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[3L])
    res <- SpectriPy:::.single_pyspec_to_rspec(p, vars)
    expect_equal(mz(res), mz(sps[3L]))
    expect_equal(intensity(res), intensity(sps[3L]))
    expect_equal(rtime(res), rtime(sps[3L]))
    expect_equal(msLevel(res), msLevel(sps[3L]))

    basiliskStop(cl)
})

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
