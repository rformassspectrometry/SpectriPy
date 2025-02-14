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

test_that(".py_matchms_spectrum_spectra_data works", {
    p <- r_to_py(sps)
    res <- .py_matchms_spectrum_spectra_data(p[[1]])
    expect_true(is.data.frame(res))
    expect_true(all(names(SpectriPy:::.SPECTRA_2_MATCHMS) %in% colnames(res)))
    expect_equal(res$msLevel, sps$msLevel[2])
    expect_equal(res$rtime, sps$rtime[2])
    expect_equal(res$precursorCharge, sps$precursorCharge[2])
    expect_equal(res$precursorMz, sps$precursorMz[2])
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

test_that("py_to_r and .single_pyspec_to_rspec work", {
    p <- r_to_py(sps)
    res <- py_to_r(p[0])
    expect_s4_class(res, "Spectra")
    expect_equal(peaksData(res), peaksData(sps[1L]))
    expect_equal(res$msLevel, sps[1L]$msLevel)
    expect_equal(res$rtime, sps[1L]$rtime)

    res <- .single_pyspec_to_rspec(
        p[1L], mapping = c(rtime = "retention_time"))
    expect_equal(res$rtime, sps[2L]$rtime)
    expect_true(is.na(res$msLevel))

    res <- .single_pyspec_to_rspec(p[1L], mapping = character())
    expect_true(is.na(res$rtime))
    expect_true(is.na(res$msLevel))
})

test_that("pyspec_to_rspec works", {
    p <- r_to_py(sps)
    res <- pyspec_to_rspec(p)
    expect_s4_class(res, "Spectra")
    expect_equal(peaksData(res), peaksData(sps))
    expect_equal(spectraData(res, names(spectraVariableMapping())),
                 spectraData(sps, names(spectraVariableMapping())))

    ## custom mapping.
    maps <- c(rtime = "retention_time", new_col = "r_new_col")
    p <- rspec_to_pyspec(sps, mapping = maps)
    res <- pyspec_to_rspec(p, mapping = maps)
    expect_equal(peaksData(res), peaksData(sps))
    expect_equal(res$rtime, sps$rtime)
    expect_equal(res$new_col, sps$new_col)
    expect_equal(spectraVariableMapping(), SpectriPy:::.SPECTRA_2_MATCHMS)

    p <- rspec_to_pyspec(sps, character())
    res <- pyspec_to_rspec(p)
    expect_true(all(is.na(res$msLevel)))
    expect_true(all(is.na(res$rtime)))
    expect_equal(spectraVariableMapping(), SpectriPy:::.SPECTRA_2_MATCHMS)
})


library(microbenchmark)
library(msdata)
fl <- system.file("TripleTOF-SWATH/PestMix1_DDA.mzML", package = "msdata")
pest_r <- filterMsLevel(Spectra(fl), 2L)
pest_r <- setBackend(pest_r, MsBackendMemory())
system.time(
    pest_py <- rspec_to_pyspec(pest_r)
) # 0.9sec

microbenchmark(
    r_to_py(pest_r),
    pyspec_to_rspec(pest_py),
    concatenateSpectra(py_to_r(pest_py)),
    times = 7
)
## Unit: milliseconds
##                                  expr       min        lq      mean    median       uq
##                       r_to_py(pest_r)  910.6022  927.3481  964.0274  940.4552 1002.636
##              pyspec_to_rspec(pest_py) 1018.6564 1025.0681 1089.1945 1033.7767 1153.097
##  concatenateSpectra(py_to_r(pest_py)) 2722.2907 2772.4679 2891.9192 2807.5928 3036.487
##       max neval cld
##  1037.166     7  a
##  1215.599     7  a
##  3095.642     7   b
