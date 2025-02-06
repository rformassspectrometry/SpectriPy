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
    res <- spectraVariableMapping(pythonLibrary = "matchms")
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    res <- spectraVariableMapping(pythonLibrary = "spectrum_utils")
    expect_true(is.character(res))
    expect_true(length(res) > 0)
})

test_that("rspec_to_pyspec works", {
    ## matchms
    cl <- basiliskStart(SpectriPy:::matchms_env)
    #basiliskRun(cl, function(x) {
        res <- rspec_to_pyspec(sps, pythonLibrary = "matchms")
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(sps)[[3]])
        expect_equal(rtime(sps)[1], py_to_r(res[0]$metadata$retention_time))
    #}, x = sps)

    #basiliskRun(cl, function(x) {
        sps_tmp <- sps
        sps_tmp$new_col <- c(9, 3, 1)
        res <- rspec_to_pyspec(sps_tmp, pythonLibrary = "matchms",
            spectraVariables = c(rtime = "retention_time", new_col = "new_col"))
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps_tmp))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(sps_tmp)[[3]])
        expect_equal(rtime(sps_tmp)[1], py_to_r(res[0]$metadata$retention_time))
        expect_equal(sps_tmp$new_col[1], py_to_r(res[0]$metadata$new_col))
    #}, x = sps)

    #basiliskRun(cl, function(x) {
        sps_tmp <- sps
        sps_tmp$rtime <- NULL
        res <- rspec_to_pyspec(sps_tmp, pythonLibrary = "matchms")
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps_tmp))
        expect_equal(as.numeric(py_to_r(res[2]$peaks$mz)), mz(sps_tmp)[[3]])
        ## Seems NA for rtime is changed to NULL in Spectrum
        expect_equal(py_to_r(res[0]$metadata$retention_time), NULL)
    #}, x = sps)

    basiliskStop(cl)
    
    ## spectrum_utils
    cl <- basiliskStart(SpectriPy:::spectrum_utils_env)
    #basiliskRun(cl, function(x) {
        res <- rspec_to_pyspec(sps, pythonLibrary = "spectrum_utils")
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps))
        expect_equal(as.numeric(py_to_r(res[2]$mz)), mz(sps)[[3]], tolerance = 1e-06)
        expect_equal(rtime(sps)[1], unname(unlist(py_to_r(res[0]$retention_time))))
    #}, x = sps)
    
    #basiliskRun(cl, function(x) {
        sps_tmp <- sps
        sps_tmp$new_col <- c(9, 3, 1)
        res <- rspec_to_pyspec(sps_tmp, pythonLibrary = "spectrum_utils", 
                            spectraVariables = c(rtime = "retention_time",
                                              new_col = "new_col"))
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps_tmp))
        expect_equal(as.numeric(py_to_r(res[2]$mz)), mz(sps_tmp)[[3]], tolerance = 1e-06)
        expect_equal(rtime(sps_tmp)[1], unname(unlist(py_to_r(res[0]$retention_time))))
        expect_equal(sps_tmp$new_col[1], unname(unlist(py_to_r(res[0]$new_col))))
    #}, x = sps)
    
    #basiliskRun(cl, function(x) {
        sps_tmp <- sps
        sps_tmp$rtime <- NULL
        res <- rspec_to_pyspec(sps_tmp, pythonLibrary = "spectrum_utils")
        expect_true(is(res, "python.builtin.list"))
        expect_equal(length(res), length(sps_tmp))
        expect_equal(as.numeric(py_to_r(res[2]$mz)), mz(sps_tmp)[[3]], tolerance = 1e-06)
        ## Seems NA for rtime is changed to NULL in Spectrum
        expect_equal(unname(unlist(py_to_r(res[0]$retention_time))), NaN)
    #}, x = sps)
    
    basiliskStop(cl)
    
})

test_that("pyspec_to_rspec works", {
    
    ## matchms
    cl <- basiliskStart(SpectriPy:::matchms_env)
    #basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(sps, pythonLibrary = "matchms")
        res <- pyspec_to_rspec(p, pythonLibrary = "matchms", 
            spectraVariables = spectraVariableMapping(pythonLibrary = "matchms"))
        expect_true(is(res, "Spectra"))
        expect_equal(spectraData(sps), spectraData(res))
    #}, x = sps)

    ## Map only selected values back.
    #basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(sps)
        res <- pyspec_to_rspec(p, pythonLibrary = "matchms", 
            spectraVariables = c(rtime = "retention_time",
                what = "not_exists"))
        expect_true(is(res, "Spectra"))
        expect_equal(mz(sps), mz(res))
        expect_equal(rtime(sps), rtime(res))
        expect_equal(intensity(sps), intensity(res))
        expect_true(all(is.na(msLevel(res))))
    #}, x = sps)

    #basiliskRun(cl, function(x) {
        sps$rtime <- NULL
        p <- rspec_to_pyspec(sps, pythonLibrary = "matchms")
        res <- pyspec_to_rspec(p, pythonLibrary = "matchms")
        expect_equal(mz(sps), mz(res))
        expect_equal(rtime(sps), rtime(res))
        expect_equal(intensity(sps), intensity(res))
    #}, x = sps)
    ## errors
    expect_error(pyspec_to_rspec(5), "Python list")

    basiliskStop(cl)
    
    ## spectrum_utils
    cl <- basiliskStart(SpectriPy:::spectrum_utils_env)
    #basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(sps, pythonLibrary = "spectrum_utils")
        res <- pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", 
            spectraVariables = spectraVariableMapping(pythonLibrary = "spectrum_utils"))
        expect_true(is(res, "Spectra"))
        expect_equal(spectraData(sps)[, -1], spectraData(res)[, -1]) ## remove msLevel
    #}, x = sps)
    
    ## Map only selected values back.
    #basiliskRun(cl, function(x) {
        p <- rspec_to_pyspec(sps, pythonLibrary = "spectrum_utils")
        res <- pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", 
            spectraVariables = c(rtime = "retention_time", what = "not_exists"))
        expect_true(is(res, "Spectra"))
        expect_equal(mz(sps), mz(res), tolerance = 1e-06)
        expect_equal(rtime(sps), rtime(res))
        expect_equal(intensity(sps), intensity(res))
        expect_true(all(is.na(msLevel(res))))
    ##}, x = sps)
    
    #basiliskRun(cl, function(x) {
        sps$rtime <- NULL
        p <- rspec_to_pyspec(sps, pythonLibrary = "spectrum_utils")
        res <- pyspec_to_rspec(p, pythonLibrary = "spectrum_utils")
        expect_equal(mz(sps), mz(res), tolerance = 1e-06)
        expect_equal(rtime(sps), rtime(res))
        expect_equal(intensity(sps), intensity(res))
    #}, x = sps)
    ## errors
    expect_error(pyspec_to_rspec(5), "Python list")
    
    basiliskStop(cl)
})

test_that(".single_rspec_to_pyspec works", {
    
    ## matchms
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping(pythonLibrary = "matchms")
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = vars, pythonLibrary = "matchms")
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), sort(unname(vars)))
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    #}, x = sps)
    ## No metadata
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = character(), pythonLibrary = "matchms")
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), character())
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    #}, x = sps)
    ## Only msLevel
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = c(msLevel = "msLevel"), pythonLibrary = "matchms")
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), "ms_level")
    #}, x = sps)
    basiliskStop(cl)
    
    ## spectrum_utils
    .names <- c("annotate_molecule_fragment", "annotate_mz_fragment",
        "annotate_peaks", "annotate_peptide_fragments", "annotation",
        "filter_intensity", "identifier", "intensity", "is_decoy", 
        "modifications", "mz", "peptide", "precursor_charge",
        "precursor_mz", "remove_precursor_peak", "retention_time",            
        "round", "scale_intensity", "set_mz_range")
    cl <- basiliskStart(SpectriPy:::spectrum_utils_env)
    vars <- spectraVariableMapping(pythonLibrary = "spectrum_utils")
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = vars, pythonLibrary = "spectrum_utils")
        expect_equal(class(res)[1L], "spectrum_utils.spectrum.MsmsSpectrum")
        expect_equal(sort(names(res)), sort(.names))
        expect_equal(as.numeric(res$intensity), intensity(sps)[[1L]])
    #}, x = sps)
    ## No metadata
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = character(), pythonLibrary = "spectrum_utils")
        expect_equal(class(res)[1L], "spectrum_utils.spectrum.MsmsSpectrum")
        expect_equal(sort(names(res)), sort(.names)) ## !!! previous: character()
        expect_equal(as.numeric(res$intensity), intensity(sps)[[1L]])
    #}, x = sps)
    ## Only msLevel
    #basiliskRun(cl, function(x) {
        res <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], spectraVariables = c(msLevel = "msLevel"), pythonLibrary = "spectrum_utils")
        expect_equal(class(res)[1L], "spectrum_utils.spectrum.MsmsSpectrum")
        expect_equal(sort(names(res)), c(sort(.names)[1:10], "msLevel", sort(.names)[11:19])) ## !!! previous: "ms_level"
    #}, x = sps)
    basiliskStop(cl)
})

test_that(".single_pyspec_to_rspec works", {
    
    ## matchms
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping(pythonLibrary = "matchms")

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "matchms", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[1L]))
    expect_equal(intensity(res), intensity(sps[1L]))
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_equal(msLevel(res), msLevel(sps[1L]))

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[2L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "matchms", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[2L]))
    expect_equal(intensity(res), intensity(sps[2L]))
    expect_equal(rtime(res), rtime(sps[2L]))
    expect_equal(msLevel(res), msLevel(sps[2L]))

    p <- SpectriPy:::.single_rspec_to_pyspec(sps[3L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "matchms", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[3L]))
    expect_equal(intensity(res), intensity(sps[3L]))
    expect_equal(rtime(res), rtime(sps[3L]))
    expect_equal(msLevel(res), msLevel(sps[3L]))

    basiliskStop(cl)
    
    ## spectrum_utils
    cl <- basiliskStart(SpectriPy:::spectrum_utils_env)
    vars <- spectraVariableMapping(pythonLibrary = "spectrum_utils")
    
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "spectrum_utils", spectraVariables = vars)
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[1L]), tolerance = 1e-06)
    expect_equal(intensity(res), intensity(sps[1L]))
    expect_equal(rtime(res), rtime(sps[1L]))
    ##expect_equal(msLevel(res), msLevel(sps[1L])) ## !!! not possible to set msLevel in spectrum_utils
    
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[2L], pythonLibrary = "spectrum_utils")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[2L]), tolerance = 1e-06)
    expect_equal(intensity(res), intensity(sps[2L]))
    expect_equal(rtime(res), rtime(sps[2L]))
    ##expect_equal(msLevel(res), msLevel(sps[2L])) ## !!! not possible to set msLevel in spectrum_utils
    
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[3L], pythonLibrary = "spectrum_utils")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_equal(mz(res), mz(sps[3L]), tolerance = 1e-06)
    expect_equal(intensity(res), intensity(sps[3L]))
    expect_equal(rtime(res), rtime(sps[3L]))
    ##expect_equal(msLevel(res), msLevel(sps[3L])) ## !!! not possible to set msLevel in spectrum_utils
    
    basiliskStop(cl)
})

test_that(".single_pyspec_to_rspec works", {
    ## matchms
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping(pythonLibrary = "matchms")

    ## Request single spectra variable
    vars <- c(rtime = "retention_time")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "matchms", spectraVariables = vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))

    ## Request spectra variables that don't exist.
    vars <- c(rtime = "retention_time", other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p,  pythonLibrary = "matchms", spectraVariables = vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))

    ## Request spectra variables that don't exist.
    vars <- c(other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "matchms")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "matchms", spectraVariables = vars)
    expect_true(is.na(rtime(res)))
    expect_true(is.na(msLevel(res)))
    basiliskStop(cl)
    
    ## spectrum_utils
    cl <- basiliskStart(SpectriPy:::spectrum_utils_env)
    vars <- spectraVariableMapping(pythonLibrary = "spectrum_utils")
    
    ## Request single spectra variable
    vars <- c(rtime = "retention_time")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "spectrum_utils")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))
    
    ## Request spectra variables that don't exist.
    vars <- c(rtime = "retention_time", other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "spectrum_utils")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_equal(rtime(res), rtime(sps[1L]))
    expect_true(is.na(msLevel(res)))
    
    ## Request spectra variables that don't exist.
    vars <- c(other_col = "other_col", b = "b")
    p <- SpectriPy:::.single_rspec_to_pyspec(sps[1L], pythonLibrary = "spectrum_utils")
    res <- SpectriPy:::.single_pyspec_to_rspec(p, pythonLibrary = "spectrum_utils", spectraVariables = vars)
    expect_true(is.na(rtime(res)))
    expect_true(is.na(msLevel(res)))
    basiliskStop(cl)
    
})
