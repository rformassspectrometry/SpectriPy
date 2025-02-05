caf <- DataFrame(
    msLevel = c(2L, 2L),
    name = "Caffeine",
    precursorMz = c(195.0877, 195.0877)
)
caf$intensity <- list(
    c(340.0, 416, 2580, 412),
    c(388.0, 3270, 85, 54, 10111))
caf$mz <- list(
    c(135.0432, 138.0632, 163.0375, 195.0880),
    c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864))
caf <- Spectra(caf)

## select_by_intensityParam
test_that("select_by_intensityParam constructor works", {
    ## errors
    expect_error(select_by_intensityParam(intensity_from = c(1.2, 1.5)), "length 1")
    expect_error(select_by_intensityParam(intensity_from = -1.2), "length 1")
    expect_error(select_by_intensityParam(intensity_to = c(1, 2)), "length 1")
    expect_error(select_by_intensityParam(intensity_to = c(1, 2)), "length 1")
    expect_error(select_by_intensityParam(foo = 1), "unused argument")
    
    res <- select_by_intensityParam(intensity_from = 5, intensity_to = 20)
    expect_s4_class(res, "select_by_intensityParam")
    expect_equal(res@intensity_from, 5)
    expect_equal(res@intensity_to, 20)
})

## select_by_mzParam
test_that("select_by_mzParam constructor works", {
    ## errors
    expect_error(select_by_mzParam(mz_from = c(1.2, 1.5)), "length 1")
    expect_error(select_by_mzParam(mz_from = -1.2), "length 1")
    expect_error(select_by_mzParam(mz_to = c(1, 2)), "length 1")
    expect_error(select_by_mzParam(mz_to = c(1, 2)), "length 1")
    expect_error(select_by_mzParam(foo = 1), "unused argument")
    
    res <- select_by_mzParam(mz_from = 0, mz_to = 400)
    expect_s4_class(res, "select_by_mzParam")
    expect_equal(res@mz_from, 0)
    expect_equal(res@mz_to, 400)
})

## remove_peaks_around_precursor_mz
test_that("remove_peaks_around_precursor_mzParam constructor works", {
    ## errors
    expect_error(remove_peaks_around_precursor_mzParam(mz_tolerance = c(1.2, 1.5)), "length 1")
    expect_error(remove_peaks_around_precursor_mzParam(mz_tolerance = -1.2), "length 1")
    expect_error(remove_peaks_around_precursor_mzParam(foo = 1), "unused argument")
    
    res <- remove_peaks_around_precursor_mzParam(mz_tolerance = 17)
    expect_s4_class(res, "remove_peaks_around_precursor_mzParam")
    expect_equal(res@mz_tolerance, 17)
})

test_that("normalize_intensitiesParam constructor works", {
    ## errors
    expect_error(normalize_intensitiesParam(foo = 1), "unused argument")
    
    res <- normalize_intensitiesParam()
    expect_s4_class(res, "normalize_intensitiesParam")
})

test_that(".fun name works parameterized", {
    a <- select_by_intensityParam()
    expect_equal(SpectriPy:::.fun_name(a), "select_by_intensity")
    a <- select_by_mzParam()
    expect_equal(SpectriPy:::.fun_name(a), "select_by_mz")
    a <- remove_peaks_around_precursor_mzParam()
    expect_equal(SpectriPy:::.fun_name(a), "remove_peaks_around_precursor_mz")
    a <- normalize_intensitiesParam()
    expect_equal(SpectriPy:::.fun_name(a), "normalize_intensities")
})

test_that("python_command and ...._param_string work", {
    ## select_by_intensity
    a <- select_by_intensityParam(intensity_from = 0, intensity_to = 100)
    res <- .select_by_intensity_param_string(a)
    expect_equal(res, "intensity_from=0, intensity_to=100")
    res <- python_command(a)
    expect_match(res, "import matchms\nfrom matchms.filtering import ")
    expect_match(res, "matchms.filtering import select_by_intensity\nres = ")
    expect_match(res, "select_by_intensity")
    expect_match(res, "s, intensity_from=0, intensity_to=100) for s in py_spectrum_in]\n")
    res <- python_command(a, "py_x")
    expect_match(res, "m=0, intensity_to=100) for s in py_x]\n")
    
    ## select_by_mz
    a <- select_by_mzParam(mz_from = 0, mz_to = 400)
    res <- .select_by_mz_param_string(a)
    expect_equal(res, "mz_from=0, mz_to=400")
    res <- python_command(a)
    expect_match(res, "import matchms\nfrom matchms.filtering import ")
    expect_match(res, "matchms.filtering import select_by_mz\nres = ")
    expect_match(res, "select_by_mz")
    expect_match(res, "s, mz_from=0, mz_to=400) for s in py_spectrum_in]\n")
    res <- python_command(a, "py_x")
    expect_match(res, "m=0, mz_to=400) for s in py_x]\n")
    
    ## remove_peaks_around_precursor_mzParam
    a <- remove_peaks_around_precursor_mzParam(mz_tolerance = 17)
    res <- .remove_peaks_around_precursor_mz_param_string(a)
    expect_equal(res, "mz_tolerance=17")
    res <- python_command(a)
    expect_match(res, "import matchms\nfrom matchms.filtering import ")
    expect_match(res, "matchms.filtering import remove_peaks_around_precursor_mz\nres = ")
    expect_match(res, "remove_peaks_around_precursor_mz")
    expect_match(res, "s, mz_tolerance=17) for s in py_spectrum_in]\n")
    res <- python_command(a, "py_x")
    expect_match(res, "mz_tolerance=17) for s in py_x]\n")
    
    ## normalize_intensitiesParam
    a <- normalize_intensitiesParam()
    res <- .normalize_intensities_param_string(a)
    expect_equal(res, character())
    res <- python_command(a)
    expect_match(res, "import matchms\nfrom matchms.filtering import ")
    expect_match(res, "matchms.filtering import normalize_intensities\nres = ")
    expect_match(res, "normalize_intensities")
    expect_match(res, "s, ) for s in py_spectrum_in]\n")
    res <- python_command(a, "py_x")
    expect_match(res, ") for s in py_x]\n")
})

test_that("python commands evaluation", {
    vars <- c(precursorMz = "precursor_mz")
    
    ## invoke basilisk
    cl <- basiliskStart(matchms_env)
    
    ## select_by_intensity
    p <- select_by_intensityParam(intensity_from = 1000, intensity_to = 20000)
    pstring <- SpectriPy:::python_command(p)
    py$py_spectrum_in <- rspec_to_pyspec(caf, reference = import("matchms"),
        mapping = vars)
    py_run_string(pstring)
    expect_is(py$res, "list")
    res <- pyspec_to_rspec(py$res, mapping = vars)
    expect_equal(intensity(caf)[[1]], c(340, 416, 2580, 412), tolerance = 1e-04)
    expect_equal(intensity(caf)[[2]], c(388, 3270, 85, 54, 10111), tolerance = 1e-04)
    expect_equal(unname(intensity(res)[[1]]), 2580, tolerance = 1e-04)
    expect_equal(intensity(res)[[2]], c(3270, 10111), tolerance = 1e-04)
    
    ## select_by_mz
    p <- select_by_mzParam(mz_from = 0, mz_to = 150)
    pstring <- SpectriPy:::python_command(p)
    py$py_spectrum_in <- rspec_to_pyspec(caf, reference = import("matchms"),
        mapping = vars)
    py_run_string(pstring)
    expect_is(py$res, "list")
    res <- pyspec_to_rspec(py$res, mapping = vars)
    expect_equal(mz(caf)[[1]], c(135.0432, 138.0632, 163.0375, 195.088), tolerance = 1e-04)
    expect_equal(mz(caf)[[2]], c(110.071, 138.0655, 138.1057, 138.1742, 195.0864), tolerance = 1e-04)
    expect_equal(mz(res)[[1]], c(135.0432, 138.0632), tolerance = 1e-04)
    expect_equal(mz(res)[[2]], c(110.071, 138.0655, 138.1057, 138.1742), tolerance = 1e-04)
    
    ## remove_peaks_around_precursor_mz
    p <- remove_peaks_around_precursor_mzParam(mz_tolerance = 20)
    pstring <- SpectriPy:::python_command(p)
    py$py_spectrum_in <- rspec_to_pyspec(caf, reference = import("matchms"),
        mapping = vars)
    py_run_string(pstring)
    expect_is(py$res, "list")
    res <- pyspec_to_rspec(py$res, mapping = vars)
    expect_equal(mz(caf)[[1]], c(135.0432, 138.0632, 163.0375, 195.088), tolerance = 1e-04)
    expect_equal(mz(caf)[[2]], c(110.071, 138.0655, 138.1057, 138.1742, 195.0864), tolerance = 1e-04)
    expect_equal(mz(res)[[1]], c(135.0432, 138.0632, 163.0375), tolerance = 1e-04)
    expect_equal(mz(res)[[2]], c(110.071, 138.0655, 138.1057, 138.1742), tolerance = 1e-04)
    
    ## normalize_intensities
    p <- normalize_intensitiesParam()
    pstring <- SpectriPy:::python_command(p)
    py$py_spectrum_in <- rspec_to_pyspec(caf, reference = import("matchms"),
        mapping = vars)
    py_run_string(pstring)
    expect_is(py$res, "list")
    res <- pyspec_to_rspec(py$res, mapping = vars)
    expect_equal(intensity(caf)[[1]], c(340, 416, 2580, 412), tolerance = 1e-04)
    expect_equal(intensity(caf)[[2]], c( 388, 3270, 85, 54, 10111))
    expect_equal(intensity(res)[[1]], c(0.1317829, 0.1612403, 1.0, 0.1596899), tolerance = 1e-04)
    expect_equal(intensity(res)[[2]], c(0.038374048, 0.323410147, 0.008406686, 0.005340718, 1.0), tolerance = 1e-04)

    basiliskStop(cl)
})

test_that(".filter_spectra_python works", {
    res <- SpectriPy:::.filter_spectra_python(caf, select_by_intensityParam())
    expect_true(length(res) == 2)
    expect_is(res, "Spectra")

    ## try with empty Spectra
    res <- SpectriPy:::.filter_spectra_python(caf[integer()], select_by_intensityParam())
    expect_true(length(res) == 0)
    expect_is(res, "Spectra")
})

test_that("filterSpectriPy works", {
    
    ## check object caf before filtering
    expect_equal(mz(caf)[[1]], c(135.0432, 138.0632, 163.0375, 195.088), tolerance = 1e-04)
    expect_equal(mz(caf)[[2]], c(110.071, 138.0655, 138.1057, 138.1742, 195.0864), tolerance = 1e-04)
    expect_equal(intensity(caf)[[1]], c(340, 416, 2580, 412), tolerance = 1e-04)
    expect_equal(intensity(caf)[[2]], c( 388, 3270, 85, 54, 10111), tolerance = 1e-04)
    
    ## select_by_intensity
    res <- filterSpectriPy(caf, 
        param = select_by_intensityParam(intensity_from = 1000, intensity_to = 20000))
    expect_true(length(res) == 2)
    expect_is(res, "Spectra")
    expect_equal(unname(intensity(res)[[1]]), 2580, tolerance = 1e-04)
    expect_equal(intensity(res)[[2]], c(3270, 10111), tolerance = 1e-04)
    
    ## select_by_mz
    res <- filterSpectriPy(caf, 
        param = select_by_mzParam(mz_from = 0, mz_to = 150))
    expect_true(length(res) == 2)
    expect_is(res, "Spectra")
    expect_equal(mz(res)[[1]], c(135.0432, 138.0632), tolerance = 1e-04)
    expect_equal(mz(res)[[2]], c(110.071, 138.0655, 138.1057, 138.1742), tolerance = 1e-04)
    
    ## remove_peaks_around_precursor_mz
    res <- filterSpectriPy(caf, 
        param = remove_peaks_around_precursor_mzParam(mz_tolerance = 20))
    expect_true(length(res) == 2)
    expect_is(res, "Spectra")
    expect_equal(mz(res)[[1]], c(135.0432, 138.0632, 163.0375), tolerance = 1e-04)
    expect_equal(mz(res)[[2]], c(110.071, 138.0655, 138.1057, 138.1742), tolerance = 1e-04)
    
    ## normalize_intensities
    res <- filterSpectriPy(caf, 
        param = normalize_intensitiesParam())
    expect_true(length(res) == 2)
    expect_is(res, "Spectra")
    expect_equal(intensity(caf)[[1]], c(340, 416, 2580, 412), tolerance = 1e-04)
    expect_equal(intensity(caf)[[2]], c(388, 3270, 85, 54, 10111), tolerance = 1e-04)
    expect_equal(intensity(res)[[1]], c(0.1317829, 0.1612403, 1.0, 0.1596899), tolerance = 1e-04)
    expect_equal(intensity(res)[[2]], c(0.038374048, 0.323410147, 0.008406686, 0.005340718, 1.0), tolerance = 1e-04)
    
    ##caf_mod <- caf
    ##caf_mod$precursorMz[2] <- NA_real_
    ##expect_error(filterSpectriPy(caf_mod,
    ##    param = remove_peaks_around_precursor_mzParam()),
    ##    "Expect precursor to be positive")
})