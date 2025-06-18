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

mhd <- DataFrame(
    msLevel = c(2L, 2L),
    precursorMz = c(170.0924, 170.0924),
    id = c("HMDB0000001", "HMDB0000001"),
    name = c("1-Methylhistidine", "1-Methylhistidine"))
mhd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
mhd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
mhd <- Spectra(mhd)
s_all <- c(caf, mhd)

## select_by_intensity
test_that("select_by_intensity constructor works", {
    ## errors
    expect_error(select_by_intensity(intensity_from = c(1.2, 1.5)), "length 1")
    expect_error(select_by_intensity(intensity_from = -1.2), "length 1")
    expect_error(select_by_intensity(intensity_to = c(1, 2)), "length 1")
    expect_error(select_by_intensity(intensity_to = c(1, 2)), "length 1")
    expect_error(select_by_intensity(foo = 1), "unused argument")

    res <- select_by_intensity(intensity_from = 5, intensity_to = 20)
    expect_s4_class(res, "select_by_intensity")
    expect_equal(res@intensity_from, 5)
    expect_equal(res@intensity_to, 20)
})

test_that("py_fun,select_by_intensity works", {
    p <- select_by_intensity(intensity_from = 3, intensity_to = 24)
    res <- py_fun(p)
    expect_true(is(res, "functools.partial"))
    expect_equal(py_to_r(res$keywords["intensity_from"]), 3.0)
    expect_equal(py_to_r(res$keywords["intensity_to"]), 24.0)
})

test_that("filterSpectriPy,Spectra,select_by_intensity works", {
    p <- select_by_intensity(intensity_from = 17, intensity_to = 300)
    res <- filterSpectriPy(s_all, p)
    expect_s4_class(res, "Spectra")
    expect_equal(intensity(res)[[1L]], numeric())
    expect_true(all(intensity(res)[[2L]] > 17 & intensity(res)[[2L]] < 300))
    expect_true(all(intensity(res)[[3L]] > 17 & intensity(res)[[3L]] < 300))
    expect_true(all(intensity(res)[[4L]] > 17 & intensity(res)[[4L]] < 300))
    expect_equal(msLevel(s_all), msLevel(res))
    expect_equal(rtime(s_all), rtime(res))
})

## select_by_mz
test_that("select_by_mz constructor works", {
    ## errors
    expect_error(select_by_mz(mz_from = c(1.2, 1.5)), "length 1")
    expect_error(select_by_mz(mz_from = -1.2), "length 1")
    expect_error(select_by_mz(mz_to = c(1, 2)), "length 1")
    expect_error(select_by_mz(mz_to = c(1, 2)), "length 1")
    expect_error(select_by_mz(foo = 1), "unused argument")

    res <- select_by_mz(mz_from = 0, mz_to = 400)
    expect_s4_class(res, "select_by_mz")
    expect_equal(res@mz_from, 0)
    expect_equal(res@mz_to, 400)
})

test_that("py_fun,select_by_mz works", {
    p <- select_by_mz(mz_from = 3, mz_to = 24)
    res <- py_fun(p)
    expect_true(is(res, "functools.partial"))
    expect_equal(py_to_r(res$keywords["mz_from"]), 3.0)
    expect_equal(py_to_r(res$keywords["mz_to"]), 24.0)
})

test_that("filterSpectriPy,Spectra,select_by_mz works", {
    p <- select_by_mz(mz_from = 100, mz_to = 180)
    res <- filterSpectriPy(s_all, p)
    expect_s4_class(res, "Spectra")
    expect_true(all(mz(res)[[1L]] > 100 & mz(res)[[1L]] < 180))
    expect_true(all(mz(res)[[2L]] > 100 & mz(res)[[2L]] < 180))
    expect_true(all(mz(res)[[3L]] > 100 & mz(res)[[3L]] < 180))
    expect_true(all(mz(res)[[4L]] > 100 & mz(res)[[4L]] < 180))
    expect_equal(msLevel(s_all), msLevel(res))
    expect_equal(rtime(s_all), rtime(res))
})

## remove_peaks_around_precursor_mz
test_that("remove_peaks_around_precursor_mz constructor works", {
    ## errors
    expect_error(remove_peaks_around_precursor_mz(
        mz_tolerance = c(1.2, 1.5)), "length 1")
    expect_error(remove_peaks_around_precursor_mz(
        mz_tolerance = -1.2), "length 1")
    expect_error(remove_peaks_around_precursor_mz(foo = 1), "unused argument")

    res <- remove_peaks_around_precursor_mz(mz_tolerance = 17)
    expect_s4_class(res, "remove_peaks_around_precursor_mz")
    expect_equal(res@mz_tolerance, 17)
})

test_that("py_fun,remove_peaks_around_precursor_mz works", {
    p <- remove_peaks_around_precursor_mz(mz_tolerance = 200)
    res <- py_fun(p)
    expect_true(is(res, "functools.partial"))
    expect_equal(py_to_r(res$keywords["mz_tolerance"]), 200)
})

test_that("filterSpectriPy,Spectra,remove_peaks_around_precursor_mz works", {
    p <- remove_peaks_around_precursor_mz(mz_tolerance = 200)
    res <- filterSpectriPy(s_all, p)
    expect_s4_class(res, "Spectra")
    expect_true(all(lengths(res) == 0))
    expect_equal(msLevel(s_all), msLevel(res))
    expect_equal(rtime(s_all), rtime(res))

    p <- remove_peaks_around_precursor_mz(mz_tolerance = 10)
    res <- filterSpectriPy(s_all, p)
    expect_true(all(lengths(res) < lengths(s_all)))
})

## normalize_intensities
test_that("normalize_intensities constructor works", {
    ## errors
    expect_error(normalize_intensities(foo = 1), "unused argument")

    res <- normalize_intensities()
    expect_s4_class(res, "normalize_intensities")
})

## test_that("filterSpectriPy,normalize_intensitites works", {
##     p <- normalize_intensities()
##     res <- filterSpectriPy(s_all, p)
## })

test_that("py_fun,select_by_intensity works", {
    p <- select_by_mz(mz_from = 3, mz_to = 24)
    res <- py_fun(p)
    expect_true(is(res, "functools.partial"))
    expect_equal(py_to_r(res$keywords["mz_from"]), 3.0)
    expect_equal(py_to_r(res$keywords["mz_to"]), 24.0)
})

test_that("filterSpectriPy,Spectra,select_by_mz works", {
    p <- select_by_mz(mz_from = 100, mz_to = 180)
    res <- filterSpectriPy(s_all, p)
    expect_s4_class(res, "Spectra")
    expect_true(all(mz(res)[[1L]] > 100 & mz(res)[[1L]] < 180))
    expect_true(all(mz(res)[[2L]] > 100 & mz(res)[[2L]] < 180))
    expect_true(all(mz(res)[[3L]] > 100 & mz(res)[[3L]] < 180))
    expect_true(all(mz(res)[[4L]] > 100 & mz(res)[[4L]] < 180))
    expect_equal(msLevel(s_all), msLevel(res))
    expect_equal(rtime(s_all), rtime(res))
})

test_that(".filter_spectra_python works", {
    res <- .filter_spectra_python(Spectra())
    expect_s4_class(res, "Spectra")
    expect_true(length(res) == 0L)
})
