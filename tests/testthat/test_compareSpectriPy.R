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

test_that("param constructors work", {
    ## Errors
    expect_error(CosineGreedy(tolerance = c(1.2, 1.5)), "length 1")
    expect_error(CosineGreedy(tolerance = -1.2), "length 1")
    expect_error(CosineGreedy(mz_power = c(1, 2)), "length 1")
    expect_error(CosineGreedy(intensity_power = c(1, 2)), "length 1")
    expect_error(NeutralLossesCosine(
        ignore_peaks_above_precursor = c(TRUE, FALSE)), "length 1")
})

test_that("CosineGreedy constructor works", {
    res <- CosineGreedy(tolerance = 5)
    expect_s4_class(res, "CosineGreedy")
    expect_equal(res@tolerance, 5)
})

test_that("CosineHungarian constructor works", {
    res <- CosineHungarian(intensity_power = 1.3)
    expect_s4_class(res, "CosineGreedy")
    expect_equal(res@intensityPower, 1.3)
})

test_that("ModifiedCosine constructor works", {
    res <- ModifiedCosine(mz_power = 4.3)
    expect_s4_class(res, "ModifiedCosine")
    expect_equal(res@mzPower, 4.3)
})

test_that("NeutralLossesCosine constructor works", {
    res <- NeutralLossesCosine(ignore_peaks_above_precursor = FALSE)
    expect_s4_class(res, "NeutralLossesCosine")
    expect_false(res@ignorePeaksAbovePrecursor)
})

test_that(".fun_name works", {
    a <- CosineGreedy()
    expect_equal(.fun_name(a), "CosineGreedy")
    a <- CosineHungarian()
    expect_equal(.fun_name(a), "CosineHungarian")
    a <- ModifiedCosine()
    expect_equal(.fun_name(a), "ModifiedCosine")
    a <- NeutralLossesCosine()
    expect_equal(.fun_name(a), "NeutralLossesCosine")
})

test_that("py_fun works", {
    res <- py_fun(CosineGreedy(tolerance = 0.5, mz_power = 0.3,
                                           intensity_power = 0.2))
    expect_equal(class(res)[1L],
                 "matchms.similarity.CosineGreedy.CosineGreedy")
    res <- py_fun(CosineHungarian(tolerance = 0.4, mz_power = 0.3,
                                           intensity_power = 0.2))
    expect_equal(
        class(res)[1L], "matchms.similarity.CosineHungarian.CosineHungarian")
    res <- py_fun(ModifiedCosine(tolerance = 0.1, mz_power = 0.3,
                                           intensity_power = 0.2))
    expect_equal(
        class(res)[1L], "matchms.similarity.ModifiedCosine.ModifiedCosine")
    res <- py_fun(NeutralLossesCosine(tolerance = 0.1, mz_power = 0.3,
                                                  intensity_power = 0.2,
                                                  ignore_peaks_above_precursor = FALSE))
    expect_equal(
        class(res)[1L],
        "matchms.similarity.NeutralLossesCosine.NeutralLossesCosine")
})

test_that(".compare_spectra_python works", {
    all <- c(caf, mhd)
    res <- .compare_spectra_python(all, caf, CosineGreedy())
    expect_true(nrow(res) == 4)
    expect_true(ncol(res) == 2)
    expect_equal(res[1, 1], 1)
    expect_equal(res[2, 2], 1)
    expect_true(all(res[3:4, ]== 0))

    res_sps <- compareSpectra(all, caf, tolerance = 0.1, ppm = 0)
    diffs <- abs(as.numeric(res - res_sps))
    expect_true(all(diffs < 0.01))

    ## only one spectra
    res <- .compare_spectra_python(all, param = CosineGreedy())
    res_2 <- .compare_spectra_python(all, all, param = CosineGreedy())
    expect_equal(res, res_2)

    ## try with empty Spectra
    res <- .compare_spectra_python(all, all[integer()], CosineGreedy())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == length(all))
    expect_true(ncol(res) == 0)

    res <- .compare_spectra_python(all[integer()], param = CosineGreedy())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == 0)

    res <- .compare_spectra_python(all[integer()], all, CosineGreedy())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == length(all))
})

test_that("compareSpectriPy works", {
    all <- c(caf, mhd)
    res_all <- compareSpectriPy(all, param = CosineGreedy())
    expect_true(is.numeric(res_all))
    expect_true(nrow(res_all) == length(all))
    expect_true(ncol(res_all) == length(all))
    expect_equal(diag(res_all), c(1, 1, 1, 1))

    ## CosineGreedy
    res <- compareSpectriPy(caf, mhd, param = CosineGreedy())
    expect_true(nrow(res) == 2)
    expect_true(ncol(res) == 2)
    expect_true(all(res == 0))
    res <- compareSpectriPy(caf, mhd, param = CosineGreedy(tolerance = 10))
    expect_true(any(res > 0))

    res <- compareSpectriPy(caf[1L], mhd[1L], param = CosineGreedy())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 1)
    expect_true(ncol(res) == 1)
    expect_true(all(res == 0))
    expect_equal(res[1, 1], res_all[1, 3])

    ## CosineHungarian
    res <- compareSpectriPy(caf, mhd, param = CosineHungarian())
    expect_true(nrow(res) == 2)
    expect_true(ncol(res) == 2)
    expect_true(all(res == 0))
    res <- compareSpectriPy(caf, mhd, param = CosineHungarian(tolerance = 10))
    expect_true(any(res > 0))

    ## ModifiedCosine
    res <- compareSpectriPy(caf, mhd, param = ModifiedCosine())
    expect_true(nrow(res) == 2)
    expect_true(ncol(res) == 2)
    expect_true(all(res > 0))
    res <- compareSpectriPy(caf, mhd, param = ModifiedCosine(tolerance = 10))
    expect_true(all(res > 0.3))

    all_mod <- all
    all_mod$precursorMz[3] <- NA_real_
    expect_error(compareSpectriPy(all_mod, all, param = ModifiedCosine()),
                 "Expect precursor to be positive")

    ## NeutralLossesCosine
    res <- compareSpectriPy(caf, mhd, param = NeutralLossesCosine())
    expect_true(nrow(res) == 2)
    expect_true(ncol(res) == 2)
    expect_true(all(res == 0))
    res <- compareSpectriPy(caf, mhd, param = NeutralLossesCosine(tolerance = 10))
    expect_true(all(res > 0))
})
