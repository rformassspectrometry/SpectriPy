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
    expect_error(CosineGreedyParam(tolerance = c(1.2, 1.5)), "length 1")
    expect_error(CosineGreedyParam(tolerance = -1.2), "length 1")
    expect_error(CosineGreedyParam(mzPower = c(1, 2)), "length 1")
    expect_error(CosineGreedyParam(intensityPower = c(1, 2)), "length 1")
    expect_error(NeutralLossesCosineParam(
        ignorePeaksAbovePrecursor = c(TRUE, FALSE)), "length 1")

    res <- CosineGreedyParam(tolerance = 5)
    expect_s4_class(res, "CosineGreedyParam")
    expect_equal(res@tolerance, 5)

    res <- CosineHungarianParam(intensityPower = 1.3)
    expect_s4_class(res, "CosineHungarianParam")
    expect_equal(res@intensityPower, 1.3)

    res <- ModifiedCosineParam(mzPower = 4.3)
    expect_s4_class(res, "ModifiedCosineParam")
    expect_equal(res@mzPower, 4.3)

    res <- NeutralLossesCosineParam(ignorePeaksAbovePrecursor = FALSE)
    expect_s4_class(res, "NeutralLossesCosineParam")
    expect_false(res@ignorePeaksAbovePrecursor)
})

test_that(".fun_name works", {
    a <- CosineHungarianParam()
    expect_equal(.fun_name(a), "CosineHungarian")
    expect_equal(.fun_name(CosineGreedyParam()), "CosineGreedy")
    expect_equal(.fun_name(ModifiedCosineParam()), "ModifiedCosine")
    expect_equal(.fun_name(NeutralLossesCosineParam()), "NeutralLossesCosine")
})

test_that("python_command and .cosine_param_string work", {
    a <- ModifiedCosineParam(tolerance = 0.9)
    res <- .cosine_param_string(a)
    expect_equal(res, "tolerance=0.9, mz_power=0, intensity_power=1")
    a <- CosineHungarianParam(mzPower = 4, intensityPower = 10)
    res <- .cosine_param_string(a)
    expect_equal(res, "tolerance=0.1, mz_power=4, intensity_power=10")
    res <- python_command(a)
    expect_match(res, "CosineHungarian")
    expect_match(res, "mz_power=4, intensity_power=10")
    expect_match(res, "py_x, py_y")
    expect_match(res, "is_symmetric=False")
    res <- python_command(a, "py_x", is_symmetric = "True")
    expect_match(res, "py_x, Co")
    expect_match(res, "is_symmetric=True")

    a <- NeutralLossesCosineParam(ignorePeaksAbovePrecursor = FALSE,
                                  tolerance = 0.9)
    res <- .cosine_param_string(a)
    expect_equal(res, "tolerance=0.9, mz_power=0, intensity_power=1")
    res <- python_command(a)
    expect_match(res, "ignore_peaks_above_precursor=False")

    a <- NeutralLossesCosineParam(ignorePeaksAbovePrecursor = TRUE,
                                  tolerance = 0.9)
    res <- .cosine_param_string(a)
    expect_equal(res, "tolerance=0.9, mz_power=0, intensity_power=1")
    res <- python_command(a)
    expect_match(res, "ignore_peaks_above_precursor=True")
})

test_that(".compare_spectra_python works", {
    all <- c(caf, mhd)
    res <- .compare_spectra_python(all, caf, CosineGreedyParam())
    expect_true(nrow(res) == 4)
    expect_true(ncol(res) == 2)

    res_sps <- compareSpectra(all, caf)
    diffs <- abs(as.numeric(res - res_sps))
    expect_true(all(diffs < 0.01))

    ## only one spectra
    res <- .compare_spectra_python(all, param = CosineGreedyParam())
    res_2 <- .compare_spectra_python(all, all, param = CosineGreedyParam())
    expect_equal(res, res_2)

    ## try with empty Spectra
    res <- .compare_spectra_python(all, all[integer()], CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == length(all))
    expect_true(ncol(res) == 0)

    res <- .compare_spectra_python(all[integer()], param = CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == 0)

    res <- .compare_spectra_python(all[integer()], all, CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == length(all))
})
