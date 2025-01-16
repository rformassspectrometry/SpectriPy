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
})

test_that("CosineGreedyParam constructor works", {
    res <- CosineGreedyParam(tolerance = 5)
    expect_s4_class(res, "CosineGreedyParam")
    expect_equal(res@tolerance, 5)
})

test_that("CosineHungarianParam constructor works", {
    res <- CosineHungarianParam(intensityPower = 1.3)
    expect_s4_class(res, "CosineGreedyParam")
    expect_equal(res@intensityPower, 1.3)
})

test_that("ModifiedCosineParam constructor works", {
    res <- ModifiedCosineParam(mzPower = 4.3)
    expect_s4_class(res, "ModifiedCosineParam")
    expect_equal(res@mzPower, 4.3)
})

test_that("NeutralLossesCosineParam constructor works", {
    res <- NeutralLossesCosineParam(ignorePeaksAbovePrecursor = FALSE)
    expect_s4_class(res, "NeutralLossesCosineParam")
    expect_false(res@ignorePeaksAbovePrecursor)
})

test_that(".fun name works parameterized", {
    a <- CosineGreedyParam()
    expect_equal(SpectriPy:::.fun_name(a), "CosineGreedy")
    a <- CosineHungarianParam()
    expect_equal(SpectriPy:::.fun_name(a), "CosineHungarian")
    a <- ModifiedCosineParam()
    expect_equal(SpectriPy:::.fun_name(a), "ModifiedCosine")
    a <- NeutralLossesCosineParam()
    expect_equal(SpectriPy:::.fun_name(a), "NeutralLossesCosine")
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

test_that("python commands evaluation", {
    p <- CosineGreedyParam(tolerance = 0.05)
    pstring <- SpectriPy:::python_command(p)
    ## Invoke basilisk
    cl <- basiliskStart(matchms_env)
    py$py_x <- rspec_to_pyspec(caf, reference = import("matchms"),
                               mapping = c(precursorMz = "precursor_mz"))
    py$py_y <- rspec_to_pyspec(caf, reference = import("matchms"),
                               mapping = c(precursorMz = "precursor_mz"))
    py_run_string(pstring)
    ## Convert to similarity array
    py_run_string("sim = res.scores.to_array()")
    py$sim["CosineGreedy_score"]

    p <- CosineHungarianParam()
    pstring <- SpectriPy:::python_command(p)
    py_run_string(pstring)
    py_run_string("sim = res.scores.to_array()")
    py$sim["CosineHungarian_score"]

    p <- ModifiedCosineParam()
    pstring <- SpectriPy:::python_command(p)
    py_run_string(pstring)
    py_run_string("sim = res.scores.to_array()")
    py$sim["ModifiedCosine_score"]

    p <- SpectriPy:::NeutralLossesCosineParam()
    pstring <- SpectriPy:::python_command(p)
    py_run_string(pstring)
    py_run_string("sim = res.scores.to_array()[\"NeutralLossesCosine_score\"]")

    basiliskStop(cl)
})

test_that(".compare_spectra_python works", {
    all <- c(caf, mhd)
    res <- SpectriPy:::.compare_spectra_python(all, caf, CosineGreedyParam())
    expect_true(nrow(res) == 4)
    expect_true(ncol(res) == 2)

    res_sps <- compareSpectra(all, caf)
    diffs <- abs(as.numeric(res - res_sps))
    expect_true(all(diffs < 0.01))

    ## only one spectra
    res <- SpectriPy:::.compare_spectra_python(all, param = CosineGreedyParam())
    res_2 <- SpectriPy:::.compare_spectra_python(all, all, param = CosineGreedyParam())
    expect_equal(res, res_2)

    ## try with empty Spectra
    res <- SpectriPy:::.compare_spectra_python(all, all[integer()], CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == length(all))
    expect_true(ncol(res) == 0)

    res <- SpectriPy:::.compare_spectra_python(all[integer()], param = CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == 0)

    res <- SpectriPy:::.compare_spectra_python(all[integer()], all, CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 0)
    expect_true(ncol(res) == length(all))
})

test_that("compareSpectriPy works", {
    all <- c(caf, mhd)
    res_all <- compareSpectriPy(all, param = CosineGreedyParam())
    expect_true(is.numeric(res_all))
    expect_true(nrow(res_all) == length(all))
    expect_true(ncol(res_all) == length(all))
    expect_equal(diag(res_all), c(1, 1, 1, 1))

    ## Test python call. LLLLLLL
    ## cl <- basiliskStart(SpectriPy:::matchms_env)
    ## p <- CosineGreedyParam()
    ## cmd <- SpectriPy:::python_command(p)
    ## py$py_x <- rspec_to_pyspec(caf, mapping = c(precursorMz = "precursor_mz"))
    ## py$py_y <- rspec_to_pyspec(mhd, mapping = c(precursorMz = "precursor_mz"))
    ## strng <- paste0("import matchms\n",
    ##                 "from matchms.similarity import CosineGreedy\n",
    ##                 "res = matchms.calculate_scores(py_x, py_y, CosineGreedy(), is_symmetric = False)\n")
    ## py_run_string(strng)
    ## basiliskStop(cl)
    ## res <- compareSpectriPy(caf, mhd, param = CosineGreedyParam())
    ## WHY is this not working??? Seems to be some python issue, maybe a bug
    ## in CosineGreedy?
    ##
    ## WORKS res <- compareSpectriPy(caf, caf, param = CosineGreedyParam())
    ## WORKS res <- compareSpectriPy(mhd, mhd, param = CosineGreedyParam())
    ## NOT res <- compareSpectriPy(mhd, caf, param = CosineGreedyParam())
    ## NOT res <- compareSpectriPy(caf, mhd, param = CosineGreedyParam())
    ## Maybe the issue is with spectra without any similarity (score = 0) and
    ## having a result matrix with more than 1 column or row.
    ## Seems to be the case:
    res <- compareSpectriPy(caf, mhd, param = CosineGreedyParam(tolerance = 10))

    res <- compareSpectriPy(caf[1L], mhd[1L], param = CosineGreedyParam())
    expect_true(is.numeric(res))
    expect_true(nrow(res) == 1)
    expect_true(ncol(res) == 1)
    expect_true(all(res == 0))
    expect_equal(res[1, 1], res_all[1, 3])

    res <- compareSpectriPy(caf, c(mhd, caf[1]), param = CosineGreedyParam())
    expect_true(nrow(res) == 2L)
    expect_true(ncol(res) == 3L)
    expect_equal(res[1, 3], 1)
    expect_equal(res[2, 3], res_all[1, 2])

    ## Add tests after matchms > 0.14.0
    ## res <- compareSpectriPy(all, all, param = NeutralLossesCosineParam())

    ## ModifiedCosine with and without precursor m/z
    res <- compareSpectriPy(all, all, param = ModifiedCosineParam())
    expect_true(nrow(res) == length(all))
    expect_true(ncol(res) == length(all))
    expect_true(all(res > 0))
    expect_equal(diag(res), c(1, 1, 1, 1))

    all_mod <- all
    all_mod$precursorMz[3] <- NA_real_
    expect_error(compareSpectriPy(all_mod, all, param = ModifiedCosineParam()),
                 "Expect precursor to be positive")
})