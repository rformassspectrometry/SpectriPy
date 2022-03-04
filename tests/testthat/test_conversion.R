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

test_that(".rspec_to_pyspec works", {
    cl <- basiliskStart(SpectriPy:::matchms_env)
    vars <- spectraVariableMapping()
    basiliskRun(cl, function(x) {
        res <- .rspec_to_pyspec(sps[1L], vars)
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), sort(unname(vars)))
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## No metadata
    basiliskRun(cl, function(x) {
        res <- .rspec_to_pyspec(sps[1L], character())
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), character())
        expect_equal(as.numeric(res$peaks$intensities),
                     intensity(sps)[[1L]])
    }, x = sps)
    ## Only msLevel
    basiliskRun(cl, function(x) {
        res <- .rspec_to_pyspec(sps[1L], c(msLevel = "msLevel"))
        expect_equal(class(res)[1L], "matchms.Spectrum.Spectrum")
        expect_equal(sort(names(res$metadata)), "mslevel")
    }, x = sps)
    basiliskStop(cl)

})
