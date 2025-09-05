test_that("spectraVariableMapping,Spectra works", {
    expect_true(all(names(spectraVariableMapping(s)) %in%
                    names(coreSpectraVariables())))
    expect_true(length(spectraVariableMapping(s)) <
                length(coreSpectraVariables()))
    ## <-
    a <- s
    spectraVariableMapping(a) <- c(precursorMz = "precursor_mz",
                                   msLevel = "ms_level")
    expect_equal(a@backend@spectraVariableMapping,
                 c(precursorMz = "precursor_mz", msLevel = "ms_level"))
    expect_true(all(is.na(precursorCharge(a))))
})
