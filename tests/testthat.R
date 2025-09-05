library(testthat)
library(Spectra)
library(SpectriPy)
library(MsBackendMgf)

fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
s <- Spectra(fl, source = MsBackendMgf())
s <- setBackend(s, MsBackendPy(), pythonVariableName = "py_be")

be <- s@backend

test_check("SpectriPy")

## TODO: add tests from the Spectra unit test suite!

## ## Run the MsBackend spectra variable test suite
## test_suite <- system.file("test_backends", "test_MsBackend",
##                           package = "Spectra")

## ## Run single test file.
## res <- test_dir(test_suite, stop_on_failure = TRUE)
