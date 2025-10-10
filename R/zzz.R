## based on https://rstudio.github.io/reticulate/articles/package.html#delay-loading-python-modules
matchms <- NULL
matchms_similarity <- NULL
matchms_filtering <- NULL
spectrum_utils <- NULL

.PY_PKGS <- c(matchms = "matchms>=0.31",
              spectrum_utils = "spectrum_utils>=0.3.2",
              numpy = "numpy>=2.2.0")

#' @importFrom reticulate py_require py_available
.onLoad <- function(libname, pkgname) {
    .check_environment()
    py_require(packages = .PY_PKGS, python_version = ">=3.12")
    .initialize_libraries2(FALSE, FALSE, asNamespace(pkgname))
}

.check_environment <- function() {
    ## Check for old setup
    if (.is_spectripy_use_system())
        warning("Ignoring environment variable 'SPECTRIPY_USE_SYSTEM'. ",
                "Please use 'RETICULATE_PYTHON' or 'RETICULATE_PYTHON_ENV'",
                " instead.")
    if (.is_spectripy_use_conda())
        warning("Environment variable 'SPECTRIPY_USE_CONDA' is no longer ",
                "supported. Please see the package vignette for updated ",
                "Python setup options.")
    se <- .spectripy_env()
    if (length(se) && se != "")
        warning("Environment variable 'SPECTRIPY_ENV' is no longer ",
                "supported. Please use 'RETICULATE_PYTHON_ENV' instead. See ",
                "the package vignette for updated Python setup options.")
    ## Potentially interfering global settings
    if ((res <- Sys.getenv("RETICULATE_PYTHON")) != "")
        packageStartupMessage("Using Python defined by 'RETICULATE_PYTHON': ",
                              res)
    if ((res <- Sys.getenv("RETICULATE_PYTHON_ENV")) != "")
        packageStartupMessage("Using Python environment defined by ",
                              "'RETICULATE_PYTHON_ENV': ", res)
}

.is_spectripy_use_system <- function() {
    res <- getOption("spectripy.use_system", default = NULL)
    if (!length(res))
        res <- Sys.getenv("SPECTRIPY_USE_SYSTEM")
    grepl("1|yes|true", res, ignore.case = TRUE)
}

.is_spectripy_use_conda <- function() {
    res <- getOption("spectripy.use_conda", default = NULL)
    if (!length(res))
        res <- Sys.getenv("SPECTRIPY_USE_CONDA")
    grepl("1|yes|true", res, ignore.case = TRUE)
}

.spectripy_env <- function() {
    if(length(res <- getOption("spectripy.env", default = NULL)))
        return(res)
    Sys.getenv("SPECTRIPY_ENV")
}

#' Load all required Python libraries and assign it to package-internal
#' variables
#'
#' @importFrom reticulate import py_config
#'
#' @noRd
.initialize_libraries2 <- function(delay_load = TRUE, convert = FALSE,
                                   envir = new.env()) {
    if (!reticulate::py_available(initialize = TRUE))
        stop("Unable to initialize the Python environment", call. = FALSE)
    tryCatch({
        assign("matchms", import("matchms", delay_load = delay_load,
                                 convert = convert), envir = envir)
        assign("matchms_similarity",
               import("matchms.similarity", delay_load = delay_load,
                      convert = convert), envir = envir)
        assign("matchms_filtering",
               import("matchms.filtering", delay_load = delay_load,
                      convert = convert), envir = envir)
        assign("spectrum_utils",
               import("spectrum_utils", delay_load = delay_load,
                      convert = convert), envir = envir)
    }, error = function(e) {
        stop("Failed to initialize Python environment and libraries!\n",
             "Original message: ", e, "\nPython configuration:\n",
             print(py_config()), "\nEnvironmental variables:\n",
             print(Sys.getenv()))
    })
}
