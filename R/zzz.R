## based on https://rstudio.github.io/reticulate/articles/package.html#delay-loading-python-modules
matchms <- NULL
matchms_similarity <- NULL
matchms_filtering <- NULL
spectrum_utils <- NULL

.PY_PKGS <- c(matchms = "matchms>=0.30.0",
              spectrum_utils = "spectrum_utils==0.3.2",
              numpy = "numpy>=2.2.0")

#' @importFrom reticulate py_require py_available
.onLoad <- function(libname, pkgname) {
    py_require(packages = .PY_PKGS, python_version = ">=3.12")
    .initialize_libraries2(FALSE, FALSE, asNamespace(pkgname))
    packageStartupMessage(print(py_config())) # temporarily
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
        stop("Failed to initialize Python environment and libraries.\n",
             "Original message:\n", e, "\nPython configuration:\n",
             print(py_config()), "\nEnvironmental variables:\n",
             print(Sys.getenv()))
    })
}
