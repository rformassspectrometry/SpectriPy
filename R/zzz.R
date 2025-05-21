## based on https://rstudio.github.io/reticulate/articles/package.html#delay-loading-python-modules
matchms <- NULL
matchms_similarity <- NULL
matchms_filtering <- NULL
spectrum_utils <- NULL

.PY_PKGS <- c(matchms = "matchms==0.28.2",
              spectrum_utils = "spectrum_utils==0.3.2",
              numpy = "numpy==2.0.2")

#' @importFrom reticulate import use_virtualenv use_condaenv
#'
#' @importFrom reticulate py_install virtualenv_exists virtualenv_create
#'
#' @importFrom reticulate conda_list conda_create install_miniconda()
.onLoad <- function(libname, pkgname) {
    envname <- .spectripy_env()
    use_conda <- .spectripy_use_conda()
    use_system <- .spectripy_use_system()
    if (use_system) {
        packageStartupMessage("Using system Python")
    } else {
        if (use_conda) {
            if (!.have_conda())
                install_miniconda()
            if (!(envname %in% conda_list()$name)) {
                packageStartupMessage("Creating conda environment '",
                                      envname, "'")
                conda_create(envname)
            }
            packageStartupMessage("Using conda environment '", envname, "'")
            use_condaenv(envname, required = TRUE)
            .py_check_install(pkgs = .PY_PKGS[c("matchms", "spectrum_utils")],
                              envname = envname, use_conda = use_conda)
        } else {
            if (!virtualenv_exists(envname)) {
                packageStartupMessage("Creating virtual environment '",
                                     envname, "'")
                virtualenv_create(envname, packages = FALSE)
            }
            packageStartupMessage("Using virtual environment '", envname, "'")
            use_virtualenv(envname, required = TRUE)
            .py_check_install(pkgs = .PY_PKGS[c("matchms", "spectrum_utils")],
                              envname = envname, use_conda = use_conda)
        }
    }
    if (!py_module_available("matchms"))
        warning("Required Python library 'matchms' not available!")
    matchms <<- import("matchms", delay_load = TRUE, convert = FALSE)
    matchms_similarity <<- import("matchms.similarity", delay_load = TRUE,
                                  convert = FALSE)
    matchms_filtering <<- import("matchms.filtering", delay_load = TRUE,
                                 convert = FALSE)
    spectrum_utils <<- import("spectrum_utils", delay_load = TRUE,
                              convert = FALSE)
}

.have_conda <- function() {
    !is.null(reticulate:::find_conda())
}

.spectripy_env <- function() {
    getOption(
        "spectripy.env", Sys.getenv("SPECTRIPY_ENV", unset = "r-spectripy"))
}

.spectripy_use_conda <- function() {
    as.logical(getOption(
        "spectripy.use_conda",
        Sys.getenv("SPECTRIPY_USE_CONDA", unset = "TRUE")))
}

.spectripy_use_system <- function() {
    as.logical(getOption(
        "spectripy.use_system",
        Sys.getenv("SPECTRIPY_USE_SYSTEM", unset = "FALSE")))
}

#' @importFrom reticulate py_module_available py_install
.py_check_install <- function(pkgs, envname = .spectripy_env(),
                              use_conda = .spectripy_use_conda()) {
    any_install <- FALSE
    for (pkg in pkgs) {
        if (!py_module_available(pkg)) {
            any_install <- TRUE
            message("Installing required library '", pkg,"'")
            if (use_conda)
                py_install(pkg, envname = envname, method = "conda",
                           pip = FALSE, channel = c("bioconda", "conda-forge"))
            else py_install(pkg, envname = envname, method = "virtualenv",
                            channel = c("bioconda", "conda-forge"))
        }
    }
    if (any_install)
        message("\nPlease restart R to load the freshly installed packages.\n")
}
