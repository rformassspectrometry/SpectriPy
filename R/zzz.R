## based on https://rstudio.github.io/reticulate/articles/package.html#delay-loading-python-modules
matchms <- NULL
matchms_similarity <- NULL
matchms_filtering <- NULL

#' @importFrom reticulate import use_virtualenv use_condaenv py_available py_install virtualenv_exists virtualenv_create virtualenv_remove conda_list conda_create
.onLoad <- function(libname, pkgname) {
    envname <- .spectripy_env()
    use_conda <- .spectripy_use_conda()
    use_system <- .spectripy_use_system()
    if (use_system) {
        packageStartupMessage("Using system Python")
    } else {
        if (use_conda) {
            if (!(envname %in% conda_list()$name)) {
                packageStartupMessage("Creating conda environment '",
                                      envname, "'")
                conda_create(envname)
            }
            packageStartupMessage("Using conda environment '", envname, "'")
            use_condaenv(envname, required = TRUE)
            .install_python_packages(envname, use_conda)
        } else {
            if (!virtualenv_exists(envname)) {
                packageStartupMessage("Creating virtual environment '",
                                     envname, "'")
                virtualenv_create(envname)
            }
            packageStartupMessage("Using virtual environment '", envname, "'")
            use_virtualenv(envname, required = TRUE)
            .install_python_packages(envname, use_conda)
        }
    }
    if (!py_module_available("matchms"))
        warning("Required Python library 'matchms' not available!")
    matchms <<- import("matchms", delay_load = TRUE, convert = FALSE)
    matchms_similarity <<- import("matchms.similarity", delay_load = TRUE,
                                  convert = FALSE)
    matchms_filtering <<- import("matchms.filtering", delay_load = TRUE,
                                 convert = FALSE)
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

#' @importFrom reticulate py_install py_module_available
.install_python_packages <- function(envname = .spectripy_env(),
                                     use_conda = .spectripy_use_conda(), ...) {
    if (!py_module_available("matchms")) {
        packageStartupMessage("Installing required libraries")
        if (use_conda) {
            py_install(c("matchms==0.28.2"),
                       envname = envname,
                       method = "conda",
                       channel = c("bioconda", "conda-forge"), ...)
        } else {
            py_install(c("matchms==0.28.2", "numpy==2.0.2"),
                       envname = envname, method = "virtualenv",
                       channel = c("conda-forge"), ...)
        }
        packageStartupMessage(
            "\nPlease restart R to load the freshly installed packages.\n")
    }
}
