#' @importFrom reticulate import use_virtualenv use_condaenv py_available py_install virtualenv_exists virtualenv_create virtualenv_remove conda_list conda_create
.onLoad <- function(libname, pkgname) {
    envname <- .spectripy_env()
    use_conda <- .spectripy_use_conda()
    use_system <- .spectripy_use_system()
    if (use_conda) {
        if (!(envname %in% conda_list()$name)) {
            conda_create(envname)
        }
        use_condaenv(envname, required = TRUE)
    } else if (!use_system) {
        if (!virtualenv_exists(envname)) {
            virtualenv_create(envname)
        }
        use_virtualenv(envname)
    }
    .install_python_packages(
        envname = envname, use_conda = use_conda, use_system = use_system
    )
    assign("matchms", import("matchms", delay_load = FALSE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
    assign("matchms_similarity",
        import("matchms.similarity", delay_load = FALSE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
    assign("matchms_filtering",
        import("matchms.filtering", delay_load = FALSE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
}

.spectripy_env <- function() {
    getOption(
        "spectripy.env", Sys.getenv("SPECTRIPY_ENV", unset = "r-spectripy"))
}

.spectripy_use_conda <- function() {
    getOption(
        "spectripy.use_conda",
        as.logical(Sys.getenv("SPECTRIPY_USE_CONDA", unset = "TRUE")))
}

.spectripy_use_system <- function() {
    getOption(
        "spectripy.use_system",
        as.logical(Sys.getenv("SPECTRIPY_USE_SYSTEM", unset = "FALSE")))
}

#' @importFrom reticulate py_install py_module_available
.install_python_packages <- function(..., envname = .spectripy_env(),
                                     use_conda = .spectripy_use_conda(),
                                     use_system = .spectripy_use_system()) {
    ## We don't want to modify the system Python, users are expected to manage
    ## dependencies themselves.
    if (use_system) {
        return()
    } else if (!py_module_available("matchms")) {
        if (use_conda) {
            py_install(c("matchms==0.28.2"), envname = envname,
                       method = "conda", pip = TRUE,
                       channels = c("conda-forge"), ...)
        } else {
            ## Somehow an old version of numpy gets installed and installation
            ## fails in the end.
            py_install(c("matchms==0.28.2"),
                       envname = envname, method = "virtualenv",
                       channels = c("conda-forge"), ...)
        }
    }
}
