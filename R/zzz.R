#' @importFrom reticulate import use_virtualenv use_condaenv py_available py_install virtualenv_exists virtualenv_create virtualenv_remove conda_list conda_create
.onLoad <- function(libname, pkgname) {
    envname <- getOption(
        "spectripy.env",
        Sys.getenv("SPECTRIPY_ENV", unset = "r-spectripy")
    )
    use_conda <- getOption(
        "spectripy.use_conda",
        as.logical(Sys.getenv("SPECTRIPY_USE_CONDA",
            unset = "FALSE"
        ))
    )
    use_system <- getOption(
        "spectripy.use_system",
        as.logical(Sys.getenv("SPECTRIPY_USE_SYSTEM",
            unset = "FALSE"
        ))
    )

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
        envname = envname, use_conda = use_conda,
        use_system = use_system
    )

    assign("matchms", import("matchms", delay_load = TRUE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
    assign("matchms_sim",
        import("matchms.similarity", delay_load = TRUE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
}

#' @importFrom reticulate py_install py_module_available
.install_python_packages <- function(
    ..., envname = getOption(
        "spectripy.env",
        Sys.getenv("SPECTRIPY_ENV",
            unset = "r-spectripy"
        )
    ),
    use_conda = getOption(
        "spectripy.use_conda",
        as.logical(Sys.getenv("SPECTRIPY_USE_CONDA",
            unset = "FALSE"
        ))
    ),
    use_system = getOption(
        "spectripy.use_system",
        as.logical(Sys.getenv("SPECTRIPY_USE_SYSTEM",
            unset = "FALSE"
        ))
    )) {
    ## We don't want to modify the system Python, users are expected to manage
    ## dependencies themselves.
    if (use_system) {
        return()
    } else if (!py_module_available("matchms")) {
        if (use_conda)
            py_install("matchms", envname = envname, method = "conda",
                       pip = TRUE, ...)
        else
            py_install("matchms", envname = envname, method = "virtualenv", ...)
    }
}
