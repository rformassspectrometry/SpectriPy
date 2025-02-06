#' @importFrom reticulate import use_virtualenv
.onLoad <- function(libname, pkgname) {
    reticulate::use_virtualenv("r-spectripy")
    assign("matchms", import("matchms", delay_load = TRUE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
    assign("matchms_sim", import("matchms.similarity",
        delay_load = TRUE,
        convert = FALSE
    ),
    envir = asNamespace(pkgname)
    )
}

#' @importFrom reticulate py_install virtualenv_exists virtualenv_remove
install_python_packages <-
  function(..., envname = "r-spectripy",
           new_env = identical(envname, "r-spectripy")) {
    if (new_env && virtualenv_exists(envname)) {
      virtualenv_remove(envname)
    }

    py_install(packages = "matchms", envname = envname, ...)
  }
