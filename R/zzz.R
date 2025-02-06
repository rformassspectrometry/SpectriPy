#' @importFrom reticulate import use_virtualenv
.onLoad <- function(libname, pkgname) {
    reticulate::use_virtualenv("r-spectripy")
    assign("matchms", import("matchms", delay_load = TRUE, convert = FALSE),
        envir = asNamespace(pkgname)
    )
    assign("matchms.similarity", import("matchms.similarity",
        delay_load = TRUE,
        convert = FALSE
    ),
    envir = asNamespace(pkgname)
    )
}
