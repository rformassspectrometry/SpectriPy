#' @importFrom basilisk BasiliskEnvironment
matchms_env <- basilisk::BasiliskEnvironment(
    envname = "matchms_env", pkgname = "SpectriPy",
    packages = c("matchms==0.14.0"),
    channels = c("bioconda", "conda-forge")
)
