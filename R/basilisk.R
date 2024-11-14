#' @importFrom basilisk BasiliskEnvironment
matchms_env <- basilisk::BasiliskEnvironment(
    envname = "matchms_env", pkgname = "SpectriPy",
    packages = c("matchms==0.27.0", "ms2deepscore==2.0.0"),
    channels = c("bioconda", "conda-forge")
)
