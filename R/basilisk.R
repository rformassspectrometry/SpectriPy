#' @importFrom basilisk BasiliskEnvironment
matchms_env <- basilisk::BasiliskEnvironment(
    envname = "matchms_env", pkgname = "SpectriPy",
    packages = c("matchms==0.28.2"),
    channels = c("bioconda", "conda-forge")
)

#' @importFrom basilisk BasiliskEnvironment
spectrum_utils_env <- basilisk::BasiliskEnvironment(
    envname = "spectrum_utils_env", pkgname = "SpectriPy",
    packages = c("spectrum_utils==0.3.2-0"),
    channels = c("bioconda", "conda-forge")
)
