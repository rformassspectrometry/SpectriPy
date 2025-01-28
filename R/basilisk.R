#' @importFrom basilisk BasiliskEnvironment
matchms_env <- basilisk::BasiliskEnvironment(
                             envname = "matchms_env", pkgname = "SpectriPy",
                             packages = c("matchms==0.28.2"),
                             channels = c("bioconda", "conda-forge")
                         )
