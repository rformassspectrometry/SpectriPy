#' @importFrom basilisk BasiliskEnvironment
matchms_env <- BasiliskEnvironment(
    envname = "matchms_env", pkgname = "SpectriPy",
    packages = c("matchms==0.14.0"),
    channels = c("bioconda", "conda-forge")
    ## , pip = c("SpatialDE==1.1.3", "NaiveDE==1.2.0")
)
