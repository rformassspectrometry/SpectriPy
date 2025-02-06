#' @title A MS data backend for data stored in Python
#'
#' @description
#'
#' The `MsBackendPython` relies on the MS data being stored/represented in a
#' Python object.
#'
#' @details
#'
#' The `MsBackendPython` object references to the Python object with the MS
#' data purely through its variable name in the special `py` variable from
#' the *reticulate* package. It needs to be evaluated whether this will work
#' in all situations, e.g. also in parallel processing setups.
#'
#' @author Johannes Rainer and the EuBIC hackathon team
#'
#' @name MsBackendPython
NULL

setClass("MsBackendPython",
         contains = "MsBackend",
         slots = c(py_var = "character",
                   py_lib = "character"),
         prototype = prototype(py_var = character(),
                               py_lib = "matchms",
                               version = "0.1"))

MsBackendPython <- function() {
    new("MsBackendPython")
}

#' Helper to get the `py` variable from *reticulate*.
#'
#' @noRd
.get_py <- function() {
    tryCatch({
        get("py")
    }, error = function(e) {
        stop("Failed to get variable 'py'. Is the 'reticulate' package loaded?")
    })
}

#' @importFrom reticulate py_has_attr
.exists_py_var <- function(x) {
    py_has_attr(.get_py(), x)
}

.check_py_var <- function(x) {
    if (!.exists_py_var(x))
        stop("No variable of name \"", x, "\" found",
             " in the Python environment")
}

#' Need:
#'
#' - code to read MS data into python.
#' - is it possible to set column names in a python array?
#' - check: can we access an in-memory backend from Python? i.e. how to
#'   access e.g. a slot from a S4 class directly (the list of peaks matrices);
#'   but question is whether we need that or if it's better to pass the list
#'   of peaks matrices directly.
#'
#' Open questions:
#'
#' - how to deal with subsetting? keep an index in R and pass that to Python?
#' - directly subset the object in Python?
#'
#' @noRd
NULL

setMethod("backendInitialize", "MsBackendPython",
          function(object, pythonVariable = character(),
                   pythonLibrary = c("matchms"), data, ...) {
              pythonLibrary <- match.arg(pythonLibrary)
              object@py_lib <- pythonLibrary
              if (!length(pythonVariable))
                  stop("'pythonVariable' has to be provided")
              if (missing(data)) {
                  .check_py_var(pythonVariable)
                  object@py_var <- pythonVariable
              } else {
                  stop("not implemented yet")
                  ## if data provided -> convert the data to Python.
                  ## Overwriting variable...
              }
              object
          })

setMethod("show", "MsBackendPython", function(object) {
    cat(class(object), "with", length(object), "spectra\n")
})

#' @importFrom reticulate py_run_string
setMethod("length", "MsBackendPython", function(x) {
    if (length(x@py_var))
        py_run_string(
            paste0("_tmp_val = len(", x@py_var,")"))[["_tmp_val"]]
    else 0L
})

#' need to check how I can extract things from the python object.
#'
#' @noRd
setMethod(
    "spectraData", "MsBackendPython",
    function(object, columns = spectraVariables(object)) {
        stop("Not yet implememted")
    })

#' re-use the MsBackendCached instead?
setMethod("spectraVariables", "MsBackendPython", function(object) {
    stop("Not yet implemented")
    unique(c(names(coreSpectraVariables), "<get metadata from python>",
             "mz", "intensity"))
})

setMethod(
    "peaksData", "MsBackendPython", function(object,
                                             columns = c("mz", "intensity")) {
        if (length(object@py_var)) {
        } else list()
    })


#' Helper function to extract the peaks data from the Python object, assuming
#' it's a `list` of `matchms.Spectrum` objects
#'
#' @noRd
.peaks_data_matchms <- function(x) {
}

#' just random tests below...
a <- backendInitialize(MsBackendPython(), pythonVariable = "pysps")
