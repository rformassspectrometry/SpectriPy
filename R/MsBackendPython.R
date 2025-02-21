#' @include conversion.R

#' @title A MS data backend for MS data stored in Python
#'
#' @description
#'
#' The `MsBackendPy` allows to access MS data stored as `matchms.Spectrum`
#' objects from the [matchms](https://github.com/matchms/matchms) Python
#' library directly from R. The MS data (peaks data or spectra variables) get
#' translated on-the-fly when accessed. Thus, the `MsBackendPy` allows a
#' seamless integration of Python MS data structures into [Spectra()] based
#' analysis workflows.
#'
#' The `MsBackendPy` object is considered *read-only*, i.e. it does not provide
#' functionality to replace the peaks data from R. However, it is possible to
#' directly change the data in the referenced Python variable.
#'
#' @details
#'
#' The `MsBackendPy` keeps only a reference to the MS data in Python (i.e. the
#' name of the variable in Python) but no other data. Any data requested from
#' the `MsBackendPy` is accessed and translated on-the-fly from the Python
#' variable. The `MsBackendPy` is thus an interface to the MS data, but does
#' not contain any data itself. Because of this also all changes done to the
#' data in Python (which inlcudes also a subset operation performed using `[`
#' on the backend in R!) would immediately affect any `MsBackendPy` instances
#' pointing to the same Python variable.
#'
#' @note
#'
#' As mentioned in the *details* section the MS data is completely stored in
#' Python and the backend only references to this data through the name of
#' the variable in Python. Thus, each time MS data is requested from the
#' backend, it is retrieved in its **current** state.
#' If for example data was transformed or metadata added or removed in the
#' Python object, it would immediately also affect the backend.
#'
#' @author Johannes Rainer and the EuBIC hackathon team
#'
#' @name MsBackendPy
#'
#'
#' @examples
#'
#' ## Loading an example MGF file provided by the SpectriPy package.
#' ## As an alternative, the data could also be imported directly in Python
#' ## using:
#' ## import matchms
#' ## from matchms.importing import load_from_mgf
#' ## s_p = list(load_from_mgf(r.fl))
#' library(Spectra)
#' library(MsBackendMgf)
#'
#' fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
#' s <- Spectra(fl, source = MsBackendMgf())
#' s
#'
#' ## Translating the MS data to Python and assigning it to a variable
#' ## named "s_p" in the (*reticulate*'s) `py` Python environment.
#' py_set_attr(py, "s_p", rspec_to_pyspec(s))
#'
#'
#' ## Create a `MsBackendPy` representing an interface to the data in the
#' ## "s_p" variable in Python:
#' be <- backendInitialize(MsBackendPy(), "s_p")
#' be
#'
#' ## Create a Spectra object which this backend:
#' s_2 <- Spectra(be)
#' s_2
#'
#' ## Available spectraVariables:
#' spectraVariables(s_2)
#'
#' ## Get the full peaks data:
#' peaksData(s_2)
#'
#' ## Get the full spectra data:
#' spectraData(s_2)
#'
#' ## Get the m/z values
#' mz(s_2)
#'
#' ## Plot the first spectrum
#' plotSpectra(s_2[1L])
setClass("MsBackendPy",
         contains = "MsBackend",
         slots = c(py_var = "character",
                   py_lib = "character",
                   is_in_py = "logical",
                   spectraVariableMapping = "character"),
         prototype = prototype(
             py_var = character(),
             py_lib = "matchms",
             is_in_py = TRUE,
             spectraVariableMapping = defaultSpectraVariableMapping(),
             version = "0.1"))

#' @export
MsBackendPy <- function() {
    new("MsBackendPy")
}

.check_spectra_variable_mapping <- function(x) {
    if (length(x) && !length(names(x)))
        stop("'spectraVariableMapping' needs to be a named character vector")
    else TRUE
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

#' Checks if the variable with the name `x` **exists** in Python (if
#' `is_in_py = TRUE`) or in R (if `is_in_py = FALSE`). It throws and error
#' if there is none.
#'
#' @noRd
.check_py_var_exists <- function(x, is_in_py = TRUE) {
    if ((is_in_py && !.exists_py_var(x)) |
        (!is_in_py && !length(get0(x))))
        stop("No variable of name \"", x, "\" found",
             " in the Python or R environment.")
    TRUE
}

.check_py_var <- function(x, is_in_py = TRUE) {
    var <- .get_py_var(x, is_in_py)
    if (!is(var, "python.builtin.list"))
        stop("Variable \"", x, "\" is supposed to be a Python list, ",
             "but it is a ", class(x)[1L])
        TRUE
}

#' Helper function to get the variable with the Python data. Can be either
#' in R or in Python.
#'
#' @param x `character(1)` with the name of the variable.
#'
#' @param is_in_py `logical(1)` whether the variable is stored in the `py`
#'     module (`is_in_py = TRUE`, the default) or in the R environment
#'     (`is_in_py = FALSE`).
#'
#' @return the variable for which then name was provided.
#'
#' @noRd
#'
#' @importFrom reticulate py_get_attr
.get_py_var <- function(x, is_in_py = TRUE) {
    if (is_in_py) {
        py_get_attr(.get_py(), x)
    } else get0(x)
}

#' @importMethodsFrom ProtGenerics backendInitialize
#'
#' @rdname MsBackendPy
setMethod("backendInitialize", "MsBackendPy",
          function(object, pythonVariableName = character(),
                   spectraVariableMapping = defaultSpectraVariableMapping(),
                   ...,
                   data) {
              .check_spectra_variable_mapping(spectraVariableMapping)
              object@spectraVariableMapping <- spectraVariableMapping
              if (!length(pythonVariableName))
                  stop("'pythonVariableName' has to be provided")
              if (!is.character(pythonVariableName))
                  stop("'pythonVariableName' is expected to be the name ",
                       "of the variable with the data.")
              if (missing(data)) {
                  ## 1) Check if the variable exists - and whether it's in py or
                  ##    in R. Set the `is_in_py` variable depending on that.
                  if (.exists_py_var(pythonVariableName))
                      object@is_in_py <- TRUE
                  else
                      object@is_in_py <- FALSE
                  .check_py_var_exists(pythonVariableName, object@is_in_py)
                  object@py_var <- pythonVariableName
                  ## 2) Check if the variable is of the correct data type.
                  .check_py_var(pythonVariableName, object@is_in_py)
                  object@py_var <- pythonVariableName
              } else {
                  stop("not yet implemented")
                  ## if data provided -> convert the data to Python.
                  ## Overwriting variable...
              }
              object
          })

#' @importMethodsFrom methods show
setMethod("show", "MsBackendPy", function(object) {
    cat(class(object), "with", length(object), "spectra\n")
    cat("Data stored in the \"", object@py_var,"\" variable ",
        ifelse(object@is_in_py, "in Python", "in R"), "\n", sep = "")
})

#' @importFrom reticulate py_run_string py_eval
#'
#' @rdname MsBackendPy
setMethod("length", "MsBackendPy", function(x) {
    if (length(x@py_var)) {
        if (x@is_in_py) py_eval(paste0("len(", x@py_var, ")"))
        else length(.get_py_var(x@py_var, FALSE))
    } else 0L
})

#' Get the names of the metadata fields from the referenced Python object.
#' Returns a named vector, elements being the metadata names in Python and
#' names the respective spectra variable names, i.e. the names that should
#' be used in R.
#'
#' @param x `MsBackendPy`
#'
#' @noRd
.py_get_metadata_names <- function(x) {
    var <- .get_py_var(x@py_var, x@is_in_py)
    if (length(var)) {
        mt <- names(var[0L]$metadata)
        names(mt) <- names(
            x@spectraVariableMapping)[match(mt, x@spectraVariableMapping)]
        nas <- is.na(names(mt))
        names(mt)[nas] <- mt[nas]
        mt
    }
    else character()
}

#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @importFrom Spectra coreSpectraVariables
#'
#' @rdname MsBackendPy
setMethod("spectraVariables", "MsBackendPy", function(object) {
    var <- .py_get_metadata_names(object)
    i <- match(var, object@spectraVariableMapping)
    if (length(i))
        var[!is.na(i)] <- names(object@spectraVariableMapping)[i[!is.na(i)]]
    union(names(coreSpectraVariables()), var)
})

#' @importMethodsFrom ProtGenerics spectraData
#'
#' @importFrom Spectra fillCoreSpectraVariables
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom methods as
#'
#' @rdname MsBackendPy
setMethod(
    "spectraData", "MsBackendPy",
    function(object, columns = spectraVariables(object), drop = FALSE) {
        if (!length(object))
            return(as(fillCoreSpectraVariables(data.frame()), "DataFrame"))
        mt <- .py_get_metadata_names(object)
        mt <- mt[names(mt) %in% columns]
        if (length(mt)) {
            FUN <- function(x) .py_matchms_spectrum_spectra_data(x, mt)
            res <- do.call(rbindFill,
                           iterate(.get_py_var(object@py_var, object@is_in_py),
                                   FUN, simplify = FALSE))
        } else
            res <- as.data.frame(matrix(ncol = 0, nrow = length(object)))
        n <- names(coreSpectraVariables())
        res <- fillCoreSpectraVariables(res, n[!n %in% c("mz", "intensity")])
        if (any(columns == "mz")) {
            res <- as(res, "DataFrame")
            res$mz <- NumericList(
                peaksData(object, "mz", TRUE), compress = FALSE)
        }
        if (any(columns == "intensity")) {
            res <- as(res, "DataFrame")
            res$intensity <- NumericList(
                peaksData(object, "intensity", TRUE), compress = FALSE)
        }
        if (any(columns == "dataStorage"))
            res$dataStorage <- object@py_var
        if (!all(columns %in% colnames(res)))
            stop("Columns ",
                 paste0("\"", columns[!columns %in% colnames(res)], "\"",
                        collapse = ", "), " not available.", call. = FALSE)
        if (drop && length(columns) == 1)
            res[, columns, drop = TRUE]
        else as(res[, columns, drop = FALSE], "DataFrame")
    })

#' @importMethodsFrom ProtGenerics peaksData
#'
#' @importFrom reticulate iterate
#'
#' @rdname MsBackendPy
setMethod(
    "peaksData", "MsBackendPy", function(object,
                                         columns = c("mz", "intensity"),
                                         drop = FALSE) {
        if (length(object@py_var)) {
            if (identical(as.vector(columns), c("mz", "intensity"))) {
                res <- iterate(.get_py_var(object@py_var, object@is_in_py),
                               .py_matchms_spectrum_peaks_data,
                               simplify = FALSE)
            } else {
                FUN <- function(x)
                    .py_matchms_spectrum_peaks_data_columns(
                        x, columns, drop = drop)
                res <- iterate(.get_py_var(object@py_var, object@is_in_py),
                               FUN, simplify = FALSE)
            }
            res
        } else list()
    })



#' TODO:
#'
#' - backendInitialize providing data -> convert to Python and store data there
#'   using py_set_attr(py, ...)
#'
#' Open questions:
#'
#' - how to deal with subsetting? keep an index in R and pass that to Python?
#' - directly subset the object in Python?
#' - -> show bold message in help that the backend keeps only a reference to the
#'   data, backend does not keep track of any changes in the data of that
#'   variable. the backend is by design only an interface to the data in Python,
#'   that reads data, whether changed or not by another process to R.
#'
#' ## Required methods:
#'
#' - [ ] `dataStorage()`
#' - [X] `length()`
#' - [X] `backendInitialize()`
#' - [X] `spectraVariables()`
#' - [X] `spectraData()`
#' - [X] `peaksData()`
#' - [ ] `extractByIndex()`
#' - [ ] `[`
#' - [ ] `backendMerge()`
#' - [ ] `$`
#' - [ ] `lengths()`
#' - [ ] `isEmpty()`
#' - [ ] `acquisitionNum()`
#' - [ ] `centroided()`
#' - [ ] `collisionEnergy()`
#' - [ ] `dataOrigin()`
#' - [ ] `intensity()`
#' - [ ] `isolationWindowLowerMz()`
#' - [ ] `isolationWindowTargetMz()`
#' - [ ] `isolationWindowUpperMz()`
#' - [ ] `msLevel()`
#' - [ ] `mz()`
#' - [ ] `polarity()`
#' - [ ] `precScanNum()`
#' - [ ] `precursorCharge()`
#' - [ ] `precursorIntensity()`
#' - [ ] `precursorMz()`
#' - [ ] `rtime()`
#' - [ ] `scanIndex()`
#' - [ ] `smoothed()`
#' - [ ] `spectraNames()`
#' - [ ] `tic()`
#' @noRd
NULL
