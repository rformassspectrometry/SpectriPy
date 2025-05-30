#' @include conversion.R

#' @title A MS data backend for MS data stored in Python
#'
#' @description
#'
#' The `MsBackendPy` allows to access MS data stored as `matchms.Spectrum`
#' or `spectrum_utils.spectrum.MsmsSpectrum` objects from the
#' [*matchms*](https://github.com/matchms/matchms) respectively
#' [*spectrum_utils*](https://github.com/bittremieux-lab/spectrum_utils) Python
#' library directly from R. The MS data (peaks data or spectra variables) are
#' translated on-the-fly when accessed. Thus, the `MsBackendPy` allows a
#' seamless integration of Python MS data structures into [Spectra::Spectra()]
#' based analysis workflows.
#'
#' The `MsBackendPy` object is considered *read-only*, i.e. it does not provide
#' functionality to replace the peaks data from R. However, it is possible to
#' directly change the data in the referenced Python variable.
#'
#' @details
#'
#' The `MsBackendPy` keeps only a reference to the MS data in Python (i.e. the
#' name of the variable in Python) as well as an index pointing to the
#' individual spectra in Python but no other data. Any data requested from
#' the `MsBackendPy` is accessed and translated on-the-fly from the Python
#' variable. The `MsBackendPy` is thus an interface to the MS data, but not
#' a data container. All changes to the MS data in the Python variable
#' (performed e.g. in Python) immediately affect any `MsBackendPy` instances
#' pointing to this variable.
#'
#' Special care must be taken if the MS data structure in Python is subset or
#' its order is changed (e.g. by another process). In that case it might be
#' needed to re-index the backend using the `reindex()` function:
#' `object <- reindex(object)`. This will update (replace) the index to the
#' individual spectra in Python which is stored within the backend.
#'
#' @note
#'
#' As mentioned in the *details* section the MS data is completely stored in
#' Python and the backend only references to this data through the name of
#' the variable in Python. Thus, each time MS data is requested from the
#' backend, it is retrieved in its **current** state.
#' If for example data was transformed or metadata added or removed in the
#' Python object, it immediately affects the `Spectra`/backend.
#'
#' @section `MsBackendPy` methods:
#'
#' The `MsBackendPy` supports all methods defined by the [Spectra::MsBackend()]
#' interface for access to MS data. Details on the invidual functions can also
#' be found in the main documentation in the *Spectra* package (i.e. for
#' [Spectra::MsBackend()]). Here we provide information for functions with
#' specific properties of the backend.
#'
#' - `backendInitialize()`: initializes the backend with information from the
#'   referenced Python variable (attribute). The name of this attribute,
#'   ideally stored in the associated Python session, is expected to be
#'   provided with the `pythonVariableName` parameter. The optional
#'   `spectraVariableMapping` parameter allows to provide additional, or
#'   alternative, mapping of `Spectra`'s *spectra variables* to metadata in the
#'   `matchms.Spectrum` objects. See [defaultSpectraVariableMapping()] (the
#'   default) for more information. Parameter `pythonLibrary` must be used
#'   to specify the Python library representing the MS data in Python. It can
#'   be either `pythonLibrary = "matchms"` (the default) or
#'   `pythonLibrary = "spectrum_utils"`. The function returns an initialized
#'   instance of `MsBackendPy`.
#'
#' - `peaksData()`: extracts the peaks data matrices from the backend. Python
#'   code is applied to the data structure in Python to
#'   extract the *m/z* and intensity values as a list of (numpy) arrays. These
#'   are then translated into an R `list` of two-column `numeric` matrices.
#'   Because Python does not allow to name columns of an array, an additional
#'   loop in R is required to set the column names to `"mz"` and `"intensity"`.
#'
#' - `spectraData()`: extracts the spectra data from the backend. Which spectra
#'   variables are translated and retrieved from the Python objects depends on
#'   the backend's `spectraVariableMapping()`. All metadata names defined are
#'   retrieved and added to the returned `DataFrame` (with eventually missing
#'   *core* spectra variables filled with `NA`).
#'
#' - `spectraVariables()`: retrieves available spectra variables, which include
#'   the names of all metadata attributes in the `matchms.Spectrum` objects
#'   and the *core* spectra variables [Spectra::coreSpectraVariables()].
#'
#' - `spectraVariableMapping<-`: replaces the `spectraVariableMapping` of the
#'   backend (see [setSpectraVariableMapping()] for details and description
#'   of the expected format).
#'
#' @section Additional helper and utility functions:
#'
#' - `reindex()`: update the internal *index* to match `1:length(object)`.
#'   This function is useful if the original data referenced by the backend was
#'   subset or re-ordered by a different process (or a function in Python).
#'
#' @param columns For `spectraData()`: `character` with the names of
#'     columns (spectra variables) to retrieve. Defaults to
#'     `spectraVariables(object)`. For `peaksData()`: `character` with the
#'     names of the peaks variables to retrieve.
#'
#' @param data For `backendInitialize()`: `DataFrame` with the full MS data
#'     (peaks data and spectra data). Currently not supported.
#'
#' @param drop For `spectraData()` and `peaksData()`: `logical(1)` whether,
#'     when a single column is requested, the data should be returned as a
#'     `vector` instead of a `data.frame` or `matrix`.
#'
#' @param name For `$`: `character(1)` with the name of the variable to
#'     retrieve.
#'
#' @param object A `MsBackendPy` object.
#'
#' @param pythonLibrary For `backendInitialize()`: `character(1)` specifying
#'     the Python library used to represent the MS data in Python. Can be
#'     either `pythonLibrary = "matchms"` (the default) or
#'     `pythonLibrary = "spectrum_utils"`.
#'
#' @param pythonVariableName For `backendInitialize()`: `character(1)` with the
#'     name of the variable/Python attribute that contains the list of
#'     `matchms.Spectrum` objects with the MS data.
#'
#' @param spectraVariableMapping For `backendInitialize()`: named `character`
#'     with the mapping between spectra variable names and (`matchms.Spectrum`)
#'     metadata names. See [defaultSpectraVariableMapping()] for more
#'     information and details.
#'
#' @param value Replacement value(s).
#'
#' @param x A `MsBackendPy` object
#'
#' @param ... Additional parameters.
#'
#' @author Johannes Rainer and the EuBIC hackathon team
#'
#' @name MsBackendPy
#'
#' @return
#'
#' See description of individual functions for their return values.
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
#' ## named "s_p" in the (*reticulate*'s) `py` Python environment. Assigning
#' ## the variable to the Python environment has performance advantages, as
#' ## any Python code applied to the MS data does not require any data
#' ## conversions.
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
#' ## Available spectra variables: these include, next to the *core* spectra
#' ## variables, also the names of all metadata stored in the `matchms.Spectrum`
#' ## objects.
#' spectraVariables(s_2)
#'
#' ## Get the full peaks data:
#' peaksData(s_2)
#'
#' ## Get the peaks from the first spectrum
#' peaksData(s_2)[[1L]]
#'
#' ## Get the full spectra data:
#' spectraData(s_2)
#'
#' ## Get the m/z values
#' mz(s_2)
#'
#' ## Plot the first spectrum
#' plotSpectra(s_2[1L])
#'
#'
#' ########
#' ## Using the spectrum_utils Python library
#'
#' ## Below we convert the data to a list of `MsmsSpectrum` object from the
#' ## spectrum_utils library.
#' py_set_attr(py, "su_p", rspec_to_pyspec(s,
#'     spectraVariableMapping("spectrum_utils"), "spectrum_utils"))
#'
#' ## Create a MsBackendPy representing this data. Importantly, we need to
#' ## specify the Python library using the `pythonLibrary` parameter and
#' ## ideally also set the `spectraVariableMapping` to the one specific for
#' ## that library.
#' be <- backendInitialize(MsBackendPy(), "su_p",
#'     spectraVariableMapping = spectraVariableMapping("spectrum_utils"),
#'     pythonLibrary = "spectrum_utils")
#' be
#'
#' ## Get the peaks data for the first 3 spectra
#' peaksData(be[1:3])
#'
#' ## Get the full spectraData
#' spectraData(be)
#'
#' ## Extract the precursor m/z
#' be$precursorMz
NULL

#' @noRd
#'
#' @slot py_var `character(1)` with the **name** of the variable. The variable
#'     could be defined in R or in Python.
#'
#' @slot py_lib `character(1)` defining the Python library representing the
#'     data. Currently not used/defaults to `matchms`.
#'
#' @slot is_in_py `logical(1)` whether the variable is in Python (py$) or in
#'     the R environment.
#'
#' @slot i `integer` with the integer/index of the spectra from the backend
#'     we should represent. This allows to subset using `[` **without**
#'     changing the original data.
setClass("MsBackendPy",
         contains = "MsBackend",
         slots = c(py_var = "character",
                   py_lib = "character",
                   is_in_py = "logical",
                   i = "integer",
                   spectraVariableMapping = "character"),
         prototype = prototype(
             py_var = character(),
             py_lib = "matchms",
             is_in_py = TRUE,
             i = integer(),
             spectraVariableMapping = defaultSpectraVariableMapping(),
             version = "0.1"))

#' @export
MsBackendPy <- function() {
    new("MsBackendPy")
}

#' @importMethodsFrom ProtGenerics backendInitialize
#'
#' @rdname MsBackendPy
setMethod("backendInitialize", "MsBackendPy",
          function(object, pythonVariableName = character(),
                   spectraVariableMapping = defaultSpectraVariableMapping(),
                   pythonLibrary = c("matchms", "spectrum_utils"),
                   ...,
                   data) {
              .check_spectra_variable_mapping(spectraVariableMapping)
              object@spectraVariableMapping <- spectraVariableMapping
              pythonLibrary <- match.arg(pythonLibrary)
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
                  object <- reindex(object)
                  object@py_lib <- pythonLibrary
              } else {
                  stop("Support for 'data' parameter is not yet implemented")
                  ## if data provided -> convert the data to Python.
              }
              object
          })

#' @importMethodsFrom methods show
setMethod("show", "MsBackendPy", function(object) {
    .check_i(object)
    l <- length(object)
    cat(class(object), "with", l, "spectra\n")
    cat("Data stored in the \"", object@py_var,"\" variable ",
        ifelse(object@is_in_py, "in Python", "in R"), "\n", sep = "")
})

#' @rdname MsBackendPy
setMethod("length", "MsBackendPy", function(x) {
    length(x@i)
})

#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @importFrom Spectra coreSpectraVariables
#'
#' @rdname MsBackendPy
setMethod("spectraVariables", "MsBackendPy", function(object) {
    if (length(object)) {
        var <- .py_get_metadata_names(object)
        i <- match(var, object@spectraVariableMapping)
        if (length(i))
            var[!is.na(i)] <- names(object@spectraVariableMapping)[i[!is.na(i)]]
    } else var <- character()
    union(names(coreSpectraVariables()), var)
})

#' @importMethodsFrom ProtGenerics spectraData
#'
#' @importFrom Spectra fillCoreSpectraVariables coreSpectraVariables
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
            return(as(fillCoreSpectraVariables(
                data.frame()), "DataFrame")[, columns, drop = drop])
        mt <- .py_get_metadata_names(object)
        mt <- mt[names(mt) %in% columns]
        if (length(mt)) {
            SFUN <- switch(
                object@py_lib,
                matchms = .py_matchms_spectrum_spectra_data,
                spectrum_utils = .py_spectrum_utils_spectrum_spectra_data)
            ## Note: creating a pandas.DataFrame in Python is not faster
            d <- .get_py_var(object@py_var, object@is_in_py)
            if (!is.list(d)) d <- py_to_r(d)
            res <- do.call(
                rbindFill,
                lapply(d[object@i], SFUN, mt))
            ## Ensure correct data type for core variables.
            csv <- coreSpectraVariables()
            for (mtc in intersect(colnames(res), names(csv)))
                res[[mtc]] <- as(res[[mtc]], csv[mtc])
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
            stop("Column(s) ",
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
            res <- switch(
                object@py_lib,
                matchms = .py_matchms_peaks_data(
                    object@py_var, (object@i - 1L)),
                spectrum_utils = .py_spectrum_utils_peaks_data(
                    object@py_var, (object@i - 1L)))
            if (identical(as.vector(columns), c("mz", "intensity"))) {
                res <- lapply(res, function(z) {
                    colnames(z) <- c("mz", "intensity")
                    z
                })
            } else {
                res <- lapply(res, function(z) {
                    colnames(z) <- c("mz", "intensity")
                    z[, columns, drop = drop]
                })
            }
            res
        } else list()
    })

#' @importMethodsFrom ProtGenerics extractByIndex
setMethod("extractByIndex", c("MsBackendPy", "ANY"), function(object, i) {
    object@i <- object@i[i]
    .check_i(object)
    object
})

#' @rdname MsBackendPy
setMethod("$", "MsBackendPy", function(x, name) {
    spectraData(x, name, drop = TRUE)
})

setMethod("lengths", "MsBackendPy", function(x, use.names = FALSE) {
    vapply(peaksData(x), nrow, NA_integer_)
})

#' @importMethodsFrom S4Vectors isEmpty
setMethod("isEmpty", "MsBackendPy", function(x) {
    lengths(x) == 0L
})

#' @importMethodsFrom ProtGenerics acquisitionNum
setMethod("acquisitionNum", "MsBackendPy", function(object) {
    spectraData(object, "acquisitionNum", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics centroided
setMethod("centroided", "MsBackendPy", function(object) {
    spectraData(object, "centroided", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics collisionEnergy
setMethod("collisionEnergy", "MsBackendPy", function(object) {
    spectraData(object, "collisionEnergy", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics dataOrigin
setMethod("dataOrigin", "MsBackendPy", function(object) {
    spectraData(object, "dataOrigin", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics intensity
setMethod("intensity", "MsBackendPy", function(object) {
    NumericList(peaksData(object, "intensity", drop = TRUE), compress = FALSE)
})

#' @importMethodsFrom ProtGenerics isolationWindowLowerMz
setMethod("isolationWindowLowerMz", "MsBackendPy", function(object) {
    spectraData(object, "isolationWindowLowerMz", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics isolationWindowTargetMz
setMethod("isolationWindowTargetMz", "MsBackendPy", function(object) {
    spectraData(object, "isolationWindowTargetMz", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics isolationWindowUpperMz
setMethod("isolationWindowUpperMz", "MsBackendPy", function(object) {
    spectraData(object, "isolationWindowUpperMz", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics msLevel
setMethod("msLevel", "MsBackendPy", function(object) {
    spectraData(object, "msLevel", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics mz
setMethod("mz", "MsBackendPy", function(object) {
    NumericList(peaksData(object, "mz", drop = TRUE), compress = FALSE)
})

#' @importMethodsFrom ProtGenerics polarity
setMethod("polarity", "MsBackendPy", function(object) {
    spectraData(object, "polarity", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics precScanNum
setMethod("precScanNum", "MsBackendPy", function(object) {
    spectraData(object, "precScanNum", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics precursorCharge
setMethod("precursorCharge", "MsBackendPy", function(object) {
    spectraData(object, "precursorCharge", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics precursorIntensity
setMethod("precursorIntensity", "MsBackendPy", function(object) {
    spectraData(object, "precursorIntensity", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics precursorMz
setMethod("precursorMz", "MsBackendPy", function(object) {
    spectraData(object, "precursorMz", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics rtime
setMethod("rtime", "MsBackendPy", function(object) {
    res <- spectraData(object, "rtime", drop = TRUE)
    l <- length(object)
    if (length(res) != l)
        rep(NA_real_, l)
    else res
})

#' @importMethodsFrom ProtGenerics scanIndex
setMethod("scanIndex", "MsBackendPy", function(object) {
    spectraData(object, "scanIndex", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics smoothed
setMethod("smoothed", "MsBackendPy", function(object) {
    spectraData(object, "smoothed", drop = TRUE)
})

#' @importMethodsFrom ProtGenerics spectraNames
setMethod("spectraNames", "MsBackendPy", function(object) {
    rownames(spectraData(object))
})

#' @importMethodsFrom ProtGenerics tic
setMethod("tic", "MsBackendPy", function(object, initial = TRUE) {
    if (initial) {
        if (any(spectraVariables(object) == "totIonCurrent"))
            spectraData(object, "totIonCurrent")[, 1L]
        else rep(NA_real_, length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

#' @importMethodsFrom Spectra spectraVariableMapping<-
#'
#' @rdname MsBackendPy
setReplaceMethod(
    "spectraVariableMapping", "MsBackendPy", function(object, value) {
        .check_spectra_variable_mapping(value)
        object@spectraVariableMapping <- value
        object
    })

#' @rdname MsBackendPy
setReplaceMethod("spectraVariableMapping", "Spectra", function(object, value) {
    spectraVariableMapping(object@backend) <- value
    object
})

#' @rdname MsBackendPy
#'
#' @export
reindex <- function(object) {
    object@i <- seq_len(.py_var_length(object))
    object
}

################################################################################
##        HELPER FUNCTIONS
################################################################################
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

#' Check that all indices are within the range of the data in Python.
#'
#' @param object `MsBackendPy` object.
#'
#' @return `TRUE` if all is OK - otherwise it throws an error
#'
#' @noRd
.check_i <- function(object) {
    l <- .py_var_length(object)
    if (l > 0 && max(object@i) > l)
        stop("Indices are out of bound. Use `reindex()` to update indices ",
             "of 'object' if \"", object@py_var, "\" was subset by ",
             "another process.", call. = FALSE)
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

#' @importFrom reticulate py_run_string py_eval
.py_var_length <- function(x) {
    if (length(x@py_var)) {
        if (x@is_in_py) py_eval(paste0("len(", x@py_var, ")"))
        else length(.get_py_var(x@py_var, FALSE))
    } else 0L
}

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
        switch(
            x@py_lib,
            matchms = {
                mt <- names(var[0L]$metadata)
                names(mt) <- names(
                    x@spectraVariableMapping)[match(
                                      mt, x@spectraVariableMapping)]
                nas <- is.na(names(mt))
                names(mt)[nas] <- mt[nas]
                mt
            },
            spectrum_utils = {
                mt <- c("identifier", "precursor_mz",
                        "precursor_charge", "retention_time")
                mt <- x@spectraVariableMapping[x@spectraVariableMapping %in% mt]
                mt
            })
    }
    else character()
}
