#' @title Converting between R and Python MS data structures
#'
#' @name conversion
#'
#' @description
#'
#' The `rspec_to_pyspec()`, `pyspec_to_rspec()` and `r_to_py.Spectra()`
#' functions allow to convert (translate) MS data structures between R and
#' Python. At present the R [Spectra::Spectra()] objects are translated into
#' a list of [matchms](https://github.com/matchms/matchms) Python
#' `matchms.Spectrum` objects.
#' The `r_to_py()` and `py_to_r()` methods are provided to integrate
#' better with the R/Python bindings of the *reticulate* package, but the
#' dedicated `rspec_to_pyspec()` and `pyspec_to_rspec()` allow more
#' customization and are more efficient (specifically for translating from
#' Python to R).
#'
#' The mapping of spectra variables (in R) to
#' (Python) spectra metadata can be configured and defined with the
#' `setSpectraVariableMapping()` and `spectraVariableMapping()`.
#'
#' @note
#'
#' The `py_ro_r()` method allows only to translate a **single**
#' `matchms.Spectrum` object into a [Spectra::Spectra] object of length 1.
#' To translate a list of `matchms.Spectrum` objects one would need to iterate
#' over the list and then concatenate the resulting `list` of `Spectra` would thus require#'
#' See the indivudual function's documentation for more details.
#'
#'
#' @section Translation of MS data objects:
#'
#' MS data structures can be translated between R and Python using the
#' `rspec_to_pyspec()`, `r_to_py()`, `pyspec_to_rspec()` and `py_to_r()`
#' functions.
#'
#' @section Mapping of spectra variables (metadata):
#'
#' Metadata for MS spectra are represented and stored as *spectra variables*
#' in the R [Spectra::Spectra()] objects. Also Python MS data structures
#' store such metadata along with the mass peak data. While spectra metadata
#' is thus supported by data structures in both programming languages,
#' different names and naming conventions are used. The
#' `spectraVariableMapping()` and `setSpectraVariableMapping()` functions allow
#' to define how the names of spectra metadata (spectra variables) should be
#' translated between R and Python. The `r_to_py()` and `py_to_r()` functions
#' will used these to name the spectra variables accordingly. Also, only
#' spectra metadata/variables in `spectraVariableMapping()` will be translated.
#' The initial mapping is based on this
#' [definition in matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv).
#'
#' - `spectraVariableMapping()`: returns the currenctly defined spectra
#'   variable mapping as a named character vector, with names representing the
#'   names of the spectra variables in R and elements the respective names
#'   of the spectra metadata in Python. Use [Spectra::spectraVariables()] on
#'   the `Spectra` object that should be converted with `r_to_py()` to list
#'   all available spectra variables. `r_to_py()` and `py_to_r()` for MS data
#'   structures will use this mapping.
#'
#' - `setSpectraVariableMapping()`: sets/replaces the currently defined mapping
#'   of spectra variable names to Python metadata names. Setting
#'   `setSpectraVariableMapping(character())` will only convert the mass peaks
#'   data (m/z and intensity values) but no spectra metadata.
#'
#' @param x For `r_to_py.Spectra()`: `Spectra` object. For `pyspec_to_rspec()`:
#'     a Python list of `matchms.Spectrum` objects.
#'
#' @return For `r_to_py.Spectra()`: Python array of Spectrum objects, same
#'     length as `x`. For `pyspec_to_rspec()`: [Spectra::Spectra()] with the
#'     converted spectra.
#'
#' @author Michael Witting, Johannes Rainer, Wout Bittremieux
#'
#' @importFrom reticulate r_to_py py_to_r
#'
#' @importMethodsFrom Spectra spectrapply
NULL

.SPECTRA_2_MATCHMS <- c(
    precursorMz = "precursor_mz",
    precursorIntensity = "precursor_intensity",
    precursorCharge = "charge",
    rtime = "retention_time",
    collisionEnergy = "collision_energy",
    isolationWindowTargetMz = "isolation_window_target_mz",
    ## polarity = "ionmode", # Disabling since matchms does not support int.
    msLevel = "ms_level"
)

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @exportMethod spectraVariableMapping
#'
#' @rdname conversion
#'
#' @export
setMethod("spectraVariableMapping", "missing", function(object, ...) {
    getOption("spectripy.spectra_variable_mapping", .SPECTRA_2_MATCHMS)
})

#' @rdname conversion
#'
#' @export
setSpectraVariableMapping <- function(x) {
    if (!is.character(x) | length(names(x)) != length(x))
        stop("'x' is expected to be a named character vector")
    options(spectripy.spectra_variable_mapping = x)
}

## -------- R TO PY ----------------------------------------------------------##

#' @rdname conversion
#'
#' @description
#'
#' Function to convert R Spectra objects into a Python list of matchms Spectrum
#' objects using the `reticulate` package.
#'
#' @param x `Spectra` object.
#'
#' @param convert Boolean; should Python objects be automatically converted to
#' their R equivalent? Defaults to `FALSE`.
#'
#' @importFrom reticulate r_to_py
#'
#' @importFrom Spectra spectrapply
#'
#' @export
r_to_py.Spectra <- function(x, convert = FALSE) {
    .rspec_to_pyspec(x, mapping = spectraVariableMapping())
}

#' @rdname conversion
#'
#' @importFrom methods is
#'
#' @export
rspec_to_pyspec <- function(x, mapping = spectraVariableMapping(),
                            .check = TRUE) {
    if (.check && !is(x, "Spectra"))
        stop("'x' should be a Spectra object.")
    ## that could be more memory efficient, but slower.
    ## r_to_py(spectrapply(x, .single_rspec_to_pyspec,
    ##                     spectraVariables = mapping,
    ##                     BPPARAM = BPPARAM))
    .rspec_to_pyspec(x, mapping = mapping)
}

#' @description
#'
#' Function to convert a **single** R Spectra object (of length 1) into a
#' Python matchms Spectrum using the `reticulate` package.
#'
#' @param x `Spectra` object **of length 1!**.
#'
#' @param spectraVariables Named `character` vector defining the spectra
#'     variables that should be stored as metadata in `matchms`' metadata. Names
#'     are expected to be the spectra variable names and values the
#'     corresponding metadata fields in `matchms`. Defaults to
#'     `.SPECTRA_2_MATCHMS`. If `spectraVariables = character()` no
#'     metadata will be stored.
#'
#' @return `Spectrum` Single Python Spectrum.
#'
#' @author Michael Witting, Johannes Rainer, Wout Bittremieux
#'
#' @importMethodsFrom Spectra spectraData
#'
#' @importMethodsFrom Spectra peaksData
#'
#' @importMethodsFrom Spectra mz
#'
#' @importMethodsFrom Spectra intensity
#'
#' @importFrom reticulate np_array r_to_py
#'
#' @noRd
.single_rspec_to_pyspec <-
    function(x, spectraVariables = spectraVariableMapping()) {
        pks <- unname(peaksData(x, c("mz", "intensity")))[[1L]]
        if (length(spectraVariables)) {
            slist <- as.list(spectraData(x, columns = names(spectraVariables)))
            names(slist) <- spectraVariables[names(slist)]
            matchms$Spectrum(
                        mz = np_array(pks[, 1L]),
                        intensities = np_array(pks[, 2L]),
                        metadata = r_to_py(slist)
                    )
        } else {
            matchms$Spectrum(
                        mz = np_array(pks[, 1L]),
                        intensities = np_array(pks[, 2L])
                    )
        }
    }

#' Converts a `Spectra::Spectra` to a `list` of `matchms.Spectrum` by first
#' extracting the peaks and spectra data and iterating over these. This is
#' a faster, but also more memory heavy implementation as the full peaks and
#' spectra data are read into memory. With the `single_rspec_to_pyspec()` only
#' the data from one spectrum at a time are read.
#'
#' @noRd
.rspec_to_pyspec <- function(x, spectraVariables = spectraVariableMapping()) {
    pks <- peaksData(x, c("mz", "intensity"))
    sv <- spectraVariables[!spectraVariables %in% c("mz", "intensity")]
    if (length(sv)) {
        spd <- as.data.frame(spectraData(x, columns = names(sv)))
        colnames(spd) <- sv[colnames(spd)]
        rm(x)
        spd <- split(spd, seq_along(pks))
        r_to_py(mapply(function(y, z) {
            matchms$Spectrum(mz = np_array(z[, 1L]),
                             intensities = np_array(z[, 2L]),
                             metadata = r_to_py(as.list(y)))
        }, spd, pks, SIMPLIFY = FALSE, USE.NAMES = FALSE))
    } else {
        rm(x)
        r_to_py(lapply(pks, function(z)
            matchms$Spectrum(mz = np_array(z[, 1L],),
                             intensities = np_array(z[, 2L]))))
    }
}

## -------- PY TO R ----------------------------------------------------------##

#' @rdname conversion
#'
#' @importFrom Spectra concatenateSpectra
#'
#' @export
pyspec_to_rspec <- function(x, mapping = spectraVariableMapping(),
        BPPARAM = SerialParam(), .check = TRUE) {
    if (!(is(x, "list") | is(x, "python.builtin.list")))
      stop("'x' is expected to be a Python list.")
    x <- py_to_r(x)
    if (.check && !all(vapply(x, function(z)
        is(z, "matchms.Spectrum.Spectrum"), logical(1))))
        stop("'x' is expected to be a Python list of matchms.Spectrum objects.")
    spectra_list <- bplapply(x, .single_pyspec_to_rspec,
                             spectraVariables = mapping, BPPARAM = BPPARAM)
    do.call(concatenateSpectra, spectra_list)
}

#' @export
pyspec_to_rspec2 <- function(x, mapping = spectraVariableMapping(), ...) {
    if (is(x, "python.builtin.list"))
        x <- py_to_r(x)
    ## Get the list of peak matrices (in python)
    be <- MsBackendMemory()
    be@peaksData <- lapply(x, function(z) {
        m <- py_to_r(z$peaks$to_numpy)
        colnames(m) <- c("mz", "intensity")
        m
    })
    ## Not very efficient and elegant... get the indivudal elements and stuff
    ## into list. pandas.DataFrame can not be easily created unfortunately.
    be@spectraData <- do.call(rbind, lapply(x, function(z) {
        pl <- z$metadata
        map <- mapping[mapping %in% names(pl)]
        if (length(map))
            as.data.frame(lapply(map, function(y) pl[[y]]))
        else data.frame(ms_level = NA_integer_)
    }))
    Spectra(be)
}

.py_matchms_spectrum_spectra_data <-
    function(x, mapping = spectraVariableMapping()) {
        pl <- x$metadata
        map <- mapping[mapping %in% names(pl)]
        if (length(map))
            as.data.frame(lapply(map, function(z) pl[[y]]))
        else data.frame(msLevel = NA_integer_)    }

#' @export
pyspec_to_rspec3 <- function(x, mapping = spectraVariableMapping()) {
    ## Assign x to a python variable that we can access... Maybe we could even
    ## try to use the name of the variable instead???
    py$tmp_ <- x
    be <- MsBackendMemory()
    be@spectraData <- .py_matchms_spectra_data("tmp_", mapping = mapping)
    be@peaksData <- .py_matchms_peaks_data("tmp_")
    py_del_attr(py, "tmp_")
    Spectra(be)
}




#' @description
#'
#' Function to convert a single Python Spectrum object into an R Spectra using
#' the `reticulate` package.
#'
#' @param x `Spectrum` Single Python Spectrum.
#'
#' @param spectraVariables named `character` vector with the names of the
#'     spectra variables that should be extracted. Names are expected to be
#'     the spectra variable names (in `Spectra`) and values the corresponding
#'     names of the variables within the Python Spectrum.
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectra` single R Spectra
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importMethodsFrom Spectra Spectra
#'
#' @importFrom Spectra MsBackendMemory
#'
#' @noRd
.single_pyspec_to_rspec <-
    function(x, spectraVariables = spectraVariableMapping()) {
        plist <- x$metadata
        vars <- spectraVariables[spectraVariables %in% names(plist)]
        if (length(vars)) {
            rlist <- as.data.frame(lapply(vars, function(z) plist[z]))
            ## ## Drop NULL variables.
            ## spd <- DataFrame(rlist[lengths(rlist) > 0])
            if (!nrow(spd)) {
                spd <- DataFrame(msLevel = NA_integer_)
            }
        } else {
            spd <- data.frame(msLevel = NA_integer_)
        }
        spd$mz <- NumericList(as.numeric(x$peaks$mz), compress = FALSE)
        spd$intensity <- NumericList(as.numeric(x$peaks$intensities),
                                     compress = FALSE)
        Spectra(spd)
    }

py_to_r.matchms.Spectrum.Spectrum <- function(x) {
    .single_pyspec_to_rspec(x)
}

#' function to create a Python command to extract the peaks data.
#'
#' @param x `character(1)` with the variable name containing the `list` of
#'     `matchms.Spectrum` data.
#'
#' @return `character(1)` with the Python command.
#'
#' @noRd
.py_matchms_cmd_peaks_data <- function(x) {
    paste0("_res_ = list()\n",
           "for i in range(len(", x, ")):\n",
           "  _res_.append(", x, "[i].peaks.to_numpy)\n")
}

#' function to get the `list` of peaks matrices from a Python `list` of
#' `matchms.Spectrum` objects. This function calls a Python command.
#'
#' @note
#'
#' The returned peak matrices don't have column names.
#'
#' @param x `character(1)` with the name of the variables containing the data.
#'
#' @param local `logical(1)` passed to [reticulate::py_run_string()]
#'
#' @return a `list` of 2-column `matrix` with the m/z and intensity values.
#'     Note that the matrices don't have column names.
#'
#' @noRd
.py_matchms_peaks_data <- function(x, local = TRUE) {
    py_to_r(
        py_run_string(
            .py_matchms_cmd_peaks_data(x), local = local,
            convert = FALSE)[["_res_"]])
}

#' function to create a Python command to loop through the variable which name
#' was provided with parameter `x` and extract and combine the metadata of all
#' `matchms.Spectrum` objects in `x` as a `pandas.DataFrame`.
#'
#' @param x `character(1)` with the name of the (Python) variable with the
#'     `list` of `matchms.Spectrum` objects.
#'
#' @return `character(1)` defining the Python command.
#'
#' @noRd
.py_matchms_cmd_spectra_data <- function(x) {
    paste0(
        "import pandas as pd\n",
        "_res_ = pd.DataFrame()\n",
        "for i in range(len(", x, ")):\n",
        "  _res_ = pd.concat([_res_, pd.DataFrame(", x, "[i].metadata, index = [0])], ignore_index = True)\n")
}

#' Retrieve the metadata of all `matchms.Spectrum` objects in a Python `list`
#' as a R `data.frame`.
#'
#' @param x `character(1)` with the name of the variable containing the MS data
#'     in Python.
#'
#' @param local `logical(1)` passed to the [reticulate::py_run_string()]
#'     function.
#'
#' @return `data.frame()` with the metadata of the MS data.
#'
#' @noRd
.py_matchms_spectra_data <-
    function(x, local = TRUE, mapping = spectraVariableMapping(), ...) {
        res <- py_to_r(
            py_run_string(
                .py_matchms_cmd_spectra_data(x),
                local = local, convert = FALSE)[["_res_"]])
        mapping <- mapping[mapping %in% colnames(res)]
        res <- res[, mapping]
        colnames(res) <- names(mapping)
        res
    }
