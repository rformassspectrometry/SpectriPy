#' @title Converting between R and Python MS data structures
#'
#' @name conversion
#'
#' @description
#'
#' The `rspec_to_pyspec()` and `pyspec_to_rspec()` functions allow to convert
#' (translate) MS data structures between R and Python. At present the
#' R [Spectra::Spectra()] objects are translated into a list of
#' [matchms](https://github.com/matchms/matchms) Python `matchms.Spectrum`
#' objects. For better integration with the *reticulate* R package also a
#' `r_to_py.Spectra()` method is available.
#'
#' The mapping of spectra variables (in R) to (Python) spectra metadata can
#' be configured and defined with the `setSpectraVariableMapping()` and
#' `spectraVariableMapping()`. These get and set the *global* (system wide)
#' setting and are thus also used by the `r_to_py()` method.
#'
#' See the indivudual function's documentation for more details.
#'
#'
#' @section Translation of MS data objects:
#'
#' MS data structures can be translated between R and Python using the
#' `rspec_to_pyspec()` and `pyspec_to_rspec()` functions, or with the
#' `r_to_py()` method.
#'
#' - `rspec_to_pyspec()` translates an R [Spectra::Spectra()] object into a
#'   list of `matchms.Spectrum` objects. Parameter `mapping` allows to specify
#'   which spectra variables from the `Spectra` object `x` should be converted
#'   in addition to the peaks data (m/z and intensity values). It defaults to
#'   `mapping = spectraVariableMapping()` (See the respective help below for
#'   more information on the variable mapping). While being fast, this function
#'   first loads all peaks and spectra data into memory before translating to
#'   Python data structures. A less memory intense operation could be to call
#'   this function in a loop to only load parts of the data at a time into
#'   memory.
#'
#' - `pyspec_to_rspec()` translates a single, or a list of `matchms.Spectrum`
#'   objects to a [Spectra::Spectra()] object. Parameter `mapping` allows to
#'   speficy the metadata variables that should be translated and mapped in
#'   addition to the peaks data.
#'
#' - `r_to_py.Spectra()` is equivalent to `rspec_to_pyspec()`. The spectra
#'   variables that should be converted can be configures with
#'   `setSpectraVariableMapping()` (see documentation below).
#'
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
#' - `defaultSpectraVariableMapping()`: returns the *default* mapping between
#'   spectra variables and *matchms* metadata names.
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
#' @param mapping named `character()` vector defining which spectra
#'     variables/metadata should be translated between R and Python and how
#'     they should be renamed. Defaults to `spectraVariableMapping()`.
#'
#' @param x For `rspec_to_pyspec` or `r_to_py.Spectra()`: the
#'     [Spectra::Spectra()]` object that should be translated. For
#'     `pyspec_to_rspec()`: a single `matchms.Spectrum` object or a Python
#'     list of `matchms.Spectrum` objects.
#'
#' @param object For `spectraVariableMapping()`: not used.
#'
#' @param ... For `spectraVariableMapping()`: not used.
#'
#' @return For `r_to_py.Spectra()` and `rspec_to_pyspec()`: Python list of
#'     `matchms.Spectrum` objects. For `pyspec_to_rspec()`:
#'     [Spectra::Spectra()] with the MS data of all `matchms.Spectrum` objects
#'     in the submitted `list`.
#'
#' @author Michael Witting, Johannes Rainer, Wout Bittremieux
#'
#' @importFrom reticulate r_to_py py_to_r
#'
#' @importMethodsFrom Spectra spectrapply
#'
#' @examples
#'
#' ## Import a MGF file as a `Spectra` object
#' library(MsBackendMgf)
#' library(SpectriPy)
#' s <- Spectra(
#'     system.file("extdata", "mgf", "spectra2.mgf", package = "SpectriPy"),
#'     source = MsBackendMgf())
#' s
#'
#' #########################
#' ## Conversion R to Python
#'
#' ## A `Spectra` can be translated to a `list` of `matchms.Spectrum` objects
#' ## using either the `r_to_py()` method or the `rspec_to_pyspec()` function:
#' s_py <- r_to_py(s)
#' s_py
#'
#' ## The `s_py` can now be used like any other Python variable within the R
#' ## *reticulate* framework. Below we extract the m/z values of the first
#' ## spectrum
#' s_py[0]$mz
#'
#' ## Extracting that information from the `Spectra` object in R
#' s[1]$mz
#'
#' ## The `spectraVariableMapping()` defines which spectra variables (metadata)
#' ## should be translated between R and Python:
#' spectraVariableMapping()
#'
#' ## The names of that character vector represent the names of the spectra
#' ## variables in R, the elements the name of the metadata variable in Python.
#' ## Below we list the available metadata information from the first
#' ## Spectrum in Python
#' s_py[0]$metadata
#'
#' ## `setSpectraVariableMapping()` allows to replace the default mapping
#' ## of variables. Below we e.g. add a new spectra variable to the `Spectra`
#' ## object.
#' s$new_col <- 1:4
#'
#' ## To translate that variable to Python we need to include it to the
#' ## `spectraVariableMapping()`. Below we define to translate only the
#' ## precursor m/z and the new spectra variable to Python. Be aware that
#' ## `setSpectraVariableMapping()` **globally** sets the default for any
#' ## spectra variable mapping between R and Python. Thus, any subsequent
#' ## calls mapping calls will use the same mapping. It is suggested to
#' ## eventually *restore* the default mapping again after the call or
#' ## use the `rspec_to_pyspec()` function instead, that allows to configure
#' ## the mapping using a parameter `mapping`.
#' setSpectraVariableMapping(
#'     c(precursorMz = "precursor_mz", new_col = "new_col"))
#' s_py <- r_to_py(s)
#'
#' s_py[0]$metadata
#'
#' ## Restoring the global spectra variable mapping configuration to
#' ## the default mapping:
#' setSpectraVariableMapping(defaultSpectraVariableMapping())
#'
#' ## As an alternative to the `r_to_py()` we can use the `rspec_to_pyspec()`
#' ## function and provide a custom mapping using the `mapping` parameter:
#' s_py <- rspec_to_pyspec(
#'     s, mapping = c(precursorMz = "precursor_mz", new_col = "new_col"))
#'
#' #########################
#' ## Conversion Python to R
#'
#' ## A `list` of `matchms.Spectrum` objects in Python can be translated into
#' ## the corresponding MS data structure in R (i.e. a `Spectra`) object using
#' ## the `pyspec_to_rspec()` function:
#' res <- pyspec_to_rspec(s_py)
#' res
#'
#' ## All spectra from Python are thus converted into a single `Spectra` object.
#'
#' ## Or providing a custom variable mapping:
#' res <- pyspec_to_rspec(
#'     s_py, mapping = c(precursorMz = "precursor_mz", new_col = "new_col"))
#' res$new_col
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

#' @rdname conversion
#'
#' @export
defaultSpectraVariableMapping <- function() {
    .SPECTRA_2_MATCHMS
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
rspec_to_pyspec <- function(x, mapping = spectraVariableMapping()) {
    if (!is(x, "Spectra"))
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
    function(x, mapping = spectraVariableMapping()) {
        pks <- unname(peaksData(x, c("mz", "intensity")))[[1L]]
        if (length(mapping)) {
            slist <- as.list(spectraData(x, columns = names(mapping)))
            names(slist) <- mapping[names(slist)]
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
.rspec_to_pyspec <- function(x, mapping = spectraVariableMapping()) {
    pks <- peaksData(x, c("mz", "intensity"))
    sv <- mapping[!mapping %in% c("mz", "intensity")]
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

#' Function to extract the metadata from a (**single**) `matchms.Spectrum`
#' object. This function can be used in a loop over a `list` of objects to
#' create a `Spectra` from it.
#'
#' @param x `matchms.Spectrum`.
#'
#' @param mapping the mapping of spectra variables to metadata names.
#'
#' @return `data.frame()` with the spectra data.
#'
#' @noRd
.py_matchms_spectrum_spectra_data <-
    function(x, mapping = spectraVariableMapping(), ...) {
        pl <- x$metadata
        map <- mapping[mapping %in% names(pl)]
        if (length(map)) {
            ## would be nice to be able to call list(pl.items())
            res <- lapply(map, function(z) py_to_r(pl[z]))
            base::as.data.frame(res[lengths(res) > 0])
        } else data.frame(msLevel = NA_integer_)
    }

#' Function to extract the peaks data as a two-column `matrix` from a
#' (**single!**) `matchms.Spectrum` object. This function can be used in a
#' loop over a `list` of objects to create a `Spectra` from it.
#'
#' @param x `matchms.Spectrum` object.
#'
#' @return `numeric` `matrix` with two columns `"mz"` and `"intensity"`.
#'
#' @noRd
.py_matchms_spectrum_peaks_data <- function(x, ...) {
    m <- py_to_r(x$peaks$to_numpy)
    colnames(m) <- c("mz", "intensity")
    m
}

.py_matchms_spectrum_peaks_data_columns <-
    function(x, columns = c("mz", "intensity"), drop = FALSE, ...) {
        m <- py_to_r(x$peaks$to_numpy)
        colnames(m) <- c("mz", "intensity")
        m[, columns, drop = drop]
    }

#' @description
#'
#' Function to convert a single Python Spectrum object into an R Spectra using
#' the `reticulate` package.
#'
#' @param x `Spectrum` Single Python Spectrum.
#'
#' @param mapping named `character` vector with the names of the
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
    function(x, mapping = spectraVariableMapping(), ...) {
        be <- MsBackendMemory()
        be@spectraData <- .py_matchms_spectrum_spectra_data(
            x, mapping = mapping)
        be@spectraData$dataStorage <- "<memory>"
        be@peaksData <- list(.py_matchms_spectrum_peaks_data(x))
        Spectra(be)
    }

## Note: disabling this now. Generally, it would be a nice thing to have,
## but it seems R/reticulate starts automatically converting whenever a
## variable from the py is accessed, either throug `py$var` or even with
## `py_get_attr(py, "var")` - the `iterate()` call below is then converting
## the data to a `Spectra` and passing that along to the R function.
## #' @importFrom reticulate py_to_r
## #' @export
## py_to_r.matchms.Spectrum.Spectrum <- function(x) {
##     .single_pyspec_to_rspec(x)
## }

#' @export
#'
#' @importFrom MsCoreUtils rbindFill
#'
#' @importFrom reticulate iterate
#'
#' @rdname conversion
pyspec_to_rspec <- function(x, mapping = spectraVariableMapping()) {
    if (is(x, "matchms.Spectrum.Spectrum"))
        return(.single_pyspec_to_rspec(x, mapping = mapping))
    be <- MsBackendMemory()
    orig_mapping <- spectraVariableMapping()
    setSpectraVariableMapping(mapping)
    if (is.list(x))
        ITER <- lapply
    else ITER <- iterate
    be@peaksData <- ITER(x, .py_matchms_spectrum_peaks_data, simplify = FALSE)
    ## Not very efficient and elegant... get the indivudal elements and stuff
    ## into list. pandas.DataFrame can not be easily created unfortunately.
    be@spectraData <- do.call(
        rbindFill, ITER(x, .py_matchms_spectrum_spectra_data, simplify = FALSE))
    be@spectraData$dataStorage <- "<memory>"
    setSpectraVariableMapping(orig_mapping)
    Spectra(be)
}

## below are functions that would use direct python calls that iterate in
## Python.

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
#' @importFrom reticulate py_run_string
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
