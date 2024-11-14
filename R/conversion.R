#' @title Low level functions to convert between Spectra and matchms Spectrum
#'
#' @name rspec_to_pyspec
#'
#' @description
#'
#' The `rspec_to_pyspec()` and `pyspec_to_rspec()` functions allow to convert
#' R [Spectra()] objects into [matchms](https://github.com/matchms/matchms)
#' Python `Spectrum` objects. These functions are designed for
#' **advanced users or developers** who want/need to integrate Python/matchms
#' functionality into R using *reticulate*. All other users should use the
#' dedicated R functions within this package that take care of running the
#' Python code in the correct Python environment.
#'
#' Parameter `mapping` allows to define which spectra variables (metadata)
#' should be copied between the R and Python spectra. Only provided spectra
#' variables will be copied to R respectively Python. `mapping` also defines
#' the mapping between the `Spectra`'s spectra variables and the Spectrum
#' metadata. The names of the character vector `mapping` are the R spectra
#' variables and the values the corresponding names in the Python's Spectrum
#' metadata. See the output of the `spectraVariableMapping()` function for the
#' default variables and the mapping of the names.
#'
#' The `spectraVariableMapping()` function provides a default mapping of some
#' core `Spectra` variables based on this [definition in matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv).
#' The function returns a named vector that can be directly used as parameter
#' `mapping` in the `rspec_to_pyspec()` and `pyspec_to_rspec()` functions.
#'
#' @param .check Optionally disable input parameter checking. Input parameter
#'     checking should only disabled for very good reasons.
#'
#' @param BPPARAM Optional parallel processing setup.
#'
#' @param mapping Named `character` providing the spectra variable names
#'     (metadata) to convert. Names are expected to be the spectra variable
#'     names and values the corresponding names of the Python Spectrum metadata
#'     fields. See description above for more details.
#'
#' @param object ignored.
#'
#' @param reference Optional reference to Python environment `matchms`.
#'
#' @param x For `rspec_to_pyspec()`: `Spectra` object. For `pyspec_to_rspec()`:
#'     a Python list of matchms Spectrum objects.
#'
#' @param ... ignored.
#'
#' @return For `rspec_to_pyspec()`: Python array of Spectrum objects, same
#'     length than `x`. For `pyspec_to_rspec()`: [Spectra()] with the converted
#'     spectra. For `spectraVariableMapping()`: named `character` vector with
#'     names being `Spectra` variable names and values the corresponding names
#'     in `matchms`.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @export
#'
#' @importFrom reticulate r_to_py import py_to_r
#'
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @importMethodsFrom Spectra spectrapply
#'
#' @examples
#'
#' ## List the default spectra variables and their mapping.
#' spectraVariableMapping()
NULL

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @exportMethod spectraVariableMapping
#'
#' @rdname rspec_to_pyspec
setMethod("spectraVariableMapping", "missing", function(object, ...) {
    .SPECTRA_2_MATCHMS
})

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

#' @rdname rspec_to_pyspec
#'
#' @importFrom methods is
#'
#' @export
rspec_to_pyspec <- function(x, mapping = spectraVariableMapping(),
                            reference = import("matchms"),
                            BPPARAM = SerialParam(), .check = TRUE) {
    if (.check && !is(x, "Spectra"))
        stop("'x' should be a Spectra object.")
    plist <- spectrapply(x, .single_rspec_to_pyspec, spectraVariables = mapping,
                         reference = reference, BPPARAM = BPPARAM)
    r_to_py(unname(plist))
}

#' @rdname rspec_to_pyspec
#'
#' @importFrom Spectra concatenateSpectra
#'
#' @export
pyspec_to_rspec <- function(x, mapping = spectraVariableMapping(),
                            BPPARAM = SerialParam(), .check = TRUE) {
    if (!is(x, "python.builtin.list"))
        stop("'x' is expected to be a Python list.")
    x <- py_to_r(x)
    if (.check && !all(vapply(x, function(z)
        is(z, "matchms.Spectrum.Spectrum"), logical(1))))
        stop("'x' is expected to be a Python list of matchms Spectrum objects.")
    spectra_list <- bplapply(x, .single_pyspec_to_rspec,
                             spectraVariables = mapping, BPPARAM = BPPARAM)
    do.call(concatenateSpectra, spectra_list)
}

#' @description
#'
#' Function to convert a **single** R Spectra object (of length 1) into a
#' Python matchms Spectrum using the `reticulate` package.
#'
#' @param x `Spectra` object **of length 1!**.
#'
#' @param spectraVariables named `character` vector defining the spectra
#'     varibles that should be stored as metadata in `matchms`' metadata. Names
#'     are expected to be the spectra variable names and values the
#'     corresponding metadata fields in `matchms`. Defaults to
#'     [spectraVariableMapping()]. If `spectraVariables = character()` no
#'     metadata will be stored.
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectrum` Single Python Spectrum
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraData
#'
#' @importMethodsFrom Spectra mz
#'
#' @importMethodsFrom Spectra intensity
#'
#' @importFrom reticulate np_array r_to_py
#'
#' @noRd
.single_rspec_to_pyspec <- function(x,
                                    spectraVariables = spectraVariableMapping(),
                                    reference = import("matchms")) {
    if (length(spectraVariables)) {
        slist <- as.list(spectraData(x, columns = names(spectraVariables)))
        ## ## Seems matchms.Spectrum does not support NA retention times?
        ## if (any(names(slist) == "rtime") && is.na(slist$rtime))
        ##     slist$rtime <- 0
        names(slist) <- spectraVariables
        reference$Spectrum(mz = np_array(mz(x)[[1L]]),
                           intensities = np_array(intensity(x)[[1L]]),
                           metadata = r_to_py(slist))
    } else reference$Spectrum(mz = np_array(mz(x)[[1L]]),
                              intensities = np_array(intensity(x)[[1L]]))
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
#' @noRd
.single_pyspec_to_rspec <-
    function(x, spectraVariables = spectraVariableMapping()) {
        plist <- x$metadata
        vars <- spectraVariables[spectraVariables %in% names(plist)]
        if (length(vars)) {
            rlist <- lapply(vars, function(z) plist[z])
            ## Drop NULL variables.
            spd <- DataFrame(rlist[lengths(rlist) > 0])
            if (!nrow(spd))
                spd <- DataFrame(msLevel = NA_integer_)
        } else
            spd <- DataFrame(msLevel = NA_integer_)
        spd$mz <- NumericList(as.numeric(x$peaks$mz), compress = FALSE)
        spd$intensity <- NumericList(as.numeric(x$peaks$intensities),
                                     compress = FALSE)
        Spectra(spd)
    }

#' Extract all spectraData and all mz and intensity values, give them to
#' Python to create an array of Spectrum. Could be faster because loop is
#' performed in Python rather than in R.
#'
#' @noRd
.multi_rspec_to_pyspec <- function() {
}

#' Extract data in Python for all elements: python function should return a
#' list with the metadata (convert to DataFrame) and the m/z and intensity
#' values. Would avoid loop in R alltogether.
#'
#' @noRd
.multi_pyspec_to_rspec <- function() {
}