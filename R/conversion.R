#' @title Low-level functions to convert between Spectra and matchms Spectrum
#'
#' @name conversion
#'
#' @description
#'
#' The `r_to_py.Spectra()` and `pyspec_to_rspec()` functions allow to convert
#' R [Spectra::Spectra()] objects into
#' [matchms](https://github.com/matchms/matchms) Python `matchms.Spectrum`
#' objects. These functions are designed for **advanced users or developers**
#' who want/need to integrate Python/matchms functionality into R using
#' *reticulate*. All other users should use the dedicated R functions within
#' this package that take care of running the Python code in the correct Python
#' environment.
#'
#' `.SPECTRA_2_MATCHMS` provides a mapping of core `Spectra` variables, based on
#' this [definition in matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv),
#' that is used in the `r_to_py.Spectra()` and `pyspec_to_rspec()` functions.
#'
#' @param x For `r_to_py.Spectra()`: `Spectra` object. For `pyspec_to_rspec()`:
#'     a Python list of matchms Spectrum objects.
#'
#' @return For `r_to_py.Spectra()`: Python array of Spectrum objects, same
#'     length as `x`. For `pyspec_to_rspec()`: [Spectra::Spectra()] with the
#'     converted spectra.
#'
#' @author Michael Witting, Johannes Rainer, Wout Bittremieux
#'
#' @importFrom reticulate r_to_py import py_to_r
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

#' @rdname r_to_py.Spectra
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
#' @importFrom Spectra spectrapply
#'
#' @export
r_to_py.Spectra <- function(x, convert = FALSE) {
    plist <- spectrapply(x, .single_rspec_to_pyspec)
    r_to_py(unname(plist))
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
.single_rspec_to_pyspec <- function(x, spectraVariables = .SPECTRA_2_MATCHMS) {
    pks <- unname(peaksData(x, c("mz", "intensity")))[[1L]]
    if (length(spectraVariables)) {
        slist <- as.list(spectraData(x, columns = names(spectraVariables)))
        names(slist) <- spectraVariables
        matchms$Spectrum(
            mz = np_array(pks[, 1L]),
            intensities = np_array(pks[, 2L]),
            metadata = r_to_py(slist)
        )
    } else {
        matchms$Spectrum(
            mz = np_array(pks[, 1L]), intensities = np_array(pks[, 2L])
        )
    }
}

#' @rdname pyspec_to_rspec
#'
#' @importFrom Spectra concatenateSpectra
#'
#' @export
pyspec_to_rspec <- function(x, mapping = .SPECTRA_2_MATCHMS,
                            BPPARAM = SerialParam(), .check = TRUE) {
    if (!(is(x, "list") | is(x, "python.builtin.list"))) {
        stop("'x' is expected to be a Python list.")
    }
    x <- py_to_r(x)
    if (.check && !all(vapply(x, function(z) {
        is(z, "matchms.Spectrum.Spectrum")
    }, logical(1)))) {
        stop("'x' is expected to be a Python list of matchms Spectrum objects.")
    }
    spectra_list <- bplapply(x, .single_pyspec_to_rspec,
        spectraVariables = mapping, BPPARAM = BPPARAM
    )
    do.call(concatenateSpectra, spectra_list)
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
    function(x, spectraVariables = .SPECTRA_2_MATCHMS) {
        plist <- x$metadata
        vars <- spectraVariables[spectraVariables %in% names(plist)]
        if (length(vars)) {
            rlist <- lapply(vars, function(z) plist[z])
            ## Drop NULL variables.
            spd <- DataFrame(rlist[lengths(rlist) > 0])
            if (!nrow(spd)) {
                spd <- DataFrame(msLevel = NA_integer_)
            }
        } else {
            spd <- DataFrame(msLevel = NA_integer_)
        }
        spd$mz <- NumericList(as.numeric(x$peaks$mz), compress = FALSE)
        spd$intensity <- NumericList(as.numeric(x$peaks$intensities),
            compress = FALSE
        )
        Spectra(spd)
    }
