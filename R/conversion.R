#' @description
#' Function to convert R Spectra to Python Spectrum
#'
#' @param x `Spectra` object
#'
#' @param metadata `character` vector with names of spectra data fields that
#'     shall be converted as well
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `array` Python array with same length as x.
#'
#' @author Michael Witting
#'
#' @export
#'
#' @importFrom reticulate r_to_py import
#'
#' @importFrom BiocParallel SerialParam
#'
#' @importMethodsFrom Spectra spectrapply
convertRSpectraToPySpectrum <- function(x,
                                        metadata = NA,
                                        reference = import("matchms"),
                                        BPPARAM = SerialParam()) {

    # convert multiple spectra to python spectrum
    pySpectrum_list <- spectrapply(x, .rspec_to_pyspec, BPPARAM = BPPARAM)
    return(r_to_py(unname(pySpectrum_list)))
}

#' @description
#' Function to convert Python Spectrum to R Spectra
#'
#' @param x `array` Python array object with Spectrum
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectra` R spectra object
#'
#' @author Michael Witting
#'
#' @export
#'
#' @importFrom reticulate import py_to_r
#'
#' @importFrom BiocParallel bplapply
convertPySpectrumToRSpectra <- function(x,
                                        reference = import("matchms"),
                                        BPPARAM = SerialParam()) {

  # convert Python array to R list
  x <- py_to_r(x)

  # differentiate if x is a single spectrum or multiple spectra
  if(length(x) == 1) {
    #convert single spectrum to spectra
    return(.pyspec_to_rspec(x))
  } else {
    spectra_list <- bplapply(x, .pyspec_to_rspec, BPPARAM = BPPARAM)
    return(do.call(c, spectra_list))
  }
}

#' @description
#'
#' Function to convert a single R Spectra object into a Python Spectrum using
#' the `reticulate` package.
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
.rspec_to_pyspec <- function(x, spectraVariables = spectraVariableMapping(),
                             reference = import("matchms")) {
    if (length(spectraVariables)) {
        slist <- as.list(spectraData(x, columns = names(spectraVariables)))
        names(slist) <- spectraVariables
        reference$Spectrum(mz = np_array(mz(x)[[1L]]),
                           intensities = np_array(intensity(x)[[1L]]),
                           metadata = r_to_py(slist))
    } else reference$Spectrum(mz = np_array(mz(x)[[1L]]),
                              intensities = np_array(intensity(x)[[1L]]))
}

#' @title Convert between Spectra and matchms variable names
#'
#' @description
#'
#' Convert `Spectra` spectra variable names to names for matchms spectrum
#' metadata. Conversion is based on this [definition in matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv)
#'
#' @param object ignored.
#'
#' @param ... ignored.
#'
#' @return named `character` vector with names being `Spectra` variable names
#'     and values the corresponding names in `matchms`.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @exportMethod spectraVariableMapping
#'
#' @examples
#'
#' spectraVariableMapping()
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



#' @description
#'
#' Function to convert a single Python Spectrum object into an R Spectra using
#' the `reticulate` package.
#'
#' @param x `Spectrum` Single Python Spectrum.
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectra` single R Spectra
#'
#' @author Michael Witting
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importMethodsFrom Spectra Spectra
#'
#' @noRd
.pyspec_to_rspec <- function(x, reference = import("matchms")) {

  # isolate all metadata and rename according to R Spectra
  spectraData_list <- x$metadata

  # Spectra variable mapping according to matchms definitions
  names(spectraData_list)[which(names(spectraData_list) == "precursor_mz")] <- "precursorMz"
  spd <- DataFrame(spectraData_list)

  # add spectrum data
  spd$mz <- NumericList(x$peaks$mz)
  spd$intensity <- NumericList(x$peaks$intensities)

  Spectra(spd)
}
