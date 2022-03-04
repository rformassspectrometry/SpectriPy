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
#' @param x `Spectra` object.
#'
#' @param metadata `character` vector with names of spectra data fields that
#'     shall be converted as well
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectrum` Single Python Spectrum
#'
#' @author Michael Witting
#'
#' @noRd
.rspec_to_pyspec <- function(x,
                             metadata = NA,
                             reference = import("matchms")) {

  # isolate all metadata and convert to list, rename list according to Python matchms
  spectraData_list <- as.list(spectraData(x))

  if(!is.na(metadata)) {
    spectraData_list <- spectraData_list[which(spectraData_list %in% metadata)]
  }

  # Spectra variable mapping according to matchms definitions
  names(spectraData_list)[which(names(spectraData_list) == "precursorMz")] <- "precursor_mz"

  # create python spectrum
  reference$Spectrum(mz = np_array(x$mz[[1]]),
                     intensities = np_array(x$intensity[[1]]),
                     metadata = r_to_py(spectraData_list))

}

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
#' @importMethodsFrom reticulate import
#'
#' @noRd
.pyspec_to_rspec <- function(x, reference = import("matchms")) {

  # isolate all metadata and rename according to R Spectra
  spectraData_list <- x$metadata

  # Spectra variable mapping according to matchms definitions
  names(spectraData_list)[which(names(spectraData_list) == "precursor_mz")] <- "precursorMz"
  spd <- DataFrame(spectraData_list)

  # add spectrum data
  spd$mz <- IRanges::NumericList(x$peaks$mz)
  spd$intensity <- IRanges::NumericList(x$peaks$intensities)

  Spectra(spd)

}
