#'
#' Prototype for conversion of R Spectra
#' @export
convertRSpectraToPySpectrum <- function(x, reference = import("matchms")) {

  # differentiate if x is a singe spectrum or multiple spectra
  if(length(x) == 1) {
    # convert single spectrum to python spectum
    return(.rspec_to_pyspec(x, reference = reference))
  } else {
    # convert multiple spectra to python spectrum
    pySpectrum_list <- spectrapply(x, .rspec_to_pyspec)
    return(unname(pySpectrum_list))
  }
}

#'
#' Prototpye for conversion of Python Spectrum
#' @export
convertPySpectrumToRSpectra <- function(x, reference = import("matchms")) {

  # differentiate if x is a single spectrum or multiple spectra
  if(length(x) == 1) {
    #convert single spectrum to spectra
    return(.pyspec_to_rspec(x))
  } else {
    spectra_list <- lapply(x, .pyspec_to_rspec)
    return(do.call(c, spectra_list))
  }
}

# function to convert an R Spectra to a Python Spectrum
.rspec_to_pyspec <- function(x, reference = import("matchms")) {

  # isolate all metadata and convert to list, rename list according to Python matchms
  spectraData_list <- as.list(spectraData(x))
  names(spectraData_list)[which(names(spectraData_list) == "precursorMz")] <- "precursor_mz"

  # create python spectrum
  reference$Spectrum(mz = np_array(x$mz[[1]]),
                     intensities = np_array(x$intensity[[1]]),
                     metadata = r_to_py(spectraData_list))

}

# function to convert a Python Spectrum to an R Spectra
.pyspec_to_rspec <- function(x, reference = import("matchms")) {

  # isolate all metadata and rename according to R Spectra
  spectraData_list <- x$metadata
  names(spectraData_list)[which(names(spectraData_list) == "precursor_mz")] <- "precursorMz"
  spd <- DataFrame(spectraData_list)

  # add spectrum data
  spd$mz <- IRanges::NumericList(x$peaks$mz)
  spd$intensity <- IRanges::NumericList(x$peaks$intensities)

  Spectra(spd)

}