# function to convert an R Spectra to a Python Spectrum
rspect_to_pyspec <- function(x, reference = import("matchms")) {
    
    reference$Spectrum(mz = np_array(x$mz[[1]]),
                       intensities = np_array(x$intensity[[1]]))
    
}

