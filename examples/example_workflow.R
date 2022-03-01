library(reticulate)
library(Spectra)
library(MsBackendMsp)
source("R/Conversion_functions.R")
source("R/compareSpectripy.R")

#load testdata to Spectra object
querry_file <- "data/test_msp.MSP"
lib_file <- "data/test_msp.MSP"
sps_q <- Spectra(querry_file, source=MsBackendMsp())
sps_l <- Spectra(lib_file, source=MsBackendMsp())

#compare Spectra with compareSpectra
compare <- compareSpectra(sps_q,sps_l, ppm=5)

#compare Spectra with matchms
comparePy <- compareSpectripy(sps_q,sps_l)

