# Install python, and all python associated packages

Install_python_spectriPy = function(PATH_miniconda = "~/miniconda3_reticulate",ENV_name = "r-reticulate"){
  print("Install reticulate prior to run this")
  library(reticulate)
  options(timeout = 600)
  install_miniconda(path = PATH_miniconda,force = FALSE)
  
  use_python(paste(PATH_miniconda,"/bin/python3",sep = ""))
  conda_create(ENV_name)
  
  py_install("matchms==0.28.2",pip = TRUE)
  py_install("spectrum_utils",pip = TRUE)
  py_install("numpy==2.0.2",pip = FALSE)
  
  conda_install(envname = ENV_name,"matchms==0.28.2",pip = TRUE,python_version = "3.12.2")
  conda_install(envname = ENV_name,"spectrum_utils",pip = TRUE)
  conda_install(envname = ENV_name,"numpy==2.0.2",pip = FALSE)
  
  print("R packages that needs to be installed for any analysis: BiocManager , Spectra , msdata , AnnotationHub , CompoundDb , GenomeInfoDbData , mzR , basilisk , SpectriPy")
  print("Use the options(timeout = 600)")
  print("WARNING: use_condaenv(paste(PATH_miniconda,ENV_name,'/bin/python3',sep ='') needs to be run after importing reticulate to tell it where to find the python environment")
  return("Python installation finished")
}
