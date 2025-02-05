library(reticulate)

# Specify Python 3.10
use_python("/your/path/to/r-miniconda-arm64/bin/python3.10", required = TRUE)

py_config()

system("pip install ms2deepscore matchms numpy torch", intern = TRUE)

# Step 1: Initialize Python Environment and Import Required Libraries
initialize_ms2deepscore <- function(python_env = NULL) {
  if (!is.null(python_env)) {
    use_virtualenv(python_env, required = TRUE)
  }
  
  # Import Python modules
  torch <- import("torch", convert = FALSE)
  ms2ds <- import("ms2deepscore", convert = FALSE)
  matchms <- import("matchms", convert = FALSE)
  
  # Step 2: Patch PyTorch to Allow Full Model Loading
  cat("Patching PyTorch to allow full model loading...\n")
  py_run_string("
import torch

# Patch `torch.load()` globally to set `weights_only=False`
original_load = torch.load
def patched_load(*args, **kwargs):
    if 'weights_only' not in kwargs:
        kwargs['weights_only'] = False  # Force weights_only=False
    return original_load(*args, **kwargs)

torch.load = patched_load
  ")
  py_run_string("
from matchms.Pipeline import Pipeline, create_workflow
from matchms.filtering.default_pipelines import DEFAULT_FILTERS
from ms2deepscore import MS2DeepScore
")
  
  return(list(ms2ds = ms2ds, matchms = matchms))#, create_workflow = create_workflow))
}

# Step 3: Load MS2DeepScore Model
load_ms2deepscore_model <- function(model_file, env) {
  ms2ds <- env$ms2ds
  cat("Loading MS2DeepScore model...\n")
  
  model <- ms2ds$models$load_model(model_file)
  
  if (is.null(model)) {
    stop("ERROR: Model failed to load.")
  } else {
    cat("Model loaded successfully!\n")
  }
  
  return(model)
}

# Step 3: Run MS2DeepScore Pipeline
run_ms2deepscore_pipeline <- function(spectrum_path, model, env) {
  matchms <- env$matchms
  create_workflow <- env$create_workflow
  ms2ds <- env$ms2ds
  
  cat("Setting up MS2DeepScore pipeline...\n")
  
  # Create the processing pipeline using `create_workflow`
  # Create the processing pipeline using `create_workflow`
  py_run_string("
pipeline = Pipeline(create_workflow(
    query_filters=DEFAULT_FILTERS,
    score_computations=[[MS2DeepScore, {'model': r.model}]]
))
report = pipeline.run(r.spectrum_path)
similarity_matrix = pipeline.scores.to_array()
")
  
  # Run the pipeline on the spectrum file
  cat("Running pipeline...\n")
  pipeline <- py$pipeline
  #report <- py$report
  # Extract similarity matrix
  similarity_matrix <- py$similarity_matrix #pipeline$score_computations
  
  cat("Pipeline run successfully!\n")
  return(similarity_matrix)
}

# Step 5: Execute the Pipeline in R
model_path <- "/your/path/to/ms2deepscore_model.pt"  # Replace with the actual path
spectrum_path <- "/your/path/to/R_python_hackathon/pesticides.mgf"  # Replace with the actual path

# Initialize the Python environment and import libraries
env <- initialize_ms2deepscore()

# Load MS2DeepScore model
model <- load_ms2deepscore_model(model_path, env)

# Run the pipeline and get similarity scores
similarity_scores <- run_ms2deepscore_pipeline(spectrum_path, model, env)

