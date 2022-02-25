# Interface for R's `Spectra` package with Python

This repository provides example code to convert/transfer MS-related data
between R and Python. Final goal is to develop an R package that allows to
efficiently integrate python packages for MS analysis with the `Spectra`
package.

## Contributors

- Carolin Huber
- Michael Witting
- Helge Hecht
- Johannes Rainer


## Pre-requisites and installation instructions

SpectriPy requires a running Python conda environment containing minimally the matchms package.
Following installation instructions are taken from https://github.com/matchms/matchms.

Prerequisites:

Python 3.7, 3.8, or 3.9

Anaconda

Install matchms from Anaconda Cloud with

```
# install matchms in a new virtual environment to avoid dependency clashes  
conda create --name matchms python=3.8  
conda activate matchms  
conda install --channel nlesc --channel bioconda --channel conda-forge matchms
```


## Concepts and examples

