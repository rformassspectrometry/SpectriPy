# Interface for R's `Spectra` package with Python

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

This repository provides example code to convert/transfer MS-related data
between R and Python. Final goal is to develop an R package that allows to
efficiently integrate python packages for MS analysis with the `Spectra`
package. The `SpectriPy` package allows to use functionality of the `matchms`
Python packages directly in R and in addition exports also the low level
functions to convert between R `Spectra` objects and Python `matchms` Spectrum
objects that enable advanced users or developers to integrate additional Python
functionality into R *via* the `reticulate` R package.

## Contributors

- Carolin Huber
- Michael Witting
- Helge Hecht
- Johannes Rainer


## Pre-requisites and installation instructions

`SpectriPy` uses [`basilisk`](https://bioconductor.org/packages/basilisk) to
ensure all required python packages are installed and available (in the correct
version) on each system. `basilisk` installs a self-contained conda environment,
thus, the `SpectriPy` package is independent of the system's Python environment.

To install the package use

```
BiocManager::install("RforMassSpectrometry/SpectriPy")
```

## Concepts and examples
