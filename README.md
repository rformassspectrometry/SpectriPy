# Interface for R's `Spectra` package with Python

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

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

`SpectriPy` uses [`basilisk`](https://bioconductor.org/packages/basilisk) to
ensure all required python packages are installed and available (in the correct
version) on each system. `basilisk` installs a self-contained conda environment,
thus, the `SpectriPy` package is independent of the system's Python envrionment.


## Concepts and examples
