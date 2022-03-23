# Integrating [`Spectra`](https://github.com/RforMassSpectrometry/Spectra) with Python's `matchms` package

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/SpectriPy/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/SpectriPy/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/SpectriPy/coverage.svg?branch=main)](http://codecov.io/github/rformassspectrometry/SpectriPy?branch=main)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

The `SpectriPy` package allows integration of Python MS packages into a
`Spectra`-based MS analysis in `R`. Python functionality is wrapped into R
functions allowing a seamless integration of the functionality of Python's
[`matchms`](https://github.com/matchms/) package into `R`. In addition,
functions to convert between R's `Spectra` objects and Python's `matchms`
spectrum objects are available to the advanced user or developer.

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

See the package's
[vignette](https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html).
