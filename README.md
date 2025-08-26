# Enhancing Cross-Language Mass Spectrometry Data Analysis with R and Python

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/SpectriPy/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/SpectriPy/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/SpectriPy/branch/main/graph/badge.svg?token=638UZM0DXP)](https://codecov.io/gh/rformassspectrometry/SpectriPy)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.08070/status.svg)](https://doi.org/10.21105/joss.08070)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/SpectriPy.svg)](https://bioconductor.org/packages/stats/bioc/SpectriPy/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/SpectriPy.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/SpectriPy/)

![SpectriPy_logo](man/figures/logo_100.png)

The *SpectriPy* package allows integration of Python MS packages into a
[*Spectra*](https://github.com/RforMassSpectrometry/Spectra)-based MS analysis
in R. By wrapping Python functionality into R functions, *SpectriPy* allows a
seamless integration of Python libraries into R. For example, *SpectriPy* can
leverage the spectral similarity, filtering, normalization etc. calculations
from the Python [*matchms*](https://github.com/matchms) library and contains
functions to convert between R's `Spectra::Spectra` objects and
`matchms.Spectrum` and `spectrum_utils.spectrum.MsmsSpectrum` objects from the
Python [*matchms*](https://github.com/matchms) and
[*spectrum_utils*](https://github.com/bittremieux-lab/spectrum_utils) libraries,
respectively. R and Python spectral objects are easily translated and available
in one workflow (i.e., a quarto document), enabling the advanced user or
developer to create custom functions or workflows on `Spectra` objects in Python
and executing them in R using the *reticulate* R package, and vice versa.

If you use *SpectriPy* in your research, please cite:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.08070/status.svg)](https://doi.org/10.21105/joss.08070)


# Installation

*SpectriPy* needs Python (version >= 3.12) to be installed on the system. All
necessary Python libraries (listed below) are automatically installed by the
[*reticulate*](https://rstudio.github.io/reticulate) R package. *SpectriPy*'s
Python library management uses the
[`py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html)
function introduced in *reticulate* version 1.41 and should hence work on most
system without problems. To install *SpectriPy*:

```r
install.packages("BiocManager")
BiocManager("SpectriPy")
```

In addition it is possible to install the Python libraries manually (e.g., for
the system Python version) and specify the version of Python (or of the local
*virtualenv* or *conda* environment) using the `RETICULATE_PYTHON` or
`RETICULATE_PYTHON_ENV` environment variables. If any of these environment
variables are defined, all Python libraries listed below **must** be installed,
since *SpectriPy* (respectively *reticulate*) will not try to install them
automatically. The required Python libraries with the suggested and tested
versions are:

- [*matchms*](https://github.com/matchms) 0.30.0
- [*spectrum_utils*](https://github.com/bittremieux-lab/spectrum_utils) 0.3.2
- *numpy* 2.2.0

See also sections [*Startup and Python
configuration*](https://rformassspectrometry.github.io/SpectriPy/articles/detailed-installation-configuration.html#sec-python)
for more details or [*Fixing package installation or loading
problems*](https://rformassspectrometry.github.io/SpectriPy/articles/detailed-installation-configuration.html#sec-fix)
if installation or loading fails.


# Documentation for users

See the extensive documentation for the use of *SpectriPy*:

- The rendered [SpectriPy package’s vignette](https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html),
- The [function reference](https://rformassspectrometry.github.io/SpectriPy/reference/index.html),
- A [tutorial](https://rformassspectrometry.github.io/Metabonaut/articles/SpectriPy_tutorial_metabonaut.html) for the annotation of LC-MS/MS
spectra using an MGF library and the ModifiedCosine algorithm from *matchms*.

TLDR:

```r
#' R session:

library(Spectra)
library(SpectriPy)
```

## Example: Spectra similarity calculations using `matchms`

The *SpectriPy* package provides the `compareSpectriPy()` function that allows
to perform spectra similarity calculations using the scoring functions from MS
Python packages. For example, the [*CosineGreedy*](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.CosineGreedy.html) parameter `CosineGreedyParam` from the
*matchms* Python package.

1) We create some simple example spectra.

```r
#' R session:

library(Spectra)
library(SpectriPy)

#' Create a Spectra object with two MS2 spectra for Caffeine.
caf <- DataFrame(
    msLevel = c(2L, 2L),
    name = "Caffeine",
    precursorMz = c(195.0877, 195.0877)
)
caf$intensity <- list(
    c(340.0, 416, 2580, 412),
    c(388.0, 3270, 85, 54, 10111))
caf$mz <- list(
    c(135.0432, 138.0632, 163.0375, 195.0880),
    c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864))
caf <- Spectra(caf)

#' Create a Spectra object with two MS2 spectra for 1-Methylhistidine
mhd <- DataFrame(
    msLevel = c(2L, 2L),
    precursorMz = c(170.0924, 170.0924),
    id = c("HMDB0000001", "HMDB0000001"),
    name = c("1-Methylhistidine", "1-Methylhistidine"))
mhd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
mhd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
mhd <- Spectra(mhd)
```

2) We calculate pairwise similarities between all spectra defined above and
those of caffeine using *Spectra*'s built-in `compareSpectra()` function.

```r
#' R session:

all <- c(caf, mhd)
res_r <- compareSpectra(all, caf)
res_r
```

3) We calculate the similarity using the *CosineGreedy* function from *matchms*,
changing the `tolerance` to a value of `0.05` (instead of the default `0.1`).

```r
#' R session:

res <- compareSpectriPy(all, caf, param = CosineGreedy(tolerance = 0.05))
res
```

As a result `compareSpectriPy()` returns also a numeric matrix of similarities.
Note also that the first `compareSpectriPy()` call takes usually a little longer
because the Python setup has to be initialized.


# Documentation for developers

See the [developer notes](devnotes.md).


## Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
General information on the package structure and some helpful pointers are given
in the [developer notes](devnotes.md) document. Also, please check the
[coding style
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#coding-style)
and importantly, follow our [code of
conduct](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct).


# License

See the [DESCRIPTION](DESCRIPTION) and [LICENSE](LICENSE) file.
