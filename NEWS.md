# SpectriPy 0.99

## Changes in 0.99.8

- Verbose startup messages and require Python 3.12.

## Changes in 0.99.7

- Replace the *msdata* package with *MsDataHub*.

## Changes in 0.99.6

- Set `delay_load = FALSE` in `import()` calls during package loading

## Changes in 0.99.5

- Replace Python library installation *via* virtualenv or conda with
  `py_require()` from *reticulate* version > 1.41.0.

## Changes in 0.99.4

- Refactor functions to initialize Python libraries.

## Changes in 0.99.3

- Address Bioconductor review comments.

## Changes in 0.99.2

- Use Python virtualenv instead of miniconda.

## Changes in 0.99.1

- Check for presence of conda and eventually install miniconda through
  *reticulate*.

## Changes in 0.99.0

- Prepare for Bioconductor submission.


# SpectriPy 0.5

## Changes in 0.5.3

- Vignette updates.

## Changes in 0.5.2

- Fix *spectrum_utils* version.

## Changes in 0.5.1

- Add additional vignettes and update/fix documentation.

## Changes in 0.5.0

- Add support for *spectrum_utils* library.

# SpectriPy 0.4

## Changes in 0.4.0

- Add a `MsBackendPython` backend referencing to MS data residing in Python.

# SpectriPy 0.3

## Changes in 0.3.1

- Add `filterSpectriPy()` function for spectra filtering using *matchms*.
- Add a new example MGF file to the package.
- Add a new quarto vignette.

## Changes in 0.3.0

- Add changes and results introduced a the EuBIC 2025 r-python hackathon in
  Neustift, Italy.
- Remove dependency from *basilisk* and base the package entirely on
  *reticulate*.

# SpectriPy 0.2

## Changes in 0.2.1

- Add `filterSpectriPy` method and related parameter objects to perform spectra
  filtering/processing via matchms in python.

## Changes in 0.2.0

- Use *matchms* version 0.28.2.


# SpectriPy 0.1

## Changes in 0.1.1

- Small updates and fixes in the package's vignette.

## Changes in 0.1.0

- Add `compareSpectriPy` method and related parameter objects to perform spectra
  similarity calculations in python.


# SpectriPy 0.0

## Changes in 0.0.2

- Add `basilisk` environment.
- Add `spectraVariableMapping`.
- Refactor functions to convert between R and python spectrum objects and add
  unit tests.
