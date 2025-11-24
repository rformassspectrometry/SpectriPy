# Filter Spectra using Python's matchms library

The `filterSpectriPy()` function allows to filter/process a `Spectra`
object using the `select_by_intensity()`, `select_by_mz()`,
`remove_peaks_around_precursor_mz()`, and `normalize_intensities()`
functions of the Python
[matchms.filtering](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html)
module.

Selection and configuration of the algorithm can be performed with one
of the parameter objects (equivalent to *matchms*' function names):

- `select_by_intensity()`: Keeps only the peaks within defined intensity
  range (keep if `intensity_from` \>= intensity \>= `intensity_to`). See
  also the respective [documentation in
  *matchms*](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.peak_processing.select_by_intensity.html).

- `select_by_mz()`: Keeps only the peaks between `mz_from` and `mz_to`
  (keep if `mz_from` \>= m/z \>= `mz_to`). See also the respective
  [documentation in
  *matchms*](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.peak_processing.select_by_mz.html).

- `remove_peaks_around_precursor_mz()`: Removes the peaks that are
  within `mz_tolerance` (in Da) of the precursor mz, excluding the
  precursor peak.

- `normalize_intensities()`: Normalizes the intensities of peaks (and
  losses) to unit height.

## Usage

``` r
select_by_intensity(intensity_from = 10, intensity_to = 200)

select_by_mz(mz_from = 0, mz_to = 1000)

remove_peaks_around_precursor_mz(mz_tolerance = 17)

normalize_intensities()

# S4 method for class 'Spectra,filter_param'
filterSpectriPy(object, param, mapping = spectraVariableMapping(), ...)
```

## Arguments

- intensity_from:

  `numeric(1)`: Set lower threshold for peak intensity. Default is 10.

- intensity_to:

  `numeric(1)`: Set upper threshold for peak intensity. Default is 200.

- mz_from:

  `numeric(1)`: Set lower threshold for m/z peak positions. Default is
  0.

- mz_to:

  `numeric(1)`: Set upper threshold for m/z peak positions. Default is
  1000.

- mz_tolerance:

  `numeric(1)`: Tolerance of m/z values that are not allowed to lie
  within the precursor mz. Default is 17 Da.

- object:

  A [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object.

- param:

  one of parameter classes listed above (such as
  `select_by_intensity()`) defining the filter/processing function in
  Python and its parameters.

- mapping:

  named [`character()`](https://rdrr.io/r/base/character.html) defining
  which spectra variables/metadata should be converted between R and
  Python and how they should be renamed. Defaults to
  [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html).
  See
  [`setSpectraVariableMapping()`](https://rformassspectrometry.github.io/SpectriPy/reference/conversion.md)
  for more information.

- ...:

  ignored.

## Value

`filterSpectriPy()` returns a `Spectra` object on which the
filtering/processing function has been applied

## Note

The first call to the `filterSpectriPy()` function can take longer
because the Python environment needs to be first set up.

`filterSpectriPy()` first translates the `Spectra` to Python, applies
the filter functions from the *matchms* Python libraries and then
translates the filtered data back to a `Spectra` object. Thus, any
spectra variables other than those that are translated between R and
Python will be lost during the processing. Use
[`setSpectraVariableMapping()`](https://rformassspectrometry.github.io/SpectriPy/reference/conversion.md)
to define which spectra variables should be transferred/converted
between R and Python. See also examples below for more information.

The [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object returned by `filterSpectriPy()` will **always** use an in-memory
backend (i.e. the
[`Spectra::MsBackendMemory()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)),
independently of the backend used by the backend used by the input
`Spectra`.

## See also

- [`Spectra::filterIntensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  [`Spectra::filterMzRange()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  [`Spectra::scalePeaks()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
  in the `Spectra` package for pure R implementations of
  filtering/processing calculations.

- [`rspec_to_pyspec()`](https://rformassspectrometry.github.io/SpectriPy/reference/conversion.md)
  or
  [`pyspec_to_rspec()`](https://rformassspectrometry.github.io/SpectriPy/reference/conversion.md)
  for the functions used to translated the MS data between R and Python.

## Author

Thomas Naake

## Examples

``` r

library(Spectra)

## create some example Spectra
DF <- DataFrame(
    msLevel = c(2L, 2L, 2L),
    name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
    precursorMz = c(195.0877, 195.0877, 170.0924)
)
DF$intensity <- list(
    c(340.0, 416, 2580, 412),
    c(388.0, 3270, 85, 54, 10111),
    c(3.407, 47.494, 3.094, 100.0, 13.240))
DF$mz <- list(
    c(135.0432, 138.0632, 163.0375, 195.0880),
    c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
    c(109.2, 124.2, 124.5, 170.16, 170.52))
sps <- Spectra(DF)

## Filter: select_by_intensity
res <- filterSpectriPy(
    sps, select_by_intensity(intensity_from = 15, intensity_to = 300))
## Only mass peaks with intensities between the specified limits are
## retained
intensity(res)
#> NumericList of length 3
#> [[1]] numeric(0)
#> [[2]] 85 54
#> [[3]] 47.494 100
## Compared to the original intensities
intensity(sps)
#> NumericList of length 3
#> [[1]] 340 416 2580 412
#> [[2]] 388 3270 85 54 10111
#> [[3]] 3.407 47.494 3.094 100 13.24

## Note that the spectra variable `"name"` was lost during conversion of
## the MS data between R and Python:
sps$name
#> [1] "Caffeine"          "Caffeine"          "1-Methylhistidine"
any(spectraVariables(res) == "name")
#> [1] FALSE

## Only spectra variables defined by `spectraVariableMapping()` are
## converted and thus retained:
spectraVariableMapping()
#>                  precursorMz           precursorIntensity 
#>               "precursor_mz"        "precursor_intensity" 
#>              precursorCharge                        rtime 
#>                     "charge"             "retention_time" 
#>              collisionEnergy      isolationWindowTargetMz 
#>           "collision_energy" "isolation_window_target_mz" 
#>                      msLevel 
#>                   "ms_level" 

## We can also pass a custom *spectra variable mapping* with the `mapping`
## parameter to the `filterSpectriPy()` function. Below we create such
## a mapping by adding the translation of a spectra variable `"name"` to
## a metadata name `"compound_name"` to the default spectra variable
## mapping `defaultSpectraVariableMapping()`.
map <- c(defaultSpectraVariableMapping(), name = "compound_name")
map
#>                  precursorMz           precursorIntensity 
#>               "precursor_mz"        "precursor_intensity" 
#>              precursorCharge                        rtime 
#>                     "charge"             "retention_time" 
#>              collisionEnergy      isolationWindowTargetMz 
#>           "collision_energy" "isolation_window_target_mz" 
#>                      msLevel                         name 
#>                   "ms_level"              "compound_name" 

## Repeat the filtering operation passing this mapping information:
res <- filterSpectriPy(
    sps, select_by_intensity(intensity_from = 15, intensity_to = 300),
    mapping = map)
res$name
#> [1] "Caffeine"          "Caffeine"          "1-Methylhistidine"
```
