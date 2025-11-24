# Converting between R and Python MS data structures

The `rspec_to_pyspec()` and `pyspec_to_rspec()` functions allow to
convert (translate) MS data structures between R and Python. At present
the R
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
objects can be either translated into a list of
[matchms](https://github.com/matchms/matchms) Python `matchms.Spectrum`
objects or
[spectrum_utils](https://github.com/bittremieux-lab/spectrum_utils)
Python `spectrum_utils.spectrum.MsmsSpectrum` objects. For better
integration with the *reticulate* R package also a `r_to_py.Spectra()`
method is available.

The mapping of spectra variables (in R) to (Python) spectra metadata can
be configured and defined with the `setSpectraVariableMapping()` and
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html).
These get and set the *global* (system wide) setting and are thus also
used by the
[`r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
method.

Properties for translation to the MS data objects of the different
Python libraries are:

- *matchms*: the `matchms.Spectrum` objects support arbitrary metadata,
  so any spectra variable can be translated and stored in these objects.

- *spectrum_utils*: the `spectrum_utils.spectrum.MsmsSpectrum` object
  supports metadata variables *identifier* (`character`), *precursor_mz*
  (`numeric`), *precursor_charge* (`integer`) and optionally also
  *retention_time* (`numeric`).

See the indivudual function's documentation for more details.

Function to convert R Spectra objects into a Python list of matchms
Spectrum objects using the `reticulate` package.

## Usage

``` r
# S4 method for class 'character'
spectraVariableMapping(object, x = character(), ...)

# S4 method for class 'missing'
spectraVariableMapping(object, ...)

setSpectraVariableMapping(x)

defaultSpectraVariableMapping()

# S3 method for class 'Spectra'
r_to_py(x, convert = FALSE)

rspec_to_pyspec(
  x,
  mapping = spectraVariableMapping(),
  pythonLibrary = c("matchms", "spectrum_utils")
)

pyspec_to_rspec(
  x,
  mapping = spectraVariableMapping(),
  pythonLibrary = c("matchms", "spectrum_utils")
)
```

## Arguments

- object:

  For
  [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html):
  not used.

- x:

  `Spectra` object.

- ...:

  For
  [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html):
  not used.

- convert:

  Boolean; should Python objects be automatically converted to their R
  equivalent? Defaults to `FALSE`.

- mapping:

  named [`character()`](https://rdrr.io/r/base/character.html) vector
  defining which spectra variables/metadata should be translated between
  R and Python and how they should be renamed. Defaults to
  [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html).

- pythonLibrary:

  For `rspec_to_pyspec()` and `pyspec_to_rspec()`: `character(1)`
  defining the Python library to which (or from which) data structures
  the data should be converted. Possible options are `"matchms"` or
  `"spectrum_utils"` with `"matchms"` being the default.

## Value

For `r_to_py.Spectra()` and `rspec_to_pyspec()`: Python list of MS data
structures, either `matchms.Spectrum` or
`spectrum_utils.spectrum.MsmsSpectrum` objects. For `pyspec_to_rspec()`:
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
with the MS data of all `matchms.Spectrum` objects in the submitted
`list`.

## Translation of MS data objects

MS data structures can be translated between R and Python using the
`rspec_to_pyspec()` and `pyspec_to_rspec()` functions, or with the
[`r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
method.

- `rspec_to_pyspec()` translates an R
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object into a list of Python MS data objects, which can be, depending
  on parameter `pythonLibrary`, `matchms.Spectrum` objects (for
  `pythonLibrary = "matchms"`, the default) or
  `spectrum_utils.spectrum.MsmsSpectrum` objects (for
  `pythonLibrary = "spectrum_utils"`). Parameter `mapping` allows to
  specify which spectra variables from the `Spectra` object `x` should
  be converted in addition to the peaks data (m/z and intensity values).
  It defaults to `mapping = spectraVariableMapping()` (See the
  respective help below for more information on the variable mapping).
  While being fast, this function first loads all peaks and spectra data
  into memory before translating to Python data structures. A less
  memory intense operation could be to call this function in a loop to
  only load parts of the data at a time into memory.

- `pyspec_to_rspec()` translates a single, or a list of
  `matchms.Spectrum` objects (with parameter
  `pythonLibrary = "matchms"`, the default) or a list of
  `spectrum_utils.spectrum.MsmsSpectrum` objects (with parameter
  `pythonLibrary = "spectrum_utils"`) to a
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object. Parameter `mapping` allows to specify the metadata variables
  that should be translated and mapped in addition to the peaks data.
  The library used to represent the MS data in Python needs to be
  specified with parameter `pythonLibrary`.

- `r_to_py.Spectra()` is equivalent to
  `rspec_to_pyspec(pythonLibrary = "matchms")`. The spectra variables
  that should be converted can be configures with
  `setSpectraVariableMapping()` (see documentation below).

## Mapping of spectra variables (metadata)

Metadata for MS spectra are represented and stored as *spectra
variables* in the R
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
objects. Also Python MS data structures store such metadata along with
the mass peak data. While spectra metadata is thus supported by data
structures in both programming languages, different names and naming
conventions are used. The
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
and `setSpectraVariableMapping()` functions allow to define how the
names of spectra metadata (spectra variables) should be translated
between R and Python. To support also the different naming conventions
used by the Python libraries *matchms* and *spectrum_utils*,
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
defines different mapping schemes for these, using by default the
mapping for *matchms*. Note also that *spectrum_utils* supports only few
selected metadata/spectra variables, so any additional spectra variables
defined by the mapping will be ignored. The
[`r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
and
[`py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
functions will use the selected naming scheme to name the spectra
variables accordingly. Also, only spectra metadata/variables in
[`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
will be translated. The initial mapping is based on this [definition in
matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv).

- `defaultSpectraVariableMapping()`: returns the *default* mapping
  between spectra variables and Python metadata names for the *matchms*
  library.

- [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html):
  returns the currently defined spectra variable mapping as a named
  character vector, with names representing the names of the spectra
  variables in R and elements the respective names of the spectra
  metadata in Python. Use
  [`Spectra::spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  on the `Spectra` object that should be converted with
  [`r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
  to list all available spectra variables.
  [`r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
  and
  [`py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
  for MS data structures will use this default mapping. Calling
  [`spectraVariableMapping()`](https://rdrr.io/pkg/Spectra/man/spectraVariableMapping.html)
  defining also the Python library (e.g.,
  `spectraVariableMapping("matchms")` or
  `spectraVariableMapping("spectrum_utils")`) will return the variable
  mapping for the specified Python library. Optional parameter `x`
  allows to specify a (potentially names) character vector with the
  names of the spectra variables that should in addition be included in
  the mapping.

- `setSpectraVariableMapping()`: sets/replaces the currently defined
  mapping of spectra variable names to Python metadata names. Setting
  `setSpectraVariableMapping(character())` will only convert the mass
  peaks data (m/z and intensity values) but no spectra metadata.

## Author

Michael Witting, Johannes Rainer, Wout Bittremieux, Thomas Naake

## Examples

``` r

## Import a MGF file as a `Spectra` object
library(MsBackendMgf)
library(SpectriPy)
s <- Spectra(
    system.file("extdata", "mgf", "spectra2.mgf", package = "SpectriPy"),
    source = MsBackendMgf())
#> Start data import from 1 files ... 
#> done
s
#> MSn data (Spectra) with 4 spectra in a MsBackendMgf backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         2   2678.94        NA
#> 2         2   2373.51        NA
#> 3         2   2511.03        NA
#> 4         2        NA        NA
#>  ... 17 more variables/columns.

#########################
## Conversion R to Python

## A `Spectra` can be translated to a `list` of `matchms.Spectrum` objects
## using either the `r_to_py()` method or the `rspec_to_pyspec()` function:
s_py <- r_to_py(s)
s_py
#> [Spectrum(precursor m/z=718.36, 203 fragments between 101.1 and 1272.7), Spectrum(precursor m/z=507.75, 755 fragments between 130.2 and 1020.5), Spectrum(precursor m/z=507.74, 85 fragments between 110.1 and 998.4), Spectrum(precursor m/z=633.27, 431 fragments between 51.0 and 753.4)]

## The `s_py` can now be used like any other Python variable within the R
## *reticulate* framework. Below we extract the m/z values of the first
## spectrum
s_py[0]$mz
#> array([ 101.07122,  109.68925,  115.86999,  120.0811 ,  129.06595,
#>         129.1023 ,  130.06534,  130.08633,  130.95578,  131.11859,
#>         132.07968,  136.07567,  141.31848,  142.21034,  147.11273,
#>         155.08185,  159.07605,  159.09169,  160.07516,  160.09506,
#>         170.06009,  171.06319,  175.11871,  186.12335,  187.08667,
#>         197.12814,  199.10913,  199.17963,  200.10286,  201.10281,
#>         201.12341,  204.11336,  204.13414,  211.14496,  213.15979,
#>         214.09769,  214.1535 ,  215.10004,  216.64253,  217.13332,
#>         224.13977,  225.12349,  226.08267,  228.13486,  229.1201 ,
#>         233.16495,  239.08228,  242.14993,  243.10875,  244.10834,
#>         245.35168,  253.09628,  256.10834,  261.1596 ,  273.13486,
#>         282.17978,  284.10266,  285.10632,  287.17343,  289.11084,
#>         301.129  ,  302.1327 ,  305.1811 ,  310.1727 ,  313.18762,
#>         318.15506,  321.1563 ,  326.17813,  328.2328 ,  334.15466,
#>         344.00247,  344.19284,  345.1581 ,  347.17093,  352.19876,
#>         353.1817 ,  356.1928 ,  361.2311 ,  362.20374,  363.207  ,
#>         367.14435,  370.21045,  384.16763,  385.17123,  395.1348 ,
#>         402.18094,  409.20935,  412.16074,  413.16446,  417.17773,
#>         424.21753,  425.2271 ,  429.1876 ,  430.19342,  441.24667,
#>         454.23114,  462.21368,  463.1971 ,  480.2241 ,  481.22632,
#>         484.24506,  488.21954,  497.25027,  498.25397,  506.25055,
#>         508.22073,  512.25824,  516.2635 ,  517.2272 ,  521.2413 ,
#>         523.74426,  525.2426 ,  526.2479 ,  530.26794,  538.2638 ,
#>         542.2722 ,  545.24817,  548.2805 ,  549.2873 ,  550.7827 ,
#>         555.2859 ,  559.2943 ,  559.7948 ,  571.2866 ,  587.2942 ,
#>         599.28723,  599.7854 ,  600.4091 ,  608.29016,  629.78845,
#>         636.28156,  639.34784,  646.328  ,  647.32837,  653.3023 ,
#>         670.3162 ,  677.3326 ,  685.323  ,  692.32416,  695.3513 ,
#>         696.3513 ,  702.3403 ,  703.35065,  707.32135,  709.8466 ,
#>         710.34625,  717.88837,  718.3871 ,  718.88586,  724.3413 ,
#>         729.36005,  730.32825,  748.37524,  757.3715 ,  763.14374,
#>         766.2276 ,  766.38745,  767.39105,  767.7366 ,  783.813  ,
#>         831.4068 ,  835.4187 ,  858.4285 ,  859.41895,  869.44836,
#>         871.4104 ,  876.4352 ,  877.4197 ,  878.41345,  886.42194,
#>         887.41785,  894.44464,  895.44763,  904.4519 ,  905.4429 ,
#>         906.4541 ,  972.48206,  989.50586,  990.4984 , 1007.53156,
#>        1008.5328 , 1017.5188 , 1018.51953, 1037.5978 , 1046.4642 ,
#>        1082.4929 , 1100.5747 , 1117.5706 , 1118.574  , 1119.5591 ,
#>        1135.5819 , 1136.5898 , 1145.5717 , 1214.5857 , 1215.5576 ,
#>        1227.598  , 1232.5897 , 1233.6329 , 1249.6306 , 1250.6298 ,
#>        1260.6073 , 1261.614  , 1272.6572 ])

## Extracting that information from the `Spectra` object in R
s[1]$mz
#> NumericList of length 1
#> [[1]] 101.07122 109.68925 115.86999 120.0811 ... 1260.6073 1261.614 1272.6572

## The `spectraVariableMapping()` defines which spectra variables (metadata)
## should be translated between R and Python:
spectraVariableMapping()
#>                  precursorMz           precursorIntensity 
#>               "precursor_mz"        "precursor_intensity" 
#>              precursorCharge                        rtime 
#>                     "charge"             "retention_time" 
#>              collisionEnergy      isolationWindowTargetMz 
#>           "collision_energy" "isolation_window_target_mz" 
#>                      msLevel 
#>                   "ms_level" 

## The names of that character vector represent the names of the spectra
## variables in R, the elements the name of the metadata variable in Python.
## Below we list the available metadata information from the first
## Spectrum in Python
s_py[0]$metadata
#> {'precursor_mz': 718.36, 'precursor_intensity': nan, 'charge': 2, 'retention_time': 2678.94, 'collision_energy': nan, 'isolation_window_target_mz': nan, 'ms_level': 2}

## `setSpectraVariableMapping()` allows to replace the default mapping
## of variables. Below we e.g. add a new spectra variable to the `Spectra`
## object.
s$new_col <- 1:4

## To translate that variable to Python we need to include it to the
## `spectraVariableMapping()`. Below we define to translate only the
## precursor m/z and the new spectra variable to Python. Be aware that
## `setSpectraVariableMapping()` **globally** sets the default for any
## spectra variable mapping between R and Python. Thus, any subsequent
## calls mapping calls will use the same mapping. It is suggested to
## eventually *restore* the default mapping again after the call or
## use the `rspec_to_pyspec()` function instead, that allows to configure
## the mapping using a parameter `mapping`.
setSpectraVariableMapping(
    c(precursorMz = "precursor_mz", new_col = "new_col"))
s_py <- r_to_py(s)

s_py[0]$metadata
#> {'precursor_mz': 718.36, 'new_col': 1}

## Restoring the global spectra variable mapping configuration to
## the default mapping:
setSpectraVariableMapping(defaultSpectraVariableMapping())

## As an alternative to the `r_to_py()` we can use the `rspec_to_pyspec()`
## function and provide a custom mapping using the `mapping` parameter:
s_py <- rspec_to_pyspec(
    s, mapping = c(precursorMz = "precursor_mz", new_col = "new_col"))

## Convert to MS data objects from the spectrum_utils Python library
s_py2 <- rspec_to_pyspec(
    s, mapping = spectraVariableMapping("spectrum_utils"),
    pythonLibrary = "spectrum_utils")

## Convert the data back to R
pyspec_to_rspec(s_py2, pythonLibrary = "spectrum_utils")
#> MSn data (Spectra) with 4 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1        NA   2678.94        NA
#> 2        NA   2373.51        NA
#> 3        NA   2511.03        NA
#> 4        NA       NaN        NA
#>  ... 16 more variables/columns.

#########################
## Conversion Python to R

## A `list` of `matchms.Spectrum` objects in Python can be translated into
## the corresponding MS data structure in R (i.e. a `Spectra`) object using
## the `pyspec_to_rspec()` function:
res <- pyspec_to_rspec(s_py)
res
#> MSn data (Spectra) with 4 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1        NA        NA        NA
#> 2        NA        NA        NA
#> 3        NA        NA        NA
#> 4        NA        NA        NA
#>  ... 16 more variables/columns.

## All spectra from Python are thus converted into a single `Spectra` object.

## Or providing a custom variable mapping:
res <- pyspec_to_rspec(
    s_py, mapping = c(precursorMz = "precursor_mz", new_col = "new_col"))
res$new_col
#> [1] 1 2 3 4
```
