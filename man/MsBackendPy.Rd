% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MsBackendPython.R
\name{MsBackendPy}
\alias{MsBackendPy}
\alias{backendInitialize,MsBackendPy-method}
\alias{length,MsBackendPy-method}
\alias{spectraVariables,MsBackendPy-method}
\alias{spectraData,MsBackendPy-method}
\alias{peaksData,MsBackendPy-method}
\alias{$,MsBackendPy-method}
\alias{spectraVariableMapping<-,MsBackendPy-method}
\alias{spectraVariableMapping<-,Spectra-method}
\alias{reindex}
\title{A MS data backend for MS data stored in Python}
\usage{
\S4method{backendInitialize}{MsBackendPy}(
  object,
  pythonVariableName = character(),
  spectraVariableMapping = defaultSpectraVariableMapping(),
  pythonLibrary = c("matchms", "spectrum_utils"),
  ...,
  data
)

\S4method{length}{MsBackendPy}(x)

\S4method{spectraVariables}{MsBackendPy}(object)

\S4method{spectraData}{MsBackendPy}(object, columns = spectraVariables(object), drop = FALSE)

\S4method{peaksData}{MsBackendPy}(object, columns = c("mz", "intensity"), drop = FALSE)

\S4method{$}{MsBackendPy}(x, name)

\S4method{spectraVariableMapping}{MsBackendPy}(object) <- value

\S4method{spectraVariableMapping}{Spectra}(object) <- value

reindex(object)
}
\arguments{
\item{object}{A \code{MsBackendPy} object.}

\item{pythonVariableName}{For \code{backendInitialize()}: \code{character(1)} with the
name of the variable/Python attribute that contains the list of
\code{matchms.Spectrum} objects with the MS data.}

\item{spectraVariableMapping}{For \code{backendInitialize()}: named \code{character}
with the mapping between spectra variable names and (\code{matchms.Spectrum})
metadata names. See \code{\link[=defaultSpectraVariableMapping]{defaultSpectraVariableMapping()}} for more
information and details.}

\item{pythonLibrary}{For \code{backendInitialize()}: \code{character(1)} specifying
the Python library used to represent the MS data in Python. Can be
either \code{pythonLibrary = "matchms"} (the default) or
\code{pythonLibrary = "spectrum_utils"}.}

\item{...}{Additional parameters.}

\item{data}{For \code{backendInitialize()}: \code{DataFrame} with the full MS data
(peaks data and spectra data). Currently not supported.}

\item{x}{A \code{MsBackendPy} object}

\item{columns}{For \code{spectraData()}: \code{character} with the names of
columns (spectra variables) to retrieve. Defaults to
\code{spectraVariables(object)}. For \code{peaksData()}: \code{character} with the
names of the peaks variables to retrieve.}

\item{drop}{For \code{spectraData()} and \code{peaksData()}: \code{logical(1)} whether,
when a single column is requested, the data should be returned as a
\code{vector} instead of a \code{data.frame} or \code{matrix}.}

\item{name}{For \code{$}: \code{character(1)} with the name of the variable to
retrieve.}

\item{value}{Replacement value(s).}
}
\value{
See description of individual functions for their return values.
}
\description{
The \code{MsBackendPy} allows to access MS data stored as \code{matchms.Spectrum}
or \code{spectrum_utils.spectrum.MsmsSpectrum} objects from the
\href{https://github.com/matchms/matchms}{\emph{matchms}} respectively
\href{https://github.com/bittremieux-lab/spectrum_utils}{\emph{spectrum_utils}} Python
library directly from R. The MS data (peaks data or spectra variables) are
translated on-the-fly when accessed. Thus, the \code{MsBackendPy} allows a
seamless integration of Python MS data structures into \code{\link[Spectra:Spectra]{Spectra::Spectra()}}
based analysis workflows.

The \code{MsBackendPy} object is considered \emph{read-only}, i.e. it does not provide
functionality to replace the peaks data from R. However, it is possible to
directly change the data in the referenced Python variable.
}
\details{
The \code{MsBackendPy} keeps only a reference to the MS data in Python (i.e. the
name of the variable in Python) as well as an index pointing to the
individual spectra in Python but no other data. Any data requested from
the \code{MsBackendPy} is accessed and translated on-the-fly from the Python
variable. The \code{MsBackendPy} is thus an interface to the MS data, but not
a data container. All changes to the MS data in the Python variable
(performed e.g. in Python) immediately affect any \code{MsBackendPy} instances
pointing to this variable.

Special care must be taken if the MS data structure in Python is subset or
its order is changed (e.g. by another process). In that case it might be
needed to re-index the backend using the \code{reindex()} function:
\code{object <- reindex(object)}. This will update (replace) the index to the
individual spectra in Python which is stored within the backend.
}
\note{
As mentioned in the \emph{details} section the MS data is completely stored in
Python and the backend only references to this data through the name of
the variable in Python. Thus, each time MS data is requested from the
backend, it is retrieved in its \strong{current} state.
If for example data was transformed or metadata added or removed in the
Python object, it immediately affects the \code{Spectra}/backend.
}
\section{\code{MsBackendPy} methods}{


The \code{MsBackendPy} supports all methods defined by the \code{\link[Spectra:MsBackend]{Spectra::MsBackend()}}
interface for access to MS data. Details on the invidual functions can also
be found in the main documentation in the \emph{Spectra} package (i.e. for
\code{\link[Spectra:MsBackend]{Spectra::MsBackend()}}). Here we provide information for functions with
specific properties of the backend.
\itemize{
\item \code{backendInitialize()}: initializes the backend with information from the
referenced Python variable (attribute). The name of this attribute,
ideally stored in the associated Python session, is expected to be
provided with the \code{pythonVariableName} parameter. The optional
\code{spectraVariableMapping} parameter allows to provide additional, or
alternative, mapping of \code{Spectra}'s \emph{spectra variables} to metadata in the
\code{matchms.Spectrum} objects. See \code{\link[=defaultSpectraVariableMapping]{defaultSpectraVariableMapping()}} (the
default) for more information. Parameter \code{pythonLibrary} must be used
to specify the Python library representing the MS data in Python. It can
be either \code{pythonLibrary = "matchms"} (the default) or
\code{pythonLibrary = "spectrum_utils"}. The function returns an initialized
instance of \code{MsBackendPy}.
\item \code{peaksData()}: extracts the peaks data matrices from the backend. Python
code is applied to the data structure in Python to
extract the \emph{m/z} and intensity values as a list of (numpy) arrays. These
are then translated into an R \code{list} of two-column \code{numeric} matrices.
Because Python does not allow to name columns of an array, an additional
loop in R is required to set the column names to \code{"mz"} and \code{"intensity"}.
\item \code{spectraData()}: extracts the spectra data from the backend. Which spectra
variables are translated and retrieved from the Python objects depends on
the backend's \code{spectraVariableMapping()}. All metadata names defined are
retrieved and added to the returned \code{DataFrame} (with eventually missing
\emph{core} spectra variables filled with \code{NA}).
\item \code{spectraVariables()}: retrieves available spectra variables, which include
the names of all metadata attributes in the \code{matchms.Spectrum} objects
and the \emph{core} spectra variables \code{\link[Spectra:spectraData]{Spectra::coreSpectraVariables()}}.
\item \verb{spectraVariableMapping<-}: replaces the \code{spectraVariableMapping} of the
backend (see \code{\link[=setSpectraVariableMapping]{setSpectraVariableMapping()}} for details and description
of the expected format).
}
}

\section{Additional helper and utility functions}{

\itemize{
\item \code{reindex()}: update the internal \emph{index} to match \code{1:length(object)}.
This function is useful if the original data referenced by the backend was
subset or re-ordered by a different process (or a function in Python).
}
}

\examples{

## Loading an example MGF file provided by the SpectriPy package.
## As an alternative, the data could also be imported directly in Python
## using:
## import matchms
## from matchms.importing import load_from_mgf
## s_p = list(load_from_mgf(r.fl))
library(Spectra)
library(MsBackendMgf)

fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
s <- Spectra(fl, source = MsBackendMgf())
s

## Translating the MS data to Python and assigning it to a variable
## named "s_p" in the (*reticulate*'s) `py` Python environment. Assigning
## the variable to the Python environment has performance advantages, as
## any Python code applied to the MS data does not require any data
## conversions.
py_set_attr(py, "s_p", rspec_to_pyspec(s))


## Create a `MsBackendPy` representing an interface to the data in the
## "s_p" variable in Python:
be <- backendInitialize(MsBackendPy(), "s_p")
be

## Create a Spectra object which this backend:
s_2 <- Spectra(be)
s_2

## Available spectra variables: these include, next to the *core* spectra
## variables, also the names of all metadata stored in the `matchms.Spectrum`
## objects.
spectraVariables(s_2)

## Get the full peaks data:
peaksData(s_2)

## Get the peaks from the first spectrum
peaksData(s_2)[[1L]]

## Get the full spectra data:
spectraData(s_2)

## Get the m/z values
mz(s_2)

## Plot the first spectrum
plotSpectra(s_2[1L])


########
## Using the spectrum_utils Python library

## Below we convert the data to a list of `MsmsSpectrum` object from the
## spectrum_utils library.
py_set_attr(py, "su_p", rspec_to_pyspec(s,
    spectraVariableMapping("spectrum_utils"), "spectrum_utils"))

## Create a MsBackendPy representing this data. Importantly, we need to
## specify the Python library using the `pythonLibrary` parameter and
## ideally also set the `spectraVariableMapping` to the one specific for
## that library.
be <- backendInitialize(MsBackendPy(), "su_p",
    spectraVariableMapping = spectraVariableMapping("spectrum_utils"),
    pythonLibrary = "spectrum_utils")
be

## Get the peaks data for the first 3 spectra
peaksData(be[1:3])

## Get the full spectraData
spectraData(be)

## Extract the precursor m/z
be$precursorMz
}
\author{
Johannes Rainer and the EuBIC hackathon team
}
