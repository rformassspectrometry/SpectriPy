---
title: 'SpectriPy: Enhancing Cross-Language Mass Spectrometry Data Analysis with R and Python'
tags:
  - Mass Spectrometry  
  - Cross-Language Integration
  - R
  - Python
  - Computational Efficiency
  - Reproducible Workflows
  - Proteomics
  - Metabolomics
authors:
  - name: Marilyn De Graeve
    orcid: 0000-0001-6916-401X
    affiliation: 1
  - name: Wout Bittremieux
    affiliation: 2
  - name: Thomas Naake
    affiliation: 3
  - name: Matthias Anagho-Mattanovich
    affiliation: 4
  - name: Nils Hoffmann
    affiliation: 5
  - name: Pierre Marchal
    affiliation: 6
  - name: Victor Chrone
    affiliation: 7
  - name: Philippine Louail
    affiliation: 1
  - name: Carolin Huber
    affiliation: 8
  - name: Helge Hecht
    affiliation: 9
  - name: Michael Witting
    affiliation: "10, 11"
  - name: Johannes Rainer
    orcid: 0000-0002-6977-7147
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: Institute for Biomedicine, Eurac Research, Bolzano, Italy.
   index: 1
 - name: Department of Computer Science, University of Antwerp, Antwerpen, Belgium.
   index: 2
 - name: Genome Biology Unit, European Molecular Biology Laboratory (EMBL), Heidelberg, Germany.
   index: 3
 - name: Novo Nordisk Foundation Center for Basic Metabolic Research, University of Copenhagen, Copenhagen, Denmark.
   index: 4
 - name: Institute for Bio- and Geosciences (IBG-5), Forschungszentrum Jülich GmbH, Jülich, Germany.
   index: 5
 - name: Department of medical oncology, University of Bern, Bern, Switzerland.
   index: 6
 - name: Alphalyse, Odense, Denmark.
   index: 7
 - name: Department of Exposure Science, Helmholtz Centre for Environmental Research - UFZ, Leipzig, Germany.
   index: 8
 - name: Faculty of Science, Masaryk University, Brno, Czech Republic.
   index: 9
 - name: Metabolomics and Proteomics Core, Helmholtz Munich, Munich, Germany.
   index: 10
 - name: Chair of Analytical Food Chemistry, TUM School of Life Sciences, Technical University of Munich, Munich, Germany.
   index: 11
date: 06 March 2025
bibliography: paper.bib
---

TODO add logo!!! 

![SpectriPy_logo](logo.png){height="150pt"}

# Summary

Mass spectrometry data analysis is a computationally intensive task, with excellent libraries available in both R and Python. However, redundant re-implementation of methods across languages is inefficient. The SpectriPy R package bridges this gap by integrating functions from Python's matchms library. SpectriPy aims to improve documentation, extend functionality, and enhance the translation of MS data structures between R and Python. Additionally, Quarto documents enable seamless integration of R and Python for reproducible workflows. Our results demonstrate improved cross-language MS data analysis, streamlined data exchange, and enhanced workflow reproducibility.

# Statement of need

Over the past decade, tremendous efforts have been made to develop powerful algorithms and excellent data analysis software for mass spectrometry (MS) and metabolomics data analysis. These include, among others, R packages from the RforMassSpectrometry initiative such as Spectra, MsCoreUtils, MetaboAnnotation and CompoundDb, as well as Python libraries like matchms, spectrum_utils, Pyteomics and pyOpenMS. Each of these softwares covers different and in part complementary aspects in the analysis of MS data, but their integration into a single workflow remains, in particular across programming languages, challenging.

Here we present [*SpectriPy*](https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html), an R package that efficiently translates MS data structures between R and Python. By leveraging R's “reticulate” system, SpectriPy enables a seamless cross-language integration allowing R and Python MS data analysis algorithms to be combined within unified analysis workflows. A set of example use cases, implemented as Quarto documents, were developed to demonstrate the advantages and power of this approach.

To summarize, SpectriPy enables and simplifies the integration of R and Python data analysis, empowering data analysts to benefit from the full power of algorithms from both programming languages. Further, software developers can now reuse algorithms across languages rather than re-implementing them, enhancing efficiency and collaboration.

# Introduction

Mass spectrometry (MS) data analysis benefits from a diverse ecosystem of computational tools in both R and Python. R provides robust statistical and bioinformatics capabilities, while Python has extensive machine learning and spectral similarity scoring tools. However, developing and maintaining equivalent functionality in both languages is inefficient.
The Reticulate R package allows execution of Python code within R, while Bioconductor’s Basilisk package facilitates the installation and management of Conda environments for R packages. 
The SpectriPy R package allows integration of Python MS packages into a Spectra-based MS analysis in R. Python functionality is wrapped into R functions allowing a seamless interoperability between the two ecosystems. For example, the integration of the spectral similarity calculations from the Python's matchms library into R. In addition, functions to convert between R's Spectra::Spectra objects and Python's matchms.Spectrum objects are available to the advanced user or developer enabling to create custom functions or workflows on Spectra objects in Python and executing them in R using the reticulate R package.
Moreover, Quarto provides a flexible framework for combining R, Python, and other languages in a single data analysis document.

## Objective

SpectriPy aims to:

1. Improve the combined R and Python MS data analysis by enhancing documentation, incorporating functionality from Python libraries, and improving the translation of MS data structures between R and Python.
2. Provide SpectriPy use-cases in combined R/Python Quarto documents for reproducible and modular MS data analysis.

# Results

## Combined R and Python MS data analysis 

A new R data structure was introduced to facilitate efficient on-the-fly translation of MS data stored in Python, improving data interoperability. 

```r
library(Spectra)
library(SpectriPy)

#... translate func

```

## Enable flexible Python environment

Furthermore, the adoption of Quarto documents facilitated the seamless integration of R and Python within a single workflow. Unlike previous implementations relying on Basilisk, the new approach provides a more flexible and reproducible framework for MS data analysis. This allows researchers to execute algorithms in their native environments while maintaining workflow consistency.
To improve the functionality of SpectriPy, extensive enhancements were made to its codebase.

In terms of code development, the SpectryPy functions were developed to allow users to specify custom Conda environments, providing greater flexibility. Additionally, a version of SpectriPy was designed to operate independently of Basilisk, reducing dependency constraints. 

## Improved MS data translation

A key enhancement was the optimization of MS spectra object handling. Previously, redundant loading of spectral data led to inefficient memory usage. In SpectriPy, spectral objects are now loaded once in the backend, significantly reducing memory overhead and improving performance.

```r
library(Spectra)
library(SpectriPy)

#... ??

```

## Defined example use-cases

To improve the usability of SpectriPy, extensive enhancements were made to its documentation. First, the installation instructions and troubleshooting guidelines were refined to assist users in setting up and utilizing the package more effectively. The improved documentation and additional examples provided clearer guidance for users integrating R and Python-based mass spectrometry analyses. 

Furthermore, new Quarto and R Markdown (Rmd) documents were created to illustrate practical use cases, including general MS data exchange between R and Python, matchMS2DeepScore-based similarity calculations, and fingerprint-based similarity scoring using matchms. The inclusion of new functionalities from Python libraries such as matchms, spectrum_utils, ms2deepscore, spec2vec, Pyteomics, and pyOpenMS enabled a more comprehensive set of tools for MS data analysis.

In terms of code development, additional spectral processing functions from the matchms.filtering module were integrated into SpectriPy, expanding its capabilities. Support for spectrum_utils.Mass.Spectrum classes was also incorporated, ensuring seamless compatibility with Python-based mass spectrometry tools. Alternative implementations of similarity calculation functions were developed too.

## Rust

Preliminary evaluations were also conducted to explore the integration of Rust-based MS data handling into SpectriPy.

```r

#...?

```

# Conclusions

Combining R and Python for MS data analysis using Quarto provides a robust and reproducible framework. By leveraging native algorithm implementations in both languages, SpectriPy minimizes redundant development efforts and enhances cross-language data sharing. Future work will focus on refining MS data translation mechanisms and encouraging users to adopt cross-language workflows rather than reiplementing functions between R and Python.

# Perspectives

Future efforts will focus on refining and consolidating the mechanisms for cross-language MS data representation and translation to ensure seamless interoperability between R and Python. This will involve further optimization of the data structures and improvement of computational efficiency when handling large MS datasets.
In addition to software enhancements, significant emphasis will be placed on the development of educational resources to support the adoption of SpectriPy and cross-language workflows. Quarto and R Markdown documents will be converted into interactive tutorials, enabling users to integrate R and Python seamlessly in their MS data analysis pipelines. These tutorials will be integrated into larger frameworks such as the Metabonaut end-to-end workflow, providing comprehensive learning resources for the computational MS community.
Furthermore, community engagement and collaborative development will be encouraged to extend SpectriPy with additional functionalities and new use cases. Open-source contributions from developers across different domains will facilitate broader adoption and continuous improvement of the package. By fostering a collaborative environment, SpectriPy aims to serve as a versatile tool that supports a diverse range of MS data analysis applications.
Ultimately, the long-term goal is to promote cross-language compatibility and reproducibility in computational mass spectrometry. By leveraging the strengths of both R and Python, SpectriPy will contribute to the advancement of flexible and efficient MS data analysis workflows, reducing redundancy and fostering innovation in the field.

# Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
General information on the package structure and some helpful pointers are given
in the [Development notes](devnotes.md) document. Also, please check the
[coding style
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#coding-style)
and importantly, follow our [code of
conduct](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct).

# References TODO

add in paper.bib !!!

* SpectriPy: https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html
https://github.com/rformassspectrometry/SpectriPy

_TODO refs for Spectra, MsCoreUtils, MetaboAnnotation and CompoundDb, as well as Python libraries like matchms, spectrum_utils, Pyteomics and pyOpenMS etc... + “reticulate” +  + OR http links maybe?_
