## Ensembles for Time Series Forecasting

This package provides methods and tools for a dynamic ensemble learning method
for time series forecasting tasks. It is versatile and easy to use framework
that calls learning methods from other R libraries, such as nnet, kernlab
or forecast. 

Besides a new dynamic methods for time series forecasting other ensemble approaches are implemented, such as stacking or blending.

The package will be submitted to CRAN after extensive testing.

### Installation

Type **devtools::install_github("vcerqueira/tsensembler")** in your R console to install the package.

### Improvements - todo list

* extensive testing;
* (consistently) updating model weights after forecasting;
* implement parallelmap;
* sanity checks;
* auto re-train base models;

## Dev Environment

* R version 3.2.3
* Platform: x86_64-pc-linux-gnu (64-bit)
* Running under: Ubuntu 14.04.3 LTS

## Contact

Any bug report or suggesting please contact me at vmac@inesctec.pt
