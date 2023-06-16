## metaX


*metaX*: An Automatic and Comprehensive Pipeline for Processing Untargeted Metabolomics Data.

[![Join the chat at https://gitter.im/metaX-metabolomics/Lobby](https://badges.gitter.im/metaX-metabolomics/Lobby.svg)](https://gitter.im/metaX-metabolomics/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[<img src="https://github.com/wenbostar/metaX/blob/master/inst/extdata/metaX_pipeline.PNG" width=800 class="center">](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1579-y)

## Installation

```r
# The latest version of metaX is no longer available at Bioconductor.
# Install the development version from GitHub:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("remotes")
BiocManager::install("wenbostar/metaX")
```

## Usage

Please follow the instruction in the *[metaX document](https://github.com/wenbostar/metaX/blob/master/vignettes/metaX.pdf)*.

## Citation

To cite the `metaX` package in publications, please use:

> Wen B, Mei Z, Zeng C, et al. `metaX`: a flexible and comprehensive software for processing metabolomics data. BMC bioinformatics, 2017, 18(1): 183. DOI: *[10.1186/s12859-017-1579-y](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1579-y)*


## Contribution

Contributions to the package are more than welcome. 


## Function request

I'm working on developing a new version of metaX. If you have any new functions which you want to be added into metaX, please feel free to open an *[issue](https://github.com/wenbostar/metaX/issues)* to describe what you want. I will discuss with you about how to implement that in the new version.
