# SMASH - Subclone Multiplicity Allocation and Somatic Heterogeneity

This package is designed to cluster somatic mutations called from a tumor sample with a matched normal sample. Each mutation is assumed to lie in a genomic segment of clonal copy number. Each mutation's inferred clonal copy number and the tumor purity estimate is required as input to successfully run the program. This information can be obtained by running algorithms such as ABSOLUTE or ASCAT that derive tumor purity and clonal copy number estimates from SNP Array intensities.

## Installation

```R
library(devtools)
install_github("Sun-lab/SMASH")
```

## Vignette

To see the vignette, follow the code below. The vignette contains background knowledge and code to perform simulation, optimization, interpretation, and visualization.

```R
library(SMASH)
vignette(package = "SMASH",topic = "intro")
```

## Citation

Little, P., Lin, D.Y., Sun, W. (2019). Associating somatic mutations to clinical outcomes: a pan-cancer study of survival time. *Genome medicine,* 11(1), 1-15. [[HTML](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0643-9), [PDF](https://genomemedicine.biomedcentral.com/track/pdf/10.1186/s13073-019-0643-9.pdf)]

