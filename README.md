<div align="left">
<a href=""><img src="https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink&label=SMASH" height="80" /></a>
</div>

<!-- badges: start -->
![C++](https://img.shields.io/badge/C++-%2300599C.svg?style=square&logo=c%2B%2B&logoColor=gold)
![R](https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink)
![CRAN status](https://www.r-pkg.org/badges/version/SMASH)
[![DOI](https://zenodo.org/badge/DOI/10.1186/s13073-019-0643-9.svg)](https://doi.org/10.1186/s13073-019-0643-9)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/github/languages/code-size/Sun-lab/SMASH.svg)](https://github.com/Sun-lab/SMASH)
[![](https://img.shields.io/github/last-commit/Sun-lab/SMASH.svg)](https://github.com/Sun-lab/SMASH/commits/master)
<!-- badges: end -->

## What is this for?

This package is designed to cluster somatic mutations called from a 
tumor sample with a matched normal sample. Each mutation is assumed 
to lie in a genomic segment of clonal copy number. Each mutation's 
inferred clonal copy number and the tumor purity estimate is required 
as input to successfully run the program. This information can be 
obtained by running algorithms such as ABSOLUTE or ASCAT that derive 
tumor purity and clonal copy number estimates from SNP Array intensities.

<p align="center">
<img src="images/ith_configs.PNG" width="50%" />
<p align="center"><em>Visualizing subclone configurations.</em></p>
</p>

## Installation

<details>

<summary>Click to expand!</summary>

Copy/paste the following code into R/RStudio for **SMASH** installation.

```R
all_packs = as.character(installed.packages()[,1])
pandoc = Sys.getenv("RSTUDIO_PANDOC")
build_vign = !is.null(pandoc) && file.exists(pandoc)

if( !("smarter" %in% all_packs) ){
	stop("Check https://github.com/pllittle/smarter for installation")
}

library(smarter)
smarter::smart_packDeps(
	cran_packs = c("Rcpp","RcppArmadillo","devtools"),
	github_packs = c("Sun-lab/SMASH"),
	pandoc = pandoc,
	build_vign = build_vign)

```

</details>

## Vignette

To see the vignette, follow the code below. The vignette contains background 
knowledge and code to perform simulation, optimization, interpretation, 
and visualization.

```R
library(SMASH)
vignette(package = "SMASH",topic = "intro")
```

## Workflow

```mermaid
flowchart LR

%% Nodes and directions
fasta{{reference.fasta}} & tbam{{tumor.bam}} & nbam{{normal.bam}} --> caller{{Variant Caller}}
caller --> vcf{{somatic.vcf}}
array{{SNP array}} --> cnaCall{{Copy Number Algorithm}}
cnaCall & vcf --> SMASH

%% Class definitions

%% Assign classes to nodes
```

## Citation

Little, P., [Lin, D.Y.](https://sph.unc.edu/adv_profile/danyu-lin-phd/), 
[Sun, W.](https://github.com/sunway1999) (2019). 
Associating somatic mutations to clinical outcomes: a pan-cancer 
study of survival time. *Genome medicine,* 11(1), 1-15. 
[[HTML](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0643-9), 
[PDF](https://genomemedicine.biomedcentral.com/track/pdf/10.1186/s13073-019-0643-9.pdf), 
[Supplement](https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-019-0643-9/MediaObjects/13073_2019_643_MOESM1_ESM.pdf)]

