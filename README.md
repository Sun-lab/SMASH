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

## Table of Contents

* [What is this for?](#what-is-this-for)
* [Installation](#installation)
* [Vignette](#vignette)
* [Workflow](#workflow)
* [Citation](#citation)
* [Contact](#contact)
* [FAQs](#faqs)

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
pandoc = Sys.getenv("RSTUDIO_PANDOC")

cran_packs = c("Rcpp","RcppArmadillo","devtools","smarter","SMASH")
req_packs = c(cran_packs)

for(pack in req_packs){
	
	chk_pack = tryCatch(find.package(pack),
		error = function(ee){NULL})
	
	if( !is.null(chk_pack) ){
		library(pack,character.only = TRUE)
		next
	}
	
	if( pack %in% cran_packs ){
		install.packages(pack,dependencies = TRUE)
	}
	
}

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
array{{SNP array data}} --> cnaCall{{Copy Number Algorithm}}
cnaCall --> cnaEst{{Purity, Allelic Copy Numbers}}
cnaEst & vcf --> SMASH{{SMASH}}

%% Class definitions
classDef myred fill:#f44336,stroke:#f3f6f4,stroke-width:2px
classDef myblue fill:#19daf8,stroke:#f3f6f4,stroke-width:2px
classDef mygreen fill:#79d50d,stroke:#f3f6f4,stroke-width:2px
classDef mymagenta fill:#fc9ffc,stroke:#f3f6f4,stroke-width:2px
classDef myyellow fill:#f6fa13,stroke:#f3f6f4,stroke-width:2px
classDef myorange fill:#f89d3e,stroke:#f3f6f4,stroke-width:2px

%% Assign classes to nodes
class tbam myred
class nbam myblue
class array myyellow
class fasta mygreen
class caller,cnaCall,SMASH myorange
class vcf,cnaEst mymagenta
```

## Citation

Little, P., [Lin, D.Y.](https://sph.unc.edu/adv_profile/danyu-lin-phd/), 
[Sun, W.](https://github.com/sunway1999) (2019). 
Associating somatic mutations to clinical outcomes: a pan-cancer 
study of survival time. *Genome medicine,* 11(1), 1-15. 
[[HTML](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0643-9), 
[PDF](https://genomemedicine.biomedcentral.com/track/pdf/10.1186/s13073-019-0643-9.pdf), 
[Supplement](https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-019-0643-9/MediaObjects/13073_2019_643_MOESM1_ESM.pdf)]

## Contact

* [Feel free to reach out](mailto:pllittle321@gmail.com?subject=SMASH:%20Q%26A&body=Dear%20Dr.%20Little,%0A%0A%0A)

## FAQs
