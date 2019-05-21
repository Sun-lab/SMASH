# SMASH - Subclone Multiplicity Allocation and Somatic Heterogeneity

This package is designed to cluster somatic mutations called from a tumor sample with a matched normal sample. Each mutation is assumed to lie in a genomic segment of clonal copy number. Each mutation's inferred clonal copy number and the tumor purity estimate is required as input to successfully run the program. This information can be obtained by running algorithms such as ABSOLUTE or ASCAT that derive tumor purity and clonal copy number estimates from SNP Array intensities.

## Installation
```
library(devtools)
install_github("Sun-lab/SMASH")
```

## Subclone configurations
We characterize subclone configurations by square matrices. The number of rows or columns denote the number of subpopulations of cancer cells or subclones. Each row of a configuration matrix represents a somatic variant's possible allocation. For instance, a tumor sample composed of three subclones could have variants arising from one of four allocations. 

A row vector (1,1,1) indicates that the variant is arose in the first subclone and is present in all three subclones. A row vector (0,1,1) indicates that the variant arose in the second subclone and is present in the third clone. A row vector (0,1,0) also indicates that the variant arose in the second sublone but was not passed on the third subclone. Lastly, a row vector (0,0,1) indicates the variant arose and is only present in the third subclone. 

Three subclones can arise from two possible configuration matrices. Let A, B, and C denote the three subclones arising in that order. One possibility is that the three subclones arose linearly, in order words A -> B -> C. A second possibility is that three subclones arose in a branching tree where A -> B and A -> C. In SMASH `eS` is an R list, where `eS[[3]]` contains the two matrices for the respective subclone configurations. The default matrices provided with `eS` reflect two of our underlying assumptions. First, we assume each tumor sample arose from one clone and second, at most two subclones arose from a parental subclone.

## Simulate and Cluster
The package is also designed to simulate somatic mutations assuming known tumor purity and clonal copy number per mutation. The code below will simulate 200 mutations from a tumor sample with three subclones evolving from a branching tree (`eS[[3]][[2]]`) with a mean sequencing depth of 500.

```
set.seed(2)
truth = gen_subj_truth(mat_eS = eS[[3]][[2]],maxLOCI = 200)
truth$q # cancer cell proportions
sum(-truth$q * log(truth$q)) # entropy
table(truth$subj_truth$true_A) # frequency of variant allocations
unique(truth$subj_truth[,c("CN_1","CN_2")])	# unique set of clonal copy numbers
dat = gen_ITH_RD(DATA = truth$subj_truth,RD = 500)
dat = data.frame(dat,truth$subj_truth,stringsAsFactors = FALSE)
dat[1:10,]
```

`truth$purity` contains the simulated tumor purity. `tAD` and `tRD` denote the tumor variant's alternate and reference read count, respectively. `CN_1` and `CN_2` denote the variant's corresponding minor and major allelic copy numbers, respectively. `tCN` denotes the total copy number, the sum of `CN_1` and `CN_2`.

With the code below, we cluster these simulated mutations by supplying all proposed subclone configurations contained in `eS` and run 50 initial randomizations of cancer cell proportions for each configuration.
```
smash_out = grid_ITH_optim(my_data = dat[,c("tAD","tRD","CN_1","CN_2","tCN")],
		my_purity = truth$purity,
		list_eS = eS,
		trials = 50,
		max_iter = 4e3,
		my_epsilon = 1e-6)
smash_out$GRID[order(-smash_out$GRID$BIC),]
```

The R data.frame `smash_out$GRID` contains the clustered results across all feasible subclone configurations provided by the user. Each row denotes a feasible model where columns `cc` and `kk` are used to index the subclone figuration `eS[[cc]][[kk]]`. `ms` denotes the model size or number of parameters estimated for a given model. `LL`, `AIC`, and `BIC` denote the log likelihood, AIC, and BIC of each model, respectively. `q` and `entropy` denote the estimated cancer cell proportions and corresponding entropy calculated as `sum(-q * log(q))`.

`smash_out$INFER` is a R list where each nested element is a data.frame containing the inferred allocation (`infer_A`) and multiplicity (`infer_M`) of each mutation for a given model. To observe the inferred results for the 3rd model or row of `smash_out$GRID`, run `smash_out$INFER[[3]]`.


