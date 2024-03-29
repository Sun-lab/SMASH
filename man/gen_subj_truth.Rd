% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMASH.R
\name{gen_subj_truth}
\alias{gen_subj_truth}
\title{gen_subj_truth}
\usage{
gen_subj_truth(mat_eS, maxLOCI, nCN = NULL)
}
\arguments{
\item{mat_eS}{A subclone configuration matrix pre-defined in R list \code{eS}}

\item{maxLOCI}{A positive integer number of simulated somatic variant calls}

\item{nCN}{A positive integer for the number of allelic copy number pairings 
to sample from. If \code{NULL}, it will be randomly sampled between 1 and 5.}
}
\value{
A list containing the following components:
\describe{
\item{\code{subj_truth}}{dataframe of each variant's simulated minor 
(\code{CN_1}) and major (\code{CN_2}) copy number states, total copy 
number (\code{tCN}), subclone allocation (\code{true_A}), multiplicity 
(\code{true_M}), mutant allele frequency (\code{true_MAF}), and cellular 
prevalence (\code{true_CP})}
\item{\code{purity}}{tumor purity}
\item{\code{eta}}{the product of tumor purity and subclone proportions}
\item{\code{q}}{vector of subclone proportions}
}
}
\description{
Simulates copy number states, multiplicities, allocations
}
