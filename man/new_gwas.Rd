% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs.R
\name{new_gwas}
\alias{new_gwas}
\title{Data Extraction of IBD and gds data}
\usage{
new_gwas(gds_data = "name.gds", caco = "name.Rda", ...)
}
\arguments{
\item{gds_data}{CoreArray Genomic Data Structure (GDS) with hierarchical structure to store
multiple scalable structure array oriented data sets.}

\item{caco}{A List of case and control subjects.}

\item{...}{argument to pass to plots (if we decided to add graphs)}

\item{ibd_data}{Data frame containing identity by descent segments in phased data. The data is output of refined_ibd
software package.}

\item{i}{chromosome number can be from 1 to 22}
}
\value{
data frames that are derived from inputs in order to overcome computational obstacle
}
\description{
Extract sample ID, IBD's for case-control subjects, snp ID, snp position, sample index, and smp.snp object to use in the pipeline.
}
\examples{
new_gwid("chr1.gds","case_control.Rda","chr1.ibd", i = 1)
}
