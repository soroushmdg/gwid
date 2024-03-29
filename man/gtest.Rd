% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesis_test.R
\name{gtest}
\alias{gtest}
\title{perform gtest}
\usage{
gtest(haplotype_structure, ...)
}
\arguments{
\item{haplotype_structure}{object of a class}

\item{...}{other variables}
}
\value{
the output will be a test_snps (data.table) object including 3 columns:
 \dQuote{snp_pos}, \dQuote{case_control}, and \dQuote{value} which is
a p-values.
}
\description{
perform gtest
}
