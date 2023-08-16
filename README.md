
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gwid

<!-- badges: start -->
<!-- badges: end -->

GWID (Genome Wide Identity by Descent) is an R-package designed for the
analysis of IBD (Identity by Descent) data, to discover rare alleles
associated with case-control phenotype. Although Genome Wide Association
Studies (GWAS) successfully reveal numerous common variants linked to
diseases, they exhibit lack of power to identify rare alleles. To
address this limitation, we have developed a pipeline that employs IBD
data (output of refined-IBD software). This methodology encompasses a
sequential process for analyzing the aforementioned data within isolated
populations. The primary objective of this approach is to enhance the
sensitivity of variant detection by utilizing information from
genetically related individuals, thereby facilitating the identification
of causal variants. An overall representation of the pipeline is
visually depicted in the following figure.

<div class="figure" style="text-align: center">

<img src="./inst/Figures/final-copy-arrow.png" alt="gwid pipeline" width="70%" />
<p class="caption">
gwid pipeline
</p>

</div>

## Usage

The `gwid` package receives four types of inputs: SNP panel information,
IBD information, haplotype data, and data concerning subjects
categorized as cases and controls. The SNP panel data is derived from
the output of the
[SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
package in the form of a **gds** file. The IBD data takes the form of
tabulated data produced by the [Refined
IBD](https://faculty.washington.edu/browning/refined-ibd.html) software.
Haplotype data comes from the output of the
[Beagle](http://faculty.washington.edu/browning/beagle/beagle.html),
while information about case and control subjects is represented using
an R list.

## Installation

You can install the development version of `gwid` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("soroushmdg/gwid")
```

## Example

The following example is for a SNP panel data from the Marshfield
Clinic. subjects in case group has Rheumatoid Arthritis (RA).

This is a basic example which shows you how to solve a common problem:

``` r
library(gwid)
#> 
#> Attaching package: 'gwid'
#> The following objects are masked from 'package:base':
#> 
#>     print, subset
caco <- gwid::case_control(case_control_rda = case_control_data)
pieces <- gwid::build_gwas(gds_data = genome_data,caco = caco,gwas_generator = TRUE)
myphase <- gwid::build_phase(phased_vcf = phase_data,caco = caco)
myregion2 <- gwid::build_gwid(ibd_data = ibd_data,gwas = pieces)
p <- plot(myregion2,ly = FALSE)
p
```

<img src="man/figures/README-example-1.png" width="100%" />

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
#summary(cars)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
