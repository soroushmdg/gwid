
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gwid

<!-- badges: start -->
<!-- badges: end -->

GWID (Genome Wide Identity by Descent) is an R-package designed for the
analysis of IBD (Identity by Descent) data, to discover rare alleles
associated with case-control phenotypes. Although Genome Wide
Association Studies (GWAS) successfully reveal numerous common variants
linked to diseases, they exhibit lack of power to identify rare alleles.
To address this limitation, we have developed a pipeline that employs
IBD data (output of refined-IBD software). This methodology encompasses
a sequential process for analyzing the aforementioned data within
isolated populations. The primary objective of this approach is to
enhance the sensitivity of variant detection by utilizing information
from genetically related individuals, thereby facilitating the
identification of causal variants. An overall representation of the
procedural workflow is visually depicted in the accompanying figure.

<div class="figure" style="text-align: right">

<img src="./inst/Figures/final-copy-arrow.png" alt="gwid pipeline" width="70%" />
<p class="caption">
gwid pipeline
</p>

</div>

## Installation

You can install the development version of gwid from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("soroushmdg/gwid")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(gwid)
#> 
#> Attaching package: 'gwid'
#> The following objects are masked from 'package:base':
#> 
#>     print, subset
## basic example code
```

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
