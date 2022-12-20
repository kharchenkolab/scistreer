<!-- badges: start -->
[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/scistreer.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/scistreer)
[![CRAN status](https://www.r-pkg.org/badges/version/scistreer)](https://cran.r-project.org/package=scistreer)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/scistreer)](https://cran.r-project.org/package=scistreer)
<!-- badges: end -->

# ScisTreeX
Fast maximum-likelihood phylogeny inference from noisy single-cell data using the 'ScisTree' algorithm [(Wu, Bioinformatics 2019)](https://academic.oup.com/bioinformatics/article/36/3/742/5555811). 'scistreex' provides an 'R' interface and improves speed via 'Rcpp' and 'RcppParallel', making the method applicable to massive single-cell datasets (>10,000 cells).

# Installation
You can install the github version via `devtools`:
```R
devtools::install_github('https://github.com/kharchenkolab/scistreer/tree/devel')
```
# Usage
Within R, you only need to supply a genotype probability matrix (cell x mutation), where each entry is the probability that the cell harbors the mutation. For example,


```R
treeML = run_scistree(P_example, ncores = 8, init = 'UPGMA', verbose = FALSE)
```
The output maximum likelihood tree is an `ape::phylo` object. You can visualize the output and the probability matrix as follows:
```R
plot_phylo_heatmap(treeML, P_example)
``` 

<p align="center">
<img src="https://user-images.githubusercontent.com/13375875/202533038-3513f6ba-454f-4bd2-9808-70e3442808cd.png" width="600">
</p>

# Citations

For the original publication, please refer to:

> Yufeng Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Volume 36, Issue 3, 1 February 2020, Pages 742–750, https://doi.org/10.1093/bioinformatics/btz676
