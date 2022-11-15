# ScisTreeR
`scistreer` is an R/Rcpp implementation of the [ScisTree](https://doi.org/10.1093/bioinformatics/btz676) algorithm. It improves scalability of the algorithm via RcppParallel and is applicable to very large single-cell datasets (>10,000 cells).

# Installation
```
devtools::install_github('https://github.com/kharchenkolab/scistreer')
```

# Citations

For the original publication, please refer to:
> Yufeng Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Volume 36, Issue 3, 1 February 2020, Pages 742–750, https://doi.org/10.1093/bioinformatics/btz676

If you would like to cite this package, please use:
> Teng Gao, Evan Biederstedt, Peter Kharchenko (2022).
ScisTreeR: Speeding up the ScisTree Algorithm via RcppParallel. R
package version 1.0.0. https://github.com/kharchenkolab/scistreer
