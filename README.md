# ScisTreeR
`scistreer` is an R/Rcpp implementation of the [ScisTree](https://doi.org/10.1093/bioinformatics/btz676) algorithm. It improves scalability of the algorithm via RcppParallel and is applicable to very large single-cell datasets (>10,000 cells).

# Installation
```
devtools::install_github('https://github.com/kharchenkolab/scistreer')
```
# Usage
Within R, you only need to supply a genotype probability matrix (cell x mutation), where each entry is the probability that the cell does *not* harbor the mutation.

```
treeML = run_scistree(P, ncores = 8, init = 'UPGMA', verbose = F)
```
The output maximum likelihood tree is a `ape::phylo` object. You can visualize the output and the probablity matrix as follows:
```
plot_phylo_heatmap(treeML, P)
``` 

<p align="center">
<img src="https://user-images.githubusercontent.com/13375875/201975144-1f45ea93-c0e3-4148-893b-438ad7102dcb.png" width="600">
</p>

# Benchmark
![image](https://user-images.githubusercontent.com/13375875/201978296-e6cbabf2-1cd9-4c92-9e70-0ca2082b53e0.png)

# Citations

For the original publication, please refer to:
> Yufeng Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Volume 36, Issue 3, 1 February 2020, Pages 742–750, https://doi.org/10.1093/bioinformatics/btz676

If you would like to cite this package, please use:
> Teng Gao, Evan Biederstedt, Yufeng Wu, Peter Kharchenko (2022).
ScisTreeR: Speeding up the ScisTree Algorithm via RcppParallel. R
package version 1.0.0. https://github.com/kharchenkolab/scistreer
