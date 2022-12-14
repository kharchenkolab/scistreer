% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scistreer.R
\name{run_scistree}
\alias{run_scistree}
\title{Run the scistree workflow}
\usage{
run_scistree(
  P,
  init = "UPGMA",
  ncores = 1,
  max_iter = 100,
  eps = 0.01,
  verbose = TRUE
)
}
\arguments{
\item{P}{matrix Genotype probability matrix (cell x mutation). Each entry is a probability (0-1) that the cell harbors the mutation}

\item{init}{character Initialization strategy; UPGMA or NJ}

\item{ncores}{integer Number of cores to use}

\item{max_iter}{integer Maximum number of iterations}

\item{eps}{numeric Tolerance threshold in likelihood difference for stopping}

\item{verbose}{logical Verbosity}
}
\value{
phylo A maximum-likelihood phylogeny
}
\description{
Run the scistree workflow
}
\examples{
tree_small = run_scistree(P_small)
}
