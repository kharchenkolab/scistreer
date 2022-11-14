#' @import dplyr
#' @importFrom phangorn upgma 
#' @importFrom ape root drop.tip reorder nj
#' @importFrom parallelDist parDist
#' @useDynLib rscistree

#' Run the scistree workflow
#' @param P matrix Genotype probability matrix (cell x mutation). Each entry is a probability (0-1) that the locus is wildtype
#' @param init character Initialization strategy; UPGMA or NJ
#' @param max_iter integer Maximum number of iterations
#' @param eps numeric Tolerance threshold in likelihood difference for stopping
#' @param verbose logical Verbosity
#' @param ncores integer Number of cores to use
#' @export
run_scistree = function(P, init = c('UPGMA', 'NJ'), max_iter = 100, eps = 0.01, ncores = 1, verbose = TRUE) {

    # contruct initial tree
    dist_mat = parDist(rbind(P, 'outgroup' = 1), threads = ncores)

    if (init == 'UPGMA') {

        treeUPGMA = upgma(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')
            
    } else if (init == 'NJ') {

        treeNJ = nj(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')

    } else {
        stop("init has be one of 'UPGMA', 'NJ'")
    }

    tree_list = perform_nni(tree_init, P, ...)
    
}


#' Maximum likelihood tree search via NNI
#' @param tree_init phylo Intial tree
#' @param P matrix Genotype probability matrix
#' @param max_iter integer Maximum number of iterations
#' @param eps numeric Tolerance threshold in likelihood difference for stopping
#' @param verbose logical Verbosity
#' @param ncores integer Number of cores to use
#' @return multiPhylo List of trees corresponding to the rearrangement steps
#' @export
perform_nni = function(tree_init, P, max_iter = 100, eps = 0.01, ncores = 1, verbose = TRUE) {
    
    P = as.matrix(P)
    
    converge = FALSE
    
    i = 1
    max_current = score_tree(tree_init, P)$l_tree
    tree_current = tree_init
    tree_list = list()
    tree_list[[1]] = tree_current
    
    while (!converge & i <= max_iter) {
        
        i = i + 1
        
        ptm = proc.time()
        
        RcppParallel::setThreadOptions(numThreads = ncores)
        
        scores = nni_cpp_parallel(tree_current, P)
        
        if (max(scores) > max_current + eps) {
            max_id = which.max(scores)
            if (max_id %% 2 == 0) {pair_id = 2} else {pair_id = 1}
            tree_current$edge = matrix(nnin_cpp(tree_current$edge, ceiling(max_id/2))[[pair_id]], ncol = 2)
            tree_list[[i]] = tree_current
            tree_list[[i]]$likelihood = max_current = max(scores)
            converge = FALSE
        } else {
            converge = TRUE
        }
        
        runtime = proc.time() - ptm
        
        if (verbose) {
            msg = glue('Iter {i} {max_current} {signif(unname(runtime[3]),2)}s')
            message(msg)
            log_info(msg)
        }
    }
    
    class(tree_list) = 'multiPhylo'
    
    return(tree_list)
}