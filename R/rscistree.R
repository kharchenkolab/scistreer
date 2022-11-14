#' @import dplyr
#' @import tidygraph
#' @import stringr
#' @importFrom igraph vcount ecount E V V<- E<- 
#' @importFrom phangorn upgma 
#' @importFrom ape root drop.tip nj ladderize
#' @importFrom parallelDist parDist
#' @useDynLib rscistree
NULL

#' Run the scistree workflow
#' @param P matrix Genotype probability matrix (cell x mutation). Each entry is a probability (0-1) that the locus is wildtype
#' @param init character Initialization strategy; UPGMA or NJ
#' @param max_iter integer Maximum number of iterations
#' @param eps numeric Tolerance threshold in likelihood difference for stopping
#' @param verbose logical Verbosity
#' @param ncores integer Number of cores to use
#' @export
run_scistree = function(P, init = 'UPGMA', ncores = 1, ...) {

    # contruct initial tree
    dist_mat = parDist(rbind(P, 'outgroup' = 1), threads = ncores)

    if (init == 'UPGMA') {

        tree_init = upgma(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')
            
    } else if (init == 'NJ') {

        tree_init = nj(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')

    } else {
        stop("init has be one of 'UPGMA', 'NJ'")
    }

    tree_list = perform_nni(tree_init, P, ncores = ncores, ...)

    tree_final = tree_list[[length(tree_list)]]

    tree_final = ladderize(tree_final)

    # get the branch lengths
    tree_final$edge.length = rep(1, nrow(tree_final$edge))

    return(tree_final)
    
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
    max_current = -Inf
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
            msg = paste0('Iter ', i, ' ', max_current, ', ', signif(unname(runtime[3]),2), 's')
            message(msg)
        }
    }
    
    class(tree_list) = 'multiPhylo'
    
    return(tree_list)
}

#' score a tree based on maximum likelihood
#' @param tree phylo object
#' @param P genotype probability matrix
#' @param get_l_matrix whether to compute the whole likelihood matrix
#' @return list Likelihood scores of a tree
#' @keywords internal
score_tree = function(tree, P, get_l_matrix = FALSE) {
 
    tree = reorder(tree, order = 'postorder')
        
    n = nrow(P)
    m = ncol(P)

    logQ = matrix(nrow = tree$Nnode * 2 + 1, ncol = m)

    logP_0 = log(P)
    logP_1 = log(1-P)
    
    node_order = c(tree$edge[,2], n+1)
    node_order = node_order[node_order > n]
    
    logQ[1:n,] = logP_1 - logP_0

    children_dict = allChildrenCPP(tree$edge)

    logQ = CgetQ(logQ, children_dict, node_order)

    if (get_l_matrix) {
        l_matrix = sweep(logQ, 2, colSums(logP_0), FUN = '+')
        l_tree = sum(apply(l_matrix, 2, max))
    } else {
        l_matrix = NULL
        l_tree = sum(apply(logQ, 2, max)) + sum(logP_0)
    }
    
    return(list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix))
    
}


#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @return list Mutation placements
#' @export
annotate_tree = function(tree, P) {
   
    sites = colnames(P)
    n = nrow(P)
    tree_stats = score_tree(tree, P, get_l_matrix = TRUE)

    l_matrix = as.data.frame(tree_stats$l_matrix)

    colnames(l_matrix) = sites
    rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

    # mutation assignment on nodes
    mut_nodes = data.frame(
            site = sites,
            node_phylo = apply(l_matrix, 2, which.max),
            l = apply(l_matrix, 2, max)
        ) %>%
        mutate(
            name = ifelse(node_phylo <= n, tree$tip.label[node_phylo], paste0('Node', node_phylo - n))
        ) %>%
        group_by(name) %>%
        summarise(
            site = paste0(sort(site), collapse = ','),
            n_mut = n(),
            l = sum(l),
            .groups = 'drop'
        )

    gtree = tree %>%
        ladderize() %>%
        as_tbl_graph() %>%
        mutate(
            leaf = node_is_leaf(),
            root = node_is_root(),
            depth = bfs_dist(root = 1),
            id = row_number()
        )

    # leaf annotation for edges
    gtree = gtree %>%
        activate(edges) %>%
        select(-any_of(c('leaf'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, leaf),
            by = c('to' = 'id')
        )

    # annotate the tree
    gtree = mut_to_tree(gtree, mut_nodes)

    return(gtree)
}

#' Transfer mutation assignment onto a single-cell phylogeny
#'
#' @param gtree tbl_graph The single-cell phylogeny
#' @param mut_nodes dataframe Mutation placements
#' @return tbl_graph A single-cell phylogeny with mutation placements
#' @keywords internal
mut_to_tree = function(gtree, mut_nodes) {
   
    # transfer mutation to tree
    gtree = gtree %>%
        activate(nodes) %>%
        select(-any_of(c('n_mut', 'l', 'site', 'clone'))) %>%
        left_join(
            mut_nodes %>%
                mutate(n_mut = unlist(purrr::map(str_split(site, ','), length))) %>%
                select(name, n_mut, site),
            by = 'name'
        ) %>%
        mutate(n_mut = ifelse(is.na(n_mut), 0, n_mut))

    # get branch length
    gtree = gtree %>% 
        activate(edges) %>%
        select(-any_of(c('length'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, length = n_mut),
            by = c('to' = 'id')
        ) %>%
        mutate(length = ifelse(leaf, pmax(length, 0.2), length))

    # label genotype on nodes
    node_to_mut = gtree %>% activate(nodes) %>% data.frame() %>% {setNames(.$site, .$id)}

    gtree = gtree %>%
        activate(nodes) %>%
        mutate(GT = unlist(
            map_bfs(node_is_root(),
            .f = function(path, ...) { paste0(na.omit(node_to_mut[path$node]), collapse = ',') })
            ),
            last_mut = unlist(
                map_bfs(node_is_root(),
                .f = function(path, ...) { 
                    past_muts = na.omit(node_to_mut[path$node])
                    if (length(past_muts) > 0) {
                        return(past_muts[length(past_muts)])
                    } else {
                        return('')
                    }
                })
            )
        ) %>%
        mutate(GT = ifelse(GT == '' & !is.na(site), site, GT))

    # preserve the clone ids
    if ('GT' %in% colnames(mut_nodes)) {
        gtree = gtree %>% activate(nodes) %>%
            left_join(
                mut_nodes %>% select(GT, clone),
                by = 'GT'
            )
    }
    
    return(gtree)
}

#' Convert a single-cell phylogeny with mutation placements into a mutation graph
#'
#' @param gtree tbl_graph The single-cell phylogeny
#' @param mut_nodes dataframe Mutation placements
#' @return igraph Mutation graph
#' @export
get_mut_graph = function(gtree) {

    mut_nodes = as.data.frame(gtree) %>%
        filter(!is.na(site)) %>%
        distinct(name, site)

    G = gtree %>%
        activate(nodes) %>%
        arrange(last_mut) %>%
        convert(to_contracted, last_mut) %>%
        mutate(label = last_mut, id = 1:n()) %>%
        as.igraph

    G = label_edges(G)

    V(G)$node = G %>%
        igraph::as_data_frame('vertices') %>%
        left_join(
            mut_nodes %>% dplyr::rename(node = name),
            by = c('label' = 'site')
        ) %>%
        pull(node)

    G = G %>% transfer_links()

    return(G)

}


#' Annotate the direct upstream or downstream mutations on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
label_edges = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(from_label = label, id),
            by = c('from' = 'id')
        ) %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(to_label = label, id),
            by = c('to' = 'id')
        ) %>%
        mutate(label = paste0(from_label, '->', to_label))
    
    E(G)$label = edge_df$label
    E(G)$from_label = edge_df$from_label
    E(G)$to_label = edge_df$to_label
    
    return(G)
}


#' Annotate the direct upstream or downstream node on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
transfer_links = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(from_node = node, id),
                by = c('from' = 'id')
            ) %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(to_node = node, id),
                by = c('to' = 'id')
            )

    E(G)$from_node = edge_df$from_node
    E(G)$to_node = edge_df$to_node
    
    return(G)
}

#' Convert the phylogeny from tidygraph to phylo object
#'
#' @param gtree tbl_graph The single-cell phylogeny
#' @return phylo The single-cell phylogeny
#' @export
to_phylo = function(gtree) {
    
    phytree = gtree %>% ape::as.phylo()
    phytree$edge.length = gtree %>% activate(edges) %>% data.frame() %>% pull(length)
    
    n_mut_root = gtree %>% activate(nodes) %>% filter(node_is_root()) %>% pull(n_mut)
    phytree$root.edge = n_mut_root
    
    return(phytree)
}