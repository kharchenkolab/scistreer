#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppParallel;

/////////////////////////////////////// NNI ////////////////////////////////////////

// Below functions are modified from R-package `ape' by Emmanuel Paradis and Klaus Schliep
// Functions are made thread-safe using RcppArmadillo. 

// https://github.com/KlausVigo/phangorn/blob/master/src/phangorn_utils.cpp
// [[Rcpp::export]]
std::vector<std::vector<int>> allChildrenCPP(const arma::Mat<int> E) {

    arma::Col<int> parent = E.col(0);
    arma::Col<int> children = E.col(1);
    int m = max(parent);

    std::vector<std::vector<int>> out(m);

    for(int i = 0; i<parent.size(); i++) {
        out[parent(i)-1L].push_back(children(i));
    }

    return out;
}

void bar_reorderRcpp(int node, int nTips, const arma::Col<int> & e1,
    const arma::Col<int> & e2, std::vector<int> & neworder, const arma::Col<int> & L,
    const arma::Col<int> & xi, const arma::Col<int> & xj, int & iii)
{
    int i = node - nTips - 1, j, k;

    for (j = xj[i] -1; j >= 0; j--)
        neworder[iii--] = L[xi[i] + j ] + 1;

    for (j = 0; j < xj[i]; j++) {
        k = e2[L[xi[i] + j ]];
        if (k > nTips)
            bar_reorderRcpp(k, nTips, e1, e2, neworder, L, xi, xj, iii);
    }
}

// [[Rcpp::export]]
arma::Mat<int> reorder_rows(arma::Mat<int> x, arma::Col<int> y) {

    // Create an output matrix
    arma::Mat<int> out = x;

    // Loop through each row and copy the data. 
    for (int i = 0; i < y.n_elem; ++i) {
        out.row(i) = x.row(y[i]-1);
    }

    return out;
}

// [[Rcpp::export]]
arma::Mat<int> reorderRcpp(arma::Mat<int> E) {

    int n = E.n_rows;
    int nTips = n/2 + 1;
    int root = nTips + 1;

    arma::Col<int> e1 = E.col(0);
    arma::Col<int> e2 = E.col(1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
    
    arma::Col<int> L(n);
    std::vector<int> neworder(n);
    arma::Col<int> pos(nnode);
    arma::Col<int> xi(nnode);
    arma::Col<int> xj(nnode);
    for (int i = 0; i < n; i++) {
        xj[e1[i] - nTips - 1]++;
    }
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i-1] + xj[i - 1];
    }
    for (int i = 0; i < n; i++) {
        k = e1[i] - nTips - 1;
        j = pos[k]; /* the current 'column' position corresponding to k */
        L[xi[k] + j] = i;
        pos[k]++;
    }

    int iii = n - 1;

    bar_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj, iii);

    E = reorder_rows(E, neworder);

    return E;
}

// Modified from R-package `phangorn' by Klaus Schliep
// [[Rcpp::export]]
std::vector<arma::Mat<int>> nnin_cpp(const arma::Mat<int> E, const int n) {

    arma::Mat<int> E1 = E;
    arma::Mat<int> E2 = E;
    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    E1(ind1, 1) = e2;
    E1(ind2[0], 1) = e1;
    E2(ind1, 1) = e3;
    E2(ind2[1], 1) = e1;

    std::vector<arma::Mat<int>> res(2);

    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);

    return res;
}

// Computing the logQ matrix for the rearranged tree
// [[Rcpp::export]]
arma::vec nnin_score_max(const arma::Mat<int> E, const int n, arma::mat logQ, double L_0) {

    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    arma::vec scores(2);

    arma::mat logQ_1 = logQ;
    arma::mat logQ_2 = logQ;
    
    logQ_1.row(p2-1) = logQ.row(e1-1) + logQ.row(e3-1);
    logQ_2.row(p2-1) = logQ.row(e1-1) + logQ.row(e2-1);

    int m = logQ.n_cols;

    scores[0] = L_0;
    scores[1] = L_0;

    for (int i = 0; i < m; ++i) {
        scores[0] += max(logQ_1.col(i));
        scores[1] += max(logQ_2.col(i));
    }

    return scores;
}

/////////////////////////////////////// Scistree ////////////////////////////////////////


// [[Rcpp::export]]
arma::mat CgetQ(arma::mat logQ, std::vector<std::vector<int>> children_dict, arma::Col<int> node_order){

    int n = node_order.n_rows;
    std::vector<int> children;

    for (int i = 0; i < n; ++i) {
        int node = node_order(i);
        children = children_dict[node-1];
        logQ.row(node-1) = logQ.row(children[0]-1) + logQ.row(children[1]-1);
    }

    return logQ;
}

// [[Rcpp::export]]
arma::mat get_logQ(const arma::Mat<int> E, const arma::mat P) {

    int n = P.n_rows;
    int m = P.n_cols;

    arma::mat logQ(n * 2 - 1, m);

    arma::mat logP_0 = log(1-P);
    arma::mat logP_1 = log(P);

    logQ.rows(0, n-1) = logP_1 - logP_0;

    arma::Col<int> node_order(E.n_rows + 1);
    node_order.rows(0,E.n_rows-1) = E.col(1);
    node_order(E.n_rows) = n+1;
    arma::uvec ids = find(node_order > n);
    node_order = node_order.elem(ids);

    std::vector<std::vector<int>> children_dict = allChildrenCPP(E);

    logQ = CgetQ(logQ, children_dict, node_order);

    return logQ;
}

// [[Rcpp::export]]
double score_tree_cpp(const arma::Mat<int> E, const arma::mat P) {

    int m = P.n_cols;

    double L_0 = accu(log(1-P));

    arma::mat logQ = get_logQ(E, P);

    double l = 0;

    for (int i = 0; i < m; ++i) {
        l += max(logQ.col(i));
    }

    l += L_0;

    return l;
}

struct score_neighbours_max: public Worker {

    // original tree
    const arma::Mat<int> E;

    const arma::mat logQ;

    const double L_0;

    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_max(const arma::Mat<int> E, const arma::mat logQ, const double L_0, NumericVector scores): 
        E(E), logQ(logQ), L_0(L_0), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            arma::vec res = nnin_score_max(E, i+1, logQ, L_0);
            scores[2*i] = res[0];
            scores[2*i+1] = res[1];
        }
    }
};

// Note that the tree has to be already in post-order
// [[Rcpp::export]]
NumericVector nni_cpp_parallel(const List tree, arma::mat P) {
    
    arma::Mat<int> E = tree["edge"];

    // E = reorderRcpp(E);

    int n = E.n_rows/2 - 1;

    arma::mat logQ = get_logQ(E, P);

    double L_0 = accu(log(1-P));

    NumericVector scores(2*n);

    score_neighbours_max score_neighbours_max(E, logQ, L_0, scores);

    parallelFor(0, n, score_neighbours_max);

    return scores;

}

////////////////////// sum instead of max ///////////////////////

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
    // https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
    unsigned int maxi = x.index_max();
    LDOUBLE maxv = x(maxi);
    if (!(maxv > -arma::datum::inf)) {
        return -arma::datum::inf;
    }
    LDOUBLE cumsum = 0.0;
    for (unsigned int i = 0; i < x.n_elem; i++) {
        if ((i != maxi) & (x(i) > -arma::datum::inf)) {
            cumsum += EXPL(x(i) - maxv);
        }
    }
  
    return maxv + log1p(cumsum);
}

// Computing the logQ matrix for the rearranged tree
// [[Rcpp::export]]
arma::vec nnin_score_sum(const arma::Mat<int> E, const int n, arma::mat logQ, arma::rowvec l_0) {

    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    arma::vec scores(2);

    arma::mat logQ_1 = logQ;
    arma::mat logQ_2 = logQ;
    
    logQ_1.row(p2-1) = logQ.row(e1-1) + logQ.row(e3-1);
    logQ_2.row(p2-1) = logQ.row(e1-1) + logQ.row(e2-1);

    int m = logQ.n_cols;

    for (int i = 0; i < m; ++i) {
        scores[0] += logSumExp(logQ_1.col(i) + l_0[i]);
        scores[1] += logSumExp(logQ_2.col(i) + l_0[i]);
    }

    return scores;
}

struct score_neighbours_sum: public Worker {

    // original tree
    const arma::Mat<int> E;

    const arma::mat logQ;

    const arma::rowvec l_0;

    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_sum(const arma::Mat<int> E, const arma::mat logQ, const arma::rowvec l_0, NumericVector scores): 
        E(E), logQ(logQ), l_0(l_0), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            arma::vec res = nnin_score_sum(E, i+1, logQ, l_0);
            scores[2*i] = res[0];
            scores[2*i+1] = res[1];
        }
    }
};


// Note that the tree has to be already in post-order
// [[Rcpp::export]]
NumericVector nni_cpp_parallel_sum(const List tree, arma::mat P) {
    
    arma::Mat<int> E = tree["edge"];

    // E = reorderRcpp(E);

    int n = E.n_rows/2 - 1;

    arma::mat logQ = get_logQ(E, P);

    arma::rowvec l_0 = sum(log(1-P), 0);

    NumericVector scores(2*n);

    score_neighbours_sum score_neighbours_sum(E, logQ, l_0, scores);

    parallelFor(0, n, score_neighbours_sum);

    return scores;

}
