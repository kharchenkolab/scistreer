library(scistreer)

test_that("Small input works", {
    tree = run_scistree(P_small, verbose = F)
    expect_equal(round(score_tree(tree, P_small)$l_tree, 3), -1606.497)
    expect_equal(1,1)
})
