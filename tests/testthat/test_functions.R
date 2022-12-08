library(scistreer2)

test_that("Small input works", {
    tree = run_scistree(P_small, verbose = F)
    expect_equal(round(tree$likelihood, 3), -1606.497)
    expect_equal(1,1)
})
