library(scistreer)

test_that("Small input works", {
    tree = run_scistree(P_small, verbose = T)
    expect_equal(tree$likelihood, -1606.5)
})
