context("C++ samplers")


test_that(
  "Simple-MPT C++ sampling",
  {
    testthat::skip_on_cran()
    testthat::skip_on_ci()
    res <- simpleMPT(
      eqnfile = system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS"),
      data = TreeBUGS::arnold2013[1:10, ],
      restrictions = list("D1 = D2 = D3", "d1 = d2"),
      n.iter = 1e3,
      n.burnin = 5e2,
      n.thin = 2
    )
    expect_s3_class(object = res, class = "simpleMPT")
  }
)

test_that(
  "Beta-MPT C++ sampling",
  {
    testthat::skip_on_cran()
    testthat::skip_on_ci()
    res <- betaMPTcpp(
      eqnfile = system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS"),
      data = TreeBUGS::arnold2013[1:10, ],
      restrictions = list("D1 = D2 = D3", "d1 = d2"),
      n.iter = 1e3,
      n.burnin = 5e2,
      n.thin = 1
    )
    expect_s3_class(object = res, class = "betaMPT")
  }
)
