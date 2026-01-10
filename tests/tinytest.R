if (requireNamespace("tinytest", quietly = TRUE)) {
  options(Ncpus = 2)
  RcppArmadillo::armadillo_set_number_of_omp_threads(2)
  tinytest::test_package("mets")
}
