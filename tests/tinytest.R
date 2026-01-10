if (requireNamespace("tinytest", quietly = TRUE)) {
  options(Ncpus = 2)
  Sys.setenv(OMP_THREAD_LIMIT = 2)
  Sys.setenv(OMP_NUM_THREADS = 2)
  Sys.setenv(OPENBLAS_NUM_THREADS = 2)
  Sys.setenv(ARMA_OPENMP_THREADS = 2)
  RcppArmadillo::armadillo_set_number_of_omp_threads(2)
  tinytest::test_package("mets")
}
