if (requireNamespace("tinytest", quietly = TRUE)) {
  options(Ncpus = 2)
  tinytest::test_package("mets")
}
