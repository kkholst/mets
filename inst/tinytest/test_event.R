## Event function

library("tinytest")

test_event <- function() {

  dat <- data.frame(
    s = mets::Event(runif(10), rbinom(10, 1, 0.5)),
    y = runif(10),
    z = letters[1:10]
  )

  # Event, length.Event
  expect_true(nrow(dat) == 10)
  # c.Evente
  expect_true(length(c(dat$s, dat$s)) == 20)
  # as.matrix.Event
  expect_true(all(dim(as.matrix(dat$s)) == c(10, 2)))

}
test_event()
