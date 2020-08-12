context("hello")

test_that("use",{
  expect_output(hello(), "Hello, world!")
})
