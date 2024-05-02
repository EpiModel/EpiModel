context("Load parameters from data.frame")

test_that("Load parameters from data.frame", {
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    "p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p3_1", "1", "numeric", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )

  expect_silent(param <- param.net_from_table(params.df))
  expect_s3_class(param, "param.net")
  expect_type(param, "list")

  expect_silent(param <- param.net(data.frame.parameters = params.df))
  expect_s3_class(param, "param.net")
  expect_type(param, "list")

  # convert back to a `long.param.df`
  param.df_back <- param.net_from_table(params.df) |> param.net_to_table()
  expect_true(all(params.df[c("param", "value", "type")] == param.df_back))

  # wrong column name
  params.df <- dplyr::tribble(~name, ~value, ~type, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))
  params.df <- dplyr::tribble(~param, ~val, ~type, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))
  params.df <- dplyr::tribble(~param, ~value, ~class, "p1", "10", "numeric")
  expect_error(param <- param.net_from_table(params.df))

  # wrong "type" value
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    "p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p3_1", "1", "factor", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )
  expect_error(param <- param.net_from_table(params.df))

  # wrong "param" format
  params.df <- dplyr::tribble(
    ~param, ~value, ~type, ~detail,
    ".p1", "10", "numeric", "foo",
    "p2", "TRUE", "logical", "bar",
    "p_3_1", "1", "numeric", "baz",
    "p3_2", "3", "numeric", "foobar",
    "p4", "tsa", "character", "foobaz"
  )
  expect_error(param <- param.net_from_table(params.df))
})
