context("Covariation schemes")
library(bnmonitor)

data("travel")

test_that("Correct class",{
  expect_that(uniform_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),is_a(c("bn.fit")))
  expect_that(proportional_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),is_a(c("bn.fit")))
  expect_that(orderp_covar(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3),is_a(c("bn.fit")))
  expect_true(is.na(orderp_covar(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.65)))
})

test_that("Correct output CD_distance",{
  expect_equal(uniform_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = 0.65)[["T"]][["prob"]][t(c("car","emp","big"))],0.65)
  expect_equal(proportional_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = 0.65)[["T"]][["prob"]][t(c("car","emp","big"))],0.65)
  expect_equal(orderp_covar(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3)[["T"]][["prob"]][t(c("train","emp","big"))],0.3)
  expect_true(is.na(orderp_covar(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.65)))
})

test_that("error when invalid node name", {
  expect_error(uniform_covar(bnfit = travel,node = "t",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),"Invalid input for node")
  expect_error(proportional_covar(bnfit = travel,node = "t",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),"Invalid input for node")
  expect_error(orderp_covar(bnfit = travel,node = "t",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),"Invalid input for node")
})

test_that("error when not well specified parents of the node", {
  expect_error(uniform_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("ep","big"),new_value = 0.65))
  expect_error(proportional_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("ep","big"),new_value = 0.65))
  expect_error(orderp_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("ep","big"),new_value = 0.65))
})

test_that("error when not correct number of parents of the node", {
  expect_error(uniform_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp"),new_value = 0.65),"Invalid length of value_parents")
  expect_error(proportional_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("ep"),new_value = 0.65),"Invalid length of value_parents")
  expect_error(orderp_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("big"),new_value = 0.65),"Invalid length of value_parents")
})

test_that("error: order preserving scheme, varying last parameter", {
  expect_error(orderp_covar(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = 0.65),"can't work by varying the last parameter")
})
