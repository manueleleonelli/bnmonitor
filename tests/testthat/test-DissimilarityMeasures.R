context("Dissimilarity measures (CD_distance and KL)")
library(bnmonitor)

data("travel")

#test_that("Correct class",{
#  expect_that(CD(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "all")$CD,is_a("data.frame"))
#  expect_that(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "all")$KL,is_a("data.frame"))
#})

#test_that("error/warning with order preserving scheme", {
#  expect_error(CD(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = "all",covariation = "orderp"))
#  expect_warning(CD(bnfit = travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = "all",covariation = "all"))
#  expect_error(KL(travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = "all",covariation = "orderp"))
#  expect_warning(KL(travel,node = "T",value_node = "car",value_parents = c("emp","big"),new_value = "all",covariation = "all"))
#})

#test_that("error/warnings for plots margins or existence", {
#  expect_warning(expect_error(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.50,0.95,0.05),covariation = "orderp")))
#  expect_warning(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.55,0.95,0.05),covariation = "orderp"),"The plot won't be showed since all the values are not possible")
#  expect_warning(expect_error(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.50,0.95,0.05),covariation = "orderp")))
#  expect_warning(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.55,0.95,0.05),covariation = "orderp"),"The plot won't be showed since all the values are not possible")
#})

#test_that("Number of outputs CD",{
#  expect_equal(nrow(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "orderp", FALSE)$CD),20)
#  expect_equal(nrow(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3,covariation = "all", FALSE)$CD),1)
#  expect_equal(nrow(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = c(0.4,0.5),covariation = "all", FALSE)$CD),2)
#  expect_equal(nrow(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.4,0.6,0.05),covariation = "all", FALSE)$CD),5)
#  expect_equal(ncol(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "uniform", FALSE)$CD),2)
#  expect_equal(ncol(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "proportional", FALSE)$CD),2)
#  expect_equal(ncol(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "orderp", FALSE)$CD),2)
#})

#test_that("Correct output CD_distance",{
#  expect_equal(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3,covariation = "orderp", FALSE)$CD[1,2],0.486,tolerance=0.01)
#  expect_true(is.na(CD(bnfit = travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = c(0.4,0.6),covariation = "orderp", FALSE)$CD[2,2]))
#})

#test_that("Number of outputs KL",{
#  expect_equal(nrow(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "orderp")$KL),20)
#  expect_equal(nrow(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3,covariation = "all")$KL),1)
#  expect_equal(nrow(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = c(0.4,0.5),covariation = "all")$KL),2)
#  expect_equal(nrow(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = seq(0.4,0.6,0.05),covariation = "all")$KL),5)
#  expect_equal(ncol(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "uniform")$KL),2)
#  expect_equal(ncol(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "proportional")$KL),2)
#  expect_equal(ncol(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = "all",covariation = "orderp")$KL),2)
#})

#test_that("Correct output KL",{
#  expect_equal(KL(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = 0.3,covariation = "orderp", FALSE)$KL[1,2],0.012,tolerance=0.01)
#  expect_true(is.na(CD(travel,node = "T",value_node = "train",value_parents = c("emp","big"),new_value = c(0.4,0.6),covariation = "orderp", FALSE)$CD[2,2]))
#})
