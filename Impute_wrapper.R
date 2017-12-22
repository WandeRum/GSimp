# Imputation Wrapper ------------------------------------------------------
require(missForest)
require(impute)
require(magrittr)
require(imputeLCMD)
source('MVI_global.R')


RF_wrapper <- function(data, ...) {
  result <- missForest(data, ...)[[1]]
  return (result)
}

kNN_wrapper <- function(data, ...) {
  result <- data %>% data.matrix %>% impute.knn(., ...) %>% extract2(1)
  return(result)
}

SVD_wrapper <- function(data, K = 5) {
  data_sc_res <- scale_recover(data, method = 'scale')
  data_sc <- data_sc_res[[1]]
  data_sc_param <- data_sc_res[[2]]
  result <- data_sc %>% impute.wrapper.SVD(., K = K) %>% 
    scale_recover(., method = 'recover', param_df = data_sc_param) %>% extract2(1)
  return(result)
}

Mean_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- mean(x, na.rm = T)
    x
  })
  return(result)
}

Median_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- median(x, na.rm = T)
    x
  })
  return(result)
}

HM_wrapper <- function(data) {
  result <- data
  result[] <- lapply(result, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/2
    x
  })
  return(result)
}

Zero_wrapper <- function(data) {
  result <- data
  result[is.na(result)] <- 0
  return(result)
}

QRILC_wrapper <- function(data, ...) {
  result <- data %>% log %>% impute.QRILC(., ...) %>% extract2(1) %>% exp
  return(result)
}