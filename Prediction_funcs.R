# Packages ----------------------------------------------------------------
require(randomForest)
require(glmnet)
require(rpart)
require(FNN)

lm_pred <- function(x, y) {
  data <- data.frame(y=y, x)
  model <- lm(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

rlm_pred <- function(x, y) {
  data <- data.frame(y=y, x)
  model <- rlm(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

rf_pred <- function(x, y, ntree=200, ...) {
  model <- randomForest(x=x, y=y, ntree=ntree, ...)
  y_hat <- predict(model, newdata=x)
  return(y_hat)
}

glmnet_pred <- function(x, y, alpha=.5, lambda=.01) {
  x_mat <- as.matrix(x)
  model <- glmnet(x=x_mat, y=y, alpha=alpha, lambda=lambda)
  y_hat <- predict(model, newx=x_mat)[, 1]
  return(y_hat)
}

rpart_pred <- function(x, y) {
  data <- data.frame(y=y, x)
  model <- rpart(y ~ ., data=data)
  y_hat <- predict(model, newdata=data)
  return(y_hat)
}

knn_pred <- function(x, y) {
  model <- knn.reg(train=x, y=y, k=5)
  y_hat <- model$pred
  return(y_hat)
}
