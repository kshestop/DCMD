glmnet_fit <- function(X_tr, Y_tr, X_ts, Y_ts, a){
  ## read information about data
  n <- nrow (X_tr) ## numbers of obs
  
  ## find number of observations in each group
  nos_g <- as.vector(tapply(rep(1,n), INDEX = Y_tr, sum))
  
  if (any(nos_g < 2)) stop ("Less than 2 cases in some group")
  
  ## choosing the best lambda
  cvfit <- cv.glmnet(x = X_tr, y = Y_tr, alpha = a, nlambda = 500, family = "multinomial", type.measure = "class")
  lambda <- cvfit$lambda[which.min(cvfit$cvm)]
  cat("The best lambda chosen by CV:", lambda, "\n")
  ## fit model with the best lambda
  fit <- glmnet (x = X_tr, y= Y_tr, alpha = a, nlambda = 500, family = "multinomial")
  
  #betas <- coef(fit, s = lambda)
  
  ## predicting for new cases
  if (is.null (X_ts)) {
    return (betas)
  } 
  else {
    pred_matrix <- predict(fit, newx = X_ts, s =lambda, type="response")

    class_pred <- predict(fit, newx = X_ts, s = lambda, type = "class")

    list(pred_matrix=pred_matrix, class_pred=class_pred, er = mean(class_pred != Y_ts))  
  }
}


##train KNN
cv.knn <- function(dis.mat, train.set.idx, training, n.gen){
  flds <- createFolds(1:size, k = 10, list = TRUE, returnTrain = FALSE)
  er <- vector()
  seq <- seq(3,25, by=2)
  for (k in seq){
    predict.lab <- vector()
    for (i in 1:10){
      train.idx <- train.set.idx[unlist(flds[-i])]
      test.idx <- train.set.idx[flds[[i]]]
      dis.matrix <- dis.mat[test.idx, train.idx] + t(dis.mat[train.idx, test.idx])
      predict.lab[flds[[i]]] <- as.numeric(knn_test_function(training[unlist(flds[-i]), 1:n.gen], training[flds[[i]], 1:n.gen], dis.matrix, training[, n.gen+1], k=k)) 
    }
    er[(k-1)/2] <- mean(predict.lab == as.numeric(training[, n.gen+1]))
  }
  k.opt <- seq[which.min(er)]
  list(predict=cbind(seq, er), opt.k=k.opt)
}

knn_fit <- function(dis.mat, training, test, n.gen){
  opt.k <- cv.knn(dis.mat, train.set.idx, training, n.gen)$opt.k
  dist.mat <- dis.mat[test.set.idx2, train.set.idx] + t(dis.mat[train.set.idx, test.set.idx2])
  predict.knn <- knn_test_function(training[, 1:n.gen], test[, 1:n.gen], distance = dist.mat, labels= training[, n.gen+1], k=opt.k)
  list(prediction=predict.knn, er=mean(predict.knn!=test[, n.gen+1]))
}


cls.med <- function(data, metric){
  pam <- pam(data, k=1, diss = F, metric = metric)$medoids
  med.pam <- as.numeric(pam)
}


