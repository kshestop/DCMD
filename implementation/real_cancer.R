# args=(commandArgs(TRUE))
# print(args)
# for(k in 1:length(args)){
#   eval(parse(text=args[[k]]))
# }
##uncomment above when you submit in the cluster


fold <- 1
##comment above when you submit in the cluster


n.dir.rep <- 300 #number of boostrap replicates


library(nloptr)
library(doMC)
registerDoMC(cores=40)
library(tools)
###paralell tools
library(snow); library(iterators); library(foreach); library(doSNOW)
library(parallel)
library(caret)
library(e1071)
library(randomForest)
library(glmnet)
library(clue)
library(class)
library(gbm)
library(cluster)
library(FastKNN)
library(pamr)

path <- "~/niagara/med826/Microbiome_classification/"

source(paste0(path, 'scinet_part3_ex_fns_v1.r'))
source(paste0(path, 'scinet_part2_new_function_v1.r'))
source(paste0(path, 'scinet_estim_functions_v1.r'))
source(paste0(path, "utilities.R"))

out.save.dir <- "~/niagara/med826/Microbiome_classification/implementation/"

boot.save.dir <- paste0(path, "bootstrap/")
if (!exists(boot.save.dir)) dir.create(boot.save.dir, showWarnings = F)
f.dir <- paste0(boot.save.dir, "fold", fold)
if (!exists(f.dir)) dir.create(f.dir, showWarnings = F)

##load simulated dataset
load(paste0(out.save.dir, "cancer.RData"))

n.obs <- dim(cancer)[1]
n.otu <- dim(cancer)[2]-1

## create training and test
test.set.idx <- 1:n.obs
test.set.idx2 <- scan(paste0(out.save.dir, sprintf("train.set.idx_fold%d.txt", fold)))
train.set.idx <- test.set.idx[-test.set.idx2]

cancer_train <- cancer[train.set.idx, ]

# feature selection
p.value <- vector()
for (i in 1:n.otu){
  p.value[i] <- wilcox.test(x=cancer_train[which(cancer_train[, n.otu+1]==1),i], y=cancer_train[which(cancer_train[, n.otu+1]==2),i])$p.value
}

fdr <- p.adjust(p.value, method = "fdr")

p.index <- which(fdr < 0.05)

n.select <- length(p.index)

##create gen.data.array to be the same as simulation
gen.data.array <- apply(cancer[, p.index], 2, as.numeric)
dim(gen.data.array)

cls.sample <- as.numeric(as.factor(cancer[, n.otu+1]))
table(cls.sample)

summary(gen.data.array[which(cls.sample==1),])
summary(gen.data.array[which(cls.sample==2),])


ncount.list <- list()
ncount.list[[1]] <- table(cls.sample)[1]
ncount.list[[2]] <- table(cls.sample)[2]

## normalized the count data;
rds.min <- 10000; rds.max <- 20000
rds.sample <- lapply(as.list(1:2), function(x){as.integer(runif(ncount.list[[x]], rds.min, rds.max))})
rds.sample.norm <- lapply(rds.sample, function(x){x/mean(x)})
cls.names <- 1:2
gen.data.array.norm <- matrix(0, nrow=dim(gen.data.array)[1], ncol=dim(gen.data.array)[2])
gen.data.array.norm[which(cls.sample==1),] <- apply(gen.data.array[which(cls.sample==1),], 2, function(x){x/rds.sample[[1]]})
gen.data.array.norm[which(cls.sample==2),] <- apply(gen.data.array[which(cls.sample==2),], 2, function(x){x/rds.sample[[2]]})

summary(gen.data.array.norm)

gen.data.trans <- apply(gen.data.array.norm, 2, function(x) log(x+1))



###############################################################
#Define the parameters for the mixture distribution for the fit
###############################################################
lc.alp.cut <- 8
alp.max.q <- 0.85

subdiv.list <- c(0.05, 0.10, 0.30, 0.55) #percentage of data where overlap is in the ovlp list
ovlp.lvl.list <- c(0.65, 0.50, 0.35, 0.20) #cutoff for overlap to skip a distribution

## Define the model type to run ##
#Default - overspecify the model
alp.l <- c(1, 1, seq(2, lc.alp.cut-1)) #low count dist. parameters
bet.l <- c(2, 1, rep(1, lc.alp.cut-2))
ct.brk <- 20 #number of bins for the data

dist.opt <- 3 #the optimal distance to use - l2-discrete CDF
tot.mod.dist <- 3 #total distances to calculate

#precalculate the conditional matrices for all mixtures (select later for each dataset)
alg.indx <- 1 #algorithm to use (BFGS)
sig.level.round <- 2 #round the scaling factor for the conditional matrix
tmat.vals <- unique(round(seq(0.01, 8.6, by=0.01), sig.level.round))


init.param.list <- list(lc=lc.alp.cut, alpmax=alp.max.q, subdiv=subdiv.list, ovlp=ovlp.lvl.list, a=alp.l, b=bet.l, ct=ct.brk, dist=dist.opt, tdist=tot.mod.dist, alg=alg.indx, sig=sig.level.round, tmat=tmat.vals, ndir=n.dir.rep)

#gen.data.array: counts, norm. reads, class
#n.dir.rep - number of bootstrap replicates
#n.gen - number of separate OTU

#lc.alp.cut - cutoff for integer intervals
#alp.max.q - quantile for the cutoff
#subdiv.list - percentage of data where overlap is in the overlap list
#ovlp.lvl.list - cutoff for overlap to skip a distribution

#alp.l, bet.l - mix. dist'n parameters
#ct.brk - number of clusters to use

###################

n.obs <- dim(gen.data.array)[1]
n.gen <- dim(gen.data.array)[2]


for (jj in 1:n.select) {
  print(c('Iteration:', jj))
  dat.iter <- cbind(gen.data.array[,jj], unlist(rds.sample.norm))
  data.train <- dat.iter[train.set.idx,] ## the training set only
  dat.run.dist <- dat.iter[test.set.idx,] ## the training + test sets to run the distance computation on
  source(paste0(path, 'scinet_part3_implementation_v1_ex_meth_run.r'))
}

####################
#Manhattan Distance:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(gen.data.array.norm)[1], dim(gen.data.array.norm)[1]))
mx.k <- ceiling(exp(log(dim(gen.data.array.norm)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
  if (k.indx == mx.k) {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
  } else {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
  }
  comp.array <- array(0, dim=c(tot.idx, dim(gen.data.array.norm)[2], 2))
  comp.array[,,1] <- gen.data.array.norm[matrix.rc.list[idx.st:idx.en,1],]
  comp.array[,,2] <- gen.data.array.norm[matrix.rc.list[idx.st:idx.en,2],]
  dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- apply(abs(comp.array[,,1] - comp.array[,,2]), 1, sum)
}
dist.mat.array.mht <- dist.calc.pw.mat
print(Sys.time() - t.t1)

####################
#Euclidean Distance:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(gen.data.array.norm)[1], dim(gen.data.array.norm)[1]))
mx.k <- ceiling(exp(log(dim(gen.data.array.norm)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
  if (k.indx == mx.k) {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
  } else {
    idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
  }
  comp.array <- array(0, dim=c(tot.idx, dim(gen.data.array.norm)[2], 2))
  comp.array[,,1] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,1],]+1)
  comp.array[,,2] <- log(gen.data.array.norm[matrix.rc.list[idx.st:idx.en,2],]+1)
  dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- sqrt(apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum))
}
dist.mat.array.eucl <- dist.calc.pw.mat
print(Sys.time() - t.t1)

## Sum the distances over the different OTU
dist.mat.array.sum <- array(0, dim=c(length(test.set.idx), length(test.set.idx),3))
dist.mat.med.sum <- array(0, dim=c(length(test.set.idx), length(cls.names),3))


for (jj in 1:n.select) {
  model.fit.init.fn <- paste0(f.dir, '/model_train_fit_otu_iteration', jj, '.RData')
  load(model.fit.init.fn)
  dist.mat.array.sum[,,1] <- dist.mat.array.sum[,,1] + dist.mat.array[,,1] #discrete pdf
  dist.mat.array.sum[,,2] <- dist.mat.array.sum[,,2] + dist.mat.array[,,2] #cont cdf
  
  dist.mat.med.sum[,,1] <- dist.mat.med.sum[,,1] + med.store.array[,,1]
  dist.mat.med.sum[,,2] <- dist.mat.med.sum[,,2] + med.store.array[,,2]
}

predict_pdf <- apply(dist.mat.med.sum[test.set.idx2,,1], 1, which.min)
predict_con <- apply(dist.mat.med.sum[test.set.idx2,,2], 1, which.min) 

gen.data.array.new <- data.frame(gen.data.array, as.factor(cls.sample))
colnames(gen.data.array.new) <- c(paste0("OTU", 1:n.gen),"label")

training <- gen.data.array.new[train.set.idx,]
test <- gen.data.array.new[test.set.idx2,]
test_label <- as.numeric(test[, n.gen+1])
test_label
cls.sample[test.set.idx2]

er_kmeans_m1 <- mean(predict_pdf != test_label)
er_kmeans_m1
er_kmeans_m2 <- mean(predict_con != test_label)
er_kmeans_m2

##find the medoid and calculate the distance of each point to the class medoid.
med.store.array <- array(0, dim=c(n.obs, length(cls.names), 2)) ##store all the dataset to the class medoid
for (cls.idx in 1:length(cls.names)) {
  cls.index.train <- train.set.idx[which(cls.sample[train.set.idx] == cls.names[cls.idx])]
  cls.index <- which(cls.sample[test.set.idx] == cls.names[cls.idx])
  
  med.eucl <- cls.med(gen.data.trans[cls.index.train, 1:n.gen], metric = "euclidean")
  
  med.mht <- cls.med(gen.data.array.norm[cls.index.train, 1:n.gen],  metric = "manhattan")
  
  med.store.array[, cls.idx, 1] <- sqrt(apply((gen.data.trans - t(replicate(dim(gen.data.trans)[1], med.eucl)))^2, 1, sum))
  med.store.array[, cls.idx, 2] <- apply(abs(gen.data.array.norm - t(replicate(dim(gen.data.array.norm)[1], med.mht))), 1, sum)
  
}
pred_kmeans_eucl <- apply(med.store.array[test.set.idx2,,1], 1, which.min)
pred_kmeans_mht <- apply(med.store.array[test.set.idx2,,2], 1, which.min)
er_kmeans_eucl <- mean(pred_kmeans_eucl != test_label)
er_kmeans_mht <- mean(pred_kmeans_mht != test_label)

size=length(train.set.idx)
### KNN
t.t1 <- Sys.time()
knn_m1 <- knn_fit(dist.mat.array.sum[,,1], training, test, n.gen)
knn_m2 <- knn_fit(dist.mat.array.sum[,,2], training, test, n.gen)
knn_eucl <- knn_fit(dist.mat.array.eucl, training, test, n.gen)
knn_mht <- knn_fit(dist.mat.array.mht, training, test, n.gen)
pred_knn_m1 <- knn_m1$prediction
pred_knn_m2 <- knn_m2$prediction
pred_knn_eucl <- knn_eucl$prediction
pred_knn_mht <- knn_mht$prediction
er_knn_m1 <- knn_m1$er
er_knn_m2 <- knn_m2$er
er_knn_eucl <- knn_eucl$er
er_knn_mht <- knn_mht$er
print(Sys.time() - t.t1)


##compare with other methods
fit_formula <- as.formula(paste0("label", "~", paste(paste0("OTU", 1:n.gen), collapse = "+")))

### random forest
trControl <- trainControl(method = "cv", number = 10, search = "grid")
tuneGrid <- expand.grid(.mtry = c(4:10))
t.t1 <- Sys.time()
modellist <- list()
for (ntree in c(500, 1000, 100)) {
  fit <- train(fit_formula, training, method = "rf", metric= "Accuracy", tuneGrid=tuneGrid, trControl=trControl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}
results <- resamples(modellist)
rf_index <- which.max(summary(results)$statistics$Accuracy[,"Max."])
rf_fit <- modellist[[rf_index]]
predict_rf <- predict(rf_fit, test[,1:n.gen])
er_rf <- mean(predict_rf != test_label)
print(Sys.time() - t.t1)
er_rf

##gradient boosting
grid <- expand.grid(n.trees = seq(500, 1000, 100), interaction.depth = 2, shrinkage = .01, n.minobsinnode = 20)
t.t1 <- Sys.time()
fit_gb <- train(fit_formula, training, method = 'gbm', trControl = trControl, tuneGrid = grid, metric = "Accuracy", verbose = F)
print(Sys.time() - t.t1)
predict_gb <- predict(fit_gb, test[,1:n.gen])
er_gb <- mean(predict_gb != test_label)
er_gb


gen.data.array.norm.new <- data.frame(gen.data.array.norm , factor(cls.sample))
colnames(gen.data.array.norm.new) <- c(paste0("OTU", 1:n.gen),"label")

training_norm <- gen.data.array.norm.new[train.set.idx,]
test_norm <- gen.data.array.norm.new[-train.set.idx,]

##lasso with log(x+1) transformed relative abundons
X_tr <- apply(model.matrix(fit_formula, data=training_norm)[, -1], 2, function(x) log(x+1))
Y_tr <- as.numeric(training[, "label"])
X_ts <- apply(model.matrix(fit_formula, data=test_norm)[, -1], 2, function(x) log(x+1))
Y_ts <- as.numeric(test[, "label"])
t.t1 <- Sys.time()
lasso_fit1 <- glmnet_fit(X_tr, Y_tr, X_ts, Y_ts, a=1)
er_lasso1 <- lasso_fit1$er
print(Sys.time() - t.t1)

##lasso with log(x+1) count data
X_tr2 <- apply(model.matrix(fit_formula, data=training)[, -1], 2, function(x) log(x+1))
X_ts2 <- apply(model.matrix(fit_formula, data=test)[, -1], 2, function(x) log(x+1))
t.t1 <- Sys.time()
lasso_fit2 <- glmnet_fit(X_tr2, Y_tr, X_ts2, Y_ts, a=1)
er_lasso2 <- lasso_fit2$er
print(Sys.time() - t.t1)

##ridge with log(x+1) transformed relative abundons
t.t1 <- Sys.time()
ridge_fit1 <- glmnet_fit(X_tr, Y_tr, X_ts, Y_ts, a=0)
er_ridge1 <- ridge_fit1$er
print(Sys.time() - t.t1)

##lasso with log(x+1) count data
t.t1 <- Sys.time()
ridge_fit2 <- glmnet_fit(X_tr2, Y_tr, X_ts2, Y_ts, a=0)
er_ridge2 <- ridge_fit2$er
print(Sys.time() - t.t1)

####pamr####
train.x <- t(gen.data.array.norm[train.set.idx,])
train.y <- factor(cls.sample[train.set.idx])
  
test.x <- t(gen.data.array.norm[test.set.idx2,])
test.y <- factor(cls.sample[test.set.idx2])
  
training <- list(x=train.x, y=train.y)
  
pamr.fit <- pamr.train(training)
pamr.pred <- pamr.predict(pamr.fit, test.x, threshold=0, type = "class")
pamr.er <- mean(pamr.pred != test.y)


train.x2 <- gen.data.array[train.set.idx,]
train.y2 <- factor(cls.sample[train.set.idx])

test.x2 <- gen.data.array[test.set.idx2,]
test.y2 <- factor(cls.sample[test.set.idx2])

training2 <- data.frame(x=train.x2, y=train.y2)
testing2 <- data.frame(x=test.x2, y=test.y2)
trControl <- trainControl(method = "cv", number = 10, search = "grid")
svmfit <- train(y~., data=training2, method="svmLinear", trControl = trControl,  preProcess = c("center","scale"))
svm.predict <- predict(svmfit, testing2)
er.svm=mean(svm.predict != test.y2)


nb <- multinomial_naive_bayes(x=train.x2, y=train.y2)
nb.predict <- predict(nb, test.x2)
er.nb=mean(nb.predict != test.y2)


save.er <- cbind(er_kmeans_pdf=er_kmeans_m1, er_kmeans_con=er_kmeans_m2, er_kmeans_eucl=er_kmeans_eucl, er_kmeans_mht=er_kmeans_mht,
                 er_knn_pdf=er_knn_m1, er_knn_cdf=er_knn_m2, er_knn_con=er_knn_m3, er_knn_eucl=er_knn_eucl, er_knn_mht=er_knn_mht,
                 er_rf=er_rf, er_gb=er_gb, er_lasso1=er_lasso1, er_lasso2=er_lasso2, er_ridge1=er_ridge1, er_ridge2=er_ridge2, pamr.er=pamr.er,er_svm=er.svm, er_nb=er.nb)

save.er

save.list <- list(predict_kmeans=cbind(predict_pdf,predict_con,pred_kmeans_eucl, pred_kmeans_mht),
                  predict_knn=cbind(pred_knn_m1,pred_knn_m2,pred_knn_m3,pred_knn_eucl,pred_knn_mht),
                  predict_rf=predict_rf, predict_gb=predict_gb, lasso_fit1=lasso_fit1, lasso_fit2=lasso_fit2, ridge_fit1=ridge_fit1, ridge_fit2=ridge_fit2,  pamr.pred=pamr.pred, predict_svm=svm.predict, predict_nb=nb.predict, test.set.idx=test.set.idx2, test_label=test_label, save.er)
save(save.list, file=paste0(out.save.dir, "predict_fold", fold, ".RData"))


