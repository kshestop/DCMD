out.save.dir <- "~/niagara/med826/Microbiome_classification/github/implementation/"

cancer_data <- read.table(paste0(out.save.dir, "otutable.txt"), header = T)

reads <- apply(cancer_data[,-1],1,sum)
summary(reads)
##calculate the relative abundance
RA <- apply(cancer_data[,-1], 2, function(x) x/reads)
mean.RA <- apply(RA, 2, mean)

cancer <- cancer_data[, c(which(mean.RA>0.001)+1, dim(cancer_data)[2])]

colnames(cancer) <- c(paste0("OTU",which(mean.RA>0.001)),"cls.sample")
cancer$cls.sample <- as.numeric(as.factor(cancer$cls.sample))

save(cancer, file = paste0(out.save.dir, "cancer.RData"))


require(caret)
flds <- createFolds(1:172, k = 10, list = TRUE, returnTrain = FALSE)

for (i in 1:10){
  write.table(flds[[i]], file = paste0(out.save.dir, sprintf("train.set.idx_fold%d.txt", i)), row.names = F, col.names = F)
}