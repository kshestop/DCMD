path <- "M:/Microbiome_classification/"

#scenario 
dat.ver.idx <- 1

library(reshape2)
source(paste0(path, 'updated_code/scinet_part3_ex_fns_v1.r'))
source(paste0(path, 'updated_code/scinet_part2_new_function_v1.r'))
source(paste0(path, 'updated_code/scinet_estim_functions_v1.r'))
out.save.dir <- paste0(path, "simulateddata/scenario", dat.ver.idx, "/")


#use the weights for the new data
#calculate for the training set - add the individual new data after estimation of weights
#calculate the models, calculate the individual weights and all
#calculate the pairwise distances and identify the KNN for each sample

#####################
## data simulation ##
#define the range
#define the distributions along the range
#sample the weights for distributions along the range
#sample the data given the distributions and the range

#####################
#the simulation setup
otu.sim <- c(25, 50, 100)

i.otu <- 1
j.sample <- 1
zp <- array(data=NA, dim = c(3,25,100))
for (rep in 1:100){
  sim.vars.sv <- paste0("_i", i.otu, "_j", j.sample, "_ver", dat.ver.idx, "_rep", rep)
  
  n.otu.sim <- otu.sim[i.otu]
  count.lb <- 100; count.ub <- 300 ## lower the min/max range
  n.dist.lb <- 5; n.dist.ub <- 15 #no. of components for the mix. distn. min/max 
  rds.min <- 10000; rds.max <- 20000
  
  source(paste0(path, 'updated_code/cls.list.R'))
  
  cls2.list <- list(); cls2.list[[1]] <- cls2.list[[2]] <- cls2.list[[3]]  <- runif(n.otu.sim, 2, 6.5) #the beta parameter for the betas dist
  cls.names <- 1:length(cls2.list)
  
  gen.data.array <- array(0, dim=c(sum(unlist(ncount.list)), n.otu.sim))
  
  pdf(paste0(out.save.dir, "temp", sim.vars.sv, '.pdf'), height=12, width=12)
  par(mfrow=c(3, 3))
  
  rds.sample <- lapply(as.list(1:length(cls2.list)), function(x){as.integer(runif(ncount.list[[x]], rds.min, rds.max))})
  rds.sample.norm <- lapply(rds.sample, function(x){x/mean(x)})
  cls.sample <- rep(1:length(cls2.list), times=unlist(ncount.list))
  
  
  for (jj in 1:n.otu.sim) {
    t1 <- as.integer(runif(1, count.lb, count.ub)) #define the max of the data
    n.dist <- sample(n.dist.lb:n.dist.ub, 1, replace=F) #define no. of distributions
    alp.vals <- unique(c(0, as.integer(exp(seq(0, log(t1), length=n.dist))))) #assign the alpha values
    data.cls.temp <- c()
    iter.loc <- sample(1:length(cls2.list))
    for (cls in iter.loc){
      rand1 <- rbeta(ncount.list[[cls]], cls1.list[[cls]][jj], cls2.list[[cls]][jj])
      x3 <- hist(rand1, plot=F, breaks=seq(0, 1, length=length(alp.vals)+1))
      wts.prob <- x3$counts/sum(x3$counts) #assign the weights within the distributions
      
      ### Generate the random data ###
      param.vec <- cbind(alp.vals, c(0, rep(1, length(alp.vals)-1)))
      rand.rate1 <- gen.rand.rates(ncount.list[[cls]], param.vec, wts.prob, sim.type='gamma', fix.wts=FALSE)
      ind.rds1 <- rds.sample.norm[[cls]] #generate the random reads (standardized) 
      dat.out1 <- gen.rand.counts(rand.rate1, ind.rds1) #dat.out = obs. counts, ind.rds = total reads
      data.cls.temp <- c(data.cls.temp, dat.out1)
      
      #		generate plots of the data
      if (cls == 1) {
        plot(alp.vals, wts.prob, type='l', col='green', ylim=c(0, 0.5))
      } else if (cls == 2) {
        lines(alp.vals, wts.prob, col='red')
      } else {
        lines(alp.vals, wts.prob, col='blue')
      }
      print(wts.prob[1])
    }
    gen.data.array[,jj] <- data.cls.temp
  }
  
  gen.data.array.norm <- apply(gen.data.array, 2, function(x){x/unlist(rds.sample)}) #normalized array
  
  dev.off()
  
  sv.list <- c('rds.sample', 'rds.sample.norm', 'cls.names','cls.sample', 'gen.data.array', 'gen.data.array.norm')
  data.save.fn <- paste0(out.save.dir, 'simulated_data', sim.vars.sv, '.RData')
  save(list=sv.list, file=data.save.fn)
}




