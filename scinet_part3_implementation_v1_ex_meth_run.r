
############################
#Define the data and cluster

#dat.out <- data.train[,1]; ind.rds <- data.train[,2]
dat.clust <- cluster.ct.dat(data.train[,1], alp.max.q, alp.min.cut, ct.brk, lc.alp.cut) #bin the data to make a reasonable no. of bins
r.y.new <- dat.clust$ry; prob.int <- dat.clust$prob; tr.val <- dat.clust$trv
mat.indxs <- match(round(data.train[,2], sig.level.round), tmat.vals) #match the reads to the precalculated matrices

#$ry - ry counts; #$prob - bins; #$trv - truncation value at max.alp.q

#########################################
#generate the aggregate model for the fit
alp.cut.min <- max(lc.alp.cut, min(data.train[,1])) #cutoff for poisson sampling
count.cut <- alp.cut.min-2 #maximum number of counts to estimate with poisson
precalc.par.list <- dist.defn.cond.mat(alp.cut.min, lc.alp.cut, tr.val, subdiv.list, ovlp.lvl.list, alp.l, bet.l) #mxt. for fit

model.params <- define.model.list(precalc.par.list, alp.cut.min, lc.alp.cut)
model.list <- model.params$mdl
all.mixt <- model.params$mix
model.mixt.loc <- model.params$loc
model.mixt.indx <- model.params$indx
n.models <- model.params$nmdl
max.mixts <- model.params$max
model.mix.all <- model.params$all

##############################################################################
#Calculate the Conditional Probability Matrix for the defined groups (r.y.new)
cond.pmat.precalc <- cond.pmat.precalc.fn(r.y.new, precalc.par.list, tmat.vals, prob.int)

model.init.prm <- c('init.param.list', 'data.train', 'dat.clust', 'model.params', 'alp.cut.min', 'lc.alp.cut', 'count.cut', 'data.train', 'precalc.par.list')

###########################################
#Load the pre-generated bootstrap iterates:
if (!exists(f.dir)) dir.create(f.dir, showWarnings = FALSE)
f.name <- paste0('bootstrap_iterate_scenario_ex_iter', jj, '.RData', sep='')
boot.f.name <- paste(f.dir, f.name, sep='/')
if (file.exists(boot.f.name)) {
	load(boot.f.name)
	samp.cts <- boot.out.list.store$indx
	samp.cts.dat <- boot.out.list.store$dat
	samp.cts.rds <- boot.out.list.store$rds
} else {
	#Perform the bootstrap resampling:
	samp.cts <- sapply(1:n.dir.rep, function(x){sample(1:dim(data.train)[1], dim(data.train)[1], replace=T)}) #sample indeces
	samp.cts.dat <- sapply(1:n.dir.rep, function(x){data.train[samp.cts[,x],1]}) #resample data
	samp.cts.rds <- sapply(1:n.dir.rep, function(x){data.train[samp.cts[,x],2]}) #resample reads
	boot.out.list.store <- list(indx=samp.cts, dat=samp.cts.dat, rds=samp.cts.rds)
  save(list=c('boot.out.list.store'), file=boot.f.name)
}

###############################
#### BEGIN ESTIMATION HERE ####

#Generate precalculation estimates for each of the models
dist.boot.vec.store <- array(0, c(n.dir.rep, tot.mod.dist, n.models))
str.boot.storage <- str.main.storage <- out.est.wts <- list() #store the bootstrap weight estimates

#initialize the fit for the full model (model.sel <- 1)
str.est.mixt.main <- gen.single.est(cond.pmat.precalc[,model.mixt.loc[[1]], mat.indxs], data.train[,1], data.train[,2], count.cut, r.y.new)

for (model.sel in 1:n.models) {
	wts.init <- gen.start.init(model.mixt.indx, model.sel, str.est.mixt.main)
	str.z.list <- gen.single.est(cond.pmat.precalc[,model.mixt.loc[[model.sel]], mat.indxs], data.train[,1], data.train[,2], count.cut, r.y.new)
	if (str.z.list[1] == -1) {str.z.list <- rep(1, length(str.z.list))/length(str.z.list)}

	out.est.wts[[model.sel]] <- matrix(0, nrow=n.dir.rep, ncol=length(str.z.list)) #store the weights estimates
	cond.p.mat.mdl <- cond.pmat.precalc[,model.mixt.loc[[model.sel]],] #the model matrix

	print("Bootstrapping the Estimates")
	tt.st.1 <- Sys.time()
	r <- foreach(kk=1:n.dir.rep) %dopar% {
		internal.loop.boot(kk, samp.cts.dat[,kk], samp.cts.rds[,kk], tmat.vals, cond.p.mat.mdl, count.cut, str.z.list, prob.int, Sys.getpid())
	}
	print(Sys.time() - tt.st.1)
	proc.list <- unlist(lapply(r, function(x){x$pid}))
#	pskill(proc.list, SIGTERM)


	for (kk in 1:n.dir.rep) {
		out.est.wts[[model.sel]][kk,] <- r[[kk]]$wts
	}

	#Calculate the distance of data vs. bootstrap results
	str.boot.storage[[model.sel]] <- calc.bt.avg(out.est.wts[[model.sel]])
	str.main.storage[[model.sel]] <- str.z.list
	bt.vld.indx <- which(out.est.wts[[model.sel]][,1] != -1)

	#Calculate the Boostrapped Distances
	tt.st.1 <- Sys.time()
	r <- foreach(kk=1:length(bt.vld.indx)) %dopar% {
		internal.loop.dist(kk, samp.cts.dat, samp.cts.rds, out.est.wts[[model.sel]], data.train[,1], data.train[,2], wts.bt, bt.vld.indx[kk], tmat.vals, cond.p.mat.mdl, tr.val, tot.mod.dist, TRUE, prob.int, Sys.getpid())
	}
	print(Sys.time() - tt.st.1)
	proc.list <- unlist(lapply(r, function(x){x$pid}))
	
#	pskill(proc.list, SIGTERM)

	for (kk in 1:length(bt.vld.indx)) {
		dist.boot.vec.store[bt.vld.indx[kk],,model.sel] <- r[[kk]]$dist
	}
	print(Sys.time() - tt.st.1)
}

store.vars <- c('out.est.wts', 'dist.boot.vec.store', 'str.boot.storage', 'str.main.storage')


##################################
#Calculate the optimap proportions

est.gen.avg <- gen.avg.ests(dist.boot.vec.store, n.dir.rep, n.models, max.mixts, model.mixt.indx, out.est.wts)
boot.mix.store.all <- est.gen.avg$all
opt.mod.sel.ineq <- est.gen.avg$prop
mix.boot.simple.mat <- est.gen.avg$mix
mdls.all.avg <- est.gen.avg$avg

dist.sv.list <- c('boot.mix.store.all', 'opt.mod.sel.ineq', 'mix.boot.simple.mat', 'mdls.all.avg')

## dat.run.dist, max.mixts, n.models
## mdls.all.avg, model.mix.all, tr.val




##################################################################
#Calculate the weights for each individual given a model and data:

wts.all.samp.calc <- ind.wts.mix(dat.run.dist, model.params$max, model.params$nmdl, mdls.all.avg, model.params$all, tr.val, opt.mod.sel.ineq, dist.opt)
wts.all.samp.mix <- wts.all.samp.calc$mix #the individual mixture probabilitites
wts.all.sub <- wts.all.samp.calc$wtd #the weighted avg. of the mixture probabilities
#### wts.all.sub: THIS MATRIX WILL BE USED TO COMPUTE THE DISTANCE MATRIX ####

##########################################
### Compute the pairwise distance here ###

#note that this is not a function because of potential memory issues

### Generate a list of matching indeces to parallelization ###
matrix.rc.list <- array(0, dim=c(dim(dat.run.dist)[1]*(dim(dat.run.dist)[1]-1)/2, 2))
counter <- 1
for (ind1 in 1:(dim(dat.run.dist)[1]-1)) {
	for (ind2 in (ind1+1):dim(dat.run.dist)[1]) {
		matrix.rc.list[counter,] <- c(ind1, ind2); counter <- counter + 1
	}
}

#Array to store the distances
dist.mat.array <- array(0, dim=c(dim(dat.run.dist)[1], dim(dat.run.dist)[1], 3)) #l2.d.pdf/l2.d.cdf/kl/l2.c.cdf

gen.mix.prob <- function(true.ct.pr, wts.all.sub) {
	wts.all.cut <- wts.all.sub[,2:(dim(wts.all.sub)[2]-1)]
	if (is.null(dim(wts.all.cut))) {
		wts.all.cut <- matrix(wts.all.cut, nrow=1, ncol=length(wts.all.cut))
	}
	true.mix.prob <- true.ct.pr%*%t(wts.all.cut) #the mixture probabilities
	true.mix.prob[1,] <- true.mix.prob[1,] + wts.all.sub[,1]
	true.mix.prob[dim(true.mix.prob)[1],] <- true.mix.prob[dim(true.mix.prob)[1],] + wts.all.sub[,dim(wts.all.sub)[2]]
	true.mix.prob <- t(true.mix.prob)
	true.mix.prob.cdf <- t(apply(true.mix.prob, 1, cumsum))
	return(list(pdf=true.mix.prob, cdf=true.mix.prob.cdf))
}

#The associated weights and probabilites
true.ct.pr <- calc.tr.ct.pr(model.mix.all, seq(-1, tr.val+1), dat.run.dist, tr.val) #the mdl. mix. mat probabilities
true.mix <- gen.mix.prob(true.ct.pr, wts.all.sub)
true.mix.prob <- true.mix$pdf
true.mix.prob.cdf <- true.mix$cdf

# wts.all.sub - test.set.idx
# #cls.sample - original
# #train.set.idx - original
# cls.sample[test.set.idx] - test.set.idx
# test.set.idx %in% train.set.idx

#distance of each point to the class medoid
cls.idx <- 1
med.store.array <- array(0, dim=c(dim(wts.all.sub)[1], length(cls.names), 3))
for (cls.idx in 1:length(cls.names)) {
	cls.index <- which(cls.sample[train.set.idx] == cls.names[cls.idx])
	class.med.est <- apply(wts.all.sub[train.set.idx[cls.index],], 2, mean)
	mix.med <- gen.mix.prob(true.ct.pr, matrix(class.med.est, nrow=1, ncol=length(class.med.est)))

	comp.array <- array(0, dim=c(dim(wts.all.sub)[1], dim(true.mix.prob)[2], 2))
	comp.array[,,1] <- true.mix.prob
	comp.array[,,2] <- mix.med$pdf[rep(1, dim(wts.all.sub)[1]),]
	l2d.pdf.med <- apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum)

	comp.array <- array(0, dim=c(dim(wts.all.sub)[1], dim(true.mix.prob.cdf)[2], 2))
	comp.array[,,1] <- true.mix.prob.cdf
	comp.array[,,2] <- mix.med$cdf[rep(1, dim(wts.all.sub)[1]),]
	l2d.cdf.med <- apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum)

	g.mat <- gen.g.mat(wts.all.sub, model.mix.all, count.cut)
	comp.array <- array(0, dim=c(dim(wts.all.sub)[1], dim(wts.all.sub)[2], 2))
	comp.array[,,1] <- wts.all.sub
	comp.array[,,2] <- (matrix(class.med.est, 1, length(class.med.est)))[rep(1, dim(wts.all.sub)[1]),]
	l2c.cdf.med <- sapply(1:dim(comp.array)[1], function(x){wts1 <- comp.array[x,,1]; wts2 <- comp.array[x,,2]; return(t(wts1-wts2)%*%g.mat%*%(wts1-wts2))})
	med.store.array[,cls.idx,1] <- l2d.pdf.med
	med.store.array[,cls.idx,2] <- l2d.cdf.med
	med.store.array[,cls.idx,3] <- l2c.cdf.med
}

#####################
#L2-Descrete PDF Norm:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(dat.run.dist)[1], dim(dat.run.dist)[1]))
mx.k <- ceiling(exp(log(dim(true.mix.prob)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
	if (k.indx == mx.k) {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
	} else {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
	}
	comp.array <- array(0, dim=c(tot.idx, dim(true.mix.prob)[2], 2))
	comp.array[,,1] <- true.mix.prob[matrix.rc.list[idx.st:idx.en,1],]
	comp.array[,,2] <- true.mix.prob[matrix.rc.list[idx.st:idx.en,2],]
	dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum)
}
dist.mat.array[,,1] <- dist.calc.pw.mat
print(Sys.time() - t.t1)

######################
#L2-Descrete CDF Norm:
t.t1 <- Sys.time()
dist.calc.pw.mat <- array(0, dim=c(dim(dat.run.dist)[1], dim(dat.run.dist)[1]))
mx.k <- ceiling(exp(log(dim(true.mix.prob.cdf)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
	if (k.indx == mx.k) {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
	} else {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
	}
	comp.array <- array(0, dim=c(tot.idx, dim(true.mix.prob.cdf)[2], 2))
	comp.array[,,1] <- true.mix.prob.cdf[matrix.rc.list[idx.st:idx.en,1],]
	comp.array[,,2] <- true.mix.prob.cdf[matrix.rc.list[idx.st:idx.en,2],]
	dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- apply((comp.array[,,1] - comp.array[,,2])^2, 1, sum)
}
dist.mat.array[,,2] <- dist.calc.pw.mat
print(Sys.time() - t.t1)

#L2-Continous CDF Distance
t.t1 <- Sys.time()
g.mat <- gen.g.mat(wts.all.sub, model.mix.all, count.cut)
dist.calc.pw.mat <- array(0, dim=c(dim(dat.run.dist)[1], dim(dat.run.dist)[1]))
mx.k <- ceiling(exp(log(dim(true.mix.prob.cdf)[2]) + log(dim(matrix.rc.list)[1]) - log(7e+7))) #about 1Gb blocks
iter.blk <- ceiling(dim(matrix.rc.list)[1]/mx.k)
for (k.indx in 1:mx.k) {
	if (k.indx == mx.k) {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- dim(matrix.rc.list)[1]; tot.idx <- idx.en - idx.st + 1
	} else {
		idx.st <- (k.indx-1)*iter.blk+1; idx.en <- k.indx*iter.blk; tot.idx <- idx.en - idx.st + 1
	}
	comp.array <- array(0, dim=c(tot.idx, dim(wts.all.sub)[2], 2))
	comp.array[,,1] <- wts.all.sub[matrix.rc.list[idx.st:idx.en,1],]
	comp.array[,,2] <- wts.all.sub[matrix.rc.list[idx.st:idx.en,2],]
	dist.calc.pw.mat[matrix.rc.list[idx.st:idx.en,]] <- sapply(1:dim(comp.array)[1], function(x){wts1 <- comp.array[x,,1]; wts2 <- comp.array[x,,2]; return(t(wts1-wts2)%*%g.mat%*%(wts1-wts2))})
}
dist.mat.array[,,3] <- dist.calc.pw.mat
print(Sys.time() - t.t1)


storage.vars <- c('dist.mat.array', 'wts.all.samp.calc', 'true.ct.pr', 'true.mix.prob', 'true.mix.prob.cdf', 'med.store.array')

### Save the results of the iterations ###
# sim.vars.sv <- paste0("_i", i.otu, "_j", j.sample, "_ver", dat.ver.idx)
model.fit.init.fn <- paste0(f.dir, '/model_train_fit_otu_iteration', jj, '.RData')

var.save.list <- c(model.init.prm, store.vars, dist.sv.list, storage.vars)
save(list=var.save.list, file=model.fit.init.fn)









