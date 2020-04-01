

#Define the list of models to run, along with the associated parameters
define.model.list <- function(precalc.par.list, alp.cut.min, lc.alp.cut) {
	#Define the model list here
	model.list <- list()
	if (alp.cut.min == lc.alp.cut) {
		model.list[[1]] <- precalc.par.list
		model.list[[2]] <- precalc.par.list[2:dim(precalc.par.list)[1],]
		model.list[[3]] <- precalc.par.list[3:dim(precalc.par.list)[1],]
		model.list[[4]] <- precalc.par.list[4:dim(precalc.par.list)[1],]
		model.list[[5]] <- precalc.par.list[5:dim(precalc.par.list)[1],]
	} else {
		model.list[[1]] <- model.list[[2]] <- model.list[[3]] <- model.list[[4]] <- model.list[[5]] <- precalc.par.list
	}

	#Match the model formulation mixtures (list of indexes):
	all.mixt <- unique(unlist(lapply(model.list, function(x){paste(x[,1], x[,2], sep="_")})))
	model.mixt.loc <- lapply(model.list, function(x){match(paste(x[,1], x[,2], sep="_"), all.mixt)}) #location in the model/precalc matrix
	model.mixt.indx <- lapply(model.mixt.loc, function(x){c(1, (x+1), length(all.mixt)+2)}) #add the strz and truncation to get weight index

	n.models <- length(model.list) #total models
	max.mixts <- dim(precalc.par.list)[1] #total unique mixtures
	model.mix.all <- precalc.par.list #joint mixture models to be used with combined weights
	return(list(mdl=model.list, mix=all.mixt, loc=model.mixt.loc, indx=model.mixt.indx, nmdl=n.models, max=max.mixts, all=model.mix.all))
}

#Pair of functions to precalculate the conditional expectation matrix of probabilities
internal.loop.cond.precalc <- function(r.y.new, precalc.par.list, tmat.val, prob.int, pid) {
	cond.pmat.precalc.k <- cond.all.precalc(r.y.new, precalc.par.list, tmat.val, prob.int)
	return(list(mat=cond.pmat.precalc.k, pid=pid))
}

cond.pmat.precalc.fn <- function(r.y.new, precalc.par.list, tmat.vals, prob.int) {
	t.st <- Sys.time()
	r <- foreach(kk=1:length(tmat.vals)) %dopar% {
		internal.loop.cond.precalc(r.y.new, precalc.par.list, tmat.vals[kk], prob.int, Sys.getpid())
	}
	proc.list <- unlist(lapply(r, function(x){x$pid}))
#	pskill(proc.list, SIGTERM)

	cond.pmat.precalc <- array(0, dim=c(length(r.y.new), dim(precalc.par.list)[1], length(tmat.vals)))
	for (kk in 1:length(tmat.vals)) {
		cond.pmat.precalc[,,kk] <- r[[kk]]$mat
	}
	print(Sys.time() - t.st)
	return(cond.pmat.precalc)
}

#the loop function to calculate bootstrap weight estimates
internal.loop.boot <- function(kk, samp.cts.dat, samp.cts.rds, tmat.vals, cond.p.mat.mdl, count.cut, opt.inits, prob.int, pid) {
	mat.indxs <- match(round(samp.cts.rds, sig.level.round), tmat.vals) #match the reads + calculate the matrix
	cond.p.mat.temp <- apply(cond.p.mat.mdl[,,mat.indxs], c(1,2), mean)
	r.y.boot <- hist(samp.cts.dat, breaks=prob.int, plot=F)$counts #bin the counts
	init.precalc <- list(opt.inits[1:(length(opt.inits)-1)])
	lsq.est <- est.strz(samp.cts.dat, cond.p.mat.temp, samp.cts.rds, count.cut, r.y=r.y.boot, wts.init=init.precalc) #estimate the weights
	out.est.wts.iter <- lsq.est$est.compare[alg.indx,] #the weights estimate
	return(list(wts=out.est.wts.iter, pid=pid))
}

#the loop function to calculate distances for each of the observations and models
internal.loop.dist <- function(kk, samp.cts.dat, samp.cts.rds, out.est.wts, dat.out, ind.rds, wts.bt, bt.vld.indx, tmat.vals, cond.p.mat.mdl, tr.val, tot.mod.dist, cluster.tag, prob.int, pid) {
	boot.dist.calc.k <- calc.boot.dist(samp.cts.dat, samp.cts.rds, out.est.wts, dat.out, ind.rds, wts.bt, bt.vld.indx, tmat.vals, cond.p.mat.mdl, tr.val, tot.mod.dist, dis.inc=c(0,0,1,0), cluster.tag=cluster.tag, prob.int=prob.int)
	return(list(dist=boot.dist.calc.k$l3[bt.vld.indx,], pid=pid))
}


#helper  functions for initialization and computing averages
gen.start.init <- function(model.mixt.indx, model.sel, str.est.mixt.main) {
	wts.indx <- model.mixt.indx[[model.sel]]; init.temp <- str.est.mixt.main[wts.indx]
	init.temp <- init.temp/sum(init.temp); init.temp <- init.temp[1:(length(init.temp)-1)]
	init.temp[which(init.temp >= 1)] <- 1-1e-3; init.temp[which(init.temp < 0)] <- 0
	wts.init <- list(init.temp)
	return(wts.init)
}

gen.single.est <- function(cond.pmat.select, dat.out, ind.rds, count.cut, r.y.new, alg.indx=1) {
	cond.p.mat <- apply(cond.pmat.select, c(1, 2), mean)
	str.est.mixt.main <- est.strz(dat.out, cond.p.mat, ind.rds, count.cut, r.y=r.y.new)$est.compare[alg.indx,] #estimate the weights
	return(str.est.mixt.main)
}

calc.bt.avg <- function(out.wts.temp) {
	bt.vld.indx <- which(out.wts.temp[,1] != -1)
	wts.bt <- apply(out.wts.temp[bt.vld.indx,], 2, mean) #bootstrap average
	return(wts.bt)
}

#calculate the model proportions based on the distances
gen.avg.ests <- function(dist.boot.vec.store, n.dir.rep, n.models, max.mixts, model.mixt.indx, out.est.wts) {
	#1. Optimal model per distance
	mod.dist.indxs <- apply(dist.boot.vec.store, c(1, 2), function(x){max(which(x == min(x)))}) #iter. opt model (pdf, kl, cdf)
	opt.mod.sel.ineq <- apply(apply(mod.dist.indxs, 2, function(x){est.strz.ry(x, n.models-1, est.s=T)/n.dir.rep}), 2, function(x){x/sum(x)})

	boot.mix.store.all <- array(0, dim=c(n.dir.rep, max.mixts+2, n.models)) #matrix of mixed weights
	for (model.sel in 1:n.models) {
		boot.mix.store.all[,model.mixt.indx[[model.sel]],model.sel] <- out.est.wts[[model.sel]] #raw results
	}
	mdls.all.avg <- sapply(1:n.models, function(x){calc.bt.avg(boot.mix.store.all[,,x])}) #all of the model estimates averaged

	#3. All Model Estimates + All Model Estimates scaled by the optimal model proportion
	mix.boot.simple.mat <- sapply(1:tot.mod.dist, function(x){mdls.all.avg%*%opt.mod.sel.ineq[,x]}) #simple mix to get the model weights

	return(list(all=boot.mix.store.all, prop=opt.mod.sel.ineq, mix=mix.boot.simple.mat, avg=mdls.all.avg))
}

#Calculate the individual weights
ind.wts.mix <- function(dat.run.dist, max.mixts, n.models, mdls.all.avg, model.mix.all, tr.val, opt.mod.sel.ineq, dist.opt) {
	wts.all.samp.mix <- array(0, dim=c(dim(dat.run.dist)[1], max.mixts+2, n.models))
	for (mdl in 1:n.models) {
		mdl.wts.calc <- mdls.all.avg[,mdl]
		mdl.wts.calc[which(mdl.wts.calc < 0)] <- 0
		mdl.wts.calc <- mdl.wts.calc/sum(mdl.wts.calc)
		wts.all.samp.mix[,,mdl] <- calc.ind.wt.prob(dat.run.dist[,1], dat.run.dist[,2], model.mix.all, mdl.wts.calc, tr.val)
	}
	if (any(is.na(wts.all.samp.mix))) {
		stop('NaN in Individual Weight Estimates')
	}
	wts.all.samp <- apply(wts.all.samp.mix, c(1, 2), function(x){sum(x*opt.mod.sel.ineq[,dist.opt])}) #distance 3, simple scale
	return(list(mix=wts.all.samp.mix, wtd=wts.all.samp))
}

#Probabilities from the model mixture matrix (descrete computations, used for pairwise distances)
calc.tr.ct.pr <- function(model.mix.all, x.vals, dat.run.dist, tr.val) {
	if (tr.val < 600) {
		#If there are not a lot of integer categores, compute each one separately
		true.ct.pr <- apply(model.mix.all, 1, function(x){dnbinom(x.vals, size=x[1], prob=x[2]/(1+x[2]))})
		true.ct.pr[dim(true.ct.pr)[1],] <- true.ct.pr[dim(true.ct.pr)[1],] + 1-apply(true.ct.pr, 2, sum)
	} else {
		#If there are a lot of integer categories aggregate some of the counts for computation
		if (min(dat.run.dist[,1]) < 50) {
			int.blks <- c(seq(-1, 29), unique(as.integer(exp(seq(log(30), log(tr.val+1), length=500)))), Inf)
			true.ct.pr <- cond.all.precalc(int.blks, model.mix.all, 1, c(-2, int.blks))[,,1]
		} else {
			int.blks <- c(unique(as.integer(exp(seq(log(min(dat.run.dist[,1])), log(tr.val+1), length=500)))), Inf)
			true.ct.pr <- cond.all.precalc(int.blks, model.mix.all, 1, c(min(int.blks)-1, int.blks))[,,1]
		}
	}
	return(true.ct.pr)
}

#Make the internal GG' matrix: (used for L2 continuous CDF computation
gam.mix.fn <- function(x.val, par1, par2) {
	return( pgamma(x.val, shape=par1[1], rate=par1[2])*pgamma(x.val, shape=par2[1], rate=par2[2]) )
}

gam.sing.fn <- function(x.val, par1) {
	return( pgamma(x.val, shape=par1[1], rate=par1[2]) )
}

gen.g.mat <- function(wts.all.samp, model.mix.all, count.cut) {
	g.mat <- array(0, dim=c(dim(wts.all.samp)[2], dim(wts.all.samp)[2]))
	for (i.mix in 1:dim(model.mix.all)[1]) {
		for (j.mix in 1:dim(model.mix.all)[1]) {
			g.mat[i.mix+1,j.mix+1] <- integrate(gam.mix.fn, lower=0, upper=count.cut, par1=model.mix.all[i.mix,], par2=model.mix.all[j.mix,])$value
		}
	}
	g.mat[1,1] <- count.cut
	g.mat[1,2:(dim(g.mat)[2]-1)] <- g.mat[2:(dim(g.mat)[1]-1),1] <- apply(model.mix.all, 1, function(x){integrate(gam.sing.fn, lower=0, upper=count.cut, par1=x)$value})
	return(g.mat)
}

