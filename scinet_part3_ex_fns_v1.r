#Estimate Cumulative probability for a range:
est.cond.p.mat.rng <- function(alp, bet, ct.range, rds.avg) {
	#rds.avg = the scaled number of reads: mean(rds.avg) == 1
	pmat.arr <- array(0, dim=c(1, length(alp), length(rds.avg)))
	for (sid in 1:length(rds.avg)) {
		l1 <- sapply(1:length(alp), function(x){nb.params <- nb.ab(alp[x], bet[x], rds.avg[sid]); pnbinom(ct.range[1]-1, size=nb.params[1], prob=nb.params[2], lower=T)})
		l2 <- sapply(1:length(alp), function(x){nb.params <- nb.ab(alp[x], bet[x], rds.avg[sid]); pnbinom(ct.range[2], size=nb.params[1], prob=nb.params[2], lower=T)})
		pmat.arr[1:dim(pmat.arr)[1],(1:length(alp)),sid] <- l2 - l1
	}
	#the unscaled conditional probability matrix to use for estimation
	cond.mat.out <- apply(pmat.arr, c(1,2), mean) #average over samples to get the probability
	return(cond.mat.out)
}

#Estimate Tail Probability for condition matrix
est.cond.p.mat.inf <- function(alp, bet, ct.range, rds.avg) {
	#rds.avg = the scaled number of reads: mean(rds.avg) == 1
	pmat.arr <- array(0, dim=c(1, length(alp), length(rds.avg)))
	for (sid in 1:length(rds.avg)) {
		pmat.arr[1:dim(pmat.arr)[1],(1:length(alp)),sid] <- sapply(1:length(alp), 
			function(x){nb.params <- nb.ab(alp[x], bet[x], rds.avg[sid]); pnbinom(ct.range[1]-1, size=nb.params[1], prob=nb.params[2], lower=F)})
	}
	cond.mat.out <- apply(pmat.arr, c(1,2), mean) #average over samples to get the probability
	return(cond.mat.out)
}

#Minimum function for gammas:
min.fn <- function(x, par1, par2) {
	return( apply(cbind(dgamma(x, shape=par1[1], rate=par1[2]), dgamma(x, shape=par2[1], rate=par2[2])), 1, min) )
}

#Define the mixture of distributions to model the data, merging those with high overlap
dist.defn.cond.mat <- function(alp.cut.min, lc.alp.cut, tr.val, subdiv.list, ovlp.lvl.list, alp.l, bet.l) {
	alp.list <- alp.cut.min:tr.val; bet.list <- rep(1, length(alp.list))
	opt.high.cut <- min.ovlp <- rep(0, length(alp.list))
	prop.count <- ceiling(subdiv.list*length(alp.list))
	ovlp.cut <- rep(ovlp.lvl.list, times=prop.count) #control for different values of overlap
	ovlp.cut <- ovlp.cut[1:length(alp.list)]

	t.st <- Sys.time()
	indx1 <- 1
	counter <- 1; counter.max <- 160
	while (indx1 < length(alp.list) && counter < counter.max) {
		for (indx2 in (indx1+1):length(bet.list)) {
			p1 <- c(alp.list[indx1], bet.list[indx1])
			p2 <- c(alp.list[indx2], bet.list[indx2])
			llim <- min(qgamma(0.001, shape=p1[1], rate=p1[2]), qgamma(0.001, shape=p2[1], rate=p2[2]))
			ulim <- max(qgamma(0.995, shape=p1[1], rate=p1[2]), qgamma(0.995, shape=p2[1], rate=p2[2]))
			ovlp <- tryCatch(integrate(min.fn, par1=p1, par2=p2, lower=llim, upper=ulim)$value, warning=function(w){-1}, error=function(e){-1})
			if (ovlp == -1) {
				print('error'); next
			}
			if (ovlp < ovlp.cut[indx1] || indx2 == length(bet.list)) {
				opt.high.cut[indx1] <- indx2; indx1 <- indx2
				min.ovlp[indx1] <- ovlp
				counter <- counter + 1
				print(c(alp.list[indx1], ovlp)); break
			}
		}
	}
	print(Sys.time() - t.st)

	#Higher counts
	alp.par.h <- alp.list[c(1, opt.high.cut[which(opt.high.cut > 0)])] #list of alpha inits
	bet.par.h <- rep(1, length(alp.par.h)) #list of beta inits
	if (length(alp.par.h) > 150) {
		alp.par.h <- cluster.ct.dat(lc.alp.cut:tr.val, alp.max.q, alp.min.cut, 150, lc.alp.cut)$prob
		alp.par.h <- alp.par.h[which(alp.par.h >= alp.cut.min)]; alp.par.h[length(alp.par.h)] <- tr.val
		bet.par.h <- rep(1, length(alp.par.h)) #list of beta inits
	}

	#Lower Counts:
	if (alp.cut.min == lc.alp.cut) {
		alp.par.l <- alp.l
		bet.par.l <- bet.l
	} else {
		alp.par.l <- numeric(0)
		bet.par.l <- numeric(0)
	}
	precalc.par.list <- rbind(cbind(alp.par.l, bet.par.l), cbind(alp.par.h, bet.par.h))
	return(precalc.par.list)
}

#Cluster the data into groups for higher counts
cluster.ct.dat <- function(dat.uid, alp.max.q, alp.min.cut, ct.brk, lc.alp.cut, min.tail.obs=6) {
	uid.ex <- dat.uid
	tr.val <- as.integer(quantile(uid.ex, alp.max.q))
	if (tr.val <= lc.alp.cut) {
		tr.val <- uid.ex[order(uid.ex, decreasing=T)][min.tail.obs] #at least this many in the truncation for outliers
		tr.val <- max(tr.val, lc.alp.cut+1)
	}
	dat.all <- trunc.data(uid.ex, tr.val)
	dat.h.c <- dat.all[which(dat.all > lc.alp.cut)]
	dat.h <- log(dat.h.c)
	dat.l <- dat.all[which(dat.all <= lc.alp.cut)]

	if (length(dat.h) > 0) {
		x1 <- seq(min(dat.h), max(dat.h), by=(max(dat.h)-min(dat.h))/ct.brk)
		x1 <- sort(c(x1, log(tr.val))); x1 <- x1[which(x1 <= log(tr.val))]

		int.sum <- c(as.integer(exp(x1)), Inf) #the probability summation breaks
		int.sum <- c(min(int.sum), unique(int.sum))
		bin.ct <- hist(dat.h.c, breaks=int.sum, plot=F) #binning is right closed
		#print(sum(bin.ct$counts) == length(dat.h)) #check for match
		r.y.h <- bin.ct$counts; 
	} else {
		int.sum <- r.y.h <- numeric(0)
	}

	if (length(dat.l) > 0) {
		r.y.l <- est.strz.ry(dat.l, lc.alp.cut-1, est.s=FALSE) #generate count frequencies
		int.sum.l <- 0:lc.alp.cut
	} else {
		r.y.l <- int.sum.l <- numeric(0)
	}

	#Intevals for calculating the binned count probabilities (rows)
	r.y.new <- c(r.y.l, r.y.h)
	if (length(int.sum.l) > 0) {
		prob.int <- c(-1, int.sum.l, unique(int.sum)) #intervals for collapsing probabilities
	} else {
		prob.int <- c(min(int.sum)-1, unique(int.sum))
	}
	return(list(ry=r.y.new, prob=prob.int, trv=tr.val))
}

#Precalculate the conditional matrix for every block of the new r.y clustering and rounded values
cond.all.precalc <- function(r.y.new, precalc.par.list, tmat.vals, prob.int) {
	cond.pmat.precalc <- array(0, dim=c(length(r.y.new), dim(precalc.par.list)[1], length(tmat.vals)))
	for (i in 1:length(tmat.vals)) {
		for (indx in 2:length(prob.int)) {
			x1 <- c(prob.int[indx-1]+1, prob.int[indx])
			if (prob.int[indx] == Inf) {
				cond.pmat.precalc[(indx-1),,i] <- est.cond.p.mat.inf(precalc.par.list[,1], precalc.par.list[,2], x1, tmat.vals[i])
			} else {
				cond.pmat.precalc[(indx-1),,i] <- est.cond.p.mat.rng(precalc.par.list[,1], precalc.par.list[,2], x1, tmat.vals[i])
			}
		}
	}
	return(cond.pmat.precalc)
}	

