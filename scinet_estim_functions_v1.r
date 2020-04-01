#Pass in the parameters + weights to sample the rates from predefined distributions
#fix.wts = set the weight proportions exactly if true, sample from a multinomial if false
gen.rand.rates <- function(n.obs, param.vec, wts.vec, sim.type='gamma', fix.wts=TRUE) {
	if (sim.type == 'gamma') {
		dat.out <- list()
		mix.cts <- rmultinom(1, size=n.obs, prob=wts.vec) #sample from the multinomial
		for (i in 1:length(wts.vec)) {
			if (fix.wts) {
				n.block <- ceiling(n.obs*wts.vec[i])
			} else {
				n.block <- mix.cts[i,1] #randomized mixtures
			}
			if (param.vec[i,2] == 0) {
				dat.out[[i]] <- rep(0, n.block) #if zeros are included
			} else {
				dat.out[[i]] <- rgamma(n.block, shape=param.vec[i,1], rate=param.vec[i,2]) #generate gamma rates (no scaling by rds)
			}
		}
		exc.val <- length(unlist(dat.out))-n.obs
		if (exc.val > 0) {
			dat.out <- trim.list(dat.out, exc.val)
		}
	}
	if (length(unlist(dat.out)) != n.obs) {
		stop('mismatch in observations')
	}
  # return(list(wts.vec=wts.vec, param.vec=param.vec, sim.type=sim.type, dat=dat.out))
	return(list(wts.vec=wts.vec, param.vec=param.vec, sim.type=sim.type, dat=sample(unlist(dat.out))))
}

#Given the rates, generate the random counts (rates are given as qi*T_bar)
#norm.sc = normalize the reads to be mean=1, otherwise use input
gen.rand.counts <- function(rand.rate, ind.rds, norm.sc=TRUE) {
	if (norm.sc) {
		rds.avg <- ind.rds/mean(ind.rds) #avg. reads scale
	} else {
		rds.avg <- ind.rds #do not normalize
	}
	rand.all <- unlist(rand.rate$dat) #extract the rates
	if (length(rds.avg == 1)) {
		rds.avg <- rep(rds.avg, length(rand.all)) #if reads are all the same, replicate
	}
	counts.out <- rep(0, length(rand.all))
	for (i in 1:length(rand.all)) {
		counts.out[i] <- rpois(1, rand.all[i]*rds.avg[i]) #generate counts
	}
	return(counts.out)
}

#Generate a set of random reads and standardize to have mu = 1
gen.rand.reads <- function(n.obs, min.rds, max.rds) {
	output <- as.integer(runif(n.obs, min.rds, max.rds))
	return( output/mean(output) )
}

#Estimate Structural Zeros - Main Function
#Given simulated data (dat.out), a cutoff for high counts (count cut) a probability matrix p(x|dist)
#optionally the true exact counts (r.ex) and the true exact rate (rand.rate)
#Estimate 1) the zeros 2) S = total # of species

est.strz <- function(dat.out, cond.p.mat, ind.rds, count.cut, rand.rate=NULL, wts.init=NULL, r.y=NULL) {
	dat.out.init <- dat.out #store the original data
	dat.out <- trunc.data(dat.out, count.cut)

	r.m <- est.strz.pmat(cond.p.mat, length(dat.out), strz=TRUE, trunc=TRUE) #modify the probability matrix
	if (is.null(r.y)) {
		r.y <- est.strz.ry(dat.out, count.cut, est.s=FALSE) #generate count frequencies (if not input)
	}

	#Estimate using the exact generated counts
	if (!is.null(rand.rate)) {
		r.ex <- calc.ex.counts(rand.rate, ind.rds, count.cut)
		r.y <- r.ex #exact counts less expected zeros
	}

	#dims: (length(r.y) == count.cut + zero + truncation, dim(r.m) <- count.cut + 0 + truncation x count.cut + zero + str.zero + trunc
	############################################

	#Generate a list of starting values
	if (is.null(wts.init)) {
		start.list <- est.strz.gen.startlist(r.y=r.y, n.wts=(dim(r.m)[2]-1)) 	#nlopt.ineq inits
	} else {
		start.list <- append(wts.init, est.strz.gen.startlist(r.y=r.y, n.wts=(dim(r.m)[2]-1)))
	}
	start.list.gen <- start.list 											#nlopt.oth inits
	start.list.mod <- lapply(start.list, function(x){c(x, 1-sum(x))}) 		#nlopt.eq inits

	#Check the derivatives
	#res <- check.derivatives(.x=x0, func=eval_f0, func_grad=eval_grad_f0, 
	#	check_derivatives_print='all', r.m = r.m, r.y = r.y, var.wt = var.wt)

	#List of three standalone algorithms for est.nlopt.oth
	ld.gen.list <- c("NLOPT_LN_COBYLA", "NLOPT_LD_MMA", "NLOPT_LD_SLSQP")
	#ld.gen.list.sel <- c("NLOPT_LN_COBYLA", "NLOPT_LD_MMA", "NLOPT_LD_SLSQP")
	ld.gen.list.sel <- c()

	#List of algorithms for est.nlopt.ineq (inequality constraints)
	ld.loc.list <- c("NLOPT_LD_LBFGS", "NLOPT_LD_VAR1", "NLOPT_LD_VAR2", "NLOPT_LD_TNEWTON", "NLOPT_LD_TNEWTON_RESTART", "NLOPT_LD_TNEWTON_PRECOND", "NLOPT_LD_TNEWTON_PRECOND_RESTART") #list of algorithms for inequality constraint
	#ld.loc.list.sel <- c('NLOPT_LD_LBFGS', 'NLOPT_LD_TNEWTON', 'NLOPT_LD_TNEWTON_RESTART', 'NLOPT_LD_TNEWTON_PRECOND', 'NLOPT_LD_TNEWTON_PRECOND_RESTART') #list of algorithms for equality constraint
	#ld.loc.list.sel <- c('NLOPT_LD_LBFGS', 'NLOPT_LD_TNEWTON_RESTART', 'NLOPT_LD_TNEWTON_PRECOND_RESTART')
	ld.loc.list.sel <- c('NLOPT_LD_LBFGS')

	#List of algorithms for the est.nlopt.eq (equality constraints)
	#ld.loc.list.eq <- c('NLOPT_LD_LBFGS', "NLOPT_LD_TNEWTON_RESTART", "NLOPT_LD_TNEWTON_PRECOND_RESTART") #list of selected algirthms
	#ld.loc.list.eq.sel <- c('NLOPT_LD_LBFGS', "NLOPT_LD_TNEWTON", "NLOPT_LD_TNEWTON_RESTART") #list of selected algirthms
	ld.loc.list.eq <- ld.loc.list
	#ld.loc.list.eq.sel <- ld.loc.list.sel
	ld.loc.list.eq.sel <- c()

	#Load the necessary functions

	#Initialize the algorithm values
	###################################################
	var.wt <- rep(1, length(r.y)) #define the variances
	###################################################

	func.list <- func.ineq()
	eval_f0=func.list[['eval_f0']]; eval_grad_f0=func.list[['eval_grad_f0']]
	eval_g0=func.list[['eval_g0']]; eval_jac_g0=func.list[['eval_jac_g0']]

	#Standalone Inequality Constraints Optimization (COBYLA, MMA + SLSQP):
	lbs <- rep(0, length(start.list.gen[[1]])); ubs <- rep(1, length(start.list.gen[[1]]))
	t1.oth <- est.nlopt.oth(r.y, r.m, var.wt, lbs, ubs, eval_f0, eval_grad_f0, eval_g0, eval_jac_g0, ld.gen.list, ld.gen.list.sel, start.list.gen)
	opt1.res <- cbind(t1.oth$est.compare, 1-apply(t1.oth$est.compare, 1, sum))
	t1.oth$est.compare <- opt1.res

	#Inequality Constraints OPTIMZATION:
	lbs <- rep(0, length(start.list[[1]])); ubs <- rep(1, length(start.list[[1]]))
	t1.ineq <- est.nlopt.ineq(r.y, r.m, var.wt, lbs, ubs, eval_f0, eval_grad_f0, eval_g0, eval_jac_g0, ld.loc.list, ld.loc.list.sel, start.list)
	opt2.res <- cbind(t1.ineq$est.compare, 1-apply(t1.ineq$est.compare, 1, sum))
	t1.ineq$est.compare <- opt2.res

	########################################
	#Linear Constraint + boundary conditions
	########################################

	#Load the necessary functions
	func.list <- func.eq()
	eval_f <- func.list[['eval_f']]; eval_g_eq <- func.list[['eval_g_eq']]

	#Equality Constraints:
	#NLOPT_LD_LBFGS*, NLOPT_LD_TNEWTON*, NLOPT_LD_TNEWTON_RESTART, NLOPT_LD_TNEWTON_PRECOND, NLOPT_LD_TNEWTON_PRECOND_RESTART
	lbs <- rep(0, length(start.list.mod[[1]])); ubs <- rep(1, length(start.list.mod[[1]]))
	t1.eq <- est.nlopt.eq(r.y, r.m, var.wt, lbs, ubs, eval_f, eval_g_eq, ld.loc.list.eq, ld.loc.list.eq.sel, start.list.mod)
	opt3.res <- t1.eq$est.compare

	return(t1.ineq)
	#return(opt2.res)
}

######################################################################
#FUNCTIONS FOR MINIMIZATION ON EQUALITY CONTRAINT ESTIMATING n WEIGHTS
######################################################################
#Returns a list of functions for equality based minimixation
func.eq <- function() {
	#Objective Function:
	eval_f <- function(x, r.m, r.y, var.wt) {
		objc <- sum( ((r.m%*%x - r.y)^2)/var.wt )
		gradc <- c(apply(matrix( rep(2*(r.m%*%x - r.y)/var.wt, dim(r.m)[2]), nrow=dim(r.m)[1], ncol=dim(r.m)[2], byrow=F )*r.m, 2, sum))
		return( list("objective"=objc, "gradient"=gradc) )
	}

	#Equality Constraint:
	eval_g_eq <- function(x, r.m, r.y, var.wt) {
		constr <- c(sum(x) - 1)
		grad <- rep(1, length(x))
		return( list( "constraints"=constr, "jacobian"=grad ) )
	}
	return(list(eval_f=eval_f,eval_g_eq=eval_g_eq))
}

#######################################################################
#FUNCTIONS FOR MINIMIZATION BASED ON RESTRICTED WEIGHTS FOR n-1 WEIGHTS
#######################################################################
#Returns a list of functions for inequality based minimization
func.ineq <- function() {
	#Define the objective function:
	#x = variable(s) of interest, r.m = probability matrix
	#r.y = freq. matrix, var.wt = variance scaling
	eval_f0 <- function(x, r.m, r.y, var.wt) {
		x.ext <- c(x, (1-sum(x)))
		return( sum( ((r.m%*%x.ext - r.y)^2)/var.wt ) )
	}

	#Define the gradient of the objective function
	#x = variable(s) of interest, r.m = probability matrix
	#r.y = freq. matrix, var.wt = variance scaling
	eval_grad_f0 <- function(x, r.m, r.y, var.wt) {
		x.ext <- c(x, (1-sum(x)))
		m1 <- r.m[,1:(dim(r.m)[2]-1)] - matrix(rep(r.m[,dim(r.m)[2]], dim(r.m)[2]-1), nrow=dim(r.m)[1], ncol=dim(r.m)[2]-1, byrow=F)
		m2 <- matrix(rep(2*(r.m%*%x.ext - r.y)/var.wt, dim(r.m)[2]-1), nrow=dim(r.m)[1], ncol=dim(r.m)[2]-1, byrow=F)
		return( c(apply(m1*m2, 2, sum)) )
	}

	#Constraint function: 1 - sum(x) >= 0; sum(x) - 1 <= 0
	eval_g0 <- function(x, r.m, r.y, var.wt) {
		return( c(sum(x) - 1) )
	}

	#Jacobian of the constraint function d/dwj of g0
	eval_jac_g0 <- function(x, r.m, r.y, var.wt) {
		return( c(rep(1, length(x))) )
	}
	return(list(eval_f0=eval_f0, eval_grad_f0=eval_grad_f0, eval_g0=eval_g0, eval_jac_g0=eval_jac_g0))
}

#Aggregate data above the cutoff
trunc.data <- function(dat.out, count.cut) {
	high.counts <- which(dat.out > count.cut)
	dat.out[high.counts] <- count.cut + 1 #place holder for the truncation
	return(dat.out)
}


#Given count data and a cutoff, generate the observed data frequencies
#If estimating S; no zeros should be included
est.strz.ry <- function(dat.out, count.cut, est.s=FALSE) {
	t1 <- data.frame(x=seq(0, count.cut+1)); t2 <- as.data.frame(table(dat.out)) #only keep counts below the cutoff
	t3 <- merge(t1, t2, by.x='x', by.y='dat.out', all.x=T); t3 <- t3[order(t3$x),]; t3[which(is.na(t3[,2])),2] <- 0
	r.y <- t3$Freq #create the output dataframe
	if (est.s) {
		return(r.y[2:length(r.y)])
	} else {
		return(r.y)
	}
}

#Given a matrix of probabilities, add the zero and truncation distributions
#cond.p.mat = matrix of probabilities, dat.len = number of observations (normalize matrix)
#strz = add the structural zero distribution, trunc = add the truncation distribution
est.strz.pmat <- function(cond.p.mat, dat.len, strz=TRUE, trunc=TRUE) {
	if (strz && trunc) {
		#add the structural zero distribution + truncation
		r.m <- cbind(0, cond.p.mat, 0); r.m[1,1] <- 1; r.m[dim(r.m)[1], dim(r.m)[2]] <- 1
	} else if (strz && !trunc) {
		#add the structural zero distribution only
		r.m <- cbind(0, cond.p.mat); r.m[1,1] <- 1
	} else if (!strz && trunc) {
		#add the truncation distribution
		r.m <- cbind(cond.p.mat, 0); r.m[dim(r.m)[1], dim(r.m)[2]] <- 1
	} else {
		r.m <- cond.p.mat
	}
	r.m <- r.m * dat.len #standardize the matrix so weights sum to one
	return(r.m)
}


#Generate a list of starting values for the optimization algorithm
est.strz.gen.startlist <- function(r.y, n.wts) {
	start.list <- list()
	w.t <- r.y[1:(length(r.y)-1)]/sum(r.y); #skip the last weight, based on obs. counts
	if (length(w.t) < (n.wts-1)) {
		w.t <- c(w.t, rep(0, n.wts-length(w.t)-1))
	} else if (length(w.t) > (n.wts-1)) {
		w.t <- w.t[1:(n.wts-1)]
	}
	
	start.list[[1]] <- c(w.t[1]/2, w.t[1]/2, w.t[2:length(w.t)]); #50/50
	start.list[[2]] <- c(0, w.t[1], w.t[2:length(w.t)]) #all to non-zero
	start.list[[3]] <- c(w.t[1], 0, w.t[2:length(w.t)]); #all to zero
	start.list[[4]] <- rep(1, n.wts)/(n.wts+1); #equal weights
	start.list[[5]] <- c(0, w.t[1], sum(w.t[2:length(w.t)]), rep(0, length(w.t)-2)) #aggregated to the 3rd distn
	start.list[[6]] <- c(w.t[1]/2, w.t[1]/2, sum(w.t[2:length(w.t)]), rep(0, length(w.t)-2)) #50/50
	start.list[[7]] <- c(w.t[1], 0, sum(w.t[2:length(w.t)]), rep(0, length(w.t)-2)) #all to zero
	start.list[[8]] <- c(1, rep(0, n.wts-1)); #all weight on 1
	return(start.list)
}


#Given reps and probs return the mean(1) and variance(2) of the NB
nb.mv <- function(r.tr, p.pr) {
	est.mu <- (1-p.pr)*r.tr/p.pr
	return( c((1-p.pr)*r.tr/p.pr, est.mu/p.pr) )
}

#Given the mean and var return the r(1) and p(2) of the NB
nb.rp <- function(mu.est, sig.est) {
	return( c((mu.est^2)/(sig.est - mu.est), mu.est/sig.est) )
}

#Given the mean and var of the NB return the shape(1) and rate(2) of the gamma
gam.ab <- function(mu.est, sig.est) {
	return( c((mu.est^2)/(sig.est - mu.est), mu.est/(sig.est - mu.est)) )
}

#Given the alpha and beta return the r(1) and p(2)
nb.ab <- function(alp, bet, rds=1) {
	return( c(alp, bet/(rds+bet)) )
}

#Given the r and p of the NB, return the alpha and beta of the gamma
gam.rp <- function(r.tr, p.pr, rds=1) {
	return( c(r.tr, p.pr*rds/(1-p.pr)) )
}








#Estimate the optimum based on inequality constraints of cobyla, mma and slsqp algorithms
est.nlopt.oth <- function(r.y, r.m, var.wt, lbs, ubs, eval_f0, eval_grad_f0, eval_g0, eval_jac_g0, ld.gen.list, ld.gen.list.sel, start.list.gen) {
	x0 <- start.list.gen[[1]]; counter <- 1
	est.compare <- matrix(0, nrow=length(ld.gen.list), ncol=length(x0)) #init the storage array
	est.obj <- rep(0, length(ld.gen.list))
	rownames(est.compare) <- names(est.obj) <- ld.gen.list

	for (loc_algeo in ld.gen.list) {
		if (loc_algeo == "NLOPT_LN_COBYLA" && loc_algeo %in% ld.gen.list.sel) {
			for (kk in 1:length(start.list.gen)) {
				x0 <- start.list.gen[[kk]]
				res4 <- nloptr( x0=x0, eval_f=eval_f0, lb = lbs, ub = ubs, eval_g_ineq = eval_g0,
					opts = list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-5, "maxeval"=10000),
					r.m = r.m, r.y = r.y, var.wt = var.wt ); #print(res4)
				if (res4$status != 1 && res4$status != 3 && res4$status != 4) {
					print(paste(loc_algeo, 'did not converge for init', kk)) #no convergence for starting point
					if (kk == length(start.list.gen)) {
						print(paste(loc_algeo, 'did not converge for ALL INITS'))
						#est.compare[counter,] <- res4$solution; #stop('algeo_stop')
						est.compare[counter,] <- c(-1, rep(0, length(x0)-1))
						est.obj[counter] <- res4$objective;
					} else {
						next
					}
				} else {
					est.compare[counter,] <- res4$solution; #print(res4$objective)
					est.obj[counter] <- res4$objective;
					break
				}
			}
		}

		if (loc_algeo == "NLOPT_LD_MMA" && loc_algeo %in% ld.gen.list.sel) {
			for (kk in 1:length(start.list.gen)) {
				x0 <- start.list.gen[[kk]]
				res4 <- nloptr( x0=x0, eval_f=eval_f0, eval_grad_f=eval_grad_f0, lb = lbs, ub = ubs, 
					eval_g_ineq=eval_g0, eval_jac_g_ineq=eval_jac_g0,
					opts = list("algorithm"="NLOPT_LD_MMA", "xtol_rel"=1.0e-5, "maxeval"=10000),
					r.m = r.m, r.y = r.y, var.wt = var.wt ); #print(res4)
				if (res4$status != 1 && res4$status != 3 && res4$status != 4) {
					print(paste(loc_algeo, 'did not converge for init', kk)) #no convergence for starting point
					if (kk == length(start.list.gen)) {
						print(paste(loc_algeo, 'did not converge for ALL INITS'))
						#est.compare[counter,] <- res4$solution; #stop('algeo_stop')
						est.compare[counter,] <- c(-1, rep(0, length(x0)-1))
						est.obj[counter] <- res4$objective;
					} else {
						next
					}
				} else {
					est.compare[counter,] <- res4$solution; #print(res4$objective)
					est.obj[counter] <- res4$objective;
					break
				}
			}
		}
		if (loc_algeo == "NLOPT_LD_SLSQP" && loc_algeo %in% ld.gen.list.sel) {
			for (kk in 1:length(start.list.gen)) {
				x0 <- start.list.gen[[kk]]
				res4 <- nloptr( x0=x0, eval_f=eval_f0, eval_grad_f=eval_grad_f0, lb = lbs, ub = ubs, 
					eval_g_ineq=eval_g0, eval_jac_g_ineq=eval_jac_g0,
					opts = list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-7, "maxeval"=20000),
					r.m = r.m, r.y = r.y, var.wt = var.wt ); #print(res4)
				if (res4$status != 1 && res4$status != 3 && res4$status != 4) {
					print(paste(loc_algeo, 'did not converge for init', kk)) #no convergence for starting point
					if (kk == length(start.list.gen)) {
						print(paste(loc_algeo, 'did not converge for ALL INITS'))
						#est.compare[counter,] <- res4$solution; #stop('algeo_stop')
						est.compare[counter,] <- c(-1, rep(0, length(x0)-1))
						est.obj[counter] <- res4$objective;
					} else {
						next
					}
				} else {
					est.compare[counter,] <- res4$solution; #print(res4$objective)
					est.obj[counter] <- res4$objective;
					break
				}
			}
		}
		counter <- counter + 1
	}
	return(list(est.compare=est.compare, est.obj=est.obj))
}

#Given the paramters, gradient functions and algorithms, estimate the optimum
#for a variety of starting points using 1-sum(x) and an INEQUALITY constraint
est.nlopt.ineq <- function(r.y, r.m, var.wt, lbs, ubs, eval_f0, eval_grad_f0, eval_g0, eval_jac_g0, ld.loc.list, ld.loc.list.sel, start.list) {
	x0 <- start.list[[1]]; counter <- 1
	est.compare <- matrix(0, nrow=length(ld.loc.list), ncol=length(x0)) #init the storage array
	est.obj <- rep(0, length(ld.loc.list))
	rownames(est.compare) <- names(est.obj) <- ld.loc.list

	for (loc_algeo in ld.loc.list) {
		#local_opts <- list( "algorithm" = loc_algeo, "xtol_rel"  = 1.0e-7, "maxeval" = 5000 )
		#opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-7, "maxeval" = 5000, "local_opts"  = local_opts, "print_level" = 0 )
		local_opts <- list( "algorithm" = loc_algeo, "xtol_rel"  = 1.0e-6, "maxeval" = 6000 )
		opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-6, "maxeval" = 6000, "local_opts"  = local_opts, "print_level" = 0 )
		# Do optimization.
		if (loc_algeo %in% ld.loc.list.sel) {
			for (kk in 1:length(start.list)) {
				x0 <- start.list[[kk]]
				res4 <- nloptr( x0=x0, eval_f=eval_f0, eval_grad_f=eval_grad_f0, lb = lbs, ub = ubs,
					eval_g_ineq=eval_g0, eval_jac_g_ineq=eval_jac_g0,
					opts = opts, 
					r.m = r.m, r.y = r.y, var.wt = var.wt); #print(res4)
				if (res4$status != 1 && res4$status != 3 && res4$status != 4 && sum(round(res4$solution, 3) != round(x0, 3)) != 0) {
					print(paste(loc_algeo, 'did not converge for init', kk)) #no convergence for starting point
					if (kk == length(start.list)) {
						print(paste(loc_algeo, 'did not converge for ALL INITS'))
						#est.compare[counter,] <- res4$solution; #stop('algeo_stop')
						est.compare[counter,] <- c(-1, rep(0, length(x0)-1))
						est.obj[counter] <- res4$objective;
					} else {
						next
					}
				} else {
					est.compare[counter,] <- res4$solution; #print(res4$objective)
					est.obj[counter] <- res4$objective;
					break
				}
			}
		}
		counter <- counter + 1
	}
	return(list(est.compare=est.compare, est.obj=est.obj))
}


#Given the paramters, gradient functions and algorithms, estimate the optimum for a variety of starting points using EQUALITY constraint
est.nlopt.eq <- function(r.y, r.m, var.wt, lbs, ubs, eval_f, eval_g_eq, ld.loc.list.eq, ld.loc.list.eq.sel, start.list.mod) {
	x0 <- start.list.mod[[1]]; counter <- 1
	est.compare <- matrix(0, nrow=length(ld.loc.list.eq), ncol=length(x0)) #init the storage array
	est.obj <- rep(0, length(ld.loc.list.eq))
	rownames(est.compare) <- names(est.obj) <- ld.loc.list.eq

	for (loc_algeo in ld.loc.list.eq) {
		local_opts <- list( "algorithm" = loc_algeo, "xtol_rel"  = 1.0e-7 )
		opts <- list( "algorithm" = "NLOPT_LD_AUGLAG_EQ", "xtol_rel" = 1.0e-7, "maxeval"=5000, "local_opts"=local_opts, "print_level"=0 )
		if (loc_algeo %in% ld.loc.list.eq.sel) {
			for (kk in 1:length(start.list.mod)) {
				x0 <- start.list.mod[[kk]] #; x0 <- c(x0, 1-sum(x0))
				res4 <- nloptr( x0=x0, eval_f=eval_f, lb = lbs, ub = ubs, eval_g_eq=eval_g_eq,
					opts = opts, 
					r.m = r.m, r.y = r.y, var.wt = var.wt ); #print(res4)
				if (res4$status != 1 && res4$status != 3 && res4$status != 4) {
					print(paste(loc_algeo, 'did not converge for init', kk)) #no convergence for starting point
					if (kk == length(start.list.mod)) {
						print(paste(loc_algeo, 'did not converge for ALL INITS'))
						#est.compare[counter,] <- res4$solution; #stop('algeo_stop')
						est.compare[counter,] <- c(-1, rep(0, length(x0)-1))
						est.obj[counter] <- res4$objective;
					} else {
						next
					}
				} else {
					est.compare[counter,] <- res4$solution; #print(res4$objective)
					est.obj[counter] <- res4$objective;
					break
				}
			}
		}
		counter <- counter + 1
	}
	return(list(est.compare=est.compare, est.obj=est.obj))
}







###################################################################################################
#Calculate the distances between the bootstrap iterates/data and weights for different combinations
#boot data/reads/estimate obs data/reads/boot avg. wts/converged boot iterates (match values + conditional p. matrix), cutoff, # of dist to comp
calc.boot.dist <- function(samp.cts.dat, samp.cts.rds, out.est.wts, dat.out, ind.rds, wts.bt, bt.vld.indx, tmat.vals, cond.pmat.precalc, count.cut, tot.mod.dist, dis.inc=c(1,1,1,1), cluster.tag=FALSE, prob.int=NULL) {
	#Compute the distances:
	#######################
	#1. Avg. Model vs Obs. Data
	if (dis.inc[1] == 1) {
		ineq.dist1 <- calc.dist.cut(dat.out, ind.rds, count.cut, wts.bt, tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
	} else {
		ineq.dist1 <- rep(Inf, 3)
	}

	#2. Avg. Model vs Bootstrap Iterate
	tot.boot.rep <- dim(samp.cts.dat)[2] #columns are the replicates
	ineq.dist2 <- ineq.dist3 <- ineq.dist4 <- array(0, dim=c(tot.boot.rep, tot.mod.dist))

	if (dis.inc[2] == 1) {
		for (bt.indx in bt.vld.indx) {
			ineq.dist2[bt.indx,] <- calc.dist.cut(samp.cts.dat[,bt.indx], samp.cts.rds[,bt.indx], count.cut, wts.bt, tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
		}
		ineq.dist2[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist2 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}

	if (dis.inc[3] == 1) {
		#3. Iterate wts on Obs. Data
		for (bt.indx in bt.vld.indx) {
			ineq.dist3[bt.indx,] <- calc.dist.cut(dat.out, ind.rds, count.cut, out.est.wts[bt.indx,], tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
		}
		ineq.dist3[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist3 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}
	if (dis.inc[4] == 1) {
		#4. Iterate wts on iterate data
		for (bt.indx in bt.vld.indx) {
			ineq.dist4[bt.indx,] <- calc.dist.cut(samp.cts.dat[,bt.indx], samp.cts.rds[,bt.indx], count.cut, out.est.wts[bt.indx,], tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
		}
		ineq.dist4[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist4 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}
	return(list(l1=ineq.dist1, l2=ineq.dist2, l3=ineq.dist3, l4=ineq.dist4)) 
}

#Parallel version of calc boot.dist
calc.boot.dist.prl <- function(samp.cts.dat, samp.cts.rds, out.est.wts, dat.out, ind.rds, wts.bt, bt.vld.indx, tmat.vals, cond.pmat.precalc, count.cut, tot.mod.dist, dis.inc=c(1,1,1,1), cluster.tag=FALSE, prob.int=NULL) {
	#Compute the distances:
	#######################
	#1. Avg. Model vs Obs. Data
	if (dis.inc[1] == 1) {
		ineq.dist1 <- calc.dist.cut(dat.out, ind.rds, count.cut, wts.bt, tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
	} else {
		ineq.dist1 <- rep(Inf, 3)
	}

	#2. Avg. Model vs Bootstrap Iterate
	tot.boot.rep <- dim(samp.cts.dat)[2] #columns are the replicates
	ineq.dist2 <- ineq.dist3 <- ineq.dist4 <- array(0, dim=c(tot.boot.rep, tot.mod.dist))

	if (dis.inc[2] == 1) {
		for (bt.indx in bt.vld.indx) {
			ineq.dist2[bt.indx,] <- calc.dist.cut(samp.cts.dat[,bt.indx], samp.cts.rds[,bt.indx], count.cut, wts.bt, tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
		}
		ineq.dist2[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist2 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}

	if (dis.inc[3] == 1) {
		#3. Iterate wts on Obs. Data
		internal.loop.dist <- function(dat.out, ind.rds, count.cut, out.est.wts, bt.indx, tmat.vals, cond.pmat.precalc, cluster.tag, prob.int, pid) {
			dist.precalc.k <- calc.dist.cut(dat.out, ind.rds, count.cut, out.est.wts[bt.indx,], tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
			return(list(mat=dist.precalc.k, pid=pid))
		}
		r.dist <- foreach(bt.indx=1:length(bt.vld.indx)) %dopar% {
			internal.loop.dist(dat.out, ind.rds, count.cut, out.est.wts, bt.vld.indx[bt.indx], tmat.vals, cond.pmat.precalc, cluster.tag, prob.int, Sys.getpid())
		}
		proc.list <- unlist(lapply(r.dist, function(x){x$pid}))
#		pskill(proc.list, SIGTERM)
		
		for (bt.indx in 1:length(bt.vld.indx)) {
			ineq.dist3[bt.vld.indx[bt.indx],] <- r.dist[[bt.indx]]$mat
		}
		ineq.dist3[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist3 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}
	if (dis.inc[4] == 1) {
		#4. Iterate wts on iterate data
		for (bt.indx in bt.vld.indx) {
			ineq.dist4[bt.indx,] <- calc.dist.cut(samp.cts.dat[,bt.indx], samp.cts.rds[,bt.indx], count.cut, out.est.wts[bt.indx,], tmat.vals, cond.pmat.precalc, cluster.tag, prob.int)
		}
		ineq.dist4[which(!(1:tot.boot.rep %in% bt.vld.indx)),] <- Inf
	} else {
		ineq.dist4 <- array(Inf, dim=c(tot.boot.rep, tot.mod.dist))
	}
	return(list(l1=ineq.dist1, l2=ineq.dist2, l3=ineq.dist3, l4=ineq.dist4)) 
}



#Calculate the distance between observed truncated counts and the estimated conditional matrix for sample
calc.dist.cut <- function(dat.obs, ind.rds, count.cut, wts.bt, tmat.vals, cond.pmat.precalc, cluster.tag=FALSE, prob.int=NULL) {
	#Define the observed data PDF and CDF:
	if (!cluster.tag) {
		y.dat.vls <- est.strz.ry(trunc.data(dat.obs, count.cut), count.cut) #0 + 1 .. 
	} else {
		y.dat.vls <- hist(dat.obs, breaks=prob.int, plot=F)$counts #bin the counts
	}
	y.dat.vls <- (y.dat.vls/sum(y.dat.vls)) #the observed values w/ str.z placeholder
	x.dat.vls <- seq(0, length(y.dat.vls)-1) #the associated count values, -1 == str.z placeholder
	y.dat.cdf.vls <- cumsum(y.dat.vls) #the computed CDF

	#Conditional Probability Matrix:
	mat.indxs <- match(round(ind.rds, sig.level.round), tmat.vals)
	cond.p.mat <- apply(cond.pmat.precalc[,,mat.indxs], c(1,2), mean)
	r.m <- est.strz.pmat(cond.p.mat, 1, strz=TRUE, trunc=TRUE) #modify the probability matrix

	exp.ct <- r.m%*%wts.bt; exp.ct.cdf <- cumsum(exp.ct) #Expected probability

	#Calculate the distances:
	inc.indx <- 1:(length(y.dat.vls)-1) #exclude the truncation value? ####
	l2.d.norm <- sum((y.dat.vls[inc.indx] - exp.ct[inc.indx])^2) #l2 PDF norm on observed data (less truncation)

	f.p <- y.dat.vls[inc.indx]; g.p <- exp.ct[inc.indx] #KL Distance (less truncation)
	valid.indx <- which(g.p > 0 & f.p > 0)
	f.p <- f.p[valid.indx]; g.p <- g.p[valid.indx]
	kl.d.norm <- sum(f.p*log(f.p/g.p))

	l2.c.norm <- sum((exp.ct.cdf[inc.indx] - y.dat.cdf.vls[inc.indx])^2) #on descretized (less truncation)
	return(c(l2.d.norm, abs(kl.d.norm), l2.c.norm))
}


#####################################################################################################
#Calculate the individual weights probability given data, reads, mixtures and a weight initialization
#dat.out <- dat.ex.mat[,jj,1]; ind.rds <- dat.ex.mat[,jj,2]
#bpar.est <- model.mix.all; wts.est <- mdl.wts.calc; count.cut <- tr.val

calc.ind.wt.prob <- function(dat.out, ind.rds, bpar.est, wts.est, count.cut) {
	rds.avg <- ind.rds/mean(ind.rds)
	n.gam <- dim(bpar.est)[1]; n.obs <- length(dat.out)
	wts.est[which(wts.est < 0)] <- 0; wts.est <- wts.est/sum(wts.est)
	lpt.wts <- matrix(0, nrow=n.obs, ncol=n.gam+2) #Storage Array
	zero.ind <- 1; trunc.ind <- n.gam+2 #point mass indeces
	gam.indx <- 2:(n.gam+1)

	#Calculate the log-lik and weights:
	for (j in 1:n.gam) {
		p.j <- bpar.est[j,2]/(rds.avg+bpar.est[j,2])
		lpt.wts[,gam.indx[j]] <- log(wts.est[gam.indx[j]]) + bpar.est[j,1]*log(p.j) + dat.out*log(1-p.j) + 
			lgamma(dat.out + bpar.est[j,1]) - lgamma(bpar.est[j,1]) - lgamma(dat.out+1)
	}
	lpt.wts[,zero.ind] <- log(wts.est[zero.ind]) + log(as.integer(dat.out == 0))
	lpt.wts[,trunc.ind] <- log(wts.est[trunc.ind]) + log(as.integer(dat.out > count.cut))

	lpt.wts.norm <- t(apply(lpt.wts, 1, function(x){ x-log(sum(exp(x))) }))
	any.na.id <- apply(lpt.wts.norm, 1, function(x){any(is.na(x))}) #check is any are NaN
	if (any(any.na.id == TRUE)) {
		lpt.wts.norm[which(any.na.id),] <- t(apply(matrix(c(lpt.wts[which(any.na.id),]), nrow=sum(any.na.id), ncol=dim(lpt.wts.norm)[2], byrow=F), 1, function(x){ x1 <- (x-max(x)); x1 - log(sum(exp(x1))) }))
	}
	pt.wts.norm <- exp(lpt.wts.norm)
	return(pt.wts.norm)
}

