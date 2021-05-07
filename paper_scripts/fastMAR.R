
#############################################################################
### fastMAR

fastMAR = function(
			data_vars, 
			log_transform = FALSE, 
			demeaned = FALSE, 
			abstol_ = 0.5, 
			maxit_ = 100,
			includeNoise = TRUE,
			estInit = FALSE,
			tSteps_ = NULL,
			nsim_ = 1
			) 
	{
    	### function to prepare data for MARSS ###
    	# data_vars # abundances of the dependent variables by column. the columns should by named.
    	# no data_time because I do not know if I can feed MARSS data with NAs 
	# log_transform # FALSE - use original values of the dep. variables. # TRUE use ln-values of the dep. variables.
	# demeaned # FALSE = use values as they are # TRUE = removes the mean before estimation, it returns the mean for the data estimates
	# abstol_ # The logLik.(iter-1)-logLik.(iter) convergence tolerance for the maximization routine. To meet convergence both the abstol and slope tests must be passed. 
	# maxit_ # Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). 
	# noise # TRUE - include noise in the  time series estimation process # FALSE - do not include noise
	# estInit # FALSE - first values of the data set used as initial conditions # TRUE - initial conditions are estimated.


	if ( log_transform == TRUE ) {data_vars = log(data_vars)}

	n = dim(data_vars)[1]

	if (is.null(n)) 
		{
		numVar = 1 
		n = length(data_vars) 
		data_vars_toUSE = data_vars

		data_vars_means = mean(data_vars)
		data_vars_means_mat = rep(data_vars_means,n)
		data_vars_demeaned = data_vars - data_vars_means_mat
		if ( demeaned == TRUE ) {data_vars_toUSE = data_vars_demeaned} 

		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal", 		# matrix(list("a1","a2"),nrow=2,ncol=1), 							## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  			## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero"	                          											## space eq. var-cov matrix
			)

		if (estInit == FALSE) { model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[1])[[1]],nrow=numVar,ncol=1) )) }

		MARvalEst1 = x_t = data_vars_toUSE[1]

		} else
		{
		numVar = dim(data_vars)[2]		
		data_vars_toUSE = data_vars

		data_vars_means = rep(NA,numVar)
		data_vars_means_mat = matrix(rep(NA,n*numVar),ncol=numVar)
		for (i in 1:numVar)
			{
			data_vars_means[i] = mean(data_vars[,i],na.rm=TRUE)
			data_vars_means_mat[,i] = rep(data_vars_means[i],n)
			}
		data_vars_demeaned = data_vars - data_vars_means_mat
		if ( demeaned == TRUE ) {data_vars_toUSE = data_vars_demeaned} 		

		data_vars_toUSE = t(data_vars_toUSE)

		### fitting the MAR model.
		model.gen = list(
	     		B = "unconstrained", 		# matrix(list("b11","b12","b21","b22"),nrow=2,ncol=2,byrow=T),  		## Specify the "Full model" matrix for the MAR community matrix
           		U = "unequal",		# matrix(list("a1","a2"),nrow=2,ncol=1), 							## state eq. intercept
           		Q = "diagonal and unequal",	# matrix(list("Q",0,0,"Q"),nrow=2,ncol=2),   					## state eq. var-cov matrix
			Z = "identity",               									  	  			## conversion matrice from space to state
           		A = "zero",     			# matrix(list(0,0),nrow=2,ncol=1), 							## space eq. intercept
           		R = "zero"	                          											## space eq. var-cov matrix
           		#x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1)							## mean of the distrib from which starting values are drawn
			)

		if (estInit == FALSE) { model.gen = c(model.gen,list( x0 = matrix(list(data_vars_toUSE[,1])[[1]],nrow=numVar,ncol=1) )) }

		MARvalEst1 = x_t = data_vars_toUSE[,1]
		}

	library(MARSS)
	library(MASS)
	MARSS_output1 = MARSS(data_vars_toUSE ,model=model.gen,control=list(conv.test.slope.tol=0.01,abstol=abstol_,maxit=maxit_,allow.degen=TRUE),form="marxss",silent=F) 	# abstol=0.0001,maxit=10000
	try(m.ci1 <- MARSSparamCIs(MARSS_output1,method="hessian",alpha=0.05,nboot=100),silent=T)
		
	MARparEst1 = cbind(matrix(MARSS_output1$par$B, ncol=numVar),MARSS_output1$par$U,MARSS_output1$par$Q)

	if (estInit == TRUE) { MARvalEst1 = x_t = as.vector(MARSS_output1$par$x0) }

	for (i in 2:n)
		{
		x_tplus1 = MARparEst1[,1:numVar] %*% x_t + MARparEst1[,numVar+1] + includeNoise*rnorm(n=numVar, mean = rep(0,numVar), sd = sqrt(MARparEst1[,(numVar+2)]))		#
		MARvalEst1 = rbind(MARvalEst1,t(x_tplus1))
		x_t = x_tplus1	
		}

	if ( is.null(tSteps_) ) { tSteps_ = n }
	MARSSsim1 = MARSSsimulate(object = MARSS_output1, tSteps = tSteps_, nsim = nsim_)$sim.data
	sim1dim = dim(MARSSsim1)

	remeaner = matrix(rep(0,numVar*n),ncol=numVar)
	remeaner_sim = array(rep(0,sim1dim[1]*sim1dim[2]*sim1dim[3]),dim=sim1dim)
	if ( demeaned == TRUE ) 
		{ 
		remeaner = data_vars_means_mat 
		for (iv in 1:nsim_) { remeaner_sim[,,iv] = data_vars_means_mat }
		}

	if ( log_transform == TRUE ) 
		{
		return( list(
			MARSSsim_data = exp(MARSSsim1 + remeaner_sim),
			MARSS_output = MARSS_output1,
			MARSSci_output = m.ci1,
			MAR_TS_Est = exp(MARvalEst1 + remeaner)
			))
		} else
		{	
		return( list(
			MARSSsim_data = MARSSsim1 + remeaner_sim,
			MARSS_output = MARSS_output1,
			MARSSci_output = m.ci1,
			MAR_TS_Est = MARvalEst1 + remeaner
			))
		}
	}


meanMaker = function(MARsim)
	{
	sim1 = MARsim$MARSSsim_data
	dim(sim1)

	if (dim(sim1)[1]==1)
		{
		meanStore = MARsim$MARSSsim_data[,,1]
		for (ii in 1:dim(sim1)[1])
			{
			for (i in 1:dim(sim1)[2])
				{
				meanStore[i] = mean(sim1[ii,i,])
				}
			}
		} else
		{
		meanStore = MARsim$MARSSsim_data[,,1]
		for (ii in 1:dim(sim1)[1])
			{
			for (i in 1:dim(sim1)[2])
				{
				meanStore[ii,i] = mean(sim1[ii,i,])
				}
			}
		}

	return(meanStore)
	}

### fastMAR
#############################################################################








# MARSSinfo("convergence")
MARSS tests both the convergence of the log-likelihood and of the
individual parameters.  If you just want the log-likelihood, say for
model selection, and don't care too much about the parameter values,
then you will be concerned mainly that the log-likelihood has
converged.  Set abstol to something fairly small like 0.0001 (in your
MARSS call pass in control=list(abstol=0.0001) ).  Then see if a
warning about logLik being converged shows up.  If it doesn't, then you
are probably fine.  The parameters are not at the MLE, but the
log-likelihood has converged.  This indicates ridges or flat spots in
the likelihood.

If you are concerned about getting the MLEs for the parameters and they
are showing up as not converged, then you'll need to run the algorithm
longer (in your MARSS call pass in control=list(maxit=10000) ). But
first think hard about whether you have created a model with ridges and
flat spots in the likelihood.

Do you have parameters that can create essentially the same pattern in
the data?  Then you may have created a model where the parameters are
confounded.  Are you trying to fit a model that cannot fit the data?
That often causes problems.  It's easy to create a MARSS model that is
logically inconsistent with your data.  Are you trying to estimate both
B and U? That is often problematic.  Try demeaning your data and
setting U to zero.  Are you trying to estimate B and you set tinitx=0?
tinit=0 is the default, so it is set to this if you did not pass in
tinitx in the model list. You should set tinitx=1 when you are trying
to estimate B.




