
#These are the required packages. 
#Poibin has the functions necessary for calculating poisson binomial probabilities, and is used in the E-step.
require("poibin")

#Foreach contains a parallelized for loop, which speeds up the E and M steps. doParallel is the backend for 
#parallelization. There are other R packages that can act as a backend (doSNOW, doMC, etc). 
#If the parallelization causes errors, an alternative backend package may work....
require("foreach")
require("doParallel")

#...which will require changing this line to the new backend function. Alternatively, removing this 
#line will force the program to run serially, which is slower but less prone to errors. 
registerDoParallel(5)

#These are the packages that contain functions for neural net training (nnet) and parameter tuning (e1071)
require("nnet")
require("e1071")
		

#ParsSNP_resources.Rdata contains the PanCan dataset and the final ParsSNP model. 
#The descriptors have already been calculated, imputed and scaled. To see how new data can be used, 
#please see the script ParsSNP_application.r

load("ParsSNP_resources.Rdata")

#Columns 1 to 24 are the descriptors. To reproduce the results in the study, we begin the descriptors
#from the training set of the data.
descriptors<-as.matrix(PanCan[PanCan$Set=="Training", 1:24])

#The learning process uses the descriptor matrix as well as the case list
samples<-PanCan[PanCan$Set=="Training", "Sample"]

#The E-step will update mutation labels within each case using Bayes' law. 
#The work-horse function bayes.case() is nested within. 
#y is the current set of labels/probabilities. 
#index is a list of indices for each case, formated by EM()
expect<-function(y, indices, ...) {
	
	#This function accepts the current labels and indices indicating the chosen case.
	#It will calculate the upper and lower bounds, given the current probabilities. 
	#Then it will update the probabilities, given the bounds. 
										
	bayes.case<-function(probs, sliding.p=0.9, fixed.base=2, fixed.constant=1, lower=1, rna=30, fixed.type="l") {
		
		#A single probability can't be updated and is immediately returned. 
		N<-length(probs)
		if (N==1) {
			return(probs)
			} 
			
		#The lower bound cannot be greater than N-1
		lower<-min(N-1, lower)
		
		
		#First we calculate the fixed upper bound.
		#The fixed upper bound is (by default) log2(N), but can be defined in a variety of ways.
		if (fixed.type=="l") {
			#The fixed upper bound cannot be larger than the number of mutations in the case. 
			fixed<-min(N,round(log(N,base=fixed.base)))
			}
		else if (fixed.type=="c") {
			#The fixed upper bound can be held constant across all cases. 
			fixed<-min(N,fixed.constant)
			}
		else {
			fixed<-min(N,round(log(N,base=fixed.base)))
			print("Invalid upperbound in E-step. Default settings used.")
			}

		#The fixed upper bound cannot be lower than the lower bound. 
		fixed<-max(fixed, lower)
		
		#The sliding upper bound is defined by the current beliefs, and also cannot be lower than the lower bound.
		sliding<-max(trunc(sum(probs)*sliding.p), lower)
		
		#Finally, the upper bound is defined as the larger of the sliding and fixed versions. 
		upper<-max(fixed, sliding)

		#The rna cutoff controls when the refined normal approximation is used by poibin, rather than an exact calculation.
		if (length(probs)>rna) {
			type<-"RNA"
		} else {
			type<-"DFT-CF"
			}
	
		#Bayes' law is applied, relying on the cumulative density function ppoibin(). 

		#The denominator is calculated. The probability of having a number of pathogenic mutations within [lower,upper].
		pmf.den<-ppoibin(upper, probs, method=type)-ppoibin(lower-1,probs,method=type)
		
		#For each mutation, we calculate the probability of OTHER mutations falling within [lower-1,upper-1]
		pmf.num<-unlist(lapply(1:N, function(i) {
			ppoibin(upper-1, probs[-i], method=type)-ppoibin(lower-2,probs[-i],method=type)
			}))
		
		#The adjusted probabilities are calculated and returned.
		adjusted<-probs/pmf.den*pmf.num
		return(adjusted)
		}
	
	#For each entry in indice (i.e. case), we update the probabilities in paralelized fashion.
	y.out<-foreach(i=1:length(indices), .inorder=F, .packages="poibin", .verbose=F) %dopar% {
		bayes.case(y[indices[[i]]], ...)
		}
	
	#Output needs to be reordered back into the original format
	y.out<-unlist(y.out)[order(unlist(indices))]

	#There should be no non-probabilistc labels. 
	if (min(y.out, na.rm=T)<=0|max(y.out, na.rm=T)>=1|sum(is.na(y.out))!=0) {
		y.out[y.out>=1]<-max(y.out[y.out<1], na.rm=T)
		y.out[y.out<=0]<-min(y.out[y.out>0], na.rm=T)
		y.out[is.na(y.out)]<-mean(y.out, na.rm=T)
		print("E-step has introduced non-probabilistic values. Vector has been coerced.")
		}
	return(y.out)
	}
	



#The maximize function simply updates labels using a neural net in cross validation. 
#x is the matrix of numeric descriptors, while y are the current labels. 
#folds is a list, each entry representing the indices for mutations belonging to the fold. 
#For a 9 mutation dataset, e.g. folds=list(fold1=c(1,2,3), fold2=c(4,5,6), fold3=c(7,8,9))
#The folds are calculated in EM() on a sample-wise basis, and do not change during learning.
maximize<-function(x, y, folds) {
	
	#nn.cv is the function that runs in parallel, one instance for each fold/testset.
	#The indices from the fold list are used to define the test set. 
	nn.cv<-function(test, x, y) {
		
		#The test set is defined.
		testset<-x[test,]
		#The training set is defined.
		x<-x[-test,]
		y<-y[-test]
		
		#The optimal parameters are found.
		param<-tune(nnet,x,y,ranges=list(size=round(ncol(x)*c(0.25,0.5,0.75)),decay=c(0.001,0.01,0.1)),
			tunecontrol=tune.control(sampling="boot",nboot=10,boot.size=min(1,10000/nrow(x))),
			maxit=200,reltol=1e-3,abstol=1e-2,trace=F)$best.parameters
		
		#A model is trained.
		model<-nnet(x,y,size=param$size,decay=param$decay,trace=F)
		
		#The predictions are calculated and returned
		y.out<-predict(model,testset)

		return(y.out)	

		}
	
	#For each set of test indices defined by folds, nn.cv is called. nncv will update the labels only at the 
	#assigned test indices, so that the updated labels can be combined by simple addition.
	y.out<-foreach (i=1:length(folds), .inorder=F, .packages=c("e1071", "nnet"), .verbose=F) %dopar% {
		nn.cv(folds[[i]], x, y)
		}
	
	#The y.out vector is reordered to the original format.
	y.out<-unlist(y.out)[order(unlist(folds))]


	#There should be no non-probabilistc labels. 
	if (min(y.out, na.rm=T)<=0|max(y.out, na.rm=T)>=1|sum(is.na(y.out))!=0) {
		y.out[y.out>=1]<-max(y.out[y.out<1], na.rm=T)
		y.out[y.out<=0]<-min(y.out[y.out>0], na.rm=T)
		y.out[is.na(y.out)]<-mean(y.out, na.rm=T)
		print("M-step has introduced non-probabilistic values. Vector has been coerced.")
		}
	return(y.out)

	}


#EM is a wrapper function that accepts the descriptor matrix and sample list, and calls the M and E functions iteratively. 
#It takes the descriptor matrix (x) and the case list (cases). 
#The maxiter, rcut, and acut variables encode the iteration, relative and absolute cutoffs, respectively. 
#A custom starting y (custom.y) can be provided.
#And the number of folds used for maximization (m.folds) can be altered.
EM<-function(x, cases, maxiter=20, rcut=0.95, acut=1e-5, custom.y=NULL,m.folds=5, ...) {

	if (length(cases) != nrow(x)) stop("The cases vector does not match the descriptor matrix.")
	
	#Here, we construct a list of indices for each case. This is useful for the E-step, which operates
	#on a case-by-case basis.
	case.index<-by(1:length(cases), cases, function(i) i)
	case.counts<-unlist(lapply(case.index, length))
	
	#Now we split the dataset into m.folds equally sized case sets, which should be roughly equal in terms of mutations. 
	#These are fixed folds that are maintained throughout the learning process. 

	repeat {
		#Each case is assigned to a fold.
		folds<-sample(rep(1:m.folds,length.out=length(case.index)))
		tab<-tapply(case.counts, folds, sum)
		if (max(tab)/min(tab)<1.05) break
	}
	
	#Finally, we construct a list of indices for each fold. This is useful for the M-step, which 
	#operates fold-by-fold.
	m.folds<-tapply(case.index, folds, unlist)

		
	#Here we initialize the labels (y). The default is to use random uniform initialization. However, 
	#the user can provide a custom initialization. This is useful for debugging, or for prioritizing
	#some mutations at the outset, e.g. recurrent mutations start with higher values on average. 
	#Of course, the custom y values must be interpretable as probabilities. 
	if (is.null(custom.y)) {
		y<-runif(nrow(x))
		print("Labels randomly initialized.")
		}
	else {
		if (mode(custom.y)=="numeric" & min(custom.y)>0 & max(custom.y)<1) {
			y<-custom.y
			print("Labels initialized with custom values.")
			}
		else {
			y<-runif(nrow(x))
			print("Custom labels were not probabilities. Labels randomly initialized.")
			}
		}	

	#The iteration count and current performance are initialized. 
	iter<-1
	msd<-1
	repeat {
		y.old<-y #Acts as a reference. 
		msd.old<-msd
		
		#The E-step is called first each iteration.
		print(paste("Iteration: ", iter, ", Expectation Calculation. ", sep=""))
		flush.console()
		y<-expect(y, indices=case.index, ...)
		
		#Then the M-step is called.
		print(paste("Iteration: ", iter, ", Maximization Calculation.", sep=""))
		flush.console()
		y<-maximize(x, y, folds=m.folds)

		#Performance metrics are printed. 
		msd<-mean((y-y.old)^2)
		print(paste("Iteration: ", iter, ", Post-Max Median: ", median(y), sep=""))
		print(paste("Iteration: ", iter, ", Post-Max MSD: ", msd, sep=""))
		flush.console()
		
		#If the absolute cut-off, the relative cut-off, or the maxiter has been hit, we stop.
		if (msd<acut) {break}
		if (msd/msd.old>rcut) {break}
		if (iter >= maxiter) {break}
		#increment the iteration.
		iter<-iter+1
		}

	return(y)
	}



label<-EM(descriptors, samples)



param<-tune(nnet,descriptors,label,ranges=list(size=round(ncol(descriptors)*c(0.25,0.5,0.75)),decay=c(0.001,0.01,0.1)),
	tunecontrol=tune.control(sampling="boot",nboot=10,boot.size=min(1,10000/nrow(descriptors))),
	maxit=200,reltol=1e-3,abstol=1e-2,trace=F)$best.parameters

new.model<-nnet(descriptors,label,size=param$size,decay=param$decay,trace=T)

plot(new.model$fitted.values, as.vector(predict(ParsSNP, descriptors)), xlab="New Model Predictions", ylab="ParsSNP Predictions", main="Reproducing ParsSNP Model", pch=20, cex=0.3)