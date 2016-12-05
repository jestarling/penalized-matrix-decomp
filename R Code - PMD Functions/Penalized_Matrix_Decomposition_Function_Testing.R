#Penalized Matrix Decomposition
#November 10 ,2016
#Jennifer Starling

#References paper:
# Witten, Tibshirani, Hastie.  2009. A penalized matrix decomposition...

source()

#################################
###   REQUIRED FUNCTIONS:     ###
#################################

#l1 penalty; Soft thresholding operator.
soft <- function(x,theta){
  return(sign(x)*pmax(0, abs(x)-theta))
}

#l1 norm of a vector.
l1norm <- function(vec){
  a <- sum(abs(vec))
  return(a)
}
 
#l2 norm of a vector.
l2norm <- function(vec){
  a <- sqrt(sum(vec^2))
  return(a)
}

#Binary search function 
#(Source: Witten, Hastie & Tibshirani R package PMA: https://cran.r-project.org/web/packages/PMA/)
#(For finding theta for soft-thresholding for each iteration)
BinarySearch <- function(argu,sumabs){

  #Define functions for the l2 norm and the soft thresholding operator.
  l2n = function(vec) {return(sqrt(sum(vec^2)))}	
  soft = function(x,theta) { return(sign(x)*pmax(0, abs(x)-theta))}
	
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

#######################################################
###   PENALIZED MATRIX DECOMPOSITION FUNCTIONS:     ###
#######################################################

#---------------------------------------------------------------------
#SPARSE MATRIX FACTORIZATION RANK 1 FUNCTION (Penalized Matrix Decomposition) 
#For a single factor, ie K=1 (rank-1 approximation to original matrix)
#Inputs:
#	X = matrix to be factorized
#	v = initial v vector.
#	lambdaU = the u penalty (c1)
# 	lambdaV = the v penalty (c2)
#	   *If lambda1 = lambda2 = 0, function returns the non-sparse Rank 1 SVD of X.
# 	maxiter = maximum number of iterations allowed
#	tol = tolerance level for convergence check
#Output: List object, including the following:
#	X.rebuilt = sparse matrix factorization of X; X = U * D * t(V)
#	U, D, V = the decomposed elements of X, where X = U * D * t(V)

sp.matrix.decomp.rank1 = function(X,v=NULL,lambdaU=1, lambdaV=1, maxiter=20, tol=1E-6){ 
	
	#Define functions for the l2 norm and the soft thresholding operator.
  	l2norm = function(vec) {return(sqrt(sum(vec^2)))}	
  	soft = function(x,theta) { return(sign(x)*pmax(0, abs(x)-theta))}
  
	#1. Housekeeping parameters.
	i=1					#Initialize iterator.
	converged <- 0		#Indicator for whether convergence met.
	p = ncol(X)			#Number of columns of X matrix.
	
	#2. Initializations
	v.old = rnorm(p)	#Initialize v.old to a random vector. (To get iterations started.)
	
	if(is.null(v)){
		v = rep(sqrt(1/p),p) #Initialize v to meet constraint l2norm(v) = 1.		
	}
	
	#3. NA Handling:
	nas <- is.na(X)
  	if(sum(nas)>0) X[nas] = mean(X[!nas])

	
	#4. Iterate until convergence.
	for (i in 1:maxiter){
		
		#1. Update u.
		
		#First, calculate theta for sign function.
		u.arg = X %*% v		#Argument to go into sign function: Xv
		u.theta = BinarySearch(u.arg,lambdaU)
		
		#Second, update u.
		u = matrix( soft(u.arg,u.theta) / l2norm(soft(u.arg,u.theta)), ncol=1)

		#------------------------
		#2. Update v.
		
		#First, calculate theta for sign function.
		v.arg = t(X) %*% u
		v.theta = BinarySearch(v.arg,lambdaV)
		
		#Second, update v.
		v = matrix( soft(v.arg,v.theta) / l2norm(soft(v.arg,v.theta)), ncol=1)
		
		#------------------------
		#3. Convergence check steps.
		
		#Exit loop if converged.
		if(sum(abs(v.old - v)) < tol){
			converged=1
			break
		}
		
		#If not converged, update v.old for next iteration.
		v.old = v	
	}
	
	#Set d value.
	d = as.numeric(t(u) %*% (X %*% v))
	
	#Reconstruct sparse X matrix.
	Xsp = d * tcrossprod(u,v)
	
	#Return	function results.
	return(list(Xsp=Xsp,u=u,d=d,v=v,lambdaU=lambdaU,lambdaV=lambdaV,converged=converged,iter=i))
}	
#---------------------------------------------------------------------

SPARSE MATRIX FACTORIZATION RANK K FUNCTION (Penalized Matrix Decomposition) 
#For a Rank K approximation of X.
#Inputs:
#	X = matrix to be factorized
#	K = rank of factorization.  Must be <= ncol(X)
#	lambdaU = the u penalty (c1)
# 	lambdaV = the v penalty (c2)
#	   *If lambda1 = lambda2 = 0, function returns the non-sparse Rank K SVD of X.
# 	maxiter = maximum number of iterations allowed
#	tol = tolerance level for convergence check
#Output: List object, including the following:
#	Xsp = sparse matrix factorization of X.
#	U, D, V = the decomposed elements of X, where X = U * D * t(V)

sparse.matrix.factorization.rankK = function(X,K=1,lambdaU=1, lambdaV=1, maxiter=20, tol=1E-6){ 
	
	#1. Housekeeping parameters.
	i=1					#Initialize iterator.
	converged <- 0		#Indicator for whether convergence met.
	
	#2. NA handling for imputing missing obs.
	nas = is.na(X)	#Identify locations of NA values in matrix X.
	X.use = X		#Initialize X values to whole matrix.
	
	#2. Initializations
	D = numeric(K)	#Initialize K-length vector to hold D values.
	U = matrix(0,nrow=nrow(X),ncol=K)	#Initialize matrix to hold U values.
	V = matrix(0,nrow=ncol(X),ncol=K)	#Initialize matrix to hold V values.
	
	#Initialize v to be a nxK matrix, with each column having an l2 norm = 1.
	p = ncol(X)
	v = rep(sqrt(1/p),p)
	v = matrix(rep(v,each=K), ncol=K, byrow=T)
	
	#3. Iterate through Rank 1 recursion to obtain Rank K factorization.
	for(k in 1:K){
		#Run Rank 1 factorization.
		output = sp.matrix.decomp.rank1(X.use,v=matrix(v[,k],ncol=1),lambdaU, lambdaV, maxiter=20, tol)
		
		#Assign values.
		U[,k] = output$u
		V[,k] = output$v
		D[k] = output$d
		
		Xnew = X.use - output$d * tcrossprod(output$u, output$v)
		X.use[!nas] = Xnew[!nas]
	}
	
	#Reconstruct the X matrix as U * D * t(V).
	if (K==1) {	X.rebuilt = as.numeric(D) * U %*% t(V) }
	if (K > 1){	X.rebuilt = U %*% diag(D) %*% t(V) }

	
	#Return function outputs.
	return(list(U=U,V=V,D=D,X.rebuilt = X.rebuilt))
}	
	
#---------------------------------------------------------------------

#######################################################
###   SIMULATIONS FOR SANITY CHECKING:              ###
#######################################################

#Simulate a matrix to chek that results behaving sensibly.

X = matrix(rnorm(20),nrow=5,ncol=4)
n = nrow(X)
p = ncol(X)

#Paper notes that if you want u and v to be equally sparse, set a constant c,
#and let lambdaU = c*sqrt(n), and let lambdaV = c * sqrt(p)
c = 2
lambdaU = c*sqrt(n)
lambdaV = c*sqrt(p)

#Test package PMA functions:

K=4
#Initialize v to be a nxK matrix, with each column having an l2 norm = 1.
p = ncol(X)
v = rep(sqrt(1/p),p)
v = matrix(rep(v,each=K), ncol=K, byrow=T)

#Built-in PMA package (1) vs my functions (2).
builtins = MultiSMD(X,lambdaU,lambdaV,K=4,niter=20,v=v,trace=F,upos=F,uneg=F,vpos=F,vneg=F)
myfctn = sparse.matrix.factorization.rankK(X,K=4,lambdaU, lambdaV, maxiter=20, tol=1E-7)

#Now test with an X matrix with missing values.
#Randomly select values to set to NA.
n.elems = nrow(X) * ncol(X)
na.locs = sample(1:n.elems,size=n.elems*.3,replace=F)
Xmiss = X
Xmiss[na.locs] = NA

#Built-in PMA package (1) vs my functions (2).
builtins = MultiSMD(Xmiss,lambdaU,lambdaV,K=4,niter=20,v=v,trace=F,upos=F,uneg=F,vpos=F,vneg=F)
myfctn = sparse.matrix.factorization.rankK(Xmiss,K=4,lambdaU, lambdaV, maxiter=20, tol=1E-7)

###########################################
###   APPLICATION TO MARKETING          ###
###########################################

#Read in marketing data.
data = read.csv('/Users/jennstarling/UTAustin/2016_Fall_SDS 385_Big_Data/Course Data/social_marketing.csv',header=T)
X.counts = as.matrix(data[,-1])	#Gets rid of ID column.

#Square Root or Anscombe-transform data, because it is count data.
#Chose square root in this case.
X = sqrt(X.counts)

#Eyeball data to see which categories might be interesting.
cbind(colmeans=sort(colMeans(X),decreasing=T))

#Try various lambda values to see which categories are interesting.
result1 = sparse.matrix.factorization.rank1(X, lambdaU=5, lambdaV=5, maxiter=20, tol=1E-6)
cbind.data.frame(colnames(X),result1$v)

result2 = sparse.matrix.factorization.rank1(X, lambdaU=2, lambdaV=2, maxiter=20, tol=1E-6)
cbind.data.frame(colnames(X),result2$v)	#Show all.
cbind.data.frame(cols=colnames(X)[which(result2$v!=0)],v=result2$v[which(result2$v!=0)])
	#show non-zero only.

result3 = sparse.matrix.factorization.rank1(X, lambdaU=1.5, lambdaV=1.5, maxiter=20, tol=1E-6)
cbind.data.frame(colnames(X),result3$v)
cbind.data.frame(cols=colnames(X)[which(result3$v!=0)],v=result3$v[which(result3$v!=0)])
	#show non-zero only.








#####################################################
### DELETE THIS                                   ###
#####################################################

SMD <- function(x, sumabsu, sumabsv, niter=20,trace=TRUE, v, upos, uneg, vpos, vneg){
	
	l2n = function(vec){ return(sqrt(sum(vec^2)))}
  		
  	soft <- function(x,theta){
  		return(sign(x)*pmax(0, abs(x)-theta))
	}	
  		
  # This gets a single factor. Do MultiSMD to get multiple factors.
  nas <- is.na(x)
  v.init <- v
  xoo <- x
  if(sum(nas)>0) xoo[nas] <- mean(x[!nas])
  oldv <- rnorm(ncol(x))
  for(iter in 1:niter){
    if(sum(abs(oldv-v))>1e-7){
      oldv <- v
      if(trace) cat(iter,fill=F)
      # update u #
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      su <- soft(argu,lamu)
      u <- matrix(su/l2n(su),ncol=1)
      # done updating u #
      # update v #
      argv <- t(u)%*%xoo
      if(vpos) argv <- pmax(argv,0)
      if(vneg) argv <- pmin(argv,0)
      lamv <- BinarySearch(argv, sumabsv)
      sv <- soft(argv,lamv)
      v <- matrix(sv/l2n(sv),ncol=1)
      # done updating v #
    } 
  }
  d <- as.numeric(t(u)%*%(xoo%*%v))
  if(trace) cat(fill=TRUE)
  return(list(d=d, u=u, v=v, v.init=v.init))
}


MultiSMD <- function(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg){
  nas <- is.na(x)
  v.init <- v
  xuse <- x
  ds <- numeric(K)
  us <- matrix(0,nrow=nrow(x),ncol=K)
  vs <- matrix(0,nrow=ncol(x),ncol=K)
  for(k in 1:K){
    out <- SMD(xuse, sumabsu=sumabsu,sumabsv=sumabsv,niter=niter,v=matrix(v[,k],ncol=1), trace=trace, upos=upos, uneg=uneg, vpos=vpos, vneg=vneg)
    us[,k] <- out$u 
    vs[,k] <- out$v 
    ds[k] <- out$d
    res <- xuse - out$d*out$u%*%t(out$v) 
    xuse[!nas] <- res[!nas] # [!nas] is new on July 24 2009
  }
  
  X.rebuilt = us %*% diag(ds) %*% t(vs)
  return(list(u=us,v=vs,d=ds, v.init=v.init,X.rebuilt = X.rebuilt))
}
