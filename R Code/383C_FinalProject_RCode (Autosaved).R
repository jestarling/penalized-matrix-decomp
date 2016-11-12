#SDS 385 Final Project
#Jesse Miller & Jennifer Starling
#Fall 2016

rm(list=ls()) #Clean workspace.

setwd("/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/penalized-matrix-decomp")

source(file="./R Code/Penalized_Matrix_Decomp_Functions.R")	#Read in Penalized Matrix Decomp functions.

#------------------------------------------------------------
#DATA LOADING:

#Read in data.
data2014 = read.csv(file="./Data/2014FA_AllData.csv",header=T)
data2015 = read.csv(file="./Data/2015FA_AllData.csv",header=T)
data = rbind(data2014,data2015)

#Alternative: Load Rdata object directly.
load("./Data/data.Rdata")

#------------------------------------------------------------
#DATA PROCESSING:

#Create Y variable as Y=1 (pass), Y=0 (fail) based on Math Grade >= C.
passing.grades = c("A+","A","A-","B+","B","B-","C+","C","C-")
data$Y = ifelse(data$Math.Grade %in% passing.grades,1,ifelse(data$Math.Grade=="",NA,0))

#Restrict data set to only cases where Y != NA, ie where pass/fail known.
data = data[!is.na(data$Y),]

#Create a subset of just the calculus classes 
data.calc = data[data$Math.Course %in% c("408C","408N","408K"),]

#Set up some initial dimensional information.
n = nrow(data)						#Total number of obs for both years combined.
n2014 = sum(data$YEAR==2014)		#Total number of obs for 2014 only.
n2015 = nrow(data$YEAR==2015)		#Total number of obs for 2015 only.
p = ncol(data)-2 					#Number of predictors (excl ID and year cols).

####################################################################
###   TOY EXAMPLES OF PENALIZED MATRIX DECOMPOSITION:            ###
####################################################################

#Example 1: A simulated matrix with no missing values.
#Illustrates how decreasing lambda penalty terms selects features.

X = matrix(rnorm(20),nrow=5,ncol=4)
n = nrow(X)
p = ncol(X)

#Paper notes that if you want u and v to be equally sparse, set a constant c,
#and let lambdaU = c*sqrt(n), and let lambdaV = c * sqrt(p)
c = 2
lambdaU = c*sqrt(n)
lambdaV = c*sqrt(p)

K=ncol(X)						#Set a K value for testing.  We'll use Rank 1 here.
K = 1
lambdas = seq(2,0,by=-.25)	#Vector of lambdaU=lambdaV values.
tests = list()				#Empty vector for holding test cases.
nonzero.x.cols = rep(0,length(lambdas))	#Empty vector for holding sparsity info.

#Loop through test cases.
for (i in 1:length(lambdas)){
	tests[[i]] = sparse.matrix.factorization.rankK(X,K,
					lambdaU= lambdas[i],
					lambdaV=lambdas[i],
					maxiter=20,tol=1E-6)
	Xnew = tests[[i]]$X.rebuilt				
	nonzero.x.cols[i] = nonzero.col.info(Xnew)$num.nonzero.cols
}

#Display results.
tests
cbind(lambdas,nonzero.x.cols)

#------------------------------------------------------------
#Example 2: Illustrate missing data imputation.

#Set up X matrix.
X = matrix(rnorm(20),nrow=5,ncol=4)
n = nrow(X)
p = ncol(X)

#Randomly select values to set to NA.
n.elems = nrow(X) * ncol(X)
na.locs = sample(1:n.elems,size=n.elems*.3,replace=F)
Xmiss = X
Xmiss[na.locs] = NA

K=ncol(X)
lambdas = 10	#Want a large lambda, not trying to induce sparsity in features here.

missing.test = sparse.matrix.factorization.rankK(X,K,
					lambdaU= lambdas,
					lambdaV=lambdas,
					maxiter=20,tol=1E-6)

Xmiss
missing.test$X.rebuilt

##################################################################################
###    IMPUTING MISSING CONTINUOUS DATA USING PENALIZED MATRIX FACTORIZATION.  ###
##################################################################################

#-----------------------------------------------------------------
#IMPUTING MISSING CONTINUOUS DATA USING PENALIZED MATRIX FACTORIZATION.

#----------------------------------
#1. Identify continuous predictors to be included in the imputation.
head(data)

#cols.continuous.subset = c(13,14,15,16,20,30,39,40,41,42,43,44)

cols.continuous.all = c(13:16,20:28,30:38,39:44)
X = data.matrix(data[,cols.continuous.all])  #Subset of just the desired continuous data cols.

#Cols.continuous.all includes 21-28 and 31-38, which subset does not.
#These cols include the breakdowns of ALEKS-1 and ALEKS-H.

	#NOTE: This technique works well when predictors have same scale.  
	#If predictors have vastly different scales, X must be scaled 
	#first to ensure logical predictions.

	#Since we are working with a few different groups of predictors, 
	#each which is related to a certain
	#type of test, it makes sense to work in groups.

#----------------------------------
#2. Try imputing all missing values for continuous columns at once.

#Scale and center data.  (By hand, so can back-transform.)
col.means = colMeans(X,na.rm=T)
col.sdevs = sqrt(apply(X,2,function(a) var(a,na.rm=T)))
X.scaled = scale(X)

#Impute missing values for all continuous predictors at once.
pmd = sparse.matrix.factorization.rankK(X.scaled,K=ncol(X),lambdaU=1000,lambdaV=1000,maxiter=20)
X.filled.scaled = pmd$X.rebuilt

head(X.scaled)
head(X.filled.scaled)

#Reverse scaling.
X.filled = (X.filled.scaled * col.sdevs) + col.means
colnames(X.filled) = colnames(X)
head(X)
head(X.filled)

#Method fails, ended up with wrong scale on many variables.
#This is a common problem, as noted by Hastie et al (1999), as scale 

#----------------------------------
#3. Impute missing variables in four related groups.
# 	These related groups are of similar scale.

#SAT VALUES:
X.sat = data.matrix(data[,39:41]) #Just the SAT variables.
pmd.sat = sparse.matrix.factorization.rankK(X.sat,K=ncol(X.sat),lambdaU=1000,lambdaV=1000,maxiter=20)
X.sat.filled = pmd.sat$X.rebuilt
colnames(X.sat.filled) = colnames(X.sat)

head(X.sat)
head(X.sat.filled)

#ACT VALUES:
X.act = data.matrix(data[,42:44]) #Just the SAT variables.
pmd.act = sparse.matrix.factorization.rankK(X.act,K=ncol(X.act),lambdaU=1000,lambdaV=1000,maxiter=20)
X.act.filled = pmd.act$X.rebuilt
colnames(X.act.filled) = colnames(X.act)

head(X.act)
head(X.act.filled)

#GPA VALUES:
X.gpa = data.matrix(data[,c(13:16,52)]) #Just the PGPA variables.
pmd.gpa = sparse.matrix.factorization.rankK(X.gpa,K=ncol(X.gpa),lambdaU=1000,lambdaV=1000,maxiter=20)
X.gpa.filled = pmd.gpa$X.rebuilt
colnames(X.gpa.filled) = colnames(X.gpa)

head(X.gpa)
head(X.gpa.filled)

#SCORE VALUES:
X.score = data.matrix(data[,c(20:28,30:38)]) #Just the PGPA variables.
pmd.score = sparse.matrix.factorization.rankK(X.score,K=ncol(X.score),lambdaU=1000,lambdaV=1000,maxiter=20)
X.score.filled = pmd.score$X.rebuilt
colnames(X.score.filled) = colnames(X.score)

head(X.score)
head(X.score.filled)

#----------------------------------

#Sanity check ranges of the output for each group.
apply(X.sat, 2, function(x) range(x,na.rm=T))
apply(X.sat.filled, 2, function(x) range(x))

apply(X.act, 2, function(x) range(x,na.rm=T))
apply(X.act.filled, 2, function(x) range(x))

apply(X.gpa, 2, function(x) range(x,na.rm=T))
apply(X.gpa.filled, 2, function(x) range(x))

apply(X.score, 2, function(x) range(x,na.rm=T)) #Some neg vals here, ask Jesse.
apply(X.score.filled, 2, function(x) range(x))

#----------------------------------

#Reconstruct data frame with imputed values.
# Construct an X version with just continuous vars,
# and construct a Data version with all vars. 
#Note: Cols not in same order as in original data set.
X.cont.new = cbind(X.gpa.filled[,1:3], X.score.filled, X.sat.filled, X.act.filled, X.gpa.filled[,4])
colnames(X.cont.new) = colnames(X)
data.cont.new = cbind(X.cont.new,Y=data$Y)	#Save new imputed data.
data.cont.old = cbind(X,Y=data$Y)			#Save old data (just continuous predictors).

#----------------------------------

#Compare the results of a logistic regression before and after data imputation.
#This quick check involves holding out 30% of the data, and obtaining
#a 'test error' for the held out 30%.  This is performed for the continuous variables 
#only, with and without imputed data.

train.rows = sample(nrow(X),nrow(X)*.7,replace=F)

train.missing = data.cont.old[train.rows,]
test.missing = data.cont.old[-train.rows,]

train.filled = data.cont.new[train.rows,]
test.filled = data.cont.new[-train.rows,]

#----------------------------------

#Test error for data with missing values.
lm.with.missing = glm(Y ~ ., family=binomial,data=as.data.frame(train.missing))
pred.missing = predict.glm(lm.with.missing,newdata=as.data.frame(test.missing[,-29]),type='response')
yhat.missing = ifelse(pred.missing >= .5,1,0)

yhat.temp = yhat.missing	#To handle values that are predicted as NA due to missing data.
yhat.temp[is.na(yhat.missing)] = 999
test.err.missing = sum(test.missing[,29] != yhat.temp) / length(yhat.filled)
test.err.missing

paste(sum(is.na(yhat.missing)),' values predicted as NA due to missing data.')

#----------------------------------
#Test error for data with imputed values.
lm.filled = glm(Y ~ ., family=binomial,data=as.data.frame(train.filled))
pred.filled = predict.glm(lm.filled,newdata=as.data.frame(test.filled[,-29]),type='response')
yhat.filled = ifelse(pred.filled >= .5,1,0)

test.err.filled = sum(test.filled[,29] != yhat.filled) / length(yhat.filled)
test.err.filled

#Conclusion: Imputing the missing data using penalized matrix decomposition drastically
#decreased the logistic regression test error.

#A few things:
#No additional model fitting or analysis has been done.  Model could of course be improved in many ways.
#Could look at using this functionality for variable selection, as well.

###################################################################
###   VARIABLE SELECTION USING PENALIZED MATRIX FACTORIZATION.  ###
###################################################################

#The following is an example of how decreasing the lambda penalty can 
#perform variable selection on the continuous variables.

#In this case, we are not worried about the scale of the values, so we
#will analyze all continuous predictors together.

cols.continuous.subset = c(13,14,15,16,20,30,39,40,41,42,43,44)
X = data.matrix(data[,cols.continuous.subset])  #Subset of just the desired continuous data cols.
X = scale(X)

#Impute missing values for all continuous predictors at once.  Vary the values of lambda.
#Will use a rank K factorization, so that we can get down to a single selected predictor.

test = list()
lambdas = c(1000,500,10,5,2.5,1,.5,0)
interesting.predictors = list()	#Empty list for holding interesting predictors for each lambda.
num.nonzero.x.cols = rep(0,length(lambdas))	#Empty vector for holding number of interesting predictors.

#Loop through test cases.
for (i in 1:length(lambdas)){
	tests[[i]] = sparse.matrix.factorization.rankK(X,K=1,
					lambdaU= lambdas[i],
					lambdaV=lambdas[i],
					maxiter=20,tol=1E-6)
	
	Xnew = tests[[i]]$X.rebuilt	
	
	Xnew.col.info = nonzero.col.info(Xnew)
			
	nonzero.x.cols[i] = Xnew.col.info$num.nonzero.cols
	interesting.predictors[[i]] = colnames(X[,Xnew.col.info$nonzero.cols.idx])
}

#Display results.
lambdas
interesting.predictors




