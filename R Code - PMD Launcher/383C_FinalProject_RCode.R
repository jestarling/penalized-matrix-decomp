#SDS 385 Final Project
#Jesse Miller & Jennifer Starling
#Fall 2016

rm(list=ls()) #Clean workspace.
library(stats)	#For KDE estimation in histograms.

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

#Create Y variable as Y=1 (pass), Y=0 (fail) based on Math.GPA >= 1.5.
#Note - this is already in data set, named Math.P.F.
data$Y = data$Math.P.F

#passing.grades = c("A+","A","A-","B+","B","B-","C+","C","C-")
#data$Y = ifelse(data$Math.Grade %in% passing.grades,1,ifelse(data$Math.Grade=="",NA,0))

#Restrict data set to only cases where Y != NA, ie where pass/fail known.
data = data[!is.na(data$Y),]

#Restrict data to just the calculus classes 408C, 408N, 408K.
data = data[data$Math.Course %in% c("408C","408N","408K"),]

#Set up some initial dimensional information.
n = nrow(data)						#Total number of obs for both years combined.
n2014 = sum(data$YEAR==2014)		#Total number of obs for 2014 only.
n2015 = nrow(data$YEAR==2015)		#Total number of obs for 2015 only.
p = ncol(data)-5 					#Number of predictors (excl ID, Year, Math.GPA, Math.P.F., Y).

####################################################################
###   TOY EXAMPLES OF PENALIZED MATRIX DECOMPOSITION:            ###
####################################################################

#Example 1: Illustrate missing data imputation.

#Set up X matrix.
X = matrix(rnorm(20),nrow=5,ncol=4)
n = nrow(X)
p = ncol(X)

#Randomly select values to set to NA.
n.elems = nrow(X) * ncol(X)
na.locs = sample(1:n.elems,size=n.elems*.5,replace=F)
Xmiss = X
Xmiss[na.locs] = NA

K=ncol(X)
lambdas = 10	#Want a large lambda, not trying to induce sparsity in features here.

missing.test = sparse.matrix.factorization.rankK(X,K,
					lambdaU= lambdas,
					lambdaV=lambdas,
					maxiter=20,tol=1E-6)

round(Xmiss,2)
round(missing.test$X.rebuilt,2)

#------------------------------------------------------------

#Example 2: A simulated matrix with no missing values.
#Illustrates how decreasing lambda penalty terms selects features.

X = matrix(rnorm(20),nrow=5,ncol=4)
n = nrow(X)
p = ncol(X)

#Paper notes that if you want u and v to be equally sparse, set a constant c,
#and let lambdaU = c*sqrt(n), and let lambdaV = c * sqrt(p).  
c = .9
lambdaU = c*sqrt(n)
lambdaV = c*sqrt(p)

K = 1	#Set a K value for testing.  We'll use Rank 1 here.

c = seq(1,.1,by=-.2)
tests = list()							#Empty vector for holding test cases.
lamU = rep(0,length(c))					#Empty vector to store lambdaU values for each c.
lamV = rep(0,length(c))					#Empty vector to store lambdaV values for each c.
nonzero.x.cols = rep(0,length(lambdas))	#Empty vector for holding sparsity info.

#Loop through test cases.
for (i in 1:length(c)){
	lambdaU = c[i]*sqrt(n)
	lambdaV = c[i]*sqrt(p)
	
	tests[[i]] = sparse.matrix.factorization.rankK(X,K,
					lambdaU,
					lambdaV,
					maxiter=20,tol=1E-6)
	Xnew = tests[[i]]$X.rebuilt				
	nonzero.x.cols[i] = nonzero.col.info(Xnew)$num.nonzero.cols
	
	#Store lambda values.
	lamU[i] = lambdaU
	lamV[i] = lambdaV	
}

#Display results.
for (i in 1:length(tests)){
	print(tests[[i]]$X.rebuilt)
}

cbind(c,lambdaU=round(lamU,3),lambdaV=round(lamV,3),nonzero.x.cols)


##################################################################################
###    IMPUTING MISSING CONTINUOUS DATA USING PENALIZED MATRIX FACTORIZATION.  ###
##################################################################################

#-----------------------------------------------------------------
#IMPUTING MISSING CONTINUOUS DATA USING PENALIZED MATRIX FACTORIZATION.

#----------------------------------
#1. Identify continuous predictors to be included in the imputation.
head(data)

cols.continuous.all = c(13:16,18,20:28,30:38,39:44,45:50)  #45-51 are AP scores. 18 is UTMA score.
X = data.matrix(data[,cols.continuous.all])  #Subset of just the desired continuous data cols.

#Cols.continuous.all includes 21-28 and 31-38, which subset does not.
#These cols include the breakdowns of ALEKS-1 and ALEKS-H.

	#NOTE: This technique works well when predictors have same scale.  
	#If predictors have vastly different scales, X must be scaled 
	#first to ensure logical predictions.

	#Since we are working with a few different groups of predictors, 
	#each which is related to a certain
	#type of test, it makes sense to work in groups.

#Some info about how much data is missing:
nrow(X) * ncol(X) 	#Total data points in X.
sum(is.na(X))		#Total points missing in X.
sum(is.na(X)) / (nrow(X) * ncol(X)) #Proportion of data missing in X.


#----------------------------------
#2. Impute all missing values for continuous columns at once.

#Scale and center data.  (By hand, so can back-transform.)
col.means = colMeans(X,na.rm=T)
col.sdevs = sqrt(apply(X,2,function(a) var(a,na.rm=T)))
X.scaled = scale(X)

#Confirm col.means and col.sdevs are the correct values used by scale(X).
attr(scale(X),"scaled:scale")
attr(scale(X),"scaled:scale")

#Impute missing values for all continuous predictors at once.
pmd = sparse.matrix.factorization.rankK(X.scaled,K=ncol(X),lambdaU=1000,lambdaV=1000,maxiter=20)
X.filled.scaled = pmd$X.rebuilt

head(X.scaled)
head(X.filled.scaled)

#Reverse scaling.
X.filled = X.filled.scaled #Placeholder to initialize X.filled matrix.

for (i in 1:ncol(X.filled.scaled)){
	X.filled[,i] = X.filled.scaled[,i] * col.sdevs[i] + col.means[i]
}
colnames(X.filled) = colnames(X)

#----------------------------------
#3. Verify Results:

#Eyeball results.
head(X)
head(X.filled)

#Sanity check ranges of the output for each column.
range_orig = apply(X, 2, function(x) round(range(x,na.rm=T),2))
range_imputed = apply(X.filled, 2, function(x) round(range(x,na.rm=T),2))

cbind.data.frame(orig.min = range_orig[1,],
	orig.max = range_orig[2,], 
	imput.min = range_imputed[1,], 
	imput.max = range_imputed[2,])
	
#Histograms to compare distributions before and after imputing data.
#SATs & ACTs:
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/LaTeX Files/SAT_ACT_hist.jpg')
par(mfrow=c(2,3))
idx = c(23:28)	#SAT and ACT variable indices.
for (i in idx){
	hist(X[,i],freq=F,main=paste(colnames(X)[i]))
	points(density(X.filled[,i]),col='blue',type='l')
}
dev.off()

#GPA, UTMA values:
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/LaTeX Files/GPA_hist.jpg')
par(mfrow=c(2,3))
idx = c(1:5)
for (i in idx){
	hist(X[,i],freq=F,main=paste(colnames(X)[i]))
	points(density(X.filled[,i]),col='blue',type='l')
}
dev.off()

#AP Score values:
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/LaTeX Files/AP_score_hist.jpg')
par(mfrow=c(2,3))
idx = c(30:35)
for (i in idx){
	hist(X[,i],freq=F,main=paste(colnames(X)[i]))
	points(density(X.filled[,i]),col='blue',type='l')
}
dev.off()

#SCORE values:
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/LaTeX Files/SCORE_hist.jpg')
par(mfrow=c(3,6))
idx = c(5:22)
for (i in idx){
	hist(X[,i],freq=F,main=paste(colnames(X)[i]))
	points(density(X.filled[,i]),col='blue',type='l')
}
dev.off()

#----------------------------------
#4. Reconstruct entire data set, and save data object.

#Data set containing only continuous variables. 
#Missing values imputed.
data.cont.new = cbind(X.filled,Y=data$Y)	

#Data set containing all predictors.  
#Missing values imputed (for continuous predictors only).
data.categorical = data[,-c(cols.continuous.all,ncol(data))]	#Save continuous vars, minus y col.
data.all.new = cbind(data.categorical,data.cont.new)			#Cbind data set back together.

#Save old data (just continuous predictors). 
#Will be used to compare regression improvements.
data.cont.old = cbind(X,Y=data$Y)			

save(data.all.new,file='./Data/data.all.new.Rdata')
save(data.cont.new,file='./Data/data.cont.new.Rdata')
save(data.cont.old,file='./Data/data.cont.old.Rdata')


######################################################################
###    LOGISTIC REGRESSION WITH IMPUTED DATA VERSUS MISSING DATA.  ###
######################################################################

#Compare the results of a logistic regression before and after data imputation.
#This quick check involves holding out 30% of the data, and obtaining
#a 'test error' for the held out 30%.  This is performed for the continuous variables 
#only, with and without imputed data.
y_idx = ncol(data.cont.new)	#Col index for y column.

#1. Set up training and test data.
train.rows = sample(nrow(X),nrow(X)*.7,replace=F)

train.missing = data.cont.old[train.rows,]
test.missing = data.cont.old[-train.rows,]

train.filled = data.cont.new[train.rows,]
test.filled = data.cont.new[-train.rows,]

#2. Perform logistic regression and calculate test error using data set with missing values.
lm.with.missing = glm(Y ~ ., family=binomial,data=as.data.frame(train.missing))
pred.missing = predict.glm(lm.with.missing,newdata=as.data.frame(test.missing[,-y_idx]),type='response')
yhat.missing = ifelse(pred.missing >= .5,1,0)

yhat.temp = yhat.missing	#To handle values that are predicted as NA due to missing data.
yhat.temp[is.na(yhat.missing)] = 999
test.err.missing = sum(test.missing[,y_idx] != yhat.temp) / length(yhat.missing)

test.err.missing	#Display results.
paste(sum(is.na(yhat.missing)),'out of ',length(yhat.missing),' values predicted as NA due to missing data.')

#----------------------------------
#Test error for data with imputed values.
lm.filled = glm(Y ~ ., family=binomial,data=as.data.frame(train.filled))
pred.filled = predict.glm(lm.filled,newdata=as.data.frame(test.filled[,-y_idx]),type='response')
yhat.filled = ifelse(pred.filled >= .5,1,0)

test.err.filled = sum(test.filled[,y_idx] != yhat.filled) / length(yhat.filled)
test.err.filled

#Conclusion: Imputing the missing data using penalized matrix decomposition drastically
#decreased the logistic regression test error from .78 to .23

#A few things:
#No additional model fitting or analysis has been done.  Model could of course be improved in many ways.
#Could look at using this functionality for variable selection, as well.

###################################################################
###   COMPARING WITH IMPUTATION USING PREDICTOR MEANS.          ###
###################################################################

#Create a data set of just continuous predictors, which will be used
#to impute the means of the missing data.
col.means = colMeans(data.cont.old[,-ncol(data.cont.old)],na.rm=T) #Col means, removing Y col.

#Replace all NA values in each column with the respective column means.
data.cont.means = data.cont.old	#Start with the data set with missing values.
for (i in 1:length(col.means)){ #-1 because excluding Y, which has no missing values.
	data.cont.means[is.na(data.cont.means[,i]), i] = col.means[i]
}

#Perform logistic regression using means-imputed data set, for comparison.

#Set up test and train data sets.
train.means = data.cont.means[train.rows,]
test.means = data.cont.means[-train.rows,]

#Test error for data with mean-imputed values.
lm.means = glm(Y ~ ., family=binomial,data=as.data.frame(train.means))
pred.means = predict.glm(lm.means,newdata=as.data.frame(test.means[,-y_idx]),type='response')
yhat.means = ifelse(pred.means >= .5,1,0)

test.err.means = sum(test.means[,y_idx] != yhat.means) / length(yhat.means)
test.err.means


#################################################################
###   COMPARISON: MATRIX METHOD vs IMPUTATION USING SVD.      ###
#################################################################
library(bcv)

#impute.svd imputes the missing entries using a low-rank SVD approximation estimated by the EM algorithm. 

data.cont.svd = impute.svd(X[,-y_idx],k=35,maxiter=20)
data.cont.svd = cbind(data.cont.svd$x,data$Y)
colnames(data.cont.svd) = c(colnames(X),"Y")

#View ranges as sanity check.
apply(data.cont.svd,2,range)

train.svd = data.cont.svd[train.rows,]
test.svd = data.cont.svd[-train.rows,]

lm.svd = glm(Y ~ ., family=binomial,data=as.data.frame(train.svd))
pred.svd = predict.glm(lm.svd,newdata=as.data.frame(test.svd[,-y_idx]),type='response')
yhat.svd = ifelse(pred.svd >= .5,1,0)

test.err.svd = sum(test.svd[,y_idx] != yhat.svd) / length(yhat.svd)
test.err.svd


###################################################################
###   VARIABLE SELECTION USING PENALIZED MATRIX FACTORIZATION.  ###
###################################################################

#The following is an example of how decreasing the lambda penalty can 
#perform variable selection on the continuous variables.

#In this case, we are not worried about the scale of the values, so we
#will analyze all continuous predictors together.

#Data set setup for variable selection.  We will use the imputed data.
Xfull = scale(data.cont.new[,-36])
n = nrow(Xfull)
p = ncol(Xfull)

#Vector of c values, since imposing equal sparisty on u and v.
#c = seq(.5,.01,-.1)
c = c(.3,.2,.1,.05,.01)

#Initialization.
tests = list()							#Empty vector for holding test cases.
models = list()							#Store each logistic regression model.
test.errors = rep(0,length(c))			#Empty vector for holding test error for each predictor subset.


lamU = rep(0,length(c))					#Empty vector to store lambdaU values for each c.
lamV = rep(0,length(c))					#Empty vector to store lambdaV values for each c.

nonzero.x.cols = rep(0,length(c))	#Empty vector for holding sparsity info.
interesting.predictors = list()			#Empty list for holding interesting predictors for each lambda.
n.nonzero.cols = rep(0,length(c))	#Empty vector for holding number of interesting predictors.
	
#Loop through test cases.
for (i in 1:length(c)){
	
	#-----------------------------------------
	#Set up penalty parameters.
	lambdaU = c[i]*sqrt(n)
	lambdaV = c[i]*sqrt(p)
	
	#-----------------------------------------
	#Run matrix factorization.
	tests[[i]] = sparse.matrix.factorization.rankK(Xfull,K=ncol(Xfull),
					lambdaU,
					lambdaV,
					maxiter=20,tol=1E-6)
	
	#Store rebuilt X matrix and number of nonzero columns.				
	Xnew = tests[[i]]$X.rebuilt	
	Xnew.col.info = nonzero.col.info(Xnew)			
	n.nonzero.cols[i] = Xnew.col.info$num.nonzero.cols
	interesting.predictors[[i]] = colnames(Xfull)[Xnew.col.info$nonzero.cols.idx]
	
	#Store lambda values.
	lamU[i] = lambdaU
	lamV[i] = lambdaV	
	
	#-----------------------------------------
	#Perform logistic regression on subset of interesting predictors.
	data.temp = cbind(Xfull[,interesting.predictors[[i]]],Y=data$Y)
	y_idx = ncol(data.temp)
	
	#Set up train/test data.
	train.rows = sample(nrow(Xfull),nrow(Xfull)*.7,replace=F)
	train.temp = data.temp[train.rows,]
	test.temp = data.temp[-train.rows,]
	
	#Set up model.
	lm.temp = glm(Y ~ ., family=binomial,data=as.data.frame(train.temp))
	pred.temp = predict.glm(lm.temp,newdata=as.data.frame(test.temp[,-y_idx]),type='response')
	yhat.temp = ifelse(pred.temp >= .5,1,0)

	test.err.temp = sum(test.temp[,y_idx] != yhat.temp) / length(yhat.temp)
	test.err.temp

	#Store model and test error.
	models[[i]] = lm.temp
	test.errors[i] = test.err.temp
}

#Pick out model with minimum test error.
test.errors
min = which(test.errors == min(test.errors))
test.errors[min]
interesting.predictors[[min]]
n.nonzero.cols[min]

#Plot test error as a function of number of predictors.
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/LaTeX Files/test.errs.jpg')
plot(n.nonzero.cols,test.errors,type='l',col='blue',main='Test Error vs Number of Predictors',
	xlab='Number of Predictors',ylab='Test Error')
abline(h=.2159,col='black')
legend(20,.25,lty=c(1,1),col=c('blue','black'),legend=c('Test error','Full model test error'))
dev.off()

npred = lapply(interesting.predictors,length)

var.sel.results = list(c=c,lambdaU=lamU,lambdaV=lamV,models=interesting.predictors,test.errors=test.errors,n.predictors=npred)
save(var.sel.results,file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Final Project/varsel.Rdata')




