mymode <- function(x) {
  x.present <- x[which(as.character(x) != "")]#exclude empty entries, in case most students are missing a value
  ux <- unique(x.present)
  as.character(ux[which.max(tabulate(match(x.present, ux)))])
}

data.source <- data.all.new[data.all.new$Math.Course %in% c("408C","408N","408K"),c("YEAR","FSE","AdmitMode","FirstSch","CurrSch","ParIncome","FatherEd","MotherEd","FirstGen","HS.Calc","PGPA_NS","PGPA_BA","PGPA_EN","PGPA_LA","UTMA","AGGR.SCORE.1","RADICALS.1","NUMBERS.1","EQUATIONS.1","FUNCTIONS.1","POLYNOMIALS.1","RATIONAL.EXP.1","LOGARITHMS.1","TRIG.1","AGGR.SCORE.H","RADICALS.H","EQUATIONS.H","FUNCTIONS.H","POLYNOMIALS.H","RATIONAL.EXP.H","LOGARITHMS.H","TRIG.H","SATQ","SATV","ACTMath","ACTENG","AP.BIO","AP.CALC.AB","AP.CALC.BC","AP.CALC.AB.SUB","AP.CHEM","AP.CS.A","Math.P.F")]
#omit AP.CS.AB (all missing)
#omit NUMBERS.H (all 100)

data.CatMiss <- data.source
data.CatMode <- data.source

#create two new data sets with modes and "missing" imputed
for (i in 3:10){
  data.mode <- mymode(data.source[,i])
  data.CatMode[,i] <- replace(data.CatMode[,i],which(as.character(data.CatMode[,i])==""),data.mode)
  data.CatMiss[,i] <- factor(data.CatMiss[,i], levels=c(levels(data.CatMiss[,i]), 'missing'))
  data.CatMiss[,i] <- replace(data.CatMiss[,i],which(as.character(data.CatMiss[,i])==""),"missing")
}

#split into 70/30 training/testing sets
train.rows.7030 = sample(nrow(data.CatMiss),nrow(data.CatMiss)*.7,replace=F)

train.7030.CatMiss = data.CatMiss[train.rows.7030,-1]
test.7030.CatMiss = data.CatMiss[-train.rows.7030,-1]

train.7030.CatMode = data.CatMode[train.rows.7030,-1]
test.7030.CatMode = data.CatMode[-train.rows.7030,-1]

#identify the number of the last column
y_idx = ncol(train.7030.CatMode)

#predictions with culled full models
lm.7030.CatMode = glm(Math.P.F ~ FirstSch+CurrSch+ParIncome+FatherEd+MotherEd+FirstGen+HS.Calc+PGPA_NS+PGPA_BA+PGPA_EN+PGPA_LA+UTMA+AGGR.SCORE.1+RADICALS.1+NUMBERS.1+EQUATIONS.1+FUNCTIONS.1+POLYNOMIALS.1+RATIONAL.EXP.1+LOGARITHMS.1+TRIG.1+AGGR.SCORE.H+RADICALS.H+EQUATIONS.H+FUNCTIONS.H+POLYNOMIALS.H+RATIONAL.EXP.H+LOGARITHMS.H+TRIG.H+SATQ+SATV+ACTMath+ACTENG+AP.BIO+AP.CALC.AB+AP.CALC.BC+AP.CALC.AB.SUB+AP.CHEM, family=binomial,data=as.data.frame(train.7030.CatMode))
lm.7030.CatMiss = glm(Math.P.F ~ FirstSch+CurrSch+ParIncome+FatherEd+MotherEd+FirstGen+HS.Calc+PGPA_NS+PGPA_BA+PGPA_EN+PGPA_LA+UTMA+AGGR.SCORE.1+RADICALS.1+NUMBERS.1+EQUATIONS.1+FUNCTIONS.1+POLYNOMIALS.1+RATIONAL.EXP.1+LOGARITHMS.1+TRIG.1+AGGR.SCORE.H+RADICALS.H+EQUATIONS.H+FUNCTIONS.H+POLYNOMIALS.H+RATIONAL.EXP.H+LOGARITHMS.H+TRIG.H+SATQ+SATV+ACTMath+ACTENG+AP.BIO+AP.CALC.AB+AP.CALC.BC+AP.CALC.AB.SUB+AP.CHEM, family=binomial,data=as.data.frame(train.7030.CatMiss))
pred.7030.CatMode = predict.glm(lm.7030.CatMode,newdata=as.data.frame(test.7030.CatMode[,-y_idx]),type='response')
pred.7030.CatMiss = predict.glm(lm.7030.CatMiss,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.7030.CatMode = ifelse(pred.7030.CatMode >= .5,1,0)
yhat.7030.CatMiss = ifelse(pred.7030.CatMiss >= .5,1,0)

#test errors with full models
test.err.7030.CatMode = sum(test.7030.CatMode[,y_idx] != yhat.7030.CatMode) / length(yhat.7030.CatMode)
test.err.7030.CatMiss = sum(test.7030.CatMiss[,y_idx] != yhat.7030.CatMiss) / length(yhat.7030.CatMiss)

#predictions with cat variables only
lm.7030.CatMode = glm(Math.P.F ~ FirstSch+CurrSch+ParIncome+FatherEd+MotherEd+FirstGen+HS.Calc, family=binomial,data=as.data.frame(train.7030.CatMode))
lm.7030.CatMiss = glm(Math.P.F ~ FirstSch+CurrSch+ParIncome+FatherEd+MotherEd+FirstGen+HS.Calc, family=binomial,data=as.data.frame(train.7030.CatMiss))
pred.7030.CatMode = predict.glm(lm.7030.CatMode,newdata=as.data.frame(test.7030.CatMode[,-y_idx]),type='response')
pred.7030.CatMiss = predict.glm(lm.7030.CatMiss,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.7030.CatMode = ifelse(pred.7030.CatMode >= .5,1,0)
yhat.7030.CatMiss = ifelse(pred.7030.CatMiss >= .5,1,0)

#test errors with cat variables only
test.err.7030.CatMode.CatOnly = sum(test.7030.CatMode[,y_idx] != yhat.7030.CatMode) / length(yhat.7030.CatMode)
test.err.7030.CatMiss.CatOnly = sum(test.7030.CatMiss[,y_idx] != yhat.7030.CatMiss) / length(yhat.7030.CatMiss)


#hereafter, use only "missing", not mode

#Test error with culled full model
lm.Full = glm(Math.P.F ~ FirstSch+CurrSch+ParIncome+FatherEd+MotherEd+FirstGen+HS.Calc+PGPA_NS+PGPA_BA+PGPA_EN+PGPA_LA+UTMA+AGGR.SCORE.1+RADICALS.1+NUMBERS.1+EQUATIONS.1+FUNCTIONS.1+POLYNOMIALS.1+RATIONAL.EXP.1+LOGARITHMS.1+TRIG.1+AGGR.SCORE.H+RADICALS.H+EQUATIONS.H+FUNCTIONS.H+POLYNOMIALS.H+RATIONAL.EXP.H+LOGARITHMS.H+TRIG.H+SATQ+SATV+ACTMath+ACTENG+AP.BIO+AP.CALC.AB+AP.CALC.BC+AP.CALC.AB.SUB+AP.CHEM, family=binomial,data=as.data.frame(train.7030.CatMiss))
pred.Full = predict.glm(lm.Full,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.Full = ifelse(pred.Full >= .5,1,0)
test.err.Full = sum(test.7030.CatMiss[,y_idx] != yhat.Full) / length(yhat.Full)

#Test error with UTMA only
lm.UTMA = glm(Math.P.F ~ UTMA, family=binomial,data=as.data.frame(train.7030.CatMiss))
pred.UTMA = predict.glm(lm.UTMA,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.UTMA = ifelse(pred.UTMA >= .5,1,0)
test.err.UTMA = sum(test.7030.CatMiss[,y_idx] != yhat.UTMA) / length(yhat.UTMA)

#Test error with SAT/ACT only
lm.SACT = glm(Math.P.F ~ SATQ + ACTMath, family=binomial,data=as.data.frame(train.7030.CatMiss))
pred.SACT = predict.glm(lm.SACT,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.SACT = ifelse(pred.SACT >= .5,1,0)
test.err.SACT = sum(test.7030.CatMiss[,y_idx] != yhat.SACT) / length(yhat.SACT)

require(MASS)

#line below does stepwise feature selection
stepAIC(lm.Full,direction="both")

#stepwise results in the model tested below
lm.step <- glm(formula = Math.P.F ~ FirstSch + ParIncome + HS.Calc + PGPA_NS + PGPA_BA + PGPA_EN + PGPA_LA + UTMA + AGGR.SCORE.1 + EQUATIONS.1 +POLYNOMIALS.1 + RATIONAL.EXP.1 + LOGARITHMS.1 + AGGR.SCORE.H +FUNCTIONS.H + LOGARITHMS.H + SATQ + SATV + ACTMath + ACTENG +AP.CALC.BC + AP.CALC.AB.SUB + AP.CHEM, family = binomial,data = as.data.frame(train.7030.CatMiss))
pred.step = predict.glm(lm.step,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.step = ifelse(pred.step >= .5,1,0)
test.err.step = sum(test.7030.CatMiss[,y_idx] != yhat.step) / length(yhat.step)

#line below does backward feature selection
stepAIC(lm.Full,direction="backward")

#backward selection results in the model tested below
lm.backward <- glm(formula = Math.P.F ~ FirstSch + ParIncome + HS.Calc + PGPA_NS +PGPA_BA + PGPA_EN + PGPA_LA + UTMA + AGGR.SCORE.1 + EQUATIONS.1 +POLYNOMIALS.1 + RATIONAL.EXP.1 + LOGARITHMS.1 + AGGR.SCORE.H +FUNCTIONS.H + LOGARITHMS.H + SATQ + SATV + ACTMath + ACTENG +AP.CALC.BC + AP.CALC.AB.SUB + AP.CHEM, family = binomial,data = as.data.frame(train.7030.CatMiss))
pred.backward = predict.glm(lm.backward,newdata=as.data.frame(test.7030.CatMiss[,-y_idx]),type='response')
yhat.backward = ifelse(pred.backward >= .5,1,0)
test.err.backward = sum(test.7030.CatMiss[,y_idx] != yhat.backward) / length(yhat.backward)


