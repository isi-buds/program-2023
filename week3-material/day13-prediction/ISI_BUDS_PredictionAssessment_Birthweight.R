##
#####  Code for examples presented in ISI BUDS - Prediction Assessment
#####	 Author: D. Gillen
##
#####
#####	Source in needed libraries and functions
#####
##
library(MASS)
library(leaps)
##
#####
#####	Read in King county BW data
#####
##
path <- "/Users/dgillen/SynologyDrive/Teaching/ISI_BUDS/2022/Lecture4_Prediction2/"
weight <- read.table( paste(path,"KingCounty2001_data.txt",sep=""), header=TRUE )
names(weight)
## All records have 'plurality' = 1; so remove from dataset
weight <- weight[,-3]

##
#####
#####	Simple descriptives
#####
##
#####	Hist of birth weight
##
  hist(weight$bwt, xlab="", ylab="", main="", freq=FALSE, nclass=25, col="lightgrey", axes=FALSE)
  axis(1, at=seq(from=0, to=5000, by=1000))
  title(main="Birth weight, grams", cex.main=2, col.main="blue", font.main=1)

##
#####	Scatterplots
##
  par(mfrow=c(2,2), mar=c(3.1, 4.1, 3.1, 2.1))
  plot(weight$age, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Mother's age, years", cex.main=2, col.main="blue", font.main=1)
  plot(weight$smokeN, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Cigarettes smoked per day", cex.main=2, col.main="blue", font.main=1)
  plot(weight$drinkN, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Alchoholic drinks per week", cex.main=2, col.main="blue", font.main=1)
  plot(weight$educ, weight$bwt, xlab="", ylab="Birth Weight")
  title(main="Highest grade completed", cex.main=2, col.main="blue", font.main=1)



##
#####
#####	Best subsets regression (need to install 'leaps' library)
#####
##
## Model with only the intercept
##
fit0 <- lm(bwt ~ 1, data=weight)

## Perform best subsets analysis
##
maxModel <- as.formula(bwt ~ 	gender + age + race + parity + 
								married + smokeN + drinkN + 
								firstep + welfare + smoker + 
								drinker + wpre + educ)
fitF <- lm(maxModel, data=weight)
summary( fitF )

bestSubRSS  <- summary(regsubsets(maxModel, data=weight, nvmax=17, nbest=1))

## What is returned?
##
names(bestSubRSS)

## 'results' contains the subset size, k, and the residual sum of squares
##
results <- c(0, sum((weight$bwt- fitted(fit0))^2))
results <- rbind(results, cbind(apply(bestSubRSS$which, 1, sum)-1, bestSubRSS$rss))

##	##	Plot of RSS (naive) vs. subset size
  par(mfrow=c(1,1), mar=c(4.1, 4.1, 1, 1), cex=2)
  plot(results[,1], results[,2], type="n", axes=F, xlab="", ylab="", ylim=range(results[,2])*c(0.95, 1.05))
  points(results[,1], results[,2], pch=5)
  title(xlab="Subset size, k")
  mtext("Residual sum of squares", 2, cex=2, line=1)
  axis(1, at=c(0:17))
  lines(results[,1], results[,2], type="b", col="red", pch=20, lwd=3)


##
#####
#####	Now let's do best subsets with Cp as the criteria
#####
##
bestSubCp  <- leaps(x=model.matrix(fitF), y=weight$bwt, int=FALSE, nbest=1, method="Cp")
## 'results' contains the subset size, k, and the Cp value
##
results <- NULL
results <- rbind(results, cbind(apply(bestSubCp$which, 1, sum)-1, bestSubCp$Cp))

##	Plot of Mallow's Cp vs. subset size
  par(mar=c(4.1, 4.1, 1, 1), cex=2)
  plot(results[,1], results[,2], type="n", axes=F, xlab="", ylab="", ylim=range(results[,2])*c(0.95, 1.05))
  points(results[,1], results[,2], pch=5)
  title(xlab="Subset size, k")
  mtext("Mallow's Cp", 2, cex=2, line=1)
  axis(1, at=c(0:17))
  lines(results[,1], results[,2], type="b", col="red", pch=20, lwd=3)


##
#####
#####	Compare which were selected in the k=10 models...
#####
##
cbind( dimnames( model.matrix(fitF) )[[2]], bestSubCp$which[11,], bestSubRSS$which[11,] )

##
#####
#####	Stepwise selection using AIC
#####
##
fitStepAIC <- stepAIC( fit0, scope=maxModel, direction="forward" )
fitStepAIC
cbind( dimnames( model.matrix(fitF) )[[2]], bestSubCp$which[11,] )


##
#####
#####	Ridge regression applied to the King County data
#####
##
lmFit  <- lm(maxModel, data=weight)
ridgeFit0  <- lm.ridge(maxModel, data=weight, lambda=0)

##	To compare coefficients use print(ridgeFit0)
print(ridgeFit0)
lmFit$coef


##
#####	Let's consider a few penalization values
##
ridgeFit  <- lm.ridge(maxModel, data=weight, lambda=c(0, 100, 1000, 10000))
ridgeFit$coef

##
#####	Calculate effective DF for ridge regression
##
calcDF <- function(Xmat, lambda)
{
  p     <- ncol(Xmat)
  Lmat  <- Xmat %*% solve(t(Xmat) %*% Xmat + diag(lambda, nrow=p)) %*% t(Xmat)
  value <- sum(diag(Lmat))
  return(value)
}


## Effective degrees of freedom
##
maxLambda <- 25000
lambdaVal <- seq(from=0, to=maxLambda, length=100)
ridgeFit  <- lm.ridge(maxModel, data=weight, lambda=lambdaVal)
ridgeCoef <- matrix(ridgeFit$coef, nrow=nrow(ridgeFit$coef))
Xmat    <- model.matrix(lm(maxModel, data=weight))
effDF   <- rep(NA, length(lambdaVal))
for(i in 1:length(lambdaVal)) effDF[i] <- calcDF(Xmat, lambda=lambdaVal[i])


##	Plot of ridge regression coefficient estimates vs. smoothing parameter
  par(mar=c(4.1, 2.1, 1, 1), cex=1.5)
  plot(x=range(lambdaVal), y=range(ridgeCoef[c(1,2,13,14),]), xlab=quote("Smoothing parameter, " * lambda), ylab="", axes=F, type="n")
  title(main="Ridge regression coefficient estimates", font.main=1, col.main="blue")
  axis(1, at=seq(from=0, to=maxLambda, length=5))
  axis(2, at=round(seq(from=min(ridgeCoef[c(1,2,13,14),]), to=max(ridgeCoef[c(1,2,13,14),]), length=5)))
  abline(h = 0, col="red", lwd=4)
  for(i in 1:4) lines(lambdaVal, ridgeCoef[c(1,2,13,14)[i],], lty=i, lwd=4)
  legend(maxLambda*0.8,60, legend=dimnames(ridgeFit$coef)[[1]][c(1,2,13,14)], lty=c(1:4), lwd=4)


## Plot of smoothing parameter vs. effective df
  par(mar=c(4.1, 4.1, 1, 1), cex=1.5)
  plot(lambdaVal, effDF, ylim=c(1, max(effDF)), xlab=quote("Smoothing parameter, " * lambda),
  ylab=quote("Effective degrees of freedom, d(" * lambda * ")"), axes=F, type="l", col="red", lwd=3)
  title(main="Smoothing parameter vs. effective degrees of freedom", font.main=1, col.main="blue")
  axis(1, at=seq(from=0, to=maxLambda, length=5))
  axis(2, at=seq(from=1, to=max(effDF), by=2))



##	Plot of ridge regression coefficient estimates vs. effective df
  par(mar=c(4.1, 2.1, 1, 1), cex=1.5)
  plot(x=c(1, max(effDF)), y=range(ridgeCoef[c(1,2,13,14),]), xlab=quote("Effective degrees of freedom, d(" * lambda * ")"), ylab="", axes=F, type="n")
  title(main="Ridge regression coefficient estimates", font.main=1, col.main="blue")
  axis(1, at=seq(from=1, to=max(effDF), by=2))
  axis(2, at=round(seq(from=min(ridgeCoef[c(1,2,13,14),]), to=max(ridgeCoef[c(1,2,13,14),]), length=5)))
  abline(h = 0, col="red", lwd=4)
  for(i in 1:4) lines(effDF, ridgeCoef[c(1,2,13,14)[i],], lty=i, lwd=4)
  legend(1, 60, legend=dimnames(ridgeFit$coef)[[1]][c(1,2,13,14)], lty=c(1:4), lwd=4)



##
#####
#####	Comparison of AIC and BIC	
#####
##
Xmat       <- model.matrix(lm(maxModel, data=weight))
lambdaVec  <- seq(from=0, to=10000, length=500)
ridgeRslt <- ridge.Err(y=weight$bwt, Xmat=Xmat, lambdaVec=lambdaVec)
##
  par(mfrow=c(2,1), cex=1.25)
  ##
  plot(ridgeRslt[,3], ridgeRslt[,4], xlab=quote("degrees of freedom, df(" * lambda * ")"), ylab="", axes=F,
  type="l", col="red", lwd=2)
  title(main="AIC", col.main="blue", font.main=1)
  axis(1, at=seq(from=5, to=19, length=5))
  axis(2, at=round(seq(from=min(ridgeRslt[,4]), to=max(ridgeRslt[,4]), length=5)))
  ##
  plot(ridgeRslt[,3], ridgeRslt[,5], xlab=quote("degrees of freedom, df(" * lambda * ")"), ylab="", axes=F,
  type="l", col="red", lwd=2)  
  title(main="BIC", col.main="blue", font.main=1)
  axis(1, at=seq(from=5, to=19, length=5))
  axis(2, at=round(seq(from=min(ridgeRslt[,5]), to=max(ridgeRslt[,5]), length=5)))




##
#####
##### How does this compare with GCV???
#####
##
maxLambda <- 25000
lambdaVal <- seq(from=0, to=maxLambda, length=100)
select(lm.ridge(maxModel, data=weight, lambda=lambdaVal))
Xmat       <- model.matrix(lm(maxModel, data=weight))
calcDF(Xmat, lambda=252.53)


##
#####
#####		Computation of all prediction criteria for a ridge.lm fit
#####
##
set.seed(12345)
maxModel <- as.formula(bwt ~ 	gender + age + race + parity + 
								married + smokeN + drinkN + 
								firstep + welfare + smoker + 
								drinker + wpre + educ)
ridgeFit  <- lm.ridge(maxModel, data=weight, lambda=252.53)
ridge.predcrit( ridgeFit, formula=maxModel, data=weight, K=10, B=500, boot=TRUE, sigmaSq="calculate" )


