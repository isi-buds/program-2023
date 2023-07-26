##
#####  Code for ISI BUDS - Prediction Assessment Functions
#####   Author: D. Gillen
##
##
#####
#####  	Below is the function lm.predcrit() which computes
#####		commonly used prediction criteria for a linear regression model
#####
##

##
#####		Computation of the CV or k-fold CV statistic for a lm fit (squared error loss)
##
cv.lm <- function( lmFit, data, K="n", GCV=FALSE ){
  y <- model.frame(lmFit)[,1]
  yhat <- lmFit$fitted
  n <- length(yhat)
  lmFormula <- formula( lmFit )
  
  Xmat <- model.matrix(lmFit)
  p <- dim( Xmat )[2]
  H <- Xmat %*% solve( t(Xmat)%*%Xmat ) %*% t(Xmat)	
  if( GCV==FALSE ) cv <- mean( ( (y-yhat) / (1-diag(H)) )^2 ) 
  else cv <- mean( ( (y-yhat) / (1-sum(diag(H))/n) )^2 ) 
  
  cv.k <- NULL
  if( K !="n" ) {
    ord <- sample(1:n, n)
    y <- y[ord]
    data <- data[ord,]
    rss.k <- rep(NA,n)
    for( i in 1:K ){
      keep <- 1:ceiling(n/K) + (i-1)*ceiling(n/K)
      if( max(keep) > n ) keep <- min(keep):n
      fit.k <- lm( lmFormula, data=data[!is.element( 1:n, keep ),] )
      yhat.k <- predict( fit.k, newdata=data[keep,] )
      rss.k[keep] <- (y[keep] - yhat.k)^2
    }
    cv.k <- mean( rss.k )
  }
  return( c(cv, cv.k) )
}


##
#####		Computation of the bootstrap estimate of MSPE for a lm() fit
##
bsMSE.lm <- function( lmFit, data, B=1000 ){
  y <- model.frame(lmFit)[,1]
  yhat <- lmFit$fitted
  n <- length(yhat)
  lmFormula <- formula( lmFit )
  
  mse.bs <- rep(NA,B)
  rss.1out <- matrix( NA, nrow=B, ncol=n )
  for( i in 1:B ){
    keep <- sample( 1:n, size=n, replace=TRUE )
    fit.b <- lm( lmFormula, data=data[keep,] )
    
    ##	Simple bootstrap estimate of MSE
    mse.bs[i] <- mean( (y - predict(fit.b, newdata=data ))^2 )
    
    ##	Leave-one-out bootstrap		
    rss.1out[i,] <- ifelse( is.element( 1:n, keep ), NA, ( y - predict(fit.b, newdata=data ) )^2 ) 
    
  }
  bsMSE <- mean(mse.bs)
  bsLve1out <- mean( apply( rss.1out, 2, mean, na.rm=TRUE ) )
  bs.632 <- .368*bsMSE + .632*bsLve1out
  return( cbind( bsMSE, bsLve1out, bs.632 ) )
}


##
#####		Computation of common prediction error estimates for a lm() fit
##
lm.predcrit <- function( lmFit, data, GCV="FALSE", K="n", B=1000, boot=FALSE, sigmaSq="calculate" ){
  if( sigmaSq=="calculate" ) sigmaSq <- summary(lmFit)$sigma^2
  rss <- sum( lmFit$residuals^2 )
  n <- length(lmFit$fitted)
  p <- summary( lmFit )$df[1]
  loglik <- (-1/2) * (n*log(2*pi*sigmaSq) + rss / sigmaSq)
  
  mse <- rss / n
  Cp <- (rss + 2*p*sigmaSq) / n
  aic <- -2*loglik + 2*p
  bic <- -2*loglik + log(n)*p
  cv <- cv.lm( lmFit, data, K=K, GCV=GCV )
  bsRslt <- NULL
  if( boot ) bsRslt <- bsMSE.lm( lmFit, data, B )
  
  
  rslt <- data.frame( t(c( p, mse, Cp, aic, bic, cv, bsRslt )) )
  rownames( rslt ) <- ""
  fullNames <- c("df", "mse", "Cp", "aic", "bic", "cv", "cv.k", "bs.mse", "bs.1out", "bs.632")
  if( K != "n" & boot==TRUE ) names( rslt ) <- fullNames
  else if( K != "n" & boot==FALSE ) names( rslt ) <- fullNames[1:7]
  else if( K == "n" & boot==TRUE ) names( rslt ) <- fullNames[c(1:6,8:10)]
  else names( rslt ) <- fullNames[1:6]
  return( rslt )
}



##
#####
#####		Below is the function ridge.predcrit() which computes
#####		commonly used prediction criteria for a ridge regression model fit
#####		with lm.ridge() in the MASS library
#####
##

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

##
#####
#####	Functions to compute AIC and BIC for ridge regression
#####
##
ridge.Err <- function(y, Xmat, lambdaVec)
{
  n   <- nrow(Xmat)
  p   <- ncol(Xmat)
  nLambda  <- length(lambdaVec)
  rss <- df <- rep(NA, nLambda)
  for(i in 1:nLambda)
  {
    Lmat   <- Xmat %*% solve(t(Xmat) %*% Xmat + diag(lambdaVec[i], nrow=p)) %*% t(Xmat)
    rss[i] <- sum((y - Lmat%*%y)^2)
    df[i]  <- sum(diag(Lmat))
  }
  sigmaSq <- rss[1] / (n - df[1])
  value <- cbind(lambdaVec, rss, df, (rss + (2*df*sigmaSq)) / n, (rss/sigmaSq) + (log(n)*df))
	dimnames(value)[[2]] <- c("lambda", "rss", "df", "aic", "bic" )
  return(value)
}


##
#####		Computation of the k-fold CV statistic for a ridge.lm fit
##
cv.ridge <- function( ridgeFit, data, formula, K=10 ){
	lmFit <- lm( formula, data=data )
	y <- lmFit$model[,1]
	dsn.matrix <- model.matrix( lmFit )
	
	lambda <- ridgeFit$lambda
	ridgeCoef <- ridgeFit$coef / ridgeFit$scales
	if( ridgeFit$Inter ) ridgeCoef <- c( ridgeFit$ym - ridgeCoef %*% ridgeFit$xm, ridgeCoef )
	yhat <- dsn.matrix %*% ridgeCoef

	n <- length( yhat )
	ord <- sample(1:n, n)
	y <- y[ord]
	data <- data[ord,]
	rss.k <- rep(NA,n)
	for( i in 1:K ){
		keep <- 1:ceiling(n/K) + (i-1)*ceiling(n/K)
		if( max(keep) > n ) keep <- min(keep):n
		ridgeFit.k <- lm.ridge( formula, data=data[!is.element( 1:n, keep ),], lambda=lambda )
		coef.k <- ridgeFit.k$coef / ridgeFit.k$scales
		if( ridgeFit.k$Inter ) coef.k <- c( ridgeFit.k$ym - coef.k %*% ridgeFit.k$xm, coef.k )
		lmFit.k <- lm( formula, data=data[keep,] )
		dsn.matrix.k <- model.matrix( lmFit.k )
		if( dim( dsn.matrix.k )[2] == length(coef.k) ){
			yhat.k <- dsn.matrix.k %*% coef.k
			rss.k[keep] <- (y[keep] - yhat.k)^2
		}
	}
	cv <- mean( rss.k, na.rm=TRUE )
	return( cv )
}


##
#####		Computation of the bootstrap estimate of MSPE for a ridge.lm() fit
##
bsMSE.ridge <- function( ridgeFit, data, formula, B=1000 ){
	lmFit <- lm( formula, data=data )
	y <- lmFit$model[,1]
	dsn.matrix <- model.matrix( lmFit )
	
	lambda <- ridgeFit$lambda
	ridgeCoef <- ridgeFit$coef / ridgeFit$scales
	if( ridgeFit$Inter ) ridgeCoef <- c( ridgeFit$ym - ridgeCoef %*% ridgeFit$xm, ridgeCoef )
	yhat <- dsn.matrix %*% ridgeCoef
	n <- length( yhat )

	mse.bs <- rep(NA,B)
	rss.1out <- matrix( NA, nrow=B, ncol=n )
	for( i in 1:B ){
		keep <- sample( 1:n, size=n, replace=TRUE )
		ridgeFit.b <- lm.ridge( formula, data=data[keep,], lambda=lambda )
		coef.b <- ridgeFit.b$coef / ridgeFit.b$scales
		if( ridgeFit.b$Inter ) coef.b <- c( ridgeFit.b$ym - coef.b %*% ridgeFit.b$xm, coef.b  )
		lmFit.b <- lm( formula, data=data )
		dsn.matrix.b <- model.matrix( lmFit.b )
		yhat.b <- dsn.matrix.b %*% coef.b

		##	Simple bootstrap estimate of MSE
		mse.bs[i] <- mean( (y - yhat.b)^2 )
		
		##	Leave-one-out bootstrap		
		rss.1out[i,] <- ifelse( is.element( 1:n, keep ), NA, ( y - yhat.b )^2 ) 
		closeAllConnections()
	}
	bsMSE <- mean(mse.bs)
	bsLve1out <- mean( apply( rss.1out, 2, mean, na.rm=TRUE ) )
	bs.632 <- .368*bsMSE + .632*bsLve1out
	return( cbind( bsMSE, bsLve1out, bs.632 ) )
}

##
#####		Computation of common prediction error estimates for a ridge.lm() fit
##
ridge.predcrit <- function( ridgeFit, data, formula, K=10, B=1000, boot=FALSE, sigmaSq="calculate" ){
	lmFit <- lm( formula, data=data )
	y <- lmFit$model[,1]
	dsn.matrix <- model.matrix( lmFit )
	
	lambda <- ridgeFit$lambda
	ridgeCoef <- ridgeFit$coef / ridgeFit$scales
	if( ridgeFit$Inter==1 ) ridgeCoef <- c( ridgeFit$ym - ridgeCoef %*% ridgeFit$xm, ridgeCoef )
	ridgeFits <- dsn.matrix %*% ridgeCoef

	n <- length( ridgeFits )
	p <- calcDF( dsn.matrix, lambda )
	rss <- sum( (y - ridgeFits)^2 )
	if( sigmaSq=="calculate" ) sigmaSq <- rss / (n-p)

	loglik <- (-1/2) * (n*log(2*pi*sigmaSq) + rss / sigmaSq)
	mse <- rss / n
	aic <- -2*loglik + 2*p
	bic <- -2*loglik + log(n)*p
	cv <- cv.ridge( ridgeFit, data, formula, K=K )
	bsRslt <- NULL
	if( boot ) bsRslt <- bsMSE.ridge( ridgeFit, data, formula, B )
	
	rslt <- data.frame( cbind( p, mse, aic, bic, cv, bsRslt ) )
	rownames( rslt ) <- ""
	fullNames <- c("df", "mse", "aic", "bic", "cv", "bs.mse", "bs.1out", "bs.632")
	if( boot==TRUE ) names( rslt ) <- fullNames
			else names( rslt ) <- fullNames[1:5]
	return( rslt )
}
