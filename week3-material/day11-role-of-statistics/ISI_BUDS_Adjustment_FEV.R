##
#####	Code for examples presented in ISI BUDS - Adjustment Lecture
#####	Author: D. Gillen
##
#####  Helper function to return vector for binary test
##
ifelse1 <- function (test, yes, no){
  if (test) yes
  else no
}

##
#####  Function to produce CIs for LM parameters (or exponentiated parameters)
##
lmCI <- function( model, expcoef=FALSE, robust=FALSE ){
  coef <- summary( model )$coef[,1]
  se <- ifelse1( robust, robust.se.lm(model)[,2], summary( model )$coef[,2] )
  tvalue <- coef / se
  pvalue <- 2*(1-pt(abs(tvalue), model$df.residual))
  if( expcoef ){
    ci95.lo <- exp( coef - qt(.975, model$df.residual) * se )
    ci95.hi <- exp( coef + qt(.975, model$df.residual) * se )
    est <- exp( coef )
  }
  else{
    ci95.lo <- coef - qt(.975, model$df.residual) * se
    ci95.hi <- coef + qt(.975, model$df.residual) * se
    est <- coef
  }
  rslt <- round( cbind( est, ci95.lo, ci95.hi, tvalue, pvalue ), 4 )
  colnames( rslt ) <- ifelse1( 	robust, 	
                                c("Est", "robust ci95.lo", "robust ci95.hi", "robust t value", "robust Pr(>|t|)"),
                                c("Est", "ci95.lo", "ci95.hi", "t value", "Pr(>|t|)") )			
  colnames( rslt )[1] <- ifelse( expcoef, "exp( Est )", "Est" )
  rslt
}

##
#####	FEV example
##
#####	Preliminary data description and management
##
path <- "/Users/dgillen/SynologyDrive/Teaching/ISI_BUDS/2022/Lecture2_DAGS_Confounding/"
fev <- read.table( paste(path,"fev.raw",sep=""), header=TRUE )
summary( fev )
fev$male <- as.numeric(fev$sex) - 1
table( fev$sex, fev$male )
fev$smoke <- as.numeric(fev$smoke) - 1
table( fev$smoke, fev$smoker )
fev <- fev[ , !is.element(names(fev), c("sex")) ]
summary( fev )
fev$logfev <- log( fev$fev )
fev$loght <- log( fev$height )

##
#####	Unadjusted comparison of log-fev by smoking status
##
fit.unadj <- lm( logfev ~ smoke, subset=age>=9, data=fev )
summary( fit.unadj )
lmCI( fit.unadj )
pdf( "/Volumes/LACIE/Teaching/Stat202/Winter2011/Lectures/Lecture4/graphs/logfevhist.pdf", width=8, height=6 )
hist( fev$logfev[ fev$age >= 9 ] )
dev.off()
lmCI( fit.unadj, expcoef=TRUE )

##
#####	Comparison of log-fev by smoking status with adjustment for age
##
fit.age <- lm( logfev ~ smoke + age, subset=age>=9, data=fev )
summary( fit.age )
lmCI( fit.age, expcoef=TRUE )
cbind( summary(fit.age)$sigma, summary(fit.unadj)$sigma )

##
#####	Additional adjustment for loght as a precision variable
##
fit.adj <- lm( logfev ~ smoke + age + loght, subset=age>=9, data=fev )
summary( fit.adj )
lmCI( fit.adj, expcoef=TRUE )
cbind( summary(fit.adj)$sigma, summary(fit.age)$sigma )



##
#####	Additional adjustment for loght as a (potential?) precision variable
##
fit.gender <- lm( logfev ~ smoke + age + loght + male, subset=age>=9, data=fev )
summary( fit.gender )
lmCI( fit.gender, expcoef=TRUE )
cbind( summary(fit.gender)$sigma, summary(fit.adj)$sigma )
