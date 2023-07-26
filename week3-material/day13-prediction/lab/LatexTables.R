##
##### Functions to aid in generating a "Table 1" and GLM estimates tables
##
ifelse1 <- function (test, yes, no){
  if (test) yes
  else no
}

##
####  summstat() is a helper function to produce summary statistics
##
summstat <- function(x, sigdigits=2, summtype="mean"){
  switch(summtype,
   mean=paste( format( mean(x,na.rm=TRUE), digits=sigdigits, nsmall=sigdigits ), 
         " (", format( sd(x,na.rm=TRUE), digits=sigdigits, nsmall=sigdigits ), ")", sep="" ),
   count=paste( "(N=", length(x), ")", sep="" ),
   freq=paste( table(x), 
               " (", format( table(x)/length(x), digits=1, nsmall=1 ), ")", sep="" ))
}


##
##### table1_latex()
#####   Produces a table 1 that can be stratified by covariate
#####     - data : dataset that variables are pulled from
#####     - vars : vector with variable names that are to be inlcuded
#####     - by.cov : variable name that table should be stratified by
#####     - var.lbls : Labels for vars that provide nice formatting
#####     - fileout : Name of the file to write the results to
#####
##### Notes:
#####   1.  For easiest implementation, make discrete variables in your
#####       dataset factors and add nice labels to them
#####   2.  The output is a .tex file that can simply be included in 
#####       a Latex document using:
#####       
              #\begin{table}
              #\caption{Your caption...}
              #\input{<path/filename}
              #\end{table}

table1_latex <- function( data, vars, by.cov, var.lbls=NULL, fileout="table1.txt" ){
  grp.names <- ifelse1( !is.null(levels(data[,by.cov])), 
                     levels(data[,by.cov]), 
                     paste("Group", 1:length(table(data[,by.cov]))) )
  n.grps <- length(table(data[,by.cov]))
  if( is.null(var.lbls) ) var.lbls <- vars
  cat("\\begin{tabular}{l", rep("c",2*n.grps-1), "}", "\n", 
      "\\hline", "\n", file=fileout )
  cat( c( "&", rbind(grp.names, c(rep("&&",n.grps-1),"\\\\ \n") ) ), 
       file=fileout, append=TRUE )
  cat( "\\multicolumn{1}{c}{Covariate}", "&", paste(unlist( lapply( split(data[,by.cov],data[,by.cov]), summstat, summtype="count" )), 
                  c(rep("&&",(n.grps-1)),"\\\\")), "\n \\hline \n", file=fileout, append=TRUE )
  for( i in 1:length(vars)){
    summtype <- ifelse(length(table(data[,vars[i]])) <= 8, "freq", "mean")
    if(summtype=="freq"){
      cat( var.lbls[i], "\\\\ \n", file=fileout, append=TRUE )
      n.lvls <- length(levels(data[,vars[i]]))
      rslt <- unlist( lapply( split(data[,vars[i]],data[,by.cov]), summstat, sigdigits=2, summtype=summtype ))
      write.table( cbind( "~~~", levels(data[,vars[i]]), "&", 
                          matrix(paste(rslt, c(rep("&&",n.lvls*(n.grps-1)),rep("\\\\",n.lvls))),nrow=n.lvls)), 
                   file=fileout, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE )
      }
    else{
      cat( c( var.lbls[i], "&", 
              rbind(unlist( lapply( split(data[,vars[i]],data[,by.cov]), summstat, sigdigits=2, summtype=summtype )), 
                    c(rep("&&",n.grps-1),"\\\\ \n") ) ), file=fileout, append=TRUE )
      }
    }
  cat( "\\hline", "\n", "\\end{tabular}", file=fileout, append=TRUE )
}


##
##### Example calls utilizing the salary data
##
##### Local path for Stat 211 files (modify for your path structure)
##
path <- "/Users/dgillen/SynologyDrive/Teaching/STAT211/Winter2022/"
#####  Read in functions defined especially for Stat 211
##
source( paste(path,"Stat211Functions.R",sep="") )

##
##### Read in salary data and create factors where needed
##
salary <- read.table( paste(path,"Data/salary95.raw",sep=""), header=TRUE )
salary$admin <- factor(salary$admin,labels=c("No","Yes"))
levels(salary$gender) <- c("Female","Male")

##
##### Table 1 stratified by sex
##
table1_latex( data=salary, 
              vars=c("deg","yrdeg","field","startyr","rank","admin"), 
              by.cov="gender", 
              var.lbls=c("Degree","Year of Degree", 
                         "Field", "Starting Year", "Rank", 
                         "Administrative Duties"), 
              fileout=paste(path,"LatexTables/Table1_Ex1.tex",sep="") )

##
##### Table 1 stratified by rank
##
table1_latex( data=salary, 
              vars=c("gender", "deg","yrdeg","field","startyr","admin"), 
              by.cov="rank", 
              var.lbls=c("Sex", "Degree", "Year of Degree", 
                         "Field", "Starting Year",  
                         "Administrative Duties"), 
              fileout=paste(path,"LatexTables/Table1_Ex2.tex",sep="")  )


##
##### glm_est_latex()
#####   Produces a regression table formatted with CI's and p-values
#####     - fit : model fit from a call to glm()
#####     - robust : use of robust variance estimator for inference?
#####     - fileout : Name of the file to write the results to
#####
##### Notes:
#####   1.  I default to robust=TRUE, and change to false as necessary
#####   2.  The labeling will identify standard LM, logistic and Poisson models
#####   3.  Just include your table in Latex as follows...
#####       
#\begin{table}
#\caption{Your caption...}
#\input{<path/filename}
#\end{table}
glm_est_latex <- function(fit, robust=TRUE, fileout){
  rslt <- glmCI(fit, robust=robust)[-c(1),]
  est_type <- switch( fit$family$family,
                      gaussian="Est. Diff in Means",
                      binomial="Est. Odds Ratio",
                      poisson="Est. Inc. Rate Ratio" )
  cat("\\begin{tabular}{l", rep("c",3), "}", "\n", 
      "\\hline", "\n", file=fileout )
  cat( paste("&", est_type, "&& \\\\ \n"), file=fileout, append=TRUE )
  cat( paste( "\\multicolumn{1}{c}{Covariate} & (95\\% CI) && p-Value \\\\ \n"), file=fileout, append=TRUE )
  cat( "\\hline \n", file=fileout, append=TRUE)
  write.table( paste(rownames(rslt), " & ", 
                     format(rslt[,1],digits=2,nsmall=2),
                     " (", format(rslt[,2],digits=2,nsmall=2), 
                     ", ", format(rslt[,3],digits=2,nsmall=2), ")",
                     " && ", ifelse(rslt[,5]<.0001,"$<$.0001", format(rslt[,5],digits=3,nsmall=3)), "\\\\", sep=""),
               file=fileout, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE )
  cat( "\\hline", "\n", "\\end{tabular}", file=fileout, append=TRUE )
  }



##
##### Example call using CHS data
##
chs <- read.table( paste(path,"Data/chsData.txt",sep="") )
chs$sbp140 <- ifelse( chs$sbp >140, 1, 0 )
chs$sbp140[ is.na(chs$sbp) ] <- NA
chs$pkyrsgrp <- cut( chs$pkyrs, c(-0.1,0,10,205) )
chs$weightgrp <- cut( chs$weight, c(0, 140, 170, 325 ) )
chs$incomegrp <- cut( chs$income, c(0,2,5,7,8) )

fit2a <- glm( sbp140 ~ age + gender + factor( pkyrsgrp ) + factor( exint0 ) + factor( weightgrp ), data=chs, family="binomial" )  
glm_est_latex(fit2a, robust=FALSE, fileout=paste(path,"LatexTables/glmEst_Ex.tex",sep=""))





