# This file contains statistical functions and a summary function
# For use in econometrics (econ104) at the University of Pennsylvania
# Created by Roman Schuster on 3/13/16

# Required dependencies:
library(lmtest) # for bptest() function
library(stats) # for logLik() function

# Creating a function that takes in a linear model and
# returns the sum of squared residuals from that model
ssr <- function(lm) {
  # Very straightforward, literally summing the squared residuals from the
  # 'residuals' attribute of the linear model
  output <- sum(lm$residuals^2)
  
  return(output)
}

# Creating a function that takes in a data set and a linear model 
# of the data and returns the standard error of regression
ser <- function(data, lm) {
  # Obtaining the sum of squared residuals from the linear model using ssr()
  ssrFromData <- ssr(lm)
  # Grabbing the sample size from the data using nrow()
  sampleSize <- nrow(data)
  # Grabbing the number of regression coefficients with length()
  numRegressors <- length(lm$coefficients)
  
  # Following the formula SER = sqrt(SSR/(N - K))
  return(sqrt(ssrFromData/(sampleSize - numRegressors)))
}

# Creating a function that takes in a data set and a linear model 
# of the data and returns the Schwarz Information Criterion:
sic <- function(data, lm) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  # Grabbing the number of regression coefficients with length()
  numRegressors <- length(lm$coefficients)
  
  # Following the formula:
  # SIC = -2*L/T + (K*Ln(T))/T
  lhs <- (((-2)*logLikFromLm)/sampleSize)
  rhs <- ((numRegressors*log(sampleSize))/sampleSize)
  output <- (lhs + rhs)
  
  # This output returns the SIC along with the number of degrees of freedom used,
  # but will only show up as the numerical SIC in the preceding summary function
  return(output)
}

# Creating a function that takes in a data set and a linear  
# model of the data and returns the Akaike Information Criterion:
aic <- function(data, lm) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  # Grabbing the number of regression coefficients with length()
  numRegressors <- length(lm$coefficients)
  
  # Following the formula:
  # AIC = -2*L/T + 2*K/T
  lhs <- (((-2)*logLikFromLm)/sampleSize)
  rhs <- ((2*numRegressors)/sampleSize)
  output <- (lhs + rhs)
  
  # This output returns the AIC along with the number of degrees of freedom used,
  # but will only show up as the numerical AIC in the preceding summary function
  return(output)
}

# Creating a function that takes in a data set and a linear 
# model of the data and returns the Hannan-Quinn Criterion
hq <- function(data, lm) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  # Grabbing the number of regression coefficients with length()
  numRegressors <- length(lm$coefficients)
  
  # Following the formula:
  # HQ = -2*(L/T) + 2*K*Ln(Ln(T))/T
  lhs <- ((-2)*(logLikFromLm/sampleSize))
  rhs = (2*(numRegressors)*log(log(sampleSize)))/sampleSize
  output <- (lhs + rhs)
  
  # This output returns the HQ along with the number of degrees of freedom used,
  # but will only show up as the numerical HQ in the preceding summary function
  return(output)
}

# Creating a function that takes in a data set and a linear 
# model of the data set and returns the Durbin-Watson Statistic
dw <- function(data, lm) {
  # Grabbing the sum of squared residuals from the linear model using ssr()
  ssrFromData <- ssr(lm)
  # Grabbing the sample size from the data using nrow()
  sampleSize <- nrow(data)
  
  # Declaring a variable named numerator that will hold the sum of squared
  # differences between adjacent disturbances
  numerator <- 0
  
  # Looping through each residual in the linear model and adding the squared
  # difference between adjacent residuals to the numerator for each iteration
  for (i in 1:(sampleSize - 1)) {
    numerator = numerator + ((lm$residuals[i + 1] - lm$residuals[i])^2)
  }
  
  # Returning the sum of squared differences divided by the sum of squared residuals
  output <- (numerator/ssrFromData)
  return(output)
}

# Creating a function that takes in a data set and 
# a linear model of the data set and returns
# a comprehensive statistical analysis
statSummary <- function(data, lm) {
  # Printing out the coefficient summary, includes:
  # Coefficient estimate
  # Std. error of coefficient estimate
  # T-Statistic of coefficient estimate
  # P-Value of T-Statistic for each coefficient estimate
  # MARK: Uses the summary function from the stats package
  print(summary(lm)$coefficients)
  
  # Grabbing the number of regression coefficients using length()
  numRegressors <- length(lm$coefficients)
  
  # Calling functions from this file, as well as the stats
  # and lmtest packages, to obtain relevant statistical data
  # and put them into variables to print out at the end
  
  # Sum of Squared Residuals
  ssrFromData <- ssr(lm)
  # Log Likelihood
  logLikFromData <- logLik(lm)
  
  # Akaike Information Criterion
  aicFromData <- aic(data, lm)
  # Schwarz Information Criterion
  sicFromData <- sic(data, lm)
  # Hannan-Quinn Criterion
  hqFromData <- hq(data, lm)
  
  # Durbin-Watson Statistic
  dwFromData <- dw(data, lm)
  
  # R Squared
  rSquaredFromData <- summary(lm)$r.squared
  # Adjusted R Squared
  adjRSquaredFromData <- summary(lm)$adj.r.squared
  
  # Std. Error of Regression
  serFromData <- ser(data, lm)
  # F-Statistic of the regression
  fStatArrayFromData <- summary(lm)$fstatistic
  # Grabbing just the F-Stat from the data to avoid printing
  # the number of degrees of freedom along with it
  fStatFromData <- fStatArrayFromData[1]
  # Calculating the p-value for the f-statistic
  pValueFStat <- 1 - pf(fStatFromData, fStatArrayFromData[2], fStatArrayFromData[3])
  
  # Breusch-Pagan-Godfrey Test for Heteroskedasticity
  breuschPaganGodfrey <- bptest(lm, data = data)
  # Grabbing just the BGP stat (not degrees of freedom)
  bpStat <- breuschPaganGodfrey$statistic
  # Grabbing the p-value of the BGP
  bpPValue <- breuschPaganGodfrey$p.value
  
  # Printing out the data in a clean, easy to read format
  cat('\nR Squared: ',substitute(rSquaredFromData),'\n')
  cat('Adjusted R Squared: ', substitute(adjRSquaredFromData), '\n')
  
  cat('\nStd. Error of Regression: ', substitute(serFromData), '\n')
  cat('Sum of Squared Residuals: ',substitute(ssrFromData),'\n')
  cat('Log Likelihood: ',substitute(logLikFromData),'\n')
  
  cat('\nF-Statistic: ', substitute(fStatFromData), '\n')
  cat('p-Value of F-Statistic: ', substitute(pValueFStat), '\n \n')
  
  cat('Akaike Info. Criterion:  ',substitute(aicFromData),'\n')
  cat('Schwarz Info. Criterion:  ',substitute(sicFromData),'\n')
  cat('Hannan-Quinn Criterion: ', substitute(hqFromData), '\n')
  cat('Durbin-Watson Statistic: ', substitute(dwFromData), '\n')
  
  cat('\nBreusch-Pagan-Godfrey Test [N*(R^2)]: ', substitute(bpStat), '\n')
  cat('BPG Test [p-value]: ', substitute(bpPValue), '\n')
}





