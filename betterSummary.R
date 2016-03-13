# This file contains statistical functions for use in econometrics (econ104)
# Created by Roman Schuster on 3/13/16

# Creating a function that takes in a linear model and returns the sum of squared
# residuals from that model
ssr <- function(lm) {
  # Very straightforward, literally summing the squared residuals from the
  # 'residuals' attribute of the linear model
  output <- sum(lm$residuals^2)
  
  return(output)
}

# Creating a function that takes in a data set, linear model of the data and number
# of regressors and returns the standard error of regression
ser <- function(lm, data, numRegressors) {
  # Obtaining the sum of squared residuals from the linear model using ssr()
  ssrFromData <- ssr(lm)
  # Grabbing the sample size from the data using nrow()
  sampleSize <- nrow(data)
  
  # Following the formula SER = sqrt(SSR/(N - K))
  return(sqrt(ssrFromData/(sampleSize - numRegressors)))
}

# R already has the logLik function that returns the log liklihood for a given
# linear model, so I'll use that to make functions for the SIC and AIC

# Creating a function that takes in the data, linear model of the data and number
# of regressors (including intercept) and returns the Schwarz Information Criterion:
sic <- function(data, lm, numRegressors) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  
  # Following the formula:
  # SIC = -2*L/T + (K*Ln(T))/T
  lhs <- (((-2)*logLikFromLm)/sampleSize)
  rhs <- ((numRegressors*log(sampleSize))/sampleSize)
  output <- (lhs + rhs)
  
  # This output returns the SIC along with the number of degrees of freedom used,
  # but will only show up as the numerical SIC in the preceding summary function
  return(output)
}

# Creating a function that takes in the data, linear model of the data and number
# of regressors (including intercept) and returns the Akaike Information Criterion
aic <- function(data, lm, numRegressors) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  
  # Following the formula:
  # AIC = -2*L/T + 2*K/T
  lhs <- (((-2)*logLikFromLm)/sampleSize)
  rhs <- ((2*numRegressors)/sampleSize)
  output <- (lhs + rhs)
  
  # This output returns the AIC along with the number of degrees of freedom used,
  # but will only show up as the numerical AIC in the preceding summary function
  return(output)
}

# Creating a function that takes in the data, a lienar model of the data, and the
# number of regressors (including intercept) and returns the Hannan-Quinn Criterion
hq <- function(data, lm, numRegressors) {
  # Using the logLik() function to get the log likelihood of the linear model
  logLikFromLm <- logLik(lm)
  # Grabbing the sample size from the data using the nrow() function
  sampleSize <- nrow(data)
  
  # Following the formula:
  # HQ = -2*(L/T) + 2*K*Ln(Ln(T))/T
  lhs <- ((-2)*(logLikFromLm/sampleSize))
  rhs = (2*(numRegressors)*log(log(sampleSize)))/sampleSize
  output <- (lhs + rhs)
  
  # This output returns the HQ along with the number of degrees of freedom used,
  # but will only show up as the numerical HQ in the preceding summary function
  return(output)
}

# Creating a function that takes in a data set and the linear model of the data set
# and returns the Durbin-Watson Statistic
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

# Rolling a summary function that takes in a data set, a linear approximation of
 # the data set, and the number of regressors (including intercept) used in the
 # linear model and returns relevant statistical data
betterSummary <- function(data, lm, numRegressors) {
  # Printing out the coefficient estimates, std. errors, t-values and p-values
  # (from the coefficients attribute of the summary function)
  print(summary(lm)$coefficients)
  
  # Calling on prior functions to obtain the relevant statistical data
  ssrFromData <- ssr(lm)
  logLikFromData <- logLik(lm)
  aicFromData <- aic(data, lm, numRegressors)
  sicFromData <- sic(data, lm, numRegressors)
  hqFromData <- hq(data, lm, numRegressors)
  dwFromData <- dw(data, lm)
  rSquaredFromData <- summary(lm)$r.squared
  adjRSquaredFromData <- summary(lm)$adj.r.squared
  serFromData <- ser(lm, data, numRegressors)
  fStatArrayFromData <- summary(lm)$fstatistic
  fStatFromData <- fStatArrayFromData[1]
  pValueFStat <- 1 - pf(fStatFromData, fStatArrayFromData[2], fStatArrayFromData[3])

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
}