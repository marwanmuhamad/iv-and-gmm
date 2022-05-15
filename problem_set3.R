library(MASS)


# Question 1 --------------------------------------------------------------

# suppose we have T = 50 for the first two questions
T. = 50

# define the parameter for betas and phi:
beta0 = 1
beta = 1
phi0 = 1

# define variable xt:
x = rep(1, T.)

# Data generating process:
set.seed(333)

#define variables z1, z2, ut and vt with multivariate normally distributed:
mu = c(0, 0, 0, 0)
omega = matrix(c(1,0,0,0,
                 0,1,0,0,
                 0,0,1,0.8,
                 0,0,0.8,1), 
               nrow = 4, byrow = TRUE)

dgp = mvrnorm(n = T., mu = mu, Sigma = omega)
z1 = dgp[, 1]
z2 = dgp[, 2]
u = dgp[, 3]
v = dgp[, 4]

df = data.frame(x = x, z1 = z1, z2 = z2)
Z = as.matrix.data.frame(df)

phi1 = 5
phi2 = 10

x <- phi0 + z1*phi1 + z2*phi2 + v

# generate y variable:
y <- beta0 + x*beta + u

df$y <- y
df$x <- x

model_ols <- lm(y ~ x)
summary(model_ols)

# 1.  is xt exogenous in (1) and (2)? why? --------------------------------

# Implementing DWH test:
# use first step of 2SLS to test the exogeneity of xt variable
model1 <- lm(x ~ z1 + z2, data = df)
# extract the residual of the model in the first step and name it res1:

df$res1 <- model1$residuals

# step 2: regress y on both x and res:
model2 <- lm(y ~ x + res1, data = df)

# inspect the F-statistic or t-statistic:
summary(model2)

# F-stat = 8855 (p-value = 2.2e-16) and t-stat = 10.57 (p-value = 5.25e-14)
# both are significant, hence indicates x is endogenous


# 2. the validity of z1 and z2 instruments --------------------------------

# take the result from first step above:
summary(model1)
# F-stats = 3323 (>10), based on Staiger and Stock (1997) rule of thumb F>10 ==> strong IVs and valid

# test for over-identifying restriction:

# chisq-table for alpha = 0.05, T = 50, and df = 1 = 3.84


# question no. 1 point 3 - 6 -----------------------------------------------

# define the function:
generate_model <- function(phi1 = 5, phi2 = 10, T. = 50) {
  #the function takes default arguments for phi1, phi2 and T are 5, 10, and 50, respectively
  if(!require("MASS")) install.packages("MASS")

  beta0 <-  1
  beta <-  1
  phi0 <-  1
  
  # define the maximum repetition R:
  R <<- 1000
  #define variable x1t:
  x <-  rep(1, T.)
  
  mu <-  c(0, 0, 0, 0)
  omega <-  matrix(c(1,0,0,0,
                   0,1,0,0,
                   0,0,1,0.8,
                   0,0,0.8,1), 
                 nrow = 4, byrow = TRUE)
  
  statistic_table = data.frame(beta_ols = rep(0, R), tstat_ols = rep(0, R), 
                               beta_2sls = rep(0, R), tstat_2sls = rep(0, R), 
                               j_stat = rep(0, R))
  
  slope_ols = rep(0, R) # slope of OLS
  std_error_ols = rep(0, R) # standard error of OLS
  t_statistic_ols = rep(0, R) # t-statistic of OLS
  
  slope_2sls = rep(0, R) # slope of 2SLS
  std_error_2sls = rep(0, R) # standard error of 2SLS
  t_statistic_2sls = rep(0, R) # t-statistic of 2SLS
  
  Jstats = rep(0, R)
  # Data generating process:
  set.seed(333)
  for (i in 1:R) {
    dgp <- mvrnorm(n = T., mu = mu, Sigma = omega)
    z1 <-  dgp[, 1]
    z2 <-  dgp[, 2]
    u <-  dgp[, 3]
    v <-  dgp[, 4]
    
    df = data.frame(x = x, z1 = z1, z2 = z2)
    Z = as.matrix.data.frame(df)
    
    #re-generate x variable:
    x <- phi0 + z1*phi1 + z2*phi2 + v
    df$x <- x
    
    
    # generate y variable:
    y <- beta0 + x*beta + u
    df$y <- y
    
    # call OLS model:
    model_ols <- lm(y ~ x, data = df)
    
    slope_ols[[i]] = model_ols$coefficients[2] # extract the estimated parameter beta1
    std_error_ols[[i]] = summary(model_ols)$coefficients[1,2] # extract standard error of the model
    t_statistic_ols[[i]] = (slope_ols[[i]] - beta)/std_error_ols[[i]] # compute t-statistics
    
    # call 2SLS model:
    # step 1:
    model_step1 <- lm(x ~ z1 + z2, data = df)
    # estimate x_het:
    df$x_het <- model_step1$fitted.values
    
    # step 2:
    model_step2 <- lm(y ~ x_het, data = df)
    
    slope_2sls[[i]] = model_step2$coefficients[1] # extract the estimated parameter beta1
    std_error_2sls[[i]] = summary(model_step2)$coefficients[1,2] # extract standard error of the model
    t_statistic_2sls[[i]] = (slope_2sls[[i]] - beta)/std_error_2sls[[i]] # compute t-statistics
    
    # Compute Sargan's (1960) J-stats:
    
    WT = solve((1/T.)*t(Z)%*%Z) # Sargan's weighting matrix
    
    # estimate the objective function:
    model <- lm(y ~ x, data = df)
    res <- model$residuals # extract the residuals
    
    # (yi - xi'beta)Zi:
    
    m = matrix(((1/T.)%*%(res)%*%Z), byrow = FALSE)
    
    QT = t(m)%*%WT%*%m
    
    # Sargan's test (1960) for over-identifying restriction:
    
    Jstats[[i]] = T.*QT[1,1]
    
    statistic_table[i, c(1,2,3,4,5)] <- c(slope_ols[[i]], t_statistic_ols[[i]], 
                                          slope_2sls[[i]], t_statistic_2sls[[i]],
                                          Jstats[[i]])
  }
  return(statistic_table) #return data.frame
}

# Plot the histogram of each result: --------------------------------------
# function to create histogram of estimated beta 1 and t-statistics of OLS Model:
plot_betas <- function(data = NULL, T. = 50, col = "lightblue") {
  # control input to only data.frame type
  if (!is.data.frame(data)) {
    stop(paste("The argument of `data` should be of `list` class, not", class(data)), call. = FALSE)
  }
  
  par(mfcol = c(1,2))
  
  hist(data$beta_ols, main = paste("Histogram of Estimated Beta (OLS) with T = ", T.), 
       xlab = paste("Estimated Beta"), ylab = "Frequency", col = col)
  hist(data$beta_2sls, main = paste("Histogram of Estimated Beta (IV) with T = ", T.), 
       xlab = paste("Estimated Beta"), ylab = "Frequency", col = col)

}

# function to create histogram of estimated beta 1 and t-statistics of IV Model:
plot_tsat_jstat <- function(data = NULL, T. = 50, col = "#FAAB56") {
  # control input to only data.frame type
  if (!is.data.frame(data)) {
    stop(paste("The argument of `data` should be of `list` class, not", class(data)), call. = FALSE)
  }
  
  par(mfcol = c(1,3))
  
  hist(data$tstat_ols, main = paste("Histogram of t-statistics (OLS) with T = ", T.), 
       xlab = paste("t-statistics"), ylab = "Frequency", col = col)
  hist(data$tstat_2sls, main = paste("Histogram of t-statistics (IV) with T = ", T.), 
       xlab = paste("t-statistics"), ylab = "Frequency", col = col)
  hist(data$j_stat, main = paste("Histogram of J-statistic with T = ", T.), 
       xlab = paste("Estimated Beta"), ylab = "Frequency", col = col)
}

bias <- function(data = NULL) {
  # control input to only data.frame type
  if (!is.data.frame(data)) {
    stop(paste("The argument of `data` should be of `list` class, not", class(data)), call. = FALSE)
  }
  # (b) Compute the estimated bias:
  
  bias_ols <- (1/R)*sum(data$beta_ols - beta)
  bias_iv <- (1/R)*sum(data$beta_2sls - beta)
  
  return(data.frame(bias_ols = bias_ols, bias_iv = bias_iv))
}

rej_freq <- function(data = NULL) {
  # (c) Compute rejection frequency at 5% level
  
  data["rej_ols"] <- ifelse(data$tstat_ols>1.96, 1, 0)
  data["rej_iv"] <- ifelse(data$tstat_2sls>1.96, 1, 0)
  data["rej_js"] <- ifelse(data$tstat_2sls>3.84, 1, 0)
  
  rej_ols <- (1/R)*sum(data$rej_ols)
  rej_iv <- (1/R)*sum(data$rej_iv)
  rej_js <- (1/R)*sum(data$rej_js)
  
  return(data.frame(rej_ols = rej_ols, rej_iv = rej_iv, rej_js = rej_js))
}

# call functions ----------------------------------------------------------


data1 <- generate_model()
debug(generate_model) # performs debugging to trace error
undebug(generate_model) # undo debugging

plot_betas(data = data1)
plot_tsat_jstat(data = data1)
rej_freq(data = data1)

# Question 2 --------------------------------------------------------------


