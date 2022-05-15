# Question 2 --------------------------------------------------------------
xdata1 <- as.matrix(read.table("cbapmewrdata.dat"))
xdata2 <- as.matrix(read.table("cbapmewrinstr.dat"))

x <- cbind(xdata1, xdata2)
#--------------------------------------------------------------------
# Numeric gradient function
# Computes numerical gradient matrix at each observation
#--------------------------------------------------------------------
numgrad <- function( fun,x, ... ) {
  f0  <-  fun(x, ... )             # n by 1    
  n   <- length( f0 )
  k   <- length( x )
  fdf <- array(0, dim=c(n,k) )
  
  # Compute step size 
  eps <- .Machine$double.eps
  dx  <- sqrt( eps )*( abs( x ) + eps )
  xh  <- x + dx
  dx  <- xh - x
  ind <- dx < sqrt(eps)
  dx[ind] <- sqrt(eps)
  
  # Compute gradient
  xdx <- diag(dx) + x  
  for (i in seq(k)) {
    fdf[,i] <- fun(xdx[,i], ...)    
  }   
  G0 <- kronecker(matrix(1, 1, k), f0 )                       # n by k 
  G1 <- kronecker(matrix(1, n, 1), t(dx) )  
  G  <- ( fdf-G0 ) / G1
  return(G)  
}

# Define the moment equations  --------------------------------------------

# cratio = x1t
# r = x2t
# zt = matrix of instrumental variables

momentEq <- function(theta,cratio,r,zt) {
  cols <- ncol(zt)
  rows <- nrow(zt)
  
  #theta[1] = Relative risk aversion parameter    
  #theta[2] = Discount parameter                  
  
  # First Order Condition based on the model
  ut <- theta[1]*cratio^(-theta[2]) * (1 + r) - 1
  
  k <- 1        
  mt <- array(0, c(rows,cols))  # create matrix with element of zeros
  for (j in seq(cols)) {
    mt[,k] <- zt[,j]*ut
    k <- k+1    
  }  
  return(mt)  
}

# Defines the mean of the moment conditions  ------------------------------

meanEq <- function(theta,cratio,r,zt) {
  ret <- colMeans(momentEq(theta,cratio,r,zt))
  return(ret)
}

# GMM objective function which also computes the optimal w  ---------------
q_twoStep <- function(theta,cratio,r,zt,lmax){
  
  d <- momentEq(theta,cratio,r,zt)
  g <- cbind(colMeans(d))
  t <- length(cratio)
  w  <- diag(5) # initialize the weighting matrix in the first step as an identity matrix 5x5
   
  w <- w/t
  ret <- t(g) %*% solve(w) %*% g
  
  return(ret)  
}

# gmm two step ------------------------------------------------------------

gmm_twoStep <- function(x) {
  options(warn = -1)
  cratio = x[, 1]
  r = x[, 2]
  zt = x[, 3:7]
  
  # step 1 --> optimization
  ret <- optim(par = c(0,0), fn = q_twoStep, cratio =cratio, r = r, zt = zt)
  pars <- ret$par
  
  # step 2 --> optimization
  ret2 <- optim(par = pars, fn = q_twoStep, cratio = cratio, r = r, zt = zt)
  pars2 <- ret2$par
  qmin <- ret2$value
  
  # Numerical gradient at optimum
  dg <- numgrad(meanEq, pars2, cratio, r, zt)

  # Compute optimal w
  d    <- momentEq(pars2, cratio, r, zt)
  g    <- cbind(colMeans(d))
  w    <- t(d) %*% d

  v <- t(dg) %*% solve(w) %*% dg
  t <- length(cratio)
  
  # Hansen Sargan J Test
  j_stat <- t*q_twoStep(theta = pars2, cratio = cratio, r = r,zt = zt,lmax = 0)
  
  cat('\nTwo Step GMM')
  cat('\nThe value of the objective function (q)  = ', qmin )
  cat('\nJ-test is                                = ', j_stat, '\n' )
  se <- sqrt(diag(solve(v)/t))
  final_result <- cbind(Estimates=pars2, se=se, t=pars2/se)
  rownames(final_result) <- c("gamma", "delta")
  print(final_result)
  
  return(list(q = qmin, j_stat = j_stat, Estimates=pars2, se=se, t=pars2/se))
  
}

res_gmm2 <- gmm_twoStep(x)
# debug(gmm_twoStep)

q_iter <- function(theta,cratio,r,zt,lmax){
  # lmax = stopping criteria
  d <- momentEq(theta,cratio,r,zt)
  g <- cbind(colMeans(d))
  
  w  <- t(d) %*% d
  k <- 1 # initialize the repetition
  while (k <= lmax) {
    w_k <- t( d[(k+1):nrow(d),] ) %*% d[1:(nrow(d)-k),]
    w    <- w + (1.0-k/(lmax+1))*(w_k + t(w_k))
    k  <- k + 1    
  }  
  t <- length(cratio)
  w <- w/t
  ret <- t(g) %*% solve(w) %*% g
  
  return(ret)  
}

# Itterative GMM ----------------------------------------------------------

gmm_iter <- function(theta = res_gmm2$Estimates, x = NULL) {
  options(warn = -1)
  
  cratio    <- x[,1]      # Ratio of real consumption at t+1 and t   
  r <-  x[, 2] # gain at period t+1
  t <- length(cratio) # number of observations
  zt <- x[, 3:7] # Instruments
  
  res <- optim(theta, q_iter, cratio = cratio, r = r, zt=zt, lmax=0, method="BFGS")
  bgmm <- res$par  
  qmin <- res$value
  
  # Numerical gradient at optimum
  dg <- numgrad(meanEq,bgmm, cratio, r, zt)
  
  # Compute optimal w
  d    <- momentEq(bgmm, cratio,r, zt)
  g    <- cbind(colMeans(d))
  w    <- t(d) %*% d
  k  <- 1 # initialize iteration
  lmax <- 0 # stopping critea for the iteration
  while (k <= lmax) {
    w_k <- t( d[(k+1):nrow(d),] ) %*% d[1:(nrow(d)-k),]
    w    <- w + (1.0-k/(lmax+1))*(w_k + t(w_k))
    k  <- k + 1    
  }
  w <- w/t
  v <- t(dg) %*% solve(w) %*% dg
  
  # Hansen Sargan J Test
  
  j_stat <- t*q(bgmm,cratio,r,zt,0)
  
  cat('\nIterative GMM')
  cat('\nThe value of the objective function (q)  = ', qmin )
  cat('\nJ-test is                                = ', j_stat, '\n' )  
  se <- sqrt(diag(solve(v)/t))
  final_result <- cbind(Estimates=bgmm, se=se, t=bgmm/se)
  rownames(final_result) <- c("gamma", "delta")
  print(final_result)
  return(list(q = qmin, j_stat = j_stat, Estimates=bgmm, se=se, t=bgmm/se))
}


gmm_iter_res <- gmm_iter(x = x)
debug(gmm_iter)
