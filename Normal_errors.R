library(doParallel)
library(MASS)  
library(DEoptim)

Cbeta=0.95
Cgamma=0.90

n=2000
#n=200,20

a=0
b=10

p=3

nu=n-p # degrees of freedom for chi-square distribution

#Set parameters
set.seed(20250716) 
alp.t <- c(1, 2, -1)  
sig.t <- 1          

#Construct design matrix X
X1 <- seq(from = a, to = b, length.out = n)
X=matrix(0,nrow = n,ncol = p)
for (i in 1:length(X1)) {
  x=X1[i]
  vx <- 1 #the vector (1, x, ..., x^(p-1))'
  for (j in 2:p) {
    vx[j] <- x^(j-1)
  }
  X[i,]=vx
}
XTX <- t(X)%*%X
#Find the inverse of XtX by using Q-R decomposition
XTX.qr <- qr(XTX)
Q<-qr.Q(XTX.qr)
R<-qr.R(XTX.qr)
invXTX <- qr.solve(R)%*%t(Q)


#Generate response variable y
eps <- rnorm(n, mean = 0, sd = sig.t)
y <- as.vector(X %*% alp.t + eps)
yy <- matrix(y, ncol = 1)

beta.hat <- invXTX %*% t(X) %*% yy
sigma.hat <- sqrt(sum((yy-X%*%beta.hat)^2)/nu)
#Store training dataset
training_data <- data.frame(x = X1, y = y)
head(training_data)

two.sti.exact=function(lam,n,p,Cbeta,a,b,B,u,invXTX,cut.p)
{
  ##  lam -- the critical constant lambda for this simulation replication
  ##  n -- sample size
  ##  p=r+1 -- 1+order of the polynomial
  ##  Cbeta -- the content
  ##  [a,b] -- the the covariate interval of interest for x
  ##  B -- the column-vector of random Z 
  ##  u -- the r.v. u 
  ##  invXTX = inv(X'X) where X is the design matrix
  ##  cut.p -- the number of cut points
  zp = qnorm((1 + Cbeta)/2) # (1+Cbeta)/2 quantile of standard normal distribution
  
  fd = function(x, lam,invXTX){
    vx <- 1 #the vector (1, x, ..., x^(p-1))'
    for (j in 2:p) {
      vx[j] <- x^(j-1)
    }
    pnorm(t(vx)%*%B + lam*(zp + sqrt(p+2)*sqrt(t(vx)%*%invXTX%*%vx))*u) -
      pnorm(t(vx)%*%B - lam*(zp + sqrt(p+2)*sqrt(t(vx)%*%invXTX%*%vx))*u)
  }
  
  #Cut the interval [a,b] into several sub-intervals
  ab.seq=seq(a,b,length.out=cut.p)
  optim.p=0
  for (j in 1:(length(ab.seq)-1)) {
    #Find the optimisation in each sub-interval
    optim.ab=optimise(fd,c(ab.seq[j],ab.seq[j+1]),lam=lam,invXTX=invXTX)$objective-Cbeta
    #Calculate the probabilities at the end points of sub-interval
    optim.a1=fd(ab.seq[j], lam,invXTX)-Cbeta
    optim.b1=fd(ab.seq[j+1], lam,invXTX)-Cbeta
    #Check the minimum value
    optim.p[j]=min(optim.ab,optim.a1,optim.b1)
  }
  #Calculate the probabilities at the end points a and b
  optim.a=fd(a, lam,invXTX)-Cbeta
  optim.b=fd(b, lam,invXTX)-Cbeta
  
  min(optim.p,optim.a,optim.b)
  
}

#Case Resampling
############################
## Bootstrap: Total time  ##
############################
t_boot_total <- system.time({
  L <- 100000
  Z.star <- matrix(NA, nrow = L, ncol = p)  
  u.star <- numeric(L)                     
  
  set.seed(20250716)  
  for (b in 1:L) {
    idx <- sample(1:n, size = n, replace = TRUE)
    X.star <- X[idx, ]
    y.star <- y[idx]
    XtX.star <- t(X.star) %*% X.star
    XtX.star.qr <- qr(XtX.star)
    Q <- qr.Q(XtX.star.qr)
    R <- qr.R(XtX.star.qr)
    invXtX.star <- qr.solve(R) %*% t(Q)
    beta.star <- invXtX.star %*% t(X.star) %*% y.star
    sigma.star <- sqrt(sum((y.star - X.star %*% beta.star)^2) / nu)
    Z.star[b, ] <- as.vector((beta.star - beta.hat) / sigma.hat)
    u.star[b] <- sigma.star / sigma.hat
  }
  
  B.sim <- Z.star
  u.sim <- u.star
  tol <- 1e-4
  
  numCores <- 24
  registerDoParallel(numCores)
  lam_boot <- foreach(i = 1:L, .combine = "c") %dopar%{
    B = B.sim[i,]
    u = u.sim[i]
    klow <- 0; kupp <- 5
    fupp <- two.sti.exact(kupp, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
    if (fupp <= 0) {
      klow <- kupp
      repeat {
        kupp <- kupp + 5
        fupp <- two.sti.exact(kupp, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
        if (fupp > 0) break
      }
    }
    while ((kupp - klow) > tol) {
      cmid <- (klow + kupp)/2
      c <- cmid
      fmid <- two.sti.exact(c, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
      if (fmid == 0) { rc = cmid; break
      } else if (fupp * fmid > 0) {
        kupp <- cmid; fupp <- fmid
      } else {
        klow <- cmid; flow <- fmid
      }
    }
    res <- cmid
  }
  stopImplicitCluster()
  
  lam_two_sti_boot <<- quantile(lam_boot, Cgamma)
})
cat("Bootstrap total elapsed (s):", t_boot_total["elapsed"], "\n")
cat("Bootstrap STI lambda (lam_two_sti_boot):", lam_two_sti_boot, "\n")

#########################
## Exact: Total time   ##
#########################
t_exact_total <- system.time({
  L <- 100000
  set.seed(10)
  # Generate the Z and u in (5)
  B.sim <- mvrnorm(L, c(rep(0,p)), invXTX)
  u.sim <- sqrt(rchisq(L, df = n - p) / (n - p))
  tol <- 1e-4
  
  numCores <- 24
  registerDoParallel(numCores)
  lam_exact <- foreach(i = 1:L, .combine = "c") %dopar%{
    B = B.sim[i,]
    u = u.sim[i]
    klow <- 0; kupp <- 5
    fupp <- two.sti.exact(kupp, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
    if (fupp <= 0) {
      klow <- kupp
      repeat {
        kupp <- kupp + 5
        fupp <- two.sti.exact(kupp, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
        if (fupp > 0) break
      }
    }
    while ((kupp - klow) > tol) {
      cmid <- (klow + kupp)/2
      c <- cmid
      fmid <- two.sti.exact(c, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
      if (fmid == 0) { rc = cmid; break
      } else if (fupp * fmid > 0) {
        kupp <- cmid; fupp <- fmid
      } else {
        klow <- cmid; flow <- fmid
      }
    }
    res <- cmid
  }
  stopImplicitCluster()
  
  lam_two_sti_exact <<- quantile(lam_exact, Cgamma)
})
cat("Exact total elapsed (s):", t_exact_total["elapsed"], "\n")
cat("Exact STI lambda (lam_two_sti_exact):", lam_two_sti_exact, "\n")


r_ratio=(lam_two_sti_boot/lam_two_sti_exact-1)*100
cat("Relative difference (r_ratio %):", r_ratio, "%\n")

zbeta = qnorm((1+Cbeta)/2) 

nsim1=10000

################################################
## Coverage A：Bootstrap Coverage Total time  ##
################################################
numCores <- 24
registerDoParallel(numCores)
t_cov_boot <- system.time({
  cover_B_vec <- foreach(j = 1:nsim1, .combine = "c") %dopar%{
    set.seed(j)
    eps <- rnorm(n, 0, sig.t)
    y <- as.vector(X %*% alp.t) + eps
    yy <- matrix(y, ncol = 1)
    alp.h <- invXTX %*% t(X) %*% yy
    sig.h <- sqrt(sum((yy - X %*% alp.h)^2) / nu)
    
    min.cov <- function(x.new, invXTX, lambda){
      vx <- 1
      for (jj in 2:p) { vx[jj] <- x.new^(jj-1) }
      d2.x <- as.numeric(t(vx) %*% invXTX %*% vx)
      L_1 <- as.numeric(t(vx) %*% alp.h) - lambda * (zbeta + sqrt((p+2) * d2.x)) * sig.h
      U_1 <- as.numeric(t(vx) %*% alp.h) + lambda * (zbeta + sqrt((p+2) * d2.x)) * sig.h
      y.mean <- as.numeric(t(vx) %*% alp.t)
      pnorm(U_1, y.mean, sig.t) - pnorm(L_1, y.mean, sig.t)
    }
    
    cut.p <- seq(a, b, length.out = 10)
    covp1 <- 0; covb1 <- 0
    for (ip in 1:(length(cut.p)-1)) {
      a1 <- cut.p[ip]; b1 <- cut.p[ip+1]
      covp1[ip] <- optimise(min.cov, c(a1, b1), invXTX = invXTX, lambda = lam_two_sti_boot)$objective
      covb1[ip] <- min.cov(a1, invXTX, lam_two_sti_boot)
    }
    covu1 <- min.cov(b, invXTX, lam_two_sti_boot)
    covmin1 <- min(c(covp1, covb1, covu1))
    covmin1
  }
})
stopImplicitCluster()
cat("Coverage elapsed (Bootstrap) (s):", t_cov_boot["elapsed"], "\n")

##############################################
## Coverage B：Exact Coverage Total time    ##
##############################################
numCores <- 24
registerDoParallel(numCores)
t_cov_exact <- system.time({
  cover_E_vec <- foreach(j = 1:nsim1, .combine = "c") %dopar%{
    set.seed(j)
    eps <- rnorm(n, 0, sig.t)
    y <- as.vector(X %*% alp.t) + eps
    yy <- matrix(y, ncol = 1)
    alp.h <- invXTX %*% t(X) %*% yy
    sig.h <- sqrt(sum((yy - X %*% alp.h)^2) / nu)
    
    min.cov <- function(x.new, invXTX, lambda){
      vx <- 1
      for (jj in 2:p) { vx[jj] <- x.new^(jj-1) }
      d2.x <- as.numeric(t(vx) %*% invXTX %*% vx)
      L_1 <- as.numeric(t(vx) %*% alp.h) - lambda * (zbeta + sqrt((p+2) * d2.x)) * sig.h
      U_1 <- as.numeric(t(vx) %*% alp.h) + lambda * (zbeta + sqrt((p+2) * d2.x)) * sig.h
      y.mean <- as.numeric(t(vx) %*% alp.t)
      pnorm(U_1, y.mean, sig.t) - pnorm(L_1, y.mean, sig.t)
    }
    
    cut.p <- seq(a, b, length.out = 10)
    covp2 <- 0; covb2 <- 0
    for (ip in 1:(length(cut.p)-1)) {
      a1 <- cut.p[ip]; b1 <- cut.p[ip+1]
      covp2[ip] <- optimise(min.cov, c(a1, b1), invXTX = invXTX, lambda = lam_two_sti_exact)$objective
      covb2[ip] <- min.cov(a1, invXTX, lam_two_sti_exact)
    }
    covu2 <- min.cov(b, invXTX, lam_two_sti_exact)
    covmin2 <- min(c(covp2, covb2, covu2))
    covmin2
  }
})
stopImplicitCluster()
cat("Coverage elapsed (Exact) (s):", t_cov_exact["elapsed"], "\n")

######################################
## Compute the coverage probability ##
######################################
con_B <- sum(cover_B_vec >= Cbeta)
con_E <- sum(cover_E_vec >= Cbeta)
cover_B <- con_B / nsim1 * 100
cover_E <- con_E / nsim1 * 100

lambda_table <- matrix(c(n,"ES",p,lam_two_sti_boot,lam_two_sti_exact,
                         (lam_two_sti_boot/lam_two_sti_exact-1)*100),
                       ncol=6, byrow=TRUE)
colnames(lambda_table) <- c('n','X','p','lambda_s^B','lambda_s','r(%)')
print(lambda_table)

coverage_table <- matrix(c(n,"ES",p,cover_B,cover_E), ncol=5, byrow=TRUE)
colnames(coverage_table) <- c('n','X','p','Bootstrap','EXACT')
print(coverage_table)


cat("Bootstrap total time (user, system, elapsed):\n")
print(t_boot_total)

cat("Exact total time (user, system, elapsed):\n")
print(t_exact_total)

cat("Bootstrap coverage time (user, system, elapsed):\n")
print(t_cov_boot)

cat("Exact coverage time (user, system, elapsed):\n")
print(t_cov_exact)