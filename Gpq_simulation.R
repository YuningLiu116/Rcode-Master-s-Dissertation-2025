library(doParallel)
library(MASS)  
library(DEoptim)
library(nlme)

Cbeta=0.95
Cgamma=0.90

beta_true <- c(10, 0.8)   
sigma_d2  <- 1.2^2       
sigma_e2  <- 0.6^2
x_min <- 0
x_max <- 10

m <- 40#Subject number
k <- 5#Number of observations per Subject
n <- m * k
p=2
nu=n-p
simulate_data <- function(x_min,x_max, m, k,
                          beta_true, sigma_d2, sigma_e2)
{
  n <- m * k                              
  subject <- rep(seq_len(m), each = k) 
  x <- rep(seq(x_min,x_max, length.out = k), times = m) 
  b_i  <- rnorm(m, 0, sqrt(sigma_d2))[subject]  
  eps  <- rnorm(n, 0, sqrt(sigma_e2)) 
  y    <- beta_true[1] + beta_true[2] * x + b_i + eps
  data.frame(y, x, subject = factor(subject))
}


set.seed(120)
dat <- simulate_data(x_min,x_max, m, k, beta_true, sigma_d2, sigma_e2)

plot(y ~ x, data = dat,
     pch = 16, 
     col = rgb(0, 0, 0, alpha = 0.5), 
     main = "Simulated Data",
     xlab = "x", ylab = "y",
     las = 1,
     frame.plot = FALSE)
grid()
for (s in levels(dat$subject)) {
  subject_data <- subset(dat, subject == s)
  lines(y ~ x, data = subject_data, col = rgb(0, 0, 0, alpha = 0.5))
}


#Ignore random intercept
fit_fix <- lm(y ~ x, dat)
beta_hat_fix <- coef(fit_fix)
sigma_e_hat_fix  <- summary(fit_fix)$sigma 
#Considering the random intercept
fit_mixed <- lme(y ~ x, random = ~1 | subject, data = dat)
beta_hat_mixed <- fixef(fit_mixed)
sigma_d_hat <- as.numeric(VarCorr(fit_mixed)["(Intercept)", "StdDev"])
sigma_e_hat <- fit_mixed$sigma
sigma_hat_x  <- sqrt(sigma_d_hat^2 + sigma_e_hat^2) 
#prepare work
X  <- model.matrix(~ x, dat)                 # n×2
Z  <- model.matrix(~ 0 + subject, dat)       # n×m
XZ <- cbind(X, Z)
y <- dat$y
N <- length(y)
r <- Matrix::rankMatrix(XZ)[]
b <- Matrix::rankMatrix(X)[]
Fmat  <- svd(XZ, nu = r, nv = 0)$u
I_r <- diag(r)
XtF  <- t(X) %*% Fmat
Fty  <- t(Fmat) %*% y
P_XZ <- Fmat %*% t(Fmat)
SS_e <- as.numeric(t(y) %*% (diag(N) - P_XZ) %*% y)
x_min <- min(dat$x)
x_max <- max(dat$x)
XTX <- t(X)%*%X
#Find the inverse of XtX by using Q-R decomposition
XTX.qr <- qr(XTX)
Q<-qr.Q(XTX.qr)
R<-qr.R(XTX.qr)
invXTX <- qr.solve(R)%*%t(Q)
#invXTX=solve(t(X)%*%X)


#prepare function
two.sti.GPQ=function(lam, U_e2, U_0_2, Zk,
                     Fmat, Z, SS_e, p,
                     XtF, Fty, I_r,
                     Cbeta, x_min, x_max, cut.p, 
                     beta_hat,sigma_hat_x)
{ 
  # Step 1: compute GPQ for sigma²_e
  G_sigma2_e <- SS_e / U_e2
  
  # Step 2: define function for root-finding to get G_sigma2_d
  find_Gsigma2_d <- function(G_guess) {
    VG <- G_guess * t(Fmat) %*% Z %*% t(Z) %*% Fmat + G_sigma2_e * I_r
    # Step 3: Find the inverse by using Q-R decomposition
    VG.qr <- qr(VG)
    Q<-qr.Q(VG.qr)
    R<-qr.R(VG.qr)
    VG_inv <- qr.solve(R)%*%t(Q)
    mid <- XtF %*% VG_inv %*% t(XtF)
    mid.qr <- qr(mid)
    Q1<-qr.Q(mid.qr)
    R1<-qr.R(mid.qr)
    mid_inv <- qr.solve(R1)%*%t(Q1)
    proj <- VG_inv - VG_inv %*% t(XtF) %*% mid_inv %*% XtF %*% VG_inv
    stat <- t(Fty) %*% proj %*% Fty
    
    as.numeric(stat - U_0_2)
  }
  
  # Step 4: solve for G_sigma2_d via uniroot
  G_sigma2_d <- uniroot(find_Gsigma2_d, lower = 1e-6, upper = 1e4, tol = 1e-6)$root
  
  # Step 5: compute VG
  VG <- G_sigma2_d * t(Fmat) %*% Z %*% t(Z) %*% Fmat + G_sigma2_e * I_r
  VG.qr <- qr(VG)
  Q2 <-qr.Q(VG.qr)
  R2 <-qr.R(VG.qr)
  VG_inv <- qr.solve(R2)%*%t(Q2)
  
  # Step 6: compute GPQ for x0'β
  A <- XtF %*% VG_inv %*% t(XtF)
  A.qr <- qr(A)
  Q3 <-qr.Q(A.qr)
  R3 <-qr.R(A.qr)
  A_inv <- qr.solve(R3)%*%t(Q3)
  mu_hat <- A_inv %*% XtF %*% VG_inv %*% Fty
  
  basis_fn <- function(x) as.numeric(x^(0:(p - 1)))  
  z_fn <- function(x) c(1, rep(0, p - 1))           
  
  # Step 7: define center and radius functions
  center_fn <- function(x) {
    xvec <- basis_fn(x)
    as.numeric(t(xvec) %*% mu_hat - Zk * sqrt(t(xvec) %*% A_inv %*% xvec))
  }
  
  radius_fn <- function(x) {
    xvec <- basis_fn(x)
    zvec <- z_fn(x)
    sqrt(as.numeric(t(zvec) %*% zvec * G_sigma2_d + G_sigma2_e + t(xvec) %*% A_inv %*% xvec))
  }
  
  fd <- function(x) {
    mu_hat_x    <- beta_hat[1] + beta_hat[2] * x
    
    cval <- center_fn(x)
    rval <- radius_fn(x)
    
    prob <-
      pnorm((cval + lam * rval - mu_hat_x) / sigma_hat_x) -
      pnorm((cval - lam * rval - mu_hat_x) / sigma_hat_x)
    
    prob
  }
  
  #Cut the interval [a,b] into several sub-intervals
  ab.seq=seq(x_min,x_max,length.out=cut.p)
  optim.p=0
  for (j in 1:(length(ab.seq)-1)) {
    #Find the optimisation in each sub-interval
    optim.ab=(optimise(fd,c(ab.seq[j],ab.seq[j+1]))$objective)-Cbeta
    #Calculate the probabilities at the end points of sub-interval
    optim.a1=fd(ab.seq[j])-Cbeta
    optim.b1=fd(ab.seq[j+1])-Cbeta
    #Check the minimum value
    optim.p[j]=min(optim.ab,optim.a1,optim.b1)
  }
  #Calculate the probabilities at the end points a and b
  optim.a=fd(x_min)-Cbeta
  optim.b=fd(x_max)-Cbeta
  min(optim.p,optim.a,optim.b)
}


two.sti.exact=function(lam,n,p,Cbeta,x_min,x_max,B,u,invXTX,cut.p)
{
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
  ab.seq=seq(x_min,x_max,length.out=cut.p)
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
  optim.a=fd(x_min, lam,invXTX)-Cbeta
  optim.b=fd(x_max, lam,invXTX)-Cbeta
  
  min(optim.p,optim.a,optim.b)
  
}

#GPQ
L=100000
set.seed(120)
U_e2  <- rchisq(L, df = N - r)
U_0_2 <- rchisq(L, df = r - b)
Zk    <- rnorm(L)

tol=0.0001
numCores = 36
registerDoParallel(numCores)
system.time({
  lam_GPQ= foreach(i = 1:L, .combine = "c",
                   .export = c("two.sti.GPQ","Fmat","Z","SS_e","p",
                               "XtF","Fty","I_r","Cbeta","x_min","x_max"))  %dopar%{
                                 this_U_e2 <- U_e2[i]
                                 this_U_02 <- U_0_2[i]
                                 this_Zk   <- Zk[i]
                                 #Find the lambda for this simulation repeat using the bisection method
                                 klow <- 0 
                                 kupp <- 5
                                 #In our numerical results, cut.p = 6 is enough to find the global optimisation.
                                 fupp <- two.sti.GPQ(kupp, this_U_e2, this_U_02, this_Zk,
                                                     Fmat, Z, SS_e, p,
                                                     XtF, Fty, I_r,
                                                     Cbeta, x_min, x_max, cut.p=10, 
                                                     beta_hat_mixed,sigma_hat_x)
                                 if (fupp <= 0) {
                                   klow <- kupp
                                   repeat {
                                     kupp <- kupp + 5
                                     fupp <- two.sti.GPQ(kupp, this_U_e2, this_U_02, this_Zk,
                                                         Fmat, Z, SS_e, p,
                                                         XtF, Fty, I_r,
                                                         Cbeta, x_min, x_max, cut.p = 10, 
                                                         beta_hat_mixed,sigma_hat_x)
                                     if (fupp > 0) break
                                   }
                                 }
                                 
                                 while ( (kupp-klow)>tol ){
                                   cmid <- (klow+kupp)/2
                                   c <- cmid
                                   fmid <- two.sti.GPQ(c, this_U_e2, this_U_02, this_Zk,
                                                       Fmat, Z, SS_e, p,
                                                       XtF, Fty, I_r,
                                                       Cbeta, x_min, x_max, cut.p=10, 
                                                       beta_hat_mixed,sigma_hat_x)
                                   if ( fmid == 0 ){
                                     rc = cmid
                                     break
                                   } else if (fupp*fmid >0){
                                     kupp <- cmid
                                     fupp <- fmid
                                   } else {
                                     klow <- cmid
                                     flow <- fmid
                                   }
                                 }
                                 res=cmid
                               }
})

stopImplicitCluster()

lam_two_sti_GPQ = quantile(lam_GPQ, Cgamma)
cat("GPQ STI lambda (lam_two_sti_GPQ):", lam_two_sti_GPQ, "\n")


#Set random seed to reproduce the results    
set.seed(120)
#Generate the Z and u in (5)
#B.sim=Z~N_p(0,invXTX)
#u.sim=sqrt(chi_{n-p}^2/(n-p))
B.sim = mvrnorm(L,c(rep(0,p)),invXTX)
u.sim = sqrt(rchisq(L,df=n-p)/(n-p))
tol=0.0001

numCores = 36
registerDoParallel(numCores)
system.time({
  lam_exact= foreach(i = 1:L,.combine = "c") %dopar%{
    B = B.sim[i,]
    u = u.sim[i]
    #Find the lambda for this simulation repeat using the bisection method
    #We use the initial interval (0, 5) as the (klow, kupp) that contains lambda
    klow <- 0 
    kupp <- 5
    #In our numerical results, cut.p = 6 is enough to find the global optimisation.
    fupp <- two.sti.exact(kupp,n,p,Cbeta,x_min,x_max,B,u,invXTX,cut.p = 6)
    
    if (fupp <= 0) {
      klow <- kupp
      repeat {
        kupp <- kupp + 5
        fupp <- two.sti.exact(kupp,n,p,Cbeta,x_min,x_max,B,u,invXTX,cut.p = 6)
        if (fupp > 0) break
      }
    }
    
    while ( (kupp-klow)>tol ){
      cmid <- (klow+kupp)/2
      c <- cmid
      fmid <- two.sti.exact(c,n,p,Cbeta,x_min,x_max,B,u,invXTX,cut.p = 6)
      if ( fmid == 0 ){
        rc = cmid
        break
      } else if (fupp*fmid >0){
        kupp <- cmid
        fupp <- fmid
      } else {
        klow <- cmid
        flow <- fmid
      }
    }
    res=cmid
  }
})


stopImplicitCluster()

lam_two_sti_exact = quantile(lam_exact, Cgamma)
cat("Exact STI lambda (lam_two_sti_exact):", lam_two_sti_exact, "\n")


#lam_two_sti_GPQ = 2.191544
#lam_two_sti_exact = 1.000595


zbeta = qnorm((1+Cbeta)/2) 

nsim1=5000

numCores = 48
registerDoParallel(numCores)
system.time({
  cover= foreach(j = 1:nsim1,.combine = "rbind") %dopar%{
    set.seed(j)
    dat <- simulate_data(x_min, x_max, m, k, beta_true, sigma_d2, sigma_e2)
    x_vec <- dat$x
    y_vec <- dat$y
    X  <- model.matrix(~ x_vec)
    beta_hat_fix <- coef(lm(y_vec ~ x_vec))
    sigma_e_hat  <- summary(lm(y_vec ~ x_vec))$sigma
    XTX <- t(X)%*%X
    XTX.qr <- qr(XTX)
    Q<-qr.Q(XTX.qr)
    R<-qr.R(XTX.qr)
    invXTX <- qr.solve(R)%*%t(Q)
    
    min.cov=function(x.new,invXTX,lambda){
      vx = c(1, x.new)
      d2.x=as.numeric(t(vx)%*%invXTX%*%vx) # The value of x'(X'X)^(-1)x
      L_1=as.numeric(t(vx)%*%beta_hat_fix)-lambda*(zbeta+sqrt((p+2)*d2.x))*sigma_e_hat # The value of l(x) exactSTI
      U_1=as.numeric(t(vx)%*%beta_hat_fix)+lambda*(zbeta+sqrt((p+2)*d2.x))*sigma_e_hat # The value of u(x) exactSTI
      y.mean=as.numeric(t(vx)%*%beta_true) # The true mean of y_x-distribution
      p.yx=pnorm(U_1,y.mean,sqrt(sigma_d2 + sigma_e2))-pnorm(L_1,y.mean,sqrt(sigma_d2 + sigma_e2))
      p.yx
    }
    
    fit_mixed <- lme(y ~ x, random = ~1 | subject, data = dat)
    beta_hat_mixed <- fixef(fit_mixed)
    sigma_d_hat <- as.numeric(VarCorr(fit_mixed)["(Intercept)", "StdDev"])
    sigma_e_hat <- fit_mixed$sigma
    sigma_hat_x  <- sqrt(sigma_d_hat^2 + sigma_e_hat^2)  
    X  <- model.matrix(~ x, dat)                 # n×2
    Z  <- model.matrix(~ 0 + subject, dat)       # n×m
    XZ <- cbind(X, Z)
    y <- dat$y
    N <- length(y)
    r <- Matrix::rankMatrix(XZ)[]
    b <- Matrix::rankMatrix(X)[]
    Fmat  <- svd(XZ, nu = r, nv = 0)$u
    I_r <- diag(r)
    XtF  <- t(X) %*% Fmat
    Fty  <- t(Fmat) %*% y
    P_XZ <- Fmat %*% t(Fmat)
    SS_e <- as.numeric(t(y) %*% (diag(N) - P_XZ) %*% y)
    x_min <- min(dat$x)
    x_max <- max(dat$x)
    
    L=1000
    set.seed(1+j)
    U_e2  <- rchisq(L, df = N - r)
    U_0_2 <- rchisq(L, df = r - b)
    Zk    <- rnorm(L)
    
    construct_center_radius <- function(U_e2, U_0_2, Zk,
                                        Fmat, Z, SS_e, p,
                                        XtF, Fty, I_r) {
      # Step 1: compute GPQ for sigma²_e
      G_sigma2_e <- SS_e / U_e2
      
      # Step 2: define function for root-finding to get G_sigma2_d
      find_Gsigma2_d <- function(G_guess) {
        VG <- G_guess * t(Fmat) %*% Z %*% t(Z) %*% Fmat + G_sigma2_e * I_r
        # Step 3: Find the inverse by using Q-R decomposition
        VG.qr <- qr(VG)
        Q<-qr.Q(VG.qr)
        R<-qr.R(VG.qr)
        VG_inv <- qr.solve(R)%*%t(Q)
        mid <- XtF %*% VG_inv %*% t(XtF)
        mid.qr <- qr(mid)
        Q1<-qr.Q(mid.qr)
        R1<-qr.R(mid.qr)
        mid_inv <- qr.solve(R1)%*%t(Q1)
        proj <- VG_inv - VG_inv %*% t(XtF) %*% mid_inv %*% XtF %*% VG_inv
        stat <- t(Fty) %*% proj %*% Fty
        
        as.numeric(stat - U_0_2)
      }
      
      # Step 4: solve for G_sigma2_d via uniroot
      G_sigma2_d <- uniroot(find_Gsigma2_d, lower = 1e-6, upper = 1e4, tol = 1e-6)$root
      
      # Step 5: compute VG
      VG <- G_sigma2_d * t(Fmat) %*% Z %*% t(Z) %*% Fmat + G_sigma2_e * I_r
      VG.qr <- qr(VG)
      Q2 <-qr.Q(VG.qr)
      R2 <-qr.R(VG.qr)
      VG_inv <- qr.solve(R2)%*%t(Q2)
      
      # Step 6: compute GPQ for x0'β
      A <- XtF %*% VG_inv %*% t(XtF)
      A.qr <- qr(A)
      Q3 <-qr.Q(A.qr)
      R3 <-qr.R(A.qr)
      A_inv <- qr.solve(R3)%*%t(Q3)
      mu_hat <- A_inv %*% XtF %*% VG_inv %*% Fty
      
      basis_fn <- function(x) as.numeric(x^(0:(p - 1)))  
      z_fn <- function(x) 1   #only random intercept        
      
      # Step 7: define center and radius functions
      center_fn <- function(x) {
        xvec <- basis_fn(x)
        as.numeric(t(xvec) %*% mu_hat - Zk * sqrt(t(xvec) %*% A_inv %*% xvec))
      }
      
      radius_fn <- function(x) {
        xvec <- basis_fn(x)
        zvec <- z_fn(x)
        sqrt(as.numeric(t(zvec) %*% zvec * G_sigma2_d + G_sigma2_e + t(xvec) %*% A_inv %*% xvec))
      }
      list(center_fn = center_fn, radius_fn = radius_fn)
    }
    
    grid.x      <- seq(x_min,x_max, length.out = 5000)
    center_mat  <- radius_mat <- matrix(NA_real_, nrow = L, ncol = length(grid.x))
    
    for (i in 1:L) {
      cr <- construct_center_radius(U_e2[i], U_0_2[i], Zk[i],
                                    Fmat, Z, SS_e, p,
                                    XtF, Fty, I_r)
      center_mat[i, ] <- vapply(grid.x, cr$center_fn, numeric(1))
      radius_mat[i, ] <- vapply(grid.x, cr$radius_fn, numeric(1))
    }
    
    center_avg <- colMeans(center_mat)
    radius_avg <- colMeans(radius_mat)
    
    cov_gpq_fun <- function(x, lambda,
                            grid.x, center_avg, radius_avg,
                            beta_true, sigma_d2, sigma_e2) {
      cval <- approx(grid.x, center_avg, xout = x, rule = 2)$y
      rval <- approx(grid.x, radius_avg, xout = x, rule = 2)$y
      
      mu_true <- beta_true[1] + beta_true[2] * x
      pnorm(cval + lambda * rval, mu_true, sqrt(sigma_d2 + sigma_e2)) -
        pnorm(cval - lambda * rval, mu_true, sqrt(sigma_d2 + sigma_e2))
    }
    
    min.cov.GPQ <- function(x.new, lambda) {
      cov_gpq_fun(x.new, lambda,
                  grid.x, center_avg, radius_avg,
                  beta_true, sigma_d2, sigma_e2)
    }
    
    #min.cov(X,invXTX,lam_two_sti_exact)
    cut.p=seq(x_min, x_max,length.out=10) # Cut [a,b] interval into small intervals
    covp2=0
    covb2=0
    for (ip in 1:(length(cut.p)-1)) {
      a1=cut.p[ip]
      b1=cut.p[ip+1]
      covp2[ip]=optimise(min.cov,c(a1,b1),invXTX=invXTX,lambda=lam_two_sti_exact)$objective
      covb2[ip]=min.cov(a1,invXTX,lam_two_sti_exact)
    }
    # Find the minimum value of P_{y_x}(l(x)<=y_x<=u(x)) among x\in[a,b]
    covu2=min.cov(x_max,invXTX,lam_two_sti_exact)
    covmin2=min(c(covp2,covb2,covu2))
    covmin1 <- optimise(min.cov.GPQ, interval = c(x_min, x_max),
                        lambda = lam_two_sti_GPQ)$objective
    res=c(covmin1,covmin2)
  }
})
stopImplicitCluster()
colnames(cover) <- c("GPQ_min_cov", "Exact_min_cov")
print(cover)
cat("γ-checking GPQ  :", mean(cover[,1] >= Cbeta), "\n")
cat("γ-checking Exact:", mean(cover[,2] >= Cbeta), "\n")

######################################
## Compute the coverage probability ##
######################################
con_G=0
con_E=0
for (i2 in 1:nsim1) {
  if(cover[i2,1]>= Cbeta){con_G=con_G+1}else{con_G=con_G}
  if(cover[i2,2]>= Cbeta){con_E=con_E+1}else{con_E=con_E}
}

cover_G=con_G/nsim1*100
cover_E=con_E/nsim1*100

#Print the values of lambda_s^M for STI in GPQ and lambda_s for exact STI and the r-ratios between them
lambda_table <- matrix(c(n,"ES",p,lam_two_sti_GPQ,lam_two_sti_exact), ncol=5, byrow=TRUE)
colnames(lambda_table) <- c('n','X','p','lambda_s^GPQ','lambda_s')
print(lambda_table)

lambda_table <- matrix(c(n,"ES",p,cover_G,cover_E), ncol=5, byrow=TRUE)
colnames(lambda_table) <- c('n','X','p','GPQ','EXACT')
print(lambda_table)

