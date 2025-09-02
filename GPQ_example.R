library(doParallel)
library(MASS)  
library(DEoptim)
library(nlme)

Cgamma  <- 0.90   # confidence level
Cbeta   <- 0.95   # coverage level

data(Orthodont, package = "nlme")
fit <- lme(distance ~ age, random = ~1 | Subject, data = Orthodont)
summary(fit)
beta_hat <- fixef(fit)               
sigma_d  <- as.numeric(VarCorr(fit)["(Intercept)", "StdDev"])
sigma_e  <- fit$sigma               

X <- model.matrix(~ age, Orthodont)
x <- Orthodont$age
x_min <- min(x)
x_max <- max(x)
sigma_hat_x  <- sqrt(sigma_d^2 + sigma_e^2)   
p <-2
Z <- model.matrix(~ Subject - 1, Orthodont) 
y <- Orthodont$distance
N <- length(y)
XZ <- cbind(X, Z)

r <- Matrix::rankMatrix(XZ)[]
b <- Matrix::rankMatrix(X)[]
I_r <- diag(r)
Fmat  <- svd(XZ, nu = r, nv = 0)$u 
XtF  <- t(X) %*% Fmat
Fty  <- t(Fmat) %*% y
P_XZ <- Fmat %*% t(Fmat)
SS_e <- as.numeric(t(y) %*% (diag(N) - P_XZ) %*% y)

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
                                 #We use the initial interval (0, 5) as the (klow, kupp) that contains lambda
                                 klow <- 0 
                                 kupp <- 5
                                 #In our numerical results, cut.p = 6 is enough to find the global optimisation.
                                 fupp <- two.sti.GPQ(kupp, this_U_e2, this_U_02, this_Zk,
                                                     Fmat, Z, SS_e, p,
                                                     XtF, Fty, I_r,
                                                     Cbeta, x_min, x_max, cut.p=10, 
                                                     beta_hat,sigma_hat_x)
                                 if (fupp <= 0) {
                                   klow <- kupp
                                   repeat {
                                     kupp <- kupp + 5
                                     fupp <- two.sti.GPQ(kupp, this_U_e2, this_U_02, this_Zk,
                                                         Fmat, Z, SS_e, p,
                                                         XtF, Fty, I_r,
                                                         Cbeta, x_min, x_max, cut.p = 10, 
                                                         beta_hat,sigma_hat_x)
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
                                                       beta_hat,sigma_hat_x)
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

lam_two_final = quantile(lam_GPQ, Cgamma)
head(lam_GPQ)

cat(sprintf("The estimated lambda at gamma = %.2f is: %.4f\n", Cgamma, lam_two_final))

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

lower_avg <- center_avg - lam_two_final * radius_avg
upper_avg <- center_avg + lam_two_final * radius_avg

pdf("lambda_histogram.pdf", width = 5, height = 4)
hist(lam_GPQ, breaks = 15, col = "skyblue",
     main = expression("Distribution of " * lambda),
     xlab = expression(lambda))
abline(v = lam_two_final, col = "red", lwd = 2)
dev.off()

lm_fit   <- lm(distance ~ age, data = Orthodont)
a <- min(Orthodont$age)
b <- max(Orthodont$age)
X_lm     <- model.matrix(lm_fit)
y_lm     <- Orthodont$distance  
p     <- 2               
nu_lm    <- nrow(X_lm) - p
n <- nrow(X_lm)
zbeta    <- qnorm((1 + Cbeta)/2)
beta_hat <- solve(t(X_lm) %*% X_lm, t(X_lm) %*% y_lm)
resid    <- y_lm - X_lm %*% beta_hat
sigma.hat <- sqrt(sum(resid^2) / nu_lm)
X_lm.qr <- qr(t(X_lm)%*%X_lm)
Q4<-qr.Q(X_lm.qr)
R4<-qr.R(X_lm.qr)
invXTX <- as.matrix( qr.solve(R4) %*% t(Q4) )

two.sti.exact=function(lam,n,p,Cbeta,a,b,B,u,invXTX,cut.p)
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

set.seed(120)
B.sim = mvrnorm(L,c(rep(0,p)),invXTX)
u.sim = sqrt(rchisq(L,df=n-p)/(n-p))
tol=0.0001
numCores = 36
registerDoParallel(numCores)
system.time({
  lam_exact= foreach(i = 1:L,.combine = "c") %dopar%{
    B = B.sim[i,]
    u = u.sim[i]
    klow <- 0 
    kupp <- 5
    fupp <- two.sti.exact(kupp,n,p,Cbeta,a,b,B,u,invXTX,cut.p = 6)
    
    if(fupp <=0){
      klow=kupp
      kupp=kupp+5
    }
    
    while ( (kupp-klow)>tol ){
      cmid <- (klow+kupp)/2
      c <- cmid
      fmid <- two.sti.exact(c,n,p,Cbeta,a,b,B,u,invXTX,cut.p = 6)
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
cat(sprintf("lambda_exact (γ-quantile) = %.4f  [L = %d]\n", lam_two_sti_exact, L))

x.new <- seq(from = a, to = b, by = 0.1)   

y.fit <- beta_hat[1] + beta_hat[2] * x.new   

dsqrt <- sapply(x.new, function(xx) {
  as.numeric(t(c(1, xx)) %*% invXTX %*% c(1, xx))
})

upper.sti <- y.fit + lam_two_sti_exact * (zbeta + sqrt(p + 2) * sqrt(dsqrt)) * sigma.hat
lower.sti <- y.fit - lam_two_sti_exact * (zbeta + sqrt(p + 2) * sqrt(dsqrt)) * sigma.hat

pdf("band_GPQ_STI.pdf", width = 6, height = 4)
plot(grid.x, center_avg, type = "n",
     ylim = range(lower_avg, upper_avg, y),
     xlab = "Age", ylab = "Distance",
     main = "GPQ-Simultaneous Tolerance Band")

polygon(c(grid.x, rev(grid.x)),
        c(lower_avg, rev(upper_avg)),
        col = rgb(0, 102/255, 204/255, 0.15), border = NA)

lines(grid.x, center_avg, lwd = 2, col = rgb(0, 102/255, 204/255))
lines(grid.x, lower_avg, lty = 2, col = rgb(0, 102/255, 204/255))
lines(grid.x, upper_avg, lty = 2, col = rgb(0, 102/255, 204/255))

points(x, y, pch = 16, col = rgb(0,0,0,0.4), cex = 0.6)

legend("topleft",
       legend = c("GPQ-STB", "Data"),
       fill   = c(rgb(0, 102/255, 204/255, 0.15), NA),
       border = c(NA, NA),
       lty    = c(1, NA),
       lwd    = c(2, NA),
       pch    = c(NA, 16),
       col    = c(rgb(0, 102/255, 204/255), rgb(0,0,0,0.4)),
       bty    = "n")
dev.off()

pdf("band_ExactSTI.pdf", width = 6, height = 4)
plot(x.new, y.fit, type = "n",
     ylim = range(lower.sti, upper.sti, Orthodont$distance),
     xlab = "Age", ylab = "Distance",
     main = "Exact Simultaneous Tolerance Band")

polygon(c(x.new, rev(x.new)),
        c(lower.sti, rev(upper.sti)),
        col = rgb(255/255, 140/255, 0, 0.15), border = NA)

lines(x.new, y.fit, lwd = 2, col = rgb(255/255, 140/255, 0))

lines(x.new, lower.sti, lty = 2, col = rgb(255/255, 140/255, 0))
lines(x.new, upper.sti, lty = 2, col = rgb(255/255, 140/255, 0))

points(Orthodont$age, Orthodont$distance, pch = 16, col = rgb(0,0,0,0.4), cex = 0.6)

legend("topleft",
       legend = c("Exact STB", "Data"),
       fill   = c(rgb(255/255, 140/255, 0, 0.15), NA),
       border = c(NA, NA),
       lty    = c(1, NA),
       lwd    = c(2, NA),
       pch    = c(NA, 16),
       col    = c(rgb(255/255, 140/255, 0), rgb(0,0,0,0.4)),
       bty    = "n")
dev.off()

pdf("band_Combined.pdf", width = 6, height = 4)
plot(grid.x, center_avg, type = "n",
     ylim = range(lower_avg, upper_avg, lower.sti, upper.sti, Orthodont$distance),
     xlab = "Age", ylab = "Distance",
     main = "Combined Bands: GPQ-STB vs Exact STB")

polygon(c(grid.x, rev(grid.x)),
        c(lower_avg, rev(upper_avg)),
        col = rgb(0, 102/255, 204/255, 0.12), border = NA)
lines(grid.x, center_avg, lwd = 2, col = rgb(0, 102/255, 204/255))
lines(grid.x, lower_avg, lty = 2, col = rgb(0, 102/255, 204/255))
lines(grid.x, upper_avg, lty = 2, col = rgb(0, 102/255, 204/255))

polygon(c(x.new, rev(x.new)),
        c(lower.sti, rev(upper.sti)),
        col = rgb(255/255, 140/255, 0, 0.12), border = NA)
lines(x.new, y.fit, lwd = 2, col = rgb(255/255, 140/255, 0))
lines(x.new, lower.sti, lty = 2, col = rgb(255/255, 140/255, 0))
lines(x.new, upper.sti, lty = 2, col = rgb(255/255, 140/255, 0))

points(Orthodont$age, Orthodont$distance, pch = 16, col = rgb(0,0,0,0.4), cex = 0.6)

legend("topleft",
       legend = c("GPQ-STB", "Exact STB", "Data"),
       fill   = c(rgb(0,102/255,204/255,0.12),
                  rgb(255/255,140/255,0,0.12),
                  NA),
       border = c(NA, NA, NA),
       lty    = c(1, 1, NA),
       lwd    = c(2, 2, NA),
       pch    = c(NA, NA, 16),
       col    = c(rgb(0,102/255,204/255),
                  rgb(255/255,140/255,0),
                  rgb(0,0,0,0.4)),
       pt.cex = c(NA, NA, 0.6),
       bty    = "n")
dev.off()

