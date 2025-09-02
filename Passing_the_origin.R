library(MASS)
library(doParallel)
library(DEoptim)

Cgamma = 0.99
Cbeta = 0.95
zbeta = qnorm((1+Cbeta)/2)

####The training data set
x <- c(1.56, 1.56, 1.27, 1, 0.46, 1.38, 1.97, 0.86, 1.66, 1.87, 0.82, 1.22, 1.43, 0.97, 1.35, 1.71, 0.72, 0.72, 0.72, 1.4, 0.83, 2, 1.45, 1.26, 1.37, 2.27, 1.5, 1.33, 1.61, 0.58, 1.12, 1.1, 1.23, 0.95, 1.43, 1.15, 1.26, 1.26, 1.39, 0.87, 1.57, 0.96, 1.51, 1.86, 0.44, 0.86, 0.57, 0.84, 0.73, 1.82, 0.4, 1.36, 2.3, 0.38, 1.41, 1.91, 1.09, 1.75, 0.3, 2.05, 0.35, 0.75, 0.19, 1.53, 0.79, 0.11, 0.22, 0.34, 0.21, 0.21, 0.51, 0.48, 0.82, 0.7, 0.18, 0.38, 0.29, 0.49, 0.49, 1.34, 0.35, 0.37, 0.28, 0.92, 0.32, 0.32, 0.96, 0.72, 0.91, 0.76, 0.4, 0.42, 0.65, 0.23, 0.5, 0.52, 0.52, 0.83, 0.29, 0.87, 1.24, 1.24, 0.35, 0.35, 2.22, 1.47, 1.06, 0.57, 0.69, 0.95, 0.93, 1.13, 0.3, 0.32, 0.38, 0.42, 0.56, 0.68, 0.68, 0.68, 0.82, 0.88, 0.88, 1.16, 1.26, 1.66, 1.1, 1.1, 1.49, 1.35, 0.97, 0.87, 1.62, 1.87, 0.69, 1.81, 1.71, 1.14, 1.08, 0.49, 0.98, 0.45, 0.45, 1.35, 1.33, 1.29, 1.25, 0.82, 0.8, 1.15, 0.72, 0.7, 0.66, 0.66, 0.66, 0.31, 0.29, 1.43, 1.04)
y <- c(0.75, 0.75, 0.61, 0.48, 0.22, 0.66, 0.94, 0.41, 0.79, 0.89, 0.39, 0.58, 0.68, 0.46, 0.64, 0.81, 0.34, 0.34, 0.34, 0.66, 0.39, 0.94, 0.68, 0.59, 0.64, 1.06, 0.7, 0.62, 0.75, 0.27, 0.52, 0.51, 0.57, 0.44, 0.66, 0.53, 0.58, 0.58, 0.64, 0.4, 0.72, 0.44, 0.69, 0.85, 0.2, 0.39, 0.27, 0.38, 0.33, 0.82, 0.18, 0.61, 1.03, 0.17, 0.63, 0.85, 0.48, 0.76, 0.13, 0.88, 0.15, 0.32, 0.08, 0.64, 0.33, 0.07, 0.14, 0.2, 0.12, 0.12, 0.29, 0.27, 0.46, 0.39, 0.1, 0.21, 0.16, 0.27, 0.27, 0.73, 0.19, 0.2, 0.15, 0.49, 0.17, 0.17, 0.51, 0.38, 0.48, 0.4, 0.21, 0.22, 0.34, 0.12, 0.26, 0.27, 0.27, 0.43, 0.15, 0.45, 0.64, 0.64, 0.18, 0.18, 1.14, 0.75, 0.54, 0.29, 0.35, 0.48, 0.47, 0.57, 0.15, 0.16, 0.19, 0.21, 0.28, 0.34, 0.34, 0.34, 0.41, 0.44, 0.44, 0.58, 0.63, 0.83, 0.5, 0.5, 0.74, 0.67, 0.48, 0.43, 0.8, 0.84, 0.34, 0.89, 0.84, 0.56, 0.53, 0.24, 0.48, 0.22, 0.22, 0.66, 0.65, 0.63, 0.61, 0.4, 0.39, 0.56, 0.35, 0.34, 0.32, 0.32, 0.32, 0.15, 0.14, 0.69, 0.5)
yy <- matrix(y,ncol=1)
n <- length(x)    
r <- 1
p <- r
nu <- n-p

##################
## [a,b]=[0,20] ##
##################

a <- min(x)
b <- max(x)  
X<-cbind(x)    #The design matrix
XtX <- t(X)%*%X

#Find the inverse of XtX by using Q-R decomposition
XtX.qr <- qr(XtX)
Q<-qr.Q(XtX.qr)
R<-qr.R(XtX.qr)
invXTX <- qr.solve(R)%*%t(Q)
invXTX <- as.matrix(invXTX)

##########################################################################
## alpha.hat: the least squares estimator of coefficient alpha          ##
## sigma.hat: the square root of the least squares estimator of sigma^2 ##
##########################################################################

alpha.hat <- invXTX %*% t(X) %*% yy
sigma.hat <- sqrt(sum((yy-X%*%alpha.hat)^2)/nu)

two.sti.exact = function(lam, n, p, Cbeta, a, b, B, u, invXTX, cut.p)
{
  zp = qnorm((1 + Cbeta)/2)
  
  fd = function(x, lam, invXTX){
    x <- as.numeric(x[1])
    vx <- matrix(x^(1:p), ncol = 1)
    mu <- drop(t(vx) %*% B)
    s  <- sqrt( drop(t(vx) %*% invXTX %*% vx) )
    pnorm(mu + lam*(zp + sqrt(p+2)*s)*u) -
      pnorm(mu - lam*(zp + sqrt(p+2)*s)*u)
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

L=100000
#Set random seed to reproduce the results
set.seed(10)
#Generate the Z and u
#B.sim=Z~N_p(0,invXTX)
#u.sim=sqrt(chi_{n-p}^2/(n-p))
B.sim = mvrnorm(L,c(rep(0,p)),invXTX)
u.sim = sqrt(rchisq(L,df=n-p)/(n-p))
tol=0.0001

numCores = 24
registerDoParallel(numCores)
system.time({
  lam_exact= foreach(i = 1:L,.combine = "c") %dopar%{
    B = as.numeric(B.sim[i, ])
    u = u.sim[i]
    #Find the lambda for this simulation repeat using the bisection method
    #We use the initial interval (0, 5) as the (klow, kupp) that contains lambda
    klow <- 0 
    kupp <- 5
    #In our numerical results, cut.p = 6 is enough to find the global optimisation.
    fupp <- two.sti.exact(kupp,n,p,Cbeta,a,b,B,u,invXTX,cut.p = 6)
    
    if (fupp <= 0) {
      klow <- kupp
      repeat {
        kupp <- kupp + 5
        fupp <- two.sti.exact(kupp, n, p, Cbeta, a, b, B, u, invXTX, cut.p = 6)
        if (fupp > 0) break
      }
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
############################################################
## The critical constant lambda_s for exact two-sided STI ##
############################################################
lam_two_sti_exact = quantile(lam_exact, Cgamma)

lambda_s_table <- matrix(c(a,b,lam_two_sti_exact), ncol=3, byrow=TRUE)
colnames(lambda_s_table) <- c('a','b','lambda_s')

###################################################
## Calculate the confidence sets for covariate x ##
###################################################

upper.bound <- function(x, y, lam) {
  # y_hat_upper(x) - y
  as.numeric(alpha.hat[1]) * x -
    lam * (zbeta + sqrt(p + 2) * sqrt(t(c(x)) %*% invXTX %*% c(x))) * sigma.hat -
    y
}
upper <- uniroot(upper.bound, c(-10, b), y = 0.4, lam = lam_two_sti_exact)$root

lower.bound <- function(x, y, lam) {
  # y_hat_lower(x) - y
  as.numeric(alpha.hat[1]) * x +
    lam * (zbeta + sqrt(p + 2) * sqrt(t(c(x)) %*% invXTX %*% c(x))) * sigma.hat -
    y
}
lower <- uniroot(lower.bound, c(-10, b), y = 0.4, lam = lam_two_sti_exact)$root

confidence_set_table <- matrix(c(a,b,0.4,lower,upper), ncol=5, byrow=TRUE)
colnames(confidence_set_table) <- c('a','b','y','lower bound','upper bound')

######################
## Compute the STIs ##
######################
x.new  <- seq(from = min(x), to = max(x), by = 0.01)
beta.h <- as.numeric(alpha.hat[1])
y.fit  <- beta.h * x.new

dsqrt <- numeric(length(x.new))
for (i in 1:length(x.new)) {
  dsqrt[i] <- as.numeric(t(c(x.new[i])) %*% invXTX %*% c(x.new[i]))
}

upper.sti <- beta.h * x.new + lam_two_sti_exact * (zbeta + sqrt(p + 2) * sqrt(dsqrt)) * sigma.hat
lower.sti <- beta.h * x.new - lam_two_sti_exact * (zbeta + sqrt(p + 2) * sqrt(dsqrt)) * sigma.hat

##############################################
#### Plot Figure for confidence sets of x ####
##############################################
upper_label <- format(round(upper, 3), nsmall = 3)
lower_label <- format(round(lower, 3), nsmall = 3)

plot(x.new, y.fit,
     xlab = 'BAC', 
     ylab = 'BrAC',
     xlim = c(0, max(x)),
     ylim = c(0, max(y)),
     type = 'l', lty = 3, col = "black",
     cex.axis = 1.3, cex.lab = 1.5)

points(x, y, pch = 8)

lines(x.new, upper.sti, col = "red", lty = 2, lwd = 2)
lines(x.new, lower.sti, col = "red", lty = 2, lwd = 2)

abline(h = 0.4, col = "blue", lty = 3, lwd = 1.5)

points(upper, 0.4, col = "red", pch = 15)
points(lower, 0.4, col = "red", pch = 15)

arrows(x0 = upper, y0 = 0.4, x1 = upper, y1 = -0.05,
       length = 0.1, col = "blue", lwd = 1.2)
arrows(x0 = lower, y0 = 0.4, x1 = lower, y1 = -0.05,
       length = 0.1, col = "blue", lwd = 1.2)

axis(1, at = c(lower, upper), labels = FALSE)

text(x = c(lower, upper),
     y = par("usr")[3] - 0.05,   
     labels = c(lower_label, upper_label),
     srt = 45,                   
     adj = 1,                   
     xpd = TRUE,                
     col = "red",
     cex = 1.2)
abline(a=0, b=beta.h, col="darkgreen", lty=1, lwd=2) 
legend("topleft", inset = .005, 
       legend = c("Fitted regression line (through origin)",
                  "Exact two-sided STI",
                  "BAC = 0.4 threshold"),
       lty = c(1, 2, 3),
       col = c("darkgreen", "red", "blue"),
       cex = 1.3,
       bty = "n")

#Print the values of lambda_s for exact STI
print(lambda_s_table)
#Print the values of confidence set for x based on exact STI
print(confidence_set_table)

