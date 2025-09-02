library(doParallel)
library(MASS)
library(DEoptim)

## ------------------ Basic Settings ------------------
Cbeta  = 0.99
Cgamma = 0.99
zbeta  = qnorm((1 + Cbeta) / 2)

n  = 20
a  = 0
b  = 10
p  = 3
nu = n - p

set.seed(20250716)
alp.t <- c(1, 2, -1)   # True model: 1 + 2x - 1x^2
sig.t <- 1

## ------------------ Build Design Matrix & Fit Model ------------------
X1 <- seq(from = a, to = b, length.out = n)
X  <- matrix(0, nrow = n, ncol = p)
for (i in 1:length(X1)) {
  x <- X1[i]; vx <- rep(1, p)
  if (p >= 2) for (j in 2:p) vx[j] <- x^(j-1)
  X[i, ] <- vx
}

XTX    <- t(X) %*% X
XTX.qr <- qr(XTX)
Q <- qr.Q(XTX.qr)
R <- qr.R(XTX.qr)
invXTX <- qr.solve(R) %*% t(Q)

# Generate training responses and fit the model
eps <- rnorm(n, mean = 0, sd = sig.t)
y   <- as.vector(X %*% alp.t + eps)
yy  <- matrix(y, ncol = 1)

alpha.hat <- invXTX %*% t(X) %*% yy
sigma.hat <- sqrt(sum((yy - X %*% alpha.hat)^2) / nu)

# Prepare lines for plotting
training_data <- data.frame(x = X1, y = y)
x.new  <- seq(from = min(X1), to = max(X1), by = 0.01)
y.true <- alp.t[1] + alp.t[2] * x.new + alp.t[3] * x.new^2
y.fit  <- alpha.hat[1] + alpha.hat[2] * x.new + alpha.hat[3] * x.new^2

## ------------------ Future points: 30 normal + additional outliers ------------------
# 30 purely normal future points
set.seed(20250801)
n_future  <- 30
x_future  <- seq(a, b, length.out = n_future)
eps_future <- rnorm(n_future, mean = 0, sd = sig.t)
y_future_norm <- alp.t[1] + alp.t[2] * x_future + alp.t[3] * x_future^2 + eps_future

data_main <- data.frame(X = x_future, Y = y_future_norm, label = 0L)  # Normal points

# Add 3 vertical outliers (large residuals, same X positions but different Y)
set.seed(110)
n_outlier_vert <- 3
idx_vert <- sample(1:n_future, n_outlier_vert)
y_future_vert <- y_future_norm[idx_vert] + rnorm(n_outlier_vert, mean = 0, sd = 8)
data_vert <- data.frame(X = x_future[idx_vert], Y = y_future_vert, label = 1L)

# Add 3 structural outliers (slightly different model parameters, eps ~ N(0,1))
set.seed(20250810)
beta_struct <- c(1.2, 1.8, -0.9)
n_outlier_struct <- 3
x_struct <- runif(n_outlier_struct, 7, 8)
y_struct <- beta_struct[1] + beta_struct[2] * x_struct + beta_struct[3] * x_struct^2 +
  rnorm(n_outlier_struct, mean = 0, sd = 1)
data_struct <- data.frame(X = x_struct, Y = y_struct, label = 1L)

# Add 3 X-measurement-error outliers (horizontal/xerror)
set.seed(20250825)
n_outlier_x <- 3
sigma_x     <- 0.1
bias_x      <- 0.5
x_true  <- runif(n_outlier_x, 6, 9)
x_obs   <- x_true + rnorm(n_outlier_x, mean = bias_x, sd = sigma_x)
y_xerr  <- alp.t[1] + alp.t[2]*x_true + alp.t[3]*x_true^2 + rnorm(n_outlier_x, 0, sig.t)
data_xerr <- data.frame(X = x_obs, Y = y_xerr, label = 1L)

# Combine: 30 normal + 3 vertical + 3 structural + 3 xerror = 39
future_data <- rbind(data_main, data_vert, data_struct, data_xerr)
row.names(future_data) <- NULL

cat("nrow(future_data) = ", nrow(future_data), "\n")
print(table(future_data$label))

lam_two_sti_exact=1.436386

x.new <- seq(from = min(X1), to = max(X1), by = 0.01)
y.fit=alpha.hat[1]+alpha.hat[2]*x.new+alpha.hat[3]*x.new^2
dsqrt=0
for (i in 1:length(x.new)) {
  dsqrt[i]=as.numeric(t(c(1,x.new[i],x.new[i]^2))%*%invXTX%*%c(1,x.new[i],x.new[i]^2))
}
upper.sti=alpha.hat[1]+alpha.hat[2]*x.new+alpha.hat[3]*x.new^2+lam_two_sti_exact*(zbeta+sqrt(p+2)*sqrt(dsqrt))*sigma.hat
lower.sti=alpha.hat[1]+alpha.hat[2]*x.new+alpha.hat[3]*x.new^2-lam_two_sti_exact*(zbeta+sqrt(p+2)*sqrt(dsqrt))*sigma.hat


## ------------------ Plot (without STI band) ------------------
plot(x.new, y.fit,
     xlab = 'X', ylab = 'Y',
     xlim = range(c(0, 10, future_data$X)),     
     ylim = range(c(y, future_data$Y)),
     type = 'l', lty = 3, col = "black",
     xaxt = 'n', cex.axis = 1.3, cex.lab = 1.5)
axis(1, at = seq(0, 10, by = 1), cex.axis = 1.1)

# Original training data (black solid points)
points(X1, y, pch = 16, col = "black", cex = 0.8)

# True regression line (green solid line)
lines(x.new, y.true, col = "green3", lwd = 2)

lines(x.new, upper.sti, col = 2, lty = 2, lwd = 2)
lines(x.new, lower.sti, col = 2, lty = 2, lwd = 2)

# Normal future points (blue hollow circles)
points(data_main$X, data_main$Y, pch = 1, col = "blue", cex = 1.2)

# Vertical outliers (red crosses)
points(data_vert$X, data_vert$Y, pch = 4, col = "red", lwd = 2, cex = 1.5)

# Structural outliers (purple squares)
points(data_struct$X, data_struct$Y, pch = 15, col = "purple", cex = 1.3)

# X-measurement-error outliers (orange triangles)
points(data_xerr$X, data_xerr$Y, pch = 17, col = "orange", cex = 1.6)

legend("bottomleft", inset = .005,
       legend = c("True regression line",
                  "Fitted regression line",
                  "STI (beta = 0.99, gamma = 0.99)",
                  "Original training points",
                  "Normal future points",
                  "Vertical outliers",
                  "Structural outliers",
                  "X-error outliers"),
       lty = c(1, 3, 2, NA, NA, NA, NA, NA),
       pch = c(NA, NA, NA, 16, 1, 4, 15, 17),
       col = c("green3", "black", "red", "black", "blue", "red", "purple", "orange"),
       cex = 1.05, bty = "n")


## -------- β*(inverse) and λ* for all 36 points --------
## If lam_fun is not defined in your session, build it here:

txt <- "     beta gamma   lambda
1   0.750  0.99 1.215900
2   0.751  0.99 1.216660
3   0.752  0.99 1.217424
4   0.753  0.99 1.218187
5   0.754  0.99 1.218796
6   0.755  0.99 1.219254
7   0.756  0.99 1.219864
8   0.757  0.99 1.220627
9   0.758  0.99 1.221393
10  0.759  0.99 1.222157
11  0.760  0.99 1.223071
12  0.761  0.99 1.223984
13  0.762  0.99 1.224290
14  0.763  0.99 1.225052
15  0.764  0.99 1.225815
16  0.765  0.99 1.226578
17  0.766  0.99 1.227039
18  0.767  0.99 1.227800
19  0.768  0.99 1.228258
20  0.769  0.99 1.228717
21  0.770  0.99 1.229477
22  0.771  0.99 1.230240
23  0.772  0.99 1.230699
24  0.773  0.99 1.231766
25  0.774  0.99 1.232529
26  0.775  0.99 1.233292
27  0.776  0.99 1.234055
28  0.777  0.99 1.234819
29  0.778  0.99 1.235582
30  0.779  0.99 1.236343
31  0.780  0.99 1.236804
32  0.781  0.99 1.237718
33  0.782  0.99 1.238634
34  0.783  0.99 1.239548
35  0.784  0.99 1.240009
36  0.785  0.99 1.240923
37  0.786  0.99 1.241380
38  0.787  0.99 1.241533
39  0.788  0.99 1.241991
40  0.789  0.99 1.242905
41  0.790  0.99 1.243365
42  0.791  0.99 1.243822
43  0.792  0.99 1.244431
44  0.793  0.99 1.245195
45  0.794  0.99 1.246109
46  0.795  0.99 1.247025
47  0.796  0.99 1.247485
48  0.797  0.99 1.248096
49  0.798  0.99 1.248856
50  0.799  0.99 1.249313
51  0.800  0.99 1.249925
52  0.801  0.99 1.250536
53  0.802  0.99 1.251144
54  0.803  0.99 1.251755
55  0.804  0.99 1.252368
56  0.805  0.99 1.253281
57  0.806  0.99 1.254045
58  0.807  0.99 1.254808
59  0.808  0.99 1.255569
60  0.809  0.99 1.256181
61  0.810  0.99 1.256793
62  0.811  0.99 1.257707
63  0.812  0.99 1.258774
64  0.813  0.99 1.259537
65  0.814  0.99 1.260300
66  0.815  0.99 1.260912
67  0.816  0.99 1.261522
68  0.817  0.99 1.262283
69  0.818  0.99 1.262743
70  0.819  0.99 1.263506
71  0.820  0.99 1.264421
72  0.821  0.99 1.265337
73  0.822  0.99 1.265945
74  0.823  0.99 1.266556
75  0.824  0.99 1.267014
76  0.825  0.99 1.267625
77  0.826  0.99 1.268539
78  0.827  0.99 1.269304
79  0.828  0.99 1.270068
80  0.829  0.99 1.270833
81  0.830  0.99 1.271745
82  0.831  0.99 1.272357
83  0.832  0.99 1.272968
84  0.833  0.99 1.273578
85  0.834  0.99 1.274187
86  0.835  0.99 1.275253
87  0.836  0.99 1.276169
88  0.837  0.99 1.276782
89  0.838  0.99 1.277695
90  0.839  0.99 1.278610
91  0.840  0.99 1.279526
92  0.841  0.99 1.279988
93  0.842  0.99 1.280600
94  0.843  0.99 1.281662
95  0.844  0.99 1.282121
96  0.845  0.99 1.282883
97  0.846  0.99 1.283951
98  0.847  0.99 1.284866
99  0.848  0.99 1.285629
100 0.849  0.99 1.286392
101 0.850  0.99 1.287157
102 0.851  0.99 1.287920
103 0.852  0.99 1.288986
104 0.853  0.99 1.289598
105 0.854  0.99 1.290208
106 0.855  0.99 1.290970
107 0.856  0.99 1.291733
108 0.857  0.99 1.292497
109 0.858  0.99 1.293411
110 0.859  0.99 1.294023
111 0.860  0.99 1.294937
112 0.861  0.99 1.295549
113 0.862  0.99 1.296313
114 0.863  0.99 1.297382
115 0.864  0.99 1.298447
116 0.865  0.99 1.299211
117 0.866  0.99 1.299823
118 0.867  0.99 1.300735
119 0.868  0.99 1.301651
120 0.869  0.99 1.301958
121 0.870  0.99 1.302570
122 0.871  0.99 1.303177
123 0.872  0.99 1.303636
124 0.873  0.99 1.304552
125 0.874  0.99 1.305315
126 0.875  0.99 1.306229
127 0.876  0.99 1.306841
128 0.877  0.99 1.307452
129 0.878  0.99 1.308517
130 0.879  0.99 1.309282
131 0.880  0.99 1.310045
132 0.881  0.99 1.310657
133 0.882  0.99 1.311722
134 0.883  0.99 1.312639
135 0.884  0.99 1.313249
136 0.885  0.99 1.314163
137 0.886  0.99 1.315079
138 0.887  0.99 1.315845
139 0.888  0.99 1.316759
140 0.889  0.99 1.317520
141 0.890  0.99 1.317979
142 0.891  0.99 1.318742
143 0.892  0.99 1.319658
144 0.893  0.99 1.320271
145 0.894  0.99 1.321184
146 0.895  0.99 1.321950
147 0.896  0.99 1.323015
148 0.897  0.99 1.323929
149 0.898  0.99 1.324846
150 0.899  0.99 1.325763
151 0.900  0.99 1.326830
152 0.901  0.99 1.327745
153 0.902  0.99 1.328659
154 0.903  0.99 1.329422
155 0.904  0.99 1.330643
156 0.905  0.99 1.331711
157 0.906  0.99 1.332477
158 0.907  0.99 1.333542
159 0.908  0.99 1.334610
160 0.909  0.99 1.335222
161 0.910  0.99 1.336288
162 0.911  0.99 1.337204
163 0.912  0.99 1.338274
164 0.913  0.99 1.339037
165 0.914  0.99 1.340103
166 0.915  0.99 1.340868
167 0.916  0.99 1.341631
168 0.917  0.99 1.342545
169 0.918  0.99 1.343614
170 0.919  0.99 1.344835
171 0.920  0.99 1.345598
172 0.921  0.99 1.346362
173 0.922  0.99 1.347429
174 0.923  0.99 1.348343
175 0.924  0.99 1.349258
176 0.925  0.99 1.349870
177 0.926  0.99 1.350632
178 0.927  0.99 1.351396
179 0.928  0.99 1.352161
180 0.929  0.99 1.353227
181 0.930  0.99 1.354143
182 0.931  0.99 1.355060
183 0.932  0.99 1.355972
184 0.933  0.99 1.356737
185 0.934  0.99 1.357657
186 0.935  0.99 1.359030
187 0.936  0.99 1.360397
188 0.937  0.99 1.361316
189 0.938  0.99 1.362534
190 0.939  0.99 1.363602
191 0.940  0.99 1.364365
192 0.941  0.99 1.365433
193 0.942  0.99 1.366197
194 0.943  0.99 1.367265
195 0.944  0.99 1.368181
196 0.945  0.99 1.369250
197 0.946  0.99 1.370164
198 0.947  0.99 1.371384
199 0.948  0.99 1.372455
200 0.949  0.99 1.373369
201 0.950  0.99 1.374132
202 0.951  0.99 1.375201
203 0.952  0.99 1.376572
204 0.953  0.99 1.377338
205 0.954  0.99 1.378102
206 0.955  0.99 1.379471
207 0.956  0.99 1.380693
208 0.957  0.99 1.382066
209 0.958  0.99 1.382985
210 0.959  0.99 1.383902
211 0.960  0.99 1.385423
212 0.961  0.99 1.386490
213 0.962  0.99 1.387711
214 0.963  0.99 1.388780
215 0.964  0.99 1.390152
216 0.965  0.99 1.391072
217 0.966  0.99 1.392747
218 0.967  0.99 1.394121
219 0.968  0.99 1.395645
220 0.969  0.99 1.397018
221 0.970  0.99 1.398546
222 0.971  0.99 1.399919
223 0.972  0.99 1.401447
224 0.973  0.99 1.402818
225 0.974  0.99 1.404190
226 0.975  0.99 1.405566
227 0.976  0.99 1.407394
228 0.977  0.99 1.409074
229 0.978  0.99 1.410753
230 0.979  0.99 1.412431
231 0.980  0.99 1.414873
232 0.981  0.99 1.416551
233 0.982  0.99 1.417924
234 0.983  0.99 1.420059
235 0.984  0.99 1.422043
236 0.985  0.99 1.424181
237 0.986  0.99 1.426317
238 0.987  0.99 1.428300
239 0.988  0.99 1.430742
240 0.989  0.99 1.433487
241 0.990  0.99 1.436386
242 0.991  0.99 1.439287
243 0.992  0.99 1.442795
244 0.993  0.99 1.446152
245 0.994  0.99 1.450273
246 0.995  0.99 1.454242
247 0.996  0.99 1.459122
248 0.997  0.99 1.465073
249 0.998  0.99 1.471178
250 0.999  0.99 1.482166"

df <- read.table(text = txt, header = TRUE)
stopifnot(all(diff(df$beta) > 0), all(diff(df$lambda) > 0))
lam_fun <- approxfun(df$beta, df$lambda, method = "linear", rule = 2)

# Helpers
yhat_one <- function(x, beta_hat) {
  p <- length(beta_hat)
  sum(beta_hat * c(1, x, x^2)[1:p])
}

beta_lambda_star_exact <- function(x0, y0,
                                   alpha.hat, sigma.hat, invXTX, p,
                                   lam_fun,
                                   beta_lo = min(df$beta),
                                   beta_hi = max(df$beta),
                                   tol = 1e-9) {
  vx <- c(1, x0, x0^2)[1:p]
  v0 <- sqrt(as.numeric(t(vx) %*% invXTX %*% vx))
  c0 <- sqrt(p + 2) * v0
  r  <- abs(y0 - sum(alpha.hat * vx))
  
  f <- function(beta) lam_fun(beta) * sigma.hat * ( qnorm((1 + beta)/2) + c0 ) - r
  
  flo <- f(beta_lo); fhi <- f(beta_hi)
  
  if (flo >= 0) {
    beta_star <- "<0.75"   
    status <- "clamped_at_low"
  } else if (fhi <= 0) {
    beta_star <- ">0.999"  
    status <- "clamped_at_high"
  } else {
    beta_star <- uniroot(f, c(beta_lo, beta_hi), tol = tol)$root
    beta_star <- round(beta_star, 6) 
    status <- "ok"
  }
  
  if (is.numeric(beta_star)) {
    lambda_star <- lam_fun(beta_star)
  } else if (beta_star == "<0.75") {
    lambda_star <- lam_fun(beta_lo)
  } else if (beta_star == ">0.999") {
    lambda_star <- lam_fun(beta_hi)
  }
  
  list(beta_star = beta_star, lambda_star = lambda_star,
       r = r, c0 = c0, status = status)
}


## === Build type/label vectors (now including xerror) ===
# We now have 30 normal + 3 vertical + 3 structural + 3 xerror = 39 points
types <- c(rep("normal",   nrow(data_main)),     # 30 normal points
           rep("vertical", nrow(data_vert)),     # 3 vertical outliers
           rep("struct",   nrow(data_struct)),   # 3 structural outliers
           rep("xerror",   nrow(data_xerr)))     # 3 X-error outliers
labels <- c(data_main$label, data_vert$label, data_struct$label, data_xerr$label)  

## === Compute beta* and lambda* for each of the 39 future points ===
res_list <- lapply(seq_len(nrow(future_data)), function(i) {
  x0 <- future_data$X[i]; y0 <- future_data$Y[i]
  tmp <- beta_lambda_star_exact(x0, y0, alpha.hat, sigma.hat, invXTX, p, lam_fun)
  yh  <- yhat_one(x0, alpha.hat)
  data.frame(
    X = x0, Y = y0, yhat = yh, r = tmp$r,
    beta_star = tmp$beta_star, lambda_star = tmp$lambda_star,
    status = tmp$status,
    label = labels[i],     # 0=normal, 1=outlier
    type  = types[i],      # normal / vertical / struct / xerror
    stringsAsFactors = FALSE
  )
})
tab39 <- do.call(rbind, res_list)

## Convert beta_star into numeric for comparisons
tab39$beta_star_num <- suppressWarnings(as.numeric(tab39$beta_star))
tab39$beta_star_num[is.na(tab39$beta_star_num) & tab39$beta_star=="<0.75"]  <- 0.75
tab39$beta_star_num[is.na(tab39$beta_star_num) & tab39$beta_star==">0.999"] <- 1.00

## Outlier decision rule based on beta*
beta_crit <- 0.99
tab39$pred_outlier <- tab39$beta_star_num >= beta_crit

## Pretty-print table
num_beta <- suppressWarnings(as.numeric(tab39$beta_star))
tab39_out <- data.frame(
  X           = tab39$X,
  Y           = tab39$Y,
  yhat        = round(tab39$yhat, 6),
  r           = round(tab39$r, 6),
  beta_star   = ifelse(is.na(num_beta),
                       as.character(tab39$beta_star),
                       sprintf("%.6f", num_beta)),
  lambda_star = round(tab39$lambda_star, 6),
  status      = tab39$status,
  stringsAsFactors = FALSE
)
print(tab39_out)

## === Confusion matrix and performance metrics for beta* method ===
cm_beta <- with(tab39, table(pred_outlier, label))
TP_b <- sum(tab39$pred_outlier & tab39$label==1)  # True positives
FP_b <- sum(tab39$pred_outlier & tab39$label==0)  # False positives
TN_b <- sum(!tab39$pred_outlier & tab39$label==0) # True negatives
FN_b <- sum(!tab39$pred_outlier & tab39$label==1) # False negatives

prec_b <- ifelse((TP_b+FP_b)==0, NA, TP_b/(TP_b+FP_b))  # Precision
rec_b  <- ifelse((TP_b+FN_b)==0, NA, TP_b/(TP_b+FN_b))  # Recall
f1_b   <- ifelse(is.na(prec_b)|(prec_b+rec_b)==0, NA, 2*prec_b*rec_b/(prec_b+rec_b)) # F1

res_summary <- list(
  beta_rule_confusion = cm_beta,
  beta_rule_metrics   = c(precision=prec_b, recall=rec_b, F1=f1_b)
)
print(res_summary)

## === Per-type hit rates (normal / vertical / struct / xerror) ===
by_type_beta <- aggregate(pred_outlier ~ type, data = tab39,
                          FUN = function(z) c(hit = sum(z), total = length(z), rate = mean(z)))

cat("\n--- Per-type hit rates (β*) ---\n")
print(by_type_beta)