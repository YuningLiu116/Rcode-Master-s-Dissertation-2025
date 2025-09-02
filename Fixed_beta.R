library(doParallel)
library(MASS)
library(DEoptim)

## ------------------ Basic Settings ------------------
Cbeta  = 0.95
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

lam_two_sti_exact=1.3741318

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
                  "STI (beta = 0.95, gamma = 0.99)",
                  "Original training points",
                  "Normal future points",
                  "Vertical outliers",
                  "Structural outliers",
                  "X-error outliers"),
       lty = c(1, 3, 2, NA, NA, NA, NA, NA),
       pch = c(NA, NA, NA, 16, 1, 4, 15, 17),
       col = c("green3", "black", "red", "black", "blue", "red", "purple", "orange"),
       cex = 1.05, bty = "n")

txt <- "beta gamma    lambda
1   0.95 0.750 0.9444427
2   0.95 0.751 0.9450531
3   0.95 0.752 0.9456635
4   0.95 0.753 0.9464264
5   0.95 0.754 0.9470367
6   0.95 0.755 0.9477997
7   0.95 0.756 0.9484100
8   0.95 0.757 0.9491730
9   0.95 0.758 0.9496307
10  0.95 0.759 0.9503937
11  0.95 0.760 0.9510040
12  0.95 0.761 0.9516144
13  0.95 0.762 0.9522247
14  0.95 0.763 0.9528351
15  0.95 0.764 0.9534454
16  0.95 0.765 0.9540558
17  0.95 0.766 0.9546661
18  0.95 0.767 0.9552765
19  0.95 0.768 0.9558868
20  0.95 0.769 0.9564972
21  0.95 0.770 0.9572601
22  0.95 0.771 0.9580231
23  0.95 0.772 0.9586334
24  0.95 0.773 0.9590912
25  0.95 0.774 0.9597015
26  0.95 0.775 0.9604645
27  0.95 0.776 0.9610748
28  0.95 0.777 0.9619904
29  0.95 0.778 0.9626007
30  0.95 0.779 0.9632111
31  0.95 0.780 0.9639740
32  0.95 0.781 0.9647369
33  0.95 0.782 0.9654999
34  0.95 0.783 0.9661102
35  0.95 0.784 0.9665680
36  0.95 0.785 0.9673309
37  0.95 0.786 0.9679413
38  0.95 0.787 0.9688568
39  0.95 0.788 0.9694672
40  0.95 0.789 0.9702301
41  0.95 0.790 0.9708405
42  0.95 0.791 0.9716034
43  0.95 0.792 0.9722137
44  0.95 0.793 0.9729767
45  0.95 0.794 0.9735870
46  0.95 0.795 0.9741974
47  0.95 0.796 0.9749603
48  0.95 0.797 0.9755707
49  0.95 0.798 0.9763336
50  0.95 0.799 0.9769440
51  0.95 0.800 0.9777069
52  0.95 0.801 0.9784698
53  0.95 0.802 0.9793854
54  0.95 0.803 0.9801483
55  0.95 0.804 0.9807587
56  0.95 0.805 0.9815216
57  0.95 0.806 0.9822845
58  0.95 0.807 0.9828949
59  0.95 0.808 0.9836871
60  0.95 0.809 0.9844208
61  0.95 0.810 0.9853363
62  0.95 0.811 0.9860992
63  0.95 0.812 0.9868622
64  0.95 0.813 0.9877777
65  0.95 0.814 0.9885406
66  0.95 0.815 0.9893036
67  0.95 0.816 0.9900946
68  0.95 0.817 0.9909821
69  0.95 0.818 0.9918976
70  0.95 0.819 0.9926605
71  0.95 0.820 0.9934235
72  0.95 0.821 0.9943663
73  0.95 0.822 0.9949765
74  0.95 0.823 0.9958649
75  0.95 0.824 0.9966278
76  0.95 0.825 0.9973907
77  0.95 0.826 0.9983063
78  0.95 0.827 0.9990692
79  0.95 0.828 0.9996796
80  0.95 0.829 1.0004425
81  0.95 0.830 1.0013580
82  0.95 0.831 1.0022736
83  0.95 0.832 1.0031891
84  0.95 0.833 1.0041046
85  0.95 0.834 1.0048676
86  0.95 0.835 1.0057831
87  0.95 0.836 1.0065460
88  0.95 0.837 1.0074615
89  0.95 0.838 1.0083771
90  0.95 0.839 1.0092926
91  0.95 0.840 1.0100555
92  0.95 0.841 1.0108185
93  0.95 0.842 1.0118866
94  0.95 0.843 1.0126495
95  0.95 0.844 1.0135651
96  0.95 0.845 1.0144806
97  0.95 0.846 1.0152435
98  0.95 0.847 1.0160065
99  0.95 0.848 1.0170746
100 0.95 0.849 1.0179901
101 0.95 0.850 1.0189056
102 0.95 0.851 1.0198212
103 0.95 0.852 1.0207367
104 0.95 0.853 1.0216522
105 0.95 0.854 1.0225900
106 0.95 0.855 1.0236359
107 0.95 0.856 1.0245514
108 0.95 0.857 1.0254887
109 0.95 0.858 1.0265350
110 0.95 0.859 1.0274721
111 0.95 0.860 1.0285400
112 0.95 0.861 1.0297394
113 0.95 0.862 1.0309601
114 0.95 0.863 1.0317230
115 0.95 0.864 1.0327911
116 0.95 0.865 1.0338593
117 0.95 0.866 1.0350800
118 0.95 0.867 1.0364532
119 0.95 0.868 1.0373688
120 0.95 0.869 1.0384369
121 0.95 0.870 1.0395050
122 0.95 0.871 1.0405731
123 0.95 0.872 1.0414886
124 0.95 0.873 1.0425568
125 0.95 0.874 1.0437775
126 0.95 0.875 1.0448456
127 0.95 0.876 1.0460663
128 0.95 0.877 1.0469818
129 0.95 0.878 1.0480499
130 0.95 0.879 1.0491180
131 0.95 0.880 1.0504913
132 0.95 0.881 1.0517120
133 0.95 0.882 1.0530853
134 0.95 0.883 1.0540009
135 0.95 0.884 1.0553741
136 0.95 0.885 1.0564423
137 0.95 0.886 1.0575104
138 0.95 0.887 1.0587311
139 0.95 0.888 1.0601044
140 0.95 0.889 1.0614777
141 0.95 0.890 1.0630035
142 0.95 0.891 1.0640717
143 0.95 0.892 1.0652924
144 0.95 0.893 1.0668182
145 0.95 0.894 1.0680551
146 0.95 0.895 1.0694122
147 0.95 0.896 1.0704803
148 0.95 0.897 1.0718694
149 0.95 0.898 1.0732269
150 0.95 0.899 1.0744630
151 0.95 0.900 1.0758209
152 0.95 0.901 1.0771942
153 0.95 0.902 1.0787350
154 0.95 0.903 1.0802460
155 0.95 0.904 1.0816339
156 0.95 0.905 1.0831451
157 0.95 0.906 1.0845184
158 0.95 0.907 1.0857533
159 0.95 0.908 1.0874176
160 0.95 0.909 1.0889435
161 0.95 0.910 1.0901779
162 0.95 0.911 1.0916901
163 0.95 0.912 1.0932294
164 0.95 0.913 1.0950470
165 0.95 0.914 1.0964203
166 0.95 0.915 1.0980988
167 0.95 0.916 1.0997772
168 0.95 0.917 1.1011505
169 0.95 0.918 1.1028290
170 0.95 0.919 1.1043549
171 0.95 0.920 1.1060455
172 0.95 0.921 1.1077118
173 0.95 0.922 1.1090970
174 0.95 0.923 1.1110687
175 0.95 0.924 1.1130524
176 0.95 0.925 1.1147308
177 0.95 0.926 1.1165619
178 0.95 0.927 1.1180878
179 0.95 0.928 1.1199188
180 0.95 0.929 1.1217499
181 0.95 0.930 1.1235809
182 0.95 0.931 1.1255646
183 0.95 0.932 1.1277008
184 0.95 0.933 1.1293793
185 0.95 0.934 1.1313629
186 0.95 0.935 1.1331940
187 0.95 0.936 1.1353302
188 0.95 0.937 1.1368657
189 0.95 0.938 1.1392975
190 0.95 0.939 1.1412811
191 0.95 0.940 1.1432648
192 0.95 0.941 1.1452484
193 0.95 0.942 1.1479950
194 0.95 0.943 1.1501312
195 0.95 0.944 1.1525726
196 0.95 0.945 1.1548615
197 0.95 0.946 1.1573029
198 0.95 0.947 1.1600494
199 0.95 0.948 1.1627960
200 0.95 0.949 1.1653900
201 0.95 0.950 1.1681366
202 0.95 0.951 1.1705780
203 0.95 0.952 1.1736371
204 0.95 0.953 1.1762238
205 0.95 0.954 1.1788177
206 0.95 0.955 1.1817169
207 0.95 0.956 1.1844702
208 0.95 0.957 1.1870641
209 0.95 0.958 1.1901093
210 0.95 0.959 1.1936188
211 0.95 0.960 1.1971283
212 0.95 0.961 1.2006378
213 0.95 0.962 1.2042999
214 0.95 0.963 1.2073573
215 0.95 0.964 1.2111664
216 0.95 0.965 1.2149811
217 0.95 0.966 1.2194113
218 0.95 0.967 1.2235310
219 0.95 0.968 1.2267303
220 0.95 0.969 1.2305450
221 0.95 0.970 1.2346649
222 0.95 0.971 1.2398529
223 0.95 0.972 1.2438245
224 0.95 0.973 1.2487030
225 0.95 0.974 1.2532846
226 0.95 0.975 1.2581635
227 0.95 0.976 1.2627411
228 0.95 0.977 1.2682378
229 0.95 0.978 1.2741852
230 0.95 0.979 1.2804413
231 0.95 0.980 1.2868500
232 0.95 0.981 1.2926483
233 0.95 0.982 1.2989072
234 0.95 0.983 1.3069942
235 0.95 0.984 1.3149310
236 0.95 0.985 1.3222527
237 0.95 0.986 1.3309500
238 0.95 0.987 1.3399525
239 0.95 0.988 1.3504810
240 0.95 0.989 1.3616197
241 0.95 0.990 1.3741318
242 0.95 0.991 1.3877133
243 0.95 0.992 1.4009857
244 0.95 0.993 1.4194510
245 0.95 0.994 1.4435577
246 0.95 0.995 1.4673615
247 0.95 0.996 1.4937634
248 0.95 0.997 1.5361809
249 0.95 0.998 1.5836334
250 0.95 0.999 1.6692384"

## === Lookup table for lambda as a function of gamma (beta fixed to Cbeta) ===
df <- read.table(text = txt, header = TRUE)
# Expect columns: beta, gamma, lambda, with beta == Cbeta for all rows
stopifnot(all(df$beta == Cbeta))
stopifnot(all(diff(df$gamma) > 0), all(diff(df$lambda) > 0))  # monotone increasing

# λ(γ) and its numeric inverse γ(λ)
lam_of_gamma   <- approxfun(df$gamma, df$lambda, method = "linear", rule = 2)
gamma_of_lambda<- approxfun(df$lambda, df$gamma, method = "linear", rule = 2)

# Helpers
yhat_one <- function(x, beta_hat) {
  p <- length(beta_hat)
  sum(beta_hat * c(1, x, x^2)[1:p])
}

## -------- gamma*-inverse (fix beta, invert for gamma) --------
gamma_lambda_star_exact <- function(x0, y0,
                                    alpha.hat, sigma.hat, invXTX, p,
                                    lam_of_gamma, gamma_grid = df$gamma,
                                    tol = 1e-9) {
  # Design vector and geometry term
  vx <- c(1, x0, x0^2)[1:p]
  v0 <- sqrt(as.numeric(t(vx) %*% invXTX %*% vx))
  c0 <- sqrt(p + 2) * v0
  
  # Residual magnitude at (x0, y0)
  r  <- abs(y0 - sum(alpha.hat * vx))
  
  # Required lambda to cover this point when beta is fixed to Cbeta
  lambda_req <- r / (sigma.hat * (zbeta + c0))
  
  # Boundaries in the table
  lam_min <- min(df$lambda); lam_max <- max(df$lambda)
  gam_min <- min(df$gamma);  gam_max <- max(df$gamma)
  
  if (lambda_req <= lam_min) {
    gamma_star <- sprintf("<%.3f", gam_min)   # needs less than the smallest gamma in table
    status <- "clamped_at_low"
    lambda_star <- lam_min
  } else if (lambda_req >= lam_max) {
    gamma_star <- sprintf(">%.3f", gam_max)   # needs more than the largest gamma in table
    status <- "clamped_at_high"
    lambda_star <- lam_max
  } else {
    # Invert λ(γ) numerically by interpolation
    gamma_star <- gamma_of_lambda(lambda_req)
    gamma_star <- round(gamma_star, 6)
    status <- "ok"
    lambda_star <- lam_of_gamma(gamma_star)
  }
  
  list(gamma_star = gamma_star, lambda_star = lambda_star,
       r = r, c0 = c0, status = status)
}

## === Build type/label vectors (now including xerror) ===
# 30 normal + 3 vertical + 3 structural + 3 xerror = 39 points
types <- c(rep("normal",   nrow(data_main)),
           rep("vertical", nrow(data_vert)),
           rep("struct",   nrow(data_struct)),
           rep("xerror",   nrow(data_xerr)))
labels <- c(data_main$label, data_vert$label, data_struct$label, data_xerr$label)

## === Compute gamma* and lambda* for each of the 39 future points ===
res_list <- lapply(seq_len(nrow(future_data)), function(i) {
  x0 <- future_data$X[i]; y0 <- future_data$Y[i]
  tmp <- gamma_lambda_star_exact(x0, y0, alpha.hat, sigma.hat, invXTX, p, lam_of_gamma)
  yh  <- yhat_one(x0, alpha.hat)
  data.frame(
    X = x0, Y = y0, yhat = yh, r = tmp$r,
    gamma_star = tmp$gamma_star, lambda_star = tmp$lambda_star,
    status = tmp$status,
    label = labels[i],      # 0=normal, 1=outlier
    type  = types[i],       # normal / vertical / struct / xerror
    stringsAsFactors = FALSE
  )
})
tab39_g <- do.call(rbind, res_list)

## Convert gamma_star to numeric for decisions
tab39_g$gamma_star_num <- suppressWarnings(as.numeric(tab39_g$gamma_star))
tab39_g$gamma_star_num[is.na(tab39_g$gamma_star_num) & grepl("^<", tab39_g$gamma_star)] <- min(df$gamma)
tab39_g$gamma_star_num[is.na(tab39_g$gamma_star_num) & grepl("^>", tab39_g$gamma_star)] <- max(df$gamma)

## Outlier decision rule based on gamma*
gamma_crit <- 0.99   # set your critical gamma here
tab39_g$pred_outlier <- tab39_g$gamma_star_num >= gamma_crit

## Pretty-print table
num_gamma <- suppressWarnings(as.numeric(tab39_g$gamma_star))
tab39_out <- data.frame(
  X           = tab39_g$X,
  Y           = tab39_g$Y,
  yhat        = round(tab39_g$yhat, 6),
  r           = round(tab39_g$r, 6),
  gamma_star  = ifelse(is.na(num_gamma),
                       as.character(tab39_g$gamma_star),
                       sprintf("%.6f", num_gamma)),
  lambda_star = round(tab39_g$lambda_star, 6),
  status      = tab39_g$status,
  stringsAsFactors = FALSE
)
print(tab39_out)

## === Confusion matrix and performance metrics (gamma*-rule) ===
cm_gamma <- with(tab39_g, table(pred_outlier, label))
TP_g <- sum(tab39_g$pred_outlier & tab39_g$label==1)   # True positives
FP_g <- sum(tab39_g$pred_outlier & tab39_g$label==0)   # False positives
TN_g <- sum(!tab39_g$pred_outlier & tab39_g$label==0)  # True negatives
FN_g <- sum(!tab39_g$pred_outlier & tab39_g$label==1)  # False negatives

prec_g <- ifelse((TP_g+FP_g)==0, NA, TP_g/(TP_g+FP_g))             # Precision
rec_g  <- ifelse((TP_g+FN_g)==0, NA, TP_g/(TP_g+FN_g))             # Recall
f1_g   <- ifelse(is.na(prec_g)|(prec_g+rec_g)==0, NA, 2*prec_g*rec_g/(prec_g+rec_g))  # F1

res_summary_gamma <- list(
  gamma_rule_confusion = cm_gamma,
  gamma_rule_metrics   = c(precision=prec_g, recall=rec_g, F1=f1_g)
)
print(res_summary_gamma)

## === Per-type hit rates for gamma*-rule (normal / vertical / struct / xerror) ===
by_type_gamma <- aggregate(pred_outlier ~ type, data = tab39_g,
                           FUN = function(z) c(hit = sum(z), total = length(z), rate = mean(z)))
cat("\n--- Per-type hit rates (γ*) ---\n")
print(by_type_gamma)