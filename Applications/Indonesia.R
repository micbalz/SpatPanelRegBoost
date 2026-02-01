### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, mlr3spatiotempcv, blockCV, sperrorest, knitr, kableExtra, Matrix, mboost, gamboostLSS, MASS, splm, plm, spdep, spData, sf, spatialreg, future, future.apply, progressr)

source("R/GSPECM.R")

set.seed(123456789)

### Load data
data(RiceFarms, package = "splm")
data(riceww)

r = length(unique(RiceFarms$time))
n = length(unique(RiceFarms$id))

RiceFarms = RiceFarms[,c("id", "time", "goutput", "size", "status", "varieties", "bimas", "seed", "urea",
                         "phosphate", "pesticide", "pseed", "purea", "pphosph", "hiredlabor", "famlabor", "totlabor",
                         "wage", "price")]

for (v in c("goutput", "size", "seed", "urea", "pseed", "purea", "pphosph", "hiredlabor", "famlabor", "totlabor", "wage", "price")) {
  RiceFarms[[v]] = log(RiceFarms[[v]])
}


for (v in c("phosphate", "pesticide")) {
  RiceFarms[[v]] = RiceFarms[[v]] / 1000
}

# for (v in c("goutput", "size", "seed", "urea", "pseed", "purea", "pphosph", "hiredlabor",
#             "famlabor", "totlabor", "wage", "price", "phosphate", "pesticide")) {
#   RiceFarms[[v]] = scale(RiceFarms[[v]], center = TRUE, scale = TRUE)
# }

RiceFarms$status     = as.factor(RiceFarms$status)
RiceFarms$varieties  = as.factor(RiceFarms$varieties)
RiceFarms$bimas      = as.factor(RiceFarms$bimas)

W = riceww

index = c("id", "time")
fml = as.formula(
  paste("goutput  ~ ", paste(c(colnames(RiceFarms[,-1:-3])), collapse = " + "))
)
X = data.frame(id = RiceFarms$id, time = RiceFarms$time, goutput = RiceFarms$goutput, model.matrix(fml, data = RiceFarms), check.names = FALSE)
X = pdata.frame(X, index)
index = X[,1]
tindex = X[,2]
ind  = index[which(names(index)%in%row.names(X))]
tind = tindex[which(names(index)%in%row.names(X))]
X = X[order(tind, ind), ]
X = as.data.frame(lapply(X, function(col) {
  if (inherits(col, "pseries")) as.numeric(col) else col
}))
Y = as.vector(X$goutput)
X = X %>% dplyr::select(!c(id, time, goutput))

WX = (diag(r) %x% W) %*% as.matrix(X[,-1])
WX = as.data.frame(WX)
colnames(WX) = paste0("W_", colnames(X[,-1]))
Z = cbind(X, WX)

# RiceFarms = RiceFarms %>%
#   arrange(id, time)
# 
# fml = as.formula(
#   paste("goutput  ~ ", paste(c(colnames(RiceFarms[,-1:-3])), collapse = " + "))
# )
# X = data.frame(model.matrix(fml, data = RiceFarms), check.names = FALSE)
# WX = as.matrix(diag(r) %x% W %*% as.matrix(X[,-1]))
# colnames(WX) = paste0("W_", colnames(X[,-1]))
# WX = as.data.frame(WX)
# Z = cbind(X,WX)
# Y = RiceFarms$goutput

cols = colnames(Z)[-1]
cols = toupper(cols)
colnames(Z) = c("(Intercept)", cols)

### Estimation
# (2) Generalized method of moments

### --- Step 1: Initial Estimation of beta ---
mod0 = lm(Y ~ ., data = Z[,-1])

# boost = glmboost(Y ~ ., data = Z[,-1], family = GStabs(stabilization = "MAD"),
#          center = TRUE, control = boost_control(mstop = 250000, trace = FALSE, nu = 0.1))
# coef(boost, off2int = TRUE)

### --- Step 2: GMM Estimation of rho2 and sigv ---
u_h = mod0$residuals
u_b = (diag(r) %x% W) %*% u_h
u_bb = (diag(r) %x% W) %*% u_b

J_b = matrix(1, nrow = r, ncol = r) / r
P = (J_b %x% diag(n))
Q = diag(n*r) - P

NT = n * (r - 1)

# Elements of G0
g11 = (2 / NT) * as.numeric(t(u_b) %*% Q %*% u_h)
g12 = -(1 / NT) * as.numeric(t(u_b) %*% Q %*% u_b)
g13 = 1

g21 = (2 / NT) * as.numeric(t(u_bb) %*% Q %*% u_b)
g22 = -(1 / NT) * as.numeric(t(u_bb) %*% Q %*% u_bb)
g23 = (1 / n) * sum(W * W)

g31 = (1 / NT) * (as.numeric(t(u_bb) %*% Q %*% u_h) + as.numeric(t(u_b) %*% Q %*% u_b))
g32 = -(1 / NT) * as.numeric(t(u_bb) %*% Q %*% u_b)
g33 = 0

# Assemble G0
G0 = matrix(c(
  g11, g12, g13,
  g21, g22, g23,
  g31, g32, g33
), nrow = 3, byrow = TRUE)

# Elements of g0
gam1 = (1 / NT) * as.numeric(t(u_h) %*% Q %*% u_h)
gam2 = (1 / NT) * as.numeric(t(u_b) %*% Q %*% u_b)
gam3 = (1 / NT) * as.numeric(t(u_b) %*% Q %*% u_h)

# Assemble g0
g0 = c(gam1, gam2, gam3)

scorr = as.numeric(crossprod((Diagonal(r) %x% W) %*% u_h, u_h) / crossprod(u_h, u_h))
scorr = scorr / (sum((Diagonal(r) %x% W)) / length(u_h))


# pars = c(scorr, crossprod(u_h, u_h) / (n*r))
pars = c(0,0.5)


obj = function(params) {
  rho2 = params[1]
  sigv = params[2]
  theta = c(rho2, rho2^2, sigv^2)
  res = G0 %*% theta - g0
  as.numeric(t(res) %*% res)
}

opt = nlminb(start = pars, objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
rho2 = opt$par[1]
sigv = opt$par[2]


### --- Step 2: GMM Estimation of rho1 and sigmu ---
S = P - (1 / (r-1) * Q)
#S = (J_b - (1 / (r-1) * (diag(r) - J_b))) %x% diag(n)

NT = n * r

# Elements of G1
g11 = (2 / NT) * as.numeric(t(u_b) %*% S %*% u_h)
g12 = -(1 / NT) * as.numeric(t(u_b) %*% S %*% u_b)
g13 = 1

g21 = (2 / NT) * as.numeric(t(u_bb) %*% S %*% u_b)
g22 = -(1 / NT) * as.numeric(t(u_bb) %*% S %*% u_bb)
g23 = (1 / n) * sum(W * W)

g31 = (1 / NT) * (as.numeric(t(u_bb) %*% S %*% u_h) + as.numeric(t(u_b) %*% S %*% u_b))
g32 = -(1 / NT) * as.numeric(t(u_bb) %*% S %*% u_b)
g33 = 0

# Assemble G1
G1 = matrix(c(
  g11, g12, g13,
  g21, g22, g23,
  g31, g32, g33
), nrow = 3, byrow = TRUE)

# Elements of g1
gam1 = (1 / NT) * as.numeric(t(u_h) %*% S %*% u_h)
gam2 = (1 / NT) * as.numeric(t(u_b) %*% S %*% u_b)
gam3 = (1 / NT) * as.numeric(t(u_b) %*% S %*% u_h)

# Assemble g1
g1 = c(gam1, gam2, gam3)

obj = function(params) {
  rho1 = params[1]
  sigmu = params[2]
  theta = c(rho1, rho1^2, sigmu^2)
  res = G1 %*% theta - g1
  as.numeric(t(res) %*% res)
}

opt = nlminb(start = pars, objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
rho1 = opt$par[1]
sigmu = opt$par[2]

sig1 = sqrt(r * sigmu**2 + sigv**2)

# (2.1) Generalized spatial panel data model (Random)
A = diag(n) - rho1 * W
B = diag(n) - rho2 * W
E = diag(r) - J_b

# phi_gls = sigmu**2 / sigv**2
# theta_gls = 1 - (sigv / sqrt(sigv**2 + r * sigmu^2))

Zgmm = as.matrix(Z)

omega_re = J_b %x% solve(r * sigmu**2 * solve(t(A) %*% A) + sigv**2 * solve(t(B) %*% B)) + 1 / sigv**2 * (E %x% (t(B) %*% B))

gls_re = as.vector(solve(t(Zgmm) %*% omega_re %*% Zgmm) %*% (t(Zgmm) %*% omega_re %*% Y))
names(gls_re) = colnames(Z)

# (2.2) Generalized spatial panel data model (Fixed)
Psi = as.matrix((E %x% (t(B) %*% B)))

gls_fix = as.vector(solve(t(Zgmm[,-1]) %*% Psi %*% Zgmm[,-1]) %*% (t(Zgmm[,-1]) %*% Psi %*% Y))
names(gls_fix) = colnames(Z[,-1])


### (3) Model-based gradient boosting (LTB)
# Generate folds 
flds = cvspat(n, r, Y, B = 10, type = "random", map = NULL)

# (3.1) Model-based gradient boosting (Random)
omega_sqrt_re = chol(omega_re)

Z_re = omega_sqrt_re %*% Zgmm
Y_re = omega_sqrt_re %*% Y
colnames(Z_re)[1] = "alpha"

mod_re = glmboost(Y_re ~ 0 + ., data = data.frame(Z_re, check.names = FALSE), family = GStabs(stabilization = "MAD"),
                      offset = 0, center = FALSE, control = boost_control(mstop = 5000, trace = FALSE, nu = 0.1))

cvr_re = cvrisk(mod_re, folds = flds)
gbm_re = coef(mod_re[mstop(cvr_re)], off2int = FALSE)

# (3.2) Model-based gradient boosting (Fixed)
Z_fix = as.matrix(Q %*% (diag(r) %x% B) %*% Zgmm[,-1])
Y_fix = as.vector(Q %*% (diag(r) %x% B) %*% Y)

mod_fix = glmboost(Y_fix ~ 0 + ., data = data.frame(Z_fix, check.names = FALSE), family = GStabs(stabilization = "MAD"),
                   center = FALSE, offset = 0, control = boost_control(mstop = 1000, trace = FALSE, nu = 0.1))

cvr_fix = cvrisk(mod_fix, folds = flds)
gbm_fix = coef(mod_fix[mstop(cvr_fix)], off2int = FALSE)


### (4) Deselection (DES)
# (4.1) Deselection (Random)
mod_des_re = DeselectBoost(mod_re, fam = Gaussian())

# (4.2) Deselection (Fixed)
mod_des_fix = DeselectBoostFE(mod_fix, fam = Gaussian())


### Results
gmm_re = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, gls_re)
gmm_fix = c(rho1 = NA, rho2 = rho2, sig2mu = NA, sig2v = sigv**2, gls_fix)

ltb_re = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, gbm_re)
ltb_fix = c(rho1 = NA, rho2 = rho2, sig2mu = NA, sig2v = sigv**2, gbm_fix)


des_re = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, coef(mod_des_re))
des_fix = c(rho1 = NA, rho2 = rho2, sig2mu = NA, sig2v = sigv**2, coef(mod_des_fix))




coefs_list = list(
  REGMM  = gmm_re,
  FESGMM = gmm_fix,
  RELTB = ltb_re,
  FELTB = ltb_fix,
  REDES = des_re,
  FEDES = des_fix
)

coefs_list = lapply(coefs_list, function(x) {
  nms = names(x)
  nms[nms == "alpha"] = "(Intercept)"
  names(x) = nms
  x
})

# All variable names used across models
all_vars = unique(unlist(lapply(coefs_list, names)))
all_vars = sort(all_vars)  # For consistent row order

# Create empty table (matrix)
coef_matrix = matrix(NA, nrow = length(all_vars), ncol = length(coefs_list),
                     dimnames = list(all_vars, names(coefs_list)))

# Fill the matrix
for (model in names(coefs_list)) {
  this_coef = coefs_list[[model]]
  coef_matrix[names(this_coef), model] = round(this_coef, 3)
}

# Convert to data.frame for nicer printing
coef_df = as.data.frame(coef_matrix)
coef_df = tibble::rownames_to_column(coef_df, var = "Variable")

# Custom variable ordering: lambda, non-W, W, sigma
non_w_vars = grep("^W_", coef_df$Variable, invert = TRUE, value = TRUE)
non_w_vars = setdiff(non_w_vars, c("rho1", "rho2", "sig2mu", "sig2v"))

w_vars = sort(grep("^W_", coef_df$Variable, value = TRUE))

ordered_vars = c("rho1", "rho2", "sig2mu", "sig2v", non_w_vars, w_vars)

coef_df = coef_df %>%
  dplyr::arrange(match(Variable, ordered_vars))


### Hausmann Test
vq = as.matrix(sigv^2* ginv(t(Zgmm[,-1]) %*% Psi %*% Zgmm[,-1]) - solve(t(Zgmm[,-1]) %*% omega_re %*% Zgmm[,-1]))
q = gls_fix - gls_re[-1]

m = as.numeric(t(q) %*% ginv(vq) %*% q)

# Degrees of freedom = length(q)
df = length(q)

# p-value
p_value = 1 - pchisq(m, df)

if (p_value < 0.05) {
  cat("Reject H0: RE is inconsistent, use FE\n")
} else {
  cat("Fail to reject H0: RE is consistent\n")
}



# Display the table nicely
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for Indonesian rice farms data set")









