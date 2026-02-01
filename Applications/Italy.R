### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, mlr3spatiotempcv, blockCV, sperrorest, knitr, kableExtra, Matrix, mboost, gamboostLSS, MASS, splm, plm, spdep, spData, sf, spatialreg, future, future.apply, progressr)

source("R/GSPECM.R")

set.seed(123456789)

### Load data
data(Insurance)
data(itaww)

itamap = st_read("Applications/Data/Prov01012002/Prov01012002_g_WGS84.shp")

r = length(unique(Insurance$year))
n = length(unique(Insurance$code))

### Keep variables in correct order
Insurance = Insurance[, c(
  "code", "year",
  "ppcd", "rgdp", "bank", "den", "rirs",
  "agen", "school", "vaagr",
  "fam", "inef", "trust"
  #,"d99", "d00", "d01", "d02"
)]

### Center and scale variables
# for (v in c("ppcd", "rgdp", "bank", "den", "agen", "fam", "inef", "trust")) {
#   Insurance[[v]] = log(Insurance[[v]])
# }
# 
# 
# for (v in c("rirs", "school", "vaagr")) {
#   Insurance[[v]] = Insurance[[v]] / 100
# }

for (v in c("ppcd", "rgdp", "bank", "den", "agen", "fam", "inef", "trust", "rirs", "school", "vaagr")) {
  Insurance[[v]] = scale(Insurance[[v]], center = TRUE, scale = TRUE)
}

W = itaww

X = data.frame(Insurance, check.names = FALSE)
WX = as.data.frame((diag(r) %x% W) %*% as.matrix(X[,-1:-3]))
colnames(WX) = paste0("W_", colnames(X[,-1:-3]))
Z = cbind(X, WX)


cols = colnames(Z)
mts = match("year", colnames(Z))
cols[(mts + 1):length(cols)] = toupper(cols[(mts + 1):length(cols)])
colnames(Z) = cols


### Estimation
# (1) Maximum likelihood
# fml = as.formula(
#   paste("PPCD ~ ", paste(c(colnames(Z[,-1:-3])), collapse = " + "))
# )
# 
# mod_mle = spreml(fml, 
#              data = Z,
#              w = itaww, 
#              errors = "semgre",
#              lag = FALSE,
#              quiet = TRUE,
#              method = "nlminb")
#
# mod_ans = spreml(fml, 
#                  data = Z,
#                  w = itaww, 
#                  errors = "semre",
#                  lag = FALSE,
#                  quiet = TRUE,
#                  method = "nlminb")
# 
# mod_kkp = spreml(fml, 
#                  data = Z,
#                  w = itaww, 
#                  errors = "sem2re",
#                  lag = FALSE,
#                  quiet = TRUE,
#                  method = "nlminb")
# 
# mod_re = spreml(fml, 
#                  data = Z,
#                  w = itaww, 
#                  errors = "re",
#                  lag = FALSE,
#                  quiet = TRUE,
#                  method = "nlminb")
# 
# mod_gmm = spgm(fml,
#              data = Z,
#              listw = itaww,
#              Durbin = FALSE,
#              model = "random",
#              method = "w2sls",
#              moments = "initial")


# (2) Generalized method of moments
Z = data.frame("(Intercept)" = 1, Z, check.names = FALSE)
Z = Z %>%
  arrange(year, code)
Y = Z$PPCD
Z = Z %>% dplyr::select(!c(code, year, PPCD))

### --- Step 1: Initial Estimation of beta ---
mod0 = lm(Y ~ ., data = Z[,-1])

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

omega_gspecm = J_b %x% solve(r * sigmu**2 * solve(t(A) %*% A) + sigv**2 * solve(t(B) %*% B)) + 1 / sigv**2 * (E %x% (t(B) %*% B))

gls_gspecm = as.vector(solve(t(Zgmm) %*% omega_gspecm %*% Zgmm) %*% (t(Zgmm) %*% omega_gspecm %*% Y))
names(gls_gspecm) = colnames(Z)

# (2.2) Anselin (ANS)
sig21_ans = as.numeric((1 / n) * t(u_h - rho2*u_b) %*% P %*% (u_h - rho2*u_b))
sig2mu_ans = (sig21_ans - sigv**2) / r

omega_ans = (1/sig21_ans * J_b + (1 / sigv**2) * E) %x% (t(B) %*% B)

gls_ans = as.vector(solve(t(Zgmm) %*% omega_ans %*% Zgmm) %*% (t(Zgmm) %*% omega_ans %*% Y))
names(gls_ans) = colnames(Z)

# (2.3) Kapoor-Keleijan-Prucha (KKP)
sig2mu_kkp = (1 / NT) * as.numeric(t(u_h) %*% S %*% u_h)

omega_kkp = J_b %x% solve(r * sig2mu_kkp * diag(n) + sigv**2 * solve(t(B) %*% B)) + 1 / sigv**2 * (E %x% (t(B) %*% B))

gls_kkp = as.vector(solve(t(Zgmm) %*% omega_kkp %*% Zgmm) %*% (t(Zgmm) %*% omega_kkp %*% Y))
names(gls_kkp) = colnames(Z)

### (3) Model-based gradient boosting (LTB)
# Generate folds 
flds = cvspat(n, r, Y, B = 10, type = "spatial", map = itamap)

# (3.1) Model-based gradient boosting (GSPECM)
omega_sqrt_gspecm = chol(omega_gspecm)

Z_gspecm = omega_sqrt_gspecm %*% Zgmm
Y_gspecm = omega_sqrt_gspecm %*% Y
colnames(Z_gspecm)[1] = "alpha"

mod_gspecm = glmboost(Y_gspecm ~ 0 + ., data = data.frame(Z_gspecm, check.names = FALSE), family = Gaussian(),
                  offset = 0, center = FALSE, control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))

cvr_gspecm = cvrisk(mod_gspecm, folds = flds)
gbm_gspecm = coef(mod_gspecm[mstop(cvr_gspecm)], off2int = FALSE)

# (3.2) Model-based gradient boosting (ANS)
omega_sqrt_ans = chol(omega_ans)

Z_ans = omega_sqrt_ans %*% Zgmm
Y_ans = omega_sqrt_ans %*% Y
colnames(Z_ans)[1] = "alpha"

mod_ans = glmboost(Y_ans ~ 0 + ., data = data.frame(Z_ans, check.names = FALSE), family = Gaussian(),
                      offset = 0, center = FALSE, control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))

cvr_ans = cvrisk(mod_ans, folds = flds)
gbm_ans = coef(mod_ans[mstop(cvr_ans)], off2int = FALSE)

# (3.3) Model-based gradient boosting (KKP)
omega_sqrt_kkp = chol(omega_kkp)

Z_kkp = omega_sqrt_kkp %*% Zgmm
Y_kkp = omega_sqrt_kkp %*% Y
colnames(Z_kkp)[1] = "alpha"

mod_kkp = glmboost(Y_kkp ~ 0 + ., data = data.frame(Z_kkp, check.names = FALSE), family = Gaussian(),
                   offset = 0, center = FALSE, control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))

cvr_kkp = cvrisk(mod_kkp, folds = flds)
gbm_kkp = coef(mod_kkp[mstop(cvr_kkp)], off2int = FALSE)


### (4) Deselection (DES)
# (4.1) Deselection (GSPECM)
mod_des_gespecm = DeselectBoost(mod_gspecm, fam = Gaussian())

# (4.2) Deselection (ANS)
mod_des_ans = DeselectBoost(mod_ans, fam = Gaussian())

# (4.3) Deselection (KKP)
mod_des_kpp = DeselectBoost(mod_kkp, fam = Gaussian())

### Results
gmm_gspecm = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, gls_gspecm)
gmm_ans = c(rho1 = 0, rho2 = rho2, sig2mu = sig2mu_ans, sig2v = sigv**2, gls_ans)
gmm_kkp = c(rho1 = rho2, rho2 = rho2, sig2mu = sig2mu_kkp, sig2v = sigv**2, gls_kkp)

ltb_gspecm = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, gbm_gspecm)
ltb_ans = c(rho1 = 0, rho2 = rho2, sig2mu = sig2mu_ans, sig2v = sigv**2, gbm_ans)
ltb_kkp = c(rho1 = rho2, rho2 = rho2, sig2mu = sig2mu_kkp, sig2v = sigv**2, gbm_kkp)

des_gspecm = c(rho1 = rho1, rho2 = rho2, sig2mu = sigmu**2, sig2v = sigv**2, coef(mod_des_gespecm))
des_ans = c(rho1 = 0, rho2 = rho2, sig2mu = sig2mu_ans, sig2v = sigv**2, coef(mod_des_ans))
des_kkp = c(rho1 = rho2, rho2 = rho2, sig2mu = sig2mu_kkp, sig2v = sigv**2, coef(mod_des_kpp))



coefs_list = list(
  GSPECMGMM  = gmm_gspecm,
  ANSGMM = gmm_ans,
  KKPGMM = gmm_kkp,
  GSPECMLTB = ltb_gspecm,
  ANSLTB = ltb_ans,
  KKPLTB = ltb_kkp,
  GSPECMDES = des_gspecm,
  ANSDES = des_ans,
  KKPDES = des_kkp
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


# Display the table nicely
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for Italian non-life insurance data set")





