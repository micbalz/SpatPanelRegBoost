### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, mlr3spatiotempcv, blockCV, sperrorest, knitr, kableExtra, Matrix, mboost, gamboostLSS, MASS, splm, plm, spdep, spData, sf, spatialreg, future, future.apply, progressr)

plan(multisession, workers = parallel::detectCores() - 1)  

source("R/GSPECM.R")

handlers(global = TRUE)
handlers("cli")

### Simulation setup
nsim = 100
n = 100
r = 5

beta_t = c(1, 3.5, -2.5, rep(0,398))
names(beta_t) = c("alpha", paste0("X", 1:((length(beta_t) - 1))))
gamma_t = c(-4, 3, rep(0,398))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sig2mu_t = 10
sig2v_t = 10
# theta_t = sig2mu_t / (sig2mu_t + sig2v_t)

### Spatial weight matrix
ncsids = st_read(system.file("shapes/sids.gpkg", package = "spData")[1])
knn = knearneigh(st_centroid(ncsids), k = 5)
nb = knn2nb(knn, row.names = ncsids$NAME)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)

### Main simulation function
run = function(n, r, rho1_t, rho2_t, beta_t, gamma_t, sig2mu_t, sig2v_t, W) {
  ### Generate data
  # Generate covariates 
  p = length(beta_t) + length(gamma_t) - 1
  p_true = sum(beta_t != 0) + sum(gamma_t != 0) - 1
  
  zeta = matrix(runif(n * (p/2), -7.5, 7.5), nrow = n, ncol = (p/2))
  zeta = do.call(rbind, replicate(r, zeta, simplify = FALSE))
  nuta = matrix(runif(n * r * (p/2), -5, 5), nrow = n * r, ncol = (p/2))
  X = zeta + nuta
  Z = as.matrix(cbind(1, X, (diag(r) %x% W) %*% X))
  colnames(Z) = c(names(beta_t), names(gamma_t))
  
  # Generate error
  mu = rnorm(n, 0, sd = sqrt(sig2mu_t))
  v = rnorm(n*r, mean = 0, sd = sqrt(sig2v_t))
  Z_mu = rep(1, times = r) %x% diag(n)
  
  u1 = solve(diag(n) - rho1_t * W) %*% mu
  u2 = solve(diag(n*r) - rho2_t * (diag(r) %x% W)) %*% v
  u = Z_mu %*% u1 + u2
  
  # Generate Y
  Y = as.vector(Z %*% c(beta_t, gamma_t) + u)
  
  ### Estimation
  # Generate folds 
  flds = cvspat(n, r, Y, B = 5, type = "spatial", map = ncsids)
  
  # (1) Estimate the beta by GBM and DES
  Zgmm = as.data.frame(Z)
  mod0 = glmboost(Y ~ 0 + ., data = Zgmm, family = Gaussian(), center = FALSE, offset = 0,
                  control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))
  cvr0 = cvrisk(mod0, folds = flds)
  mod0 = mod0[mstop(cvr0)]
  mod0 = DeselectBoost(mod0, fam = Gaussian())
  
  # (2) GMM estimator of rho2 and sigv
  u_h = as.vector(mod0$resid())
  u_b = (diag(r) %x% W) %*% u_h
  u_bb = (diag(r) %x% W) %*% u_b
  
  J_b = matrix(1, nrow = r, ncol = r) / r
  
  P = (J_b %x% diag(n))
  Q = diag(n*r) - P
  
  NT = n * (r - 1)
  
  # Elements of G0
  g11 = (2 / NT) * t(u_b) %*% Q %*% u_h
  g12 = -(1 / NT) * t(u_b) %*% Q %*% u_b
  g13 = 1
  
  g21 = (2 / NT) * t(u_bb) %*% Q %*% u_b
  g22 = -(1 / NT) * t(u_bb) %*% Q %*% u_bb
  g23 = (1 / n) * sum(W * W)
  
  g31 = (1 / NT) * (t(u_bb) %*% Q %*% u_h + t(u_b) %*% Q %*% u_b)
  g32 = -(1 / NT) * t(u_bb) %*% Q %*% u_b
  g33 = 0
  
  # Assemble G0
  G0 = matrix(c(
    g11, g12, g13,
    g21, g22, g23,
    g31, g32, g33
  ), nrow = 3, byrow = TRUE)
  
  # Elements of g0
  gam1 = (1 / NT) * t(u_h) %*% Q %*% u_h
  gam2 = (1 / NT) * t(u_b) %*% Q %*% u_b
  gam3 = (1 / NT) * t(u_b) %*% Q %*% u_h
  
  # Assemble g0
  g0 = c(gam1, gam2, gam3)
  
  obj = function(params) {
    rho2 = params[1]
    sigv = params[2]
    theta = c(rho2, rho2^2, sigv^2)
    res = G0 %*% theta - g0
    t(res) %*% res
  }
  
  # Starting values
  # scorr = as.numeric(crossprod((Diagonal(r) %x% W) %*% u_h, u_h) / crossprod(u_h, u_h))
  # scorr = scorr / (sum((Diagonal(r) %x% W)) / length(u_h))
  # 
  # 
  # pars = c(scorr, crossprod(u_h, u_h) / (n*r))
  
  pars = c(0,0.5)
  
  opt = nlminb(start = pars, objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
  rho2 = opt$par[1]
  sigv = opt$par[2]
  
  
  # (3) GMM estimator of rho1 and sigmu
  S = P - (1 / (r-1) * Q)
  # S = (J_b - (1 / (r-1) * (Diagonal(r) - J_b))) %x% Diagonal(n)
  
  NT = n * r
  
  # Elements of G1
  g11 = (2 / NT) * t(u_b) %*% S %*% u_h
  g12 = -(1 / NT) * t(u_b) %*% S %*% u_b
  g13 = 1
  
  g21 = (2 / NT) * t(u_bb) %*% S %*% u_b
  g22 = -(1 / NT) * t(u_bb) %*% S %*% u_bb
  g23 = (1 / n) * sum(W * W)
  
  g31 = (1 / NT) * (t(u_bb) %*% S %*% u_h + t(u_b) %*% S %*% u_b)
  g32 = -(1 / NT) * t(u_bb) %*% S %*% u_b
  g33 = 0
  
  # Assemble G1
  G1 = matrix(c(
    g11, g12, g13,
    g21, g22, g23,
    g31, g32, g33
  ), nrow = 3, byrow = TRUE)
  
  # Elements of g1
  gam1 = (1 / NT) * t(u_h) %*% S %*% u_h
  gam2 = (1 / NT) * t(u_b) %*% S %*% u_b
  gam3 = (1 / NT) * t(u_b) %*% S %*% u_h
  
  # Assemble g1
  g1 = c(gam1, gam2, gam3)
  
  obj = function(params) {
    rho1 = params[1]
    sigmu = params[2]
    theta = c(rho1, rho1^2, sigmu^2)
    res = G1 %*% theta - g1
    t(res) %*% res
  }
  
  # Starting values
  # scorr = as.numeric(crossprod((diag(r) %x% W) %*% S %*% u_h, S %*% u_h) / crossprod(S %*% u_h, S %*% u_h))
  # 
  # scorr = scorr / (sum((diag(r) %x% W)) / length(S %*% u_h))
  # pars = c(scorr, scale(crossprod(S %*% u_h, S %*% u_h)) / (n*r))
  
  pars = c(0,0.5)
  
  opt = nlminb(start = pars, objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
  rho1 = opt$par[1]
  sigmu = opt$par[2]
  
  # Compute omega 
  A = diag(n) - rho1 * W
  B = diag(n) - rho2 * W
  E = diag(r) - J_b
  
  omega = J_b %x% solve(r * sigmu**2 * solve(t(A) %*% A) + sigv**2 * solve(t(B) %*% B)) + 1 / sigv**2 * (E %x% (t(B) %*% B))
  
  # (1) GBM (Random)
  omega_sqrt = chol(omega)
  
  Z_gspecm = omega_sqrt %*% Z
  Y_gspecm = omega_sqrt %*% Y
  
  gspecm = glmboost(Y_gspecm ~ 0 + ., data = data.frame(Z_gspecm, check.names = FALSE), family = Gaussian(),
                    offset = 0, center = FALSE, control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))
  
  cvr = cvrisk(gspecm, folds = flds)
  gspecmB = coef(gspecm[mstop(cvr)], off2int = FALSE)
  
  # (3.2) GBM (Fixed)
  Z_fix = Q %*% (diag(r) %x% B) %*% Z[,-1]
  Y_fix = Q %*% (diag(r) %x% B) %*% Y
  
  gspecm_fix = glmboost(Y_fix ~ 0 + ., data = data.frame(Z_fix, check.names = FALSE), family = Gaussian(), offset = 0, center = FALSE,
                        control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))
  
  cvrfix = cvrisk(gspecm_fix, folds = flds)
  fixB = coef(gspecm_fix[mstop(cvrfix)], off2int = FALSE)
  
  ### (4) Deselection
  # (4.1) Deselection (Random)
  mod_des = DeselectBoost(gspecm, fam = Gaussian())
  des = coef(mod_des)
  
  # (4.2) Deselection (Fixed)
  mod_desfix = DeselectBoostFE(gspecm_fix, fam = Gaussian())
  desfix = coef(mod_desfix)
  
  # Allocate all coefficients into a list
  mods = list(
    REGBM = gspecmB,
    REDES = des,
    FEGBM = fixB,
    FEDES = desfix
  )
  
  ### Performance criteria
  # (5) Variable Selection
  nameVar = colnames(Z[,-1])[1:p]
  trueVar = nameVar[c(1:(p_true/2), length(beta_t):(length(beta_t) + p_true/2 - 1))]
  falseVar = nameVar[!nameVar %in% trueVar]
  
  # extract selected variable names from each model
  selected = lapply(mods, \(v) setdiff(names(v), c("alpha", "(Intercept)")))
  
  # compute metrics in one vectorized call
  metrics = t(sapply(selected, function(sel) c(
    TPR = sum(trueVar %in% sel) / length(trueVar),
    TNR = 1 - sum(falseVar %in% sel) / length(falseVar)
  )))
  
  # (6) Estimation
  SE = sapply(mods, function(coef) {
    mts = intersect(names(coef), names(c(beta_t, gamma_t)))
    mts = mts[!mts %in% c("alpha", "(Intercept)")]
    
    sum((c(beta_t, gamma_t)[mts] - coef[mts])^2)
  })
  
  results = data.frame(
    Model = names(mods),
    TPR   = metrics[names(mods), "TPR"],
    TNR   = metrics[names(mods), "TNR"],
    MSE   = SE[names(mods)]
  )
  rownames(results) = NULL
  
  
  return(results)
}


sims = function(param_pairs, nsim) {
  set.seed(123456789)
  
  result = list()
  
  for (i in 1:nrow(param_pairs)) {
    rho1 = param_pairs$rho1[i]
    rho2 = param_pairs$rho2[i]
    
    cat("Running simulation study for rho1 =", rho1, "and rho2 =", rho2, "\n")
    
    pb = progressor(along = 1:nsim)
    
    res_list = future_lapply(1:nsim, function(j) {
      res = run(n, r, rho1, rho2, beta_t, gamma_t, sig2mu_t, sig2v_t, W)
      pb(sprintf("rho1=%.2f, rho2=%.2f, replication %d", rho1, rho2, j))
      res
    }, future.seed = TRUE)
    
    result[[paste0("rho1=", rho1, "_rho2=", rho2)]] = res_list
  }
  
  return(result)
}

param_pairs = data.frame(
  rho1 = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8),
  rho2 = c(0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8)
)

results = sims(param_pairs, nsim)

agg = lapply(results, function(simlist) {
  rbindlist(simlist)                    
}) |> lapply(function(dt) {
  setDT(dt)[, lapply(.SD, mean), by = Model]   
})

agg = lapply(agg, function(dt) {
  dt[] = lapply(dt, function(col) if(is.numeric(col)) round(col, 3) else col)
  dt
})

agg

long = rbindlist(
  lapply(names(agg), function(nm) {
    tmp <- as.data.table(agg[[nm]])
    tmp[, c("rho1", "rho2") := tstrsplit(nm, "_", fixed = TRUE)]
    tmp[, rho1 := sub("rho1=", "", rho1)]
    tmp[, rho2 := sub("rho2=", "", rho2)]
    tmp
  })
)

long = long %>%
  mutate(
    group = case_when(
      Model %in% c("REGBM", "REDES") ~ "Random",
      Model %in% c("FEGBM", "FEDES") ~ "Fixed"
    ),
    method = case_when(
      Model %in% c("REGBM", "FEGBM") ~ "LTB",
      Model %in% c("REDES", "FEDES") ~ "DES"
    )
  )

long = long %>%
  pivot_longer(
    cols = c(TPR, TNR, MSE),
    names_to = "Metric",
    values_to = "Value"
  )

fin = long %>%
  mutate(col = paste(group, method)) %>%
  dplyr::select(rho1, rho2, Metric, col, Value) %>%
  pivot_wider(
    names_from = col,
    values_from = Value
  ) %>%
  arrange(rho1, rho2, factor(Metric, levels = c("TPR", "TNR", "MSE")))


kable(
  fin,
  digits  = 3,
  align   = "cccccccccc",
  caption = "Performance metrics for all combinations of spatial parameters (rho1, rho2)"
)
