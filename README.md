# Gradient Boosting for Spatial Panel Models with Random and Fixed Effects

This repository provides access to three-step feasible model-based gradient boosting for spatial regression models with autoregressive disturbances.
It includes:  

- Helper functions for model-based gradient boosting   
- Monte Carlo experiments under varying spatial autocorrelation parametes 
- Estimates for application settings including life expectancy in German districts.

The repository serves as a foundation for replication.  

## Technical Details  

For in-depth derivations and explanations of model-based gradient boosting for spatial panel models with random and fixed effects, refer to: 

tba.

## Example 
```
require(Matrix)
require(mboost)
require(spData)

set.seed(123456789)

# Simulate artificial data
nsim = 100
n = 100
r = 5

rho1_t = 0.4
rho2_t = -0.4
beta_t = c(1, 3.5, -2.5, rep(0,18))
names(beta_t) = c("alpha", paste0("X", 1:((length(beta_t) - 1))))
gamma_t = c(-4, 3, rep(0,18))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sig2mu_t = 10
sig2v_t = 10

# Spatial weight matrix
ncsids = st_read(system.file("shapes/sids.gpkg", package = "spData")[1])
knn = knearneigh(st_centroid(ncsids), k = 5)
nb = knn2nb(knn, row.names = ncsids$NAME)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)

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

Zgmm = as.data.frame(Z)

# (1) Estimate the beta by OLS
mod0 = lm(Y ~ ., data = Zgmm[-1])

# (2) GMM estimator of rho2 and sigv
u_h = mod0$residuals
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

opt = nlminb(start = c(0,0.5), objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
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

opt = nlminb(start = c(0,0.5), objective = obj, lower = c(-0.999, 0), upper = c(0.999, Inf))
rho1 = opt$par[1]
sigmu = opt$par[2]

# Model-based gradient boosting
A = diag(n) - rho1 * W
B = diag(n) - rho2 * W
E = diag(r) - J_b
omega = J_b %x% solve(r * sigmu**2 * solve(t(A) %*% A) + sigv**2 * solve(t(B) %*% B)) + 1 / sigv**2 * (E %x% (t(B) %*% B))

omega_sqrt = chol(omega)
  
Z_gspecm = omega_sqrt %*% Z
Y_gspecm = omega_sqrt %*% Y
  
gspecm = glmboost(Y_gspecm ~ 0 + ., data = data.frame(Z_gspecm, check.names = FALSE), family = Gaussian(),
                  offset = 0, center = FALSE, control = boost_control(mstop = 500, trace = FALSE, nu = 0.1))
  
coef(gspecm[200], off2int = FALSE)

par(mar = c(5, 4, 4, 6))  
plot(gspecm)
```
