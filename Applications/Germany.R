### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, tmap, mlr3spatiotempcv, blockCV, sperrorest, knitr, kableExtra, Matrix, mboost, gamboostLSS, MASS, splm, plm, spdep, spData, sf, spatialreg, future, future.apply, progressr)

source("R/GSPECM.R")

set.seed(123456789)

## Load Data (INKAR version 2021 (inkar_2021.csv) is freely available at https://www.inkar.de/ (DL-DE BY 2.0))
inkar = fread("Applications/Data/inkar_2021.csv")
krs_spdf = st_read(dsn = "Applications/Data/vg5000_ebenen_0101", layer = "VG5000_KRS")

### Filter for Kreise and years of interest 
years = c(2014:2019)

krs = inkar %>%
  filter(Raumbezug == "Kreise",
         Zeitbezug %in% years) %>%      
  distinct(Name, Kennziffer, Zeitbezug, Indikator, .keep_all = TRUE) %>%  
  pivot_wider(
    id_cols = c(Name, Kennziffer, Zeitbezug),
    names_from = Indikator,
    values_from = Wert
  )


### Remove variables that are completely NA
krs = dplyr::select(krs, where(~ !all(is.na(.))))

### Clean IDs
krs = krs %>% 
  filter(Name != "Eisenach, Stadt") %>%
  mutate(Kennziffer = sprintf("%05d", Kennziffer))


### Spatial weight matrix (static)
knn = knearneigh(st_centroid(krs_spdf), k = 10)
nb  = knn2nb(knn, row.names = krs_spdf$AGS)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)

### Build the map
# tm_shape(st_as_sf(krs_spdf) %>%
#            left_join(krs, by = c("AGS" = "Kennziffer"))) +
#   tm_polygons("Lebenserwartung",
#               palette = viridis::viridis(100, direction = -1, option = "G"),
#               style = "kmeans") +
#   tm_facets(by = "Zeitbezug", ncol = 3) +
#   tm_layout(title = "")

### Select variables of interest
vrs = c("Name",
        "Zeitbezug",
        "Lebenserwartung",
        "Durchschnittsalter der Bevölkerung",
        "Arbeitslosenquote",
        "Beschäftigtenquote",
        "Erwerbsquote",
        "Selbständigenquote",
        "Schuldnerquote",
        "SGB II - Quote",
        "Beschäftigte mit akademischem Berufsabschluss", 
        "Eheschließungen",
        "Ehescheidungen",
        "Ausländeranteil",
        "Medianeinkommen",
        "Haushaltseinkommen",
        "Verbraucherinsolvenzverfahren",
        "Arbeitsvolumen",
        "Bodenfläche gesamt",
        "Wohnfläche",
        "Siedlungs- und Verkehrsfläche",
        "Waldfläche",
        "Erholungsfläche",
        "Wasserfläche",
        "Mietpreise",
        "Steuereinnahmen",
        "Bruttoinlandsprodukt je Einwohner",
        "Krankenhausversorgung",
        "Ärzte",
        "Pflegebedürftige",
        "Einwohnerdichte",
        "Pkw-Dichte",
        "Pendlersaldo",
        "Straßenverkehrsunfälle",
        "Getötete im Straßenverkehr")

Life = krs[, vrs]

### Rename columns
colnames(Life) = c( "NAME", "TIME", "LIFE","AGE","UNEMPLOYMENT","EMPLOYMENT","PART","SELF",
                    "DEBT","WELFARE","ACADEMICS","MARRIAGES","DIVORCES",
                    "FOREIGN","MEDINC","HHINC","INS","LABOR","LAND","LIVE",
                    "URBAN","FOREST","RECR","WATER","RENT","TAX","GDP",
                    "HOSP","DR","CARE","POP","CAR","COM","ACC","TRF")

Life = Life %>% dplyr::select(where(~ !any(is.na(.x))))

r = length(unique(Life$TIME))
n = length(unique(Life$NAME))

### Standardization
for (i in colnames(Life)[-1:-2]) {
  Life[[i]] = scale(Life[[i]], center = TRUE, scale = TRUE)
  # Life[[i]] = asinh(Life[[i]])
}

X = data.frame(Life, check.names = FALSE)
X = pdata.frame(X, NULL)
index = X[,1]
tindex = X[,2]
ind  = index[which(names(index) %in% row.names(X))]
tind = tindex[which(names(index) %in% row.names(X))]
X = X[order(tind, ind), ]
X = data.frame("(Intercept)" = 1, X, check.names = FALSE)
Y = X$LIFE
X = X %>% dplyr::select(!c(NAME, TIME, LIFE))

WX = as.matrix(diag(r) %x% W %*% as.matrix(X[,-1]))
colnames(WX) = paste0("W_", colnames(X[,-1]))
WX = as.data.frame(WX)
Z = cbind(X,WX)

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
flds = cvspat(n, r, Y, B = 10, type = "spatial", map = krs_spdf)

# (3.1) Model-based gradient boosting (GSPECM)
omega_sqrt_gspecm = chol(omega_gspecm)

Z_gspecm = omega_sqrt_gspecm %*% Zgmm
Y_gspecm = omega_sqrt_gspecm %*% Y
colnames(Z_gspecm)[1] = "alpha"

mod_gspecm = glmboost(Y_gspecm ~ 0 + ., data = data.frame(Z_gspecm, check.names = FALSE), family = GStabs(stabilization = "MAD"),
                      offset = 0, center = FALSE, control = boost_control(mstop = 2000, trace = FALSE, nu = 0.1))

cvr_gspecm = cvrisk(mod_gspecm, folds = flds)
gbm_gspecm = coef(mod_gspecm[mstop(cvr_gspecm)], off2int = FALSE)

# (3.2) Model-based gradient boosting (ANS)
omega_sqrt_ans = chol(omega_ans)

Z_ans = omega_sqrt_ans %*% Zgmm
Y_ans = omega_sqrt_ans %*% Y
colnames(Z_ans)[1] = "alpha"

mod_ans = glmboost(Y_ans ~ 0 + ., data = data.frame(Z_ans, check.names = FALSE), family = GStabs(stabilization = "MAD"),
                   offset = 0, center = FALSE, control = boost_control(mstop = 2000, trace = FALSE, nu = 0.1))

cvr_ans = cvrisk(mod_ans, folds = flds)
gbm_ans = coef(mod_ans[mstop(cvr_ans)], off2int = FALSE)

# (3.3) Model-based gradient boosting (KKP)
omega_sqrt_kkp = chol(omega_kkp)

Z_kkp = omega_sqrt_kkp %*% Zgmm
Y_kkp = omega_sqrt_kkp %*% Y
colnames(Z_kkp)[1] = "alpha"

mod_kkp = glmboost(Y_kkp ~ 0 + ., data = data.frame(Z_kkp, check.names = FALSE), family = GStabs(stabilization = "MAD"),
                   offset = 0, center = FALSE, control = boost_control(mstop = 2000, trace = FALSE, nu = 0.1))

cvr_kkp = cvrisk(mod_kkp, folds = flds)
gbm_kkp = coef(mod_kkp[mstop(cvr_kkp)], off2int = FALSE)


### (4) Deselection (DES)
# (4.1) Deselection (GSPECM)
mod_des_gespecm = DeselectBoost(mod_gspecm, fam = GStabs(stabilization = "MAD"))

# (4.2) Deselection (ANS)
mod_des_ans = DeselectBoost(mod_ans, fam = GStabs(stabilization = "MAD"))

# (4.3) Deselection (KKP)
mod_des_kpp = DeselectBoost(mod_kkp, fam = GStabs(stabilization = "MAD"))

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
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for Germany life expectancy data set")









