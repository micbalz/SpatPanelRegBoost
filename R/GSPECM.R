GStabs = function(stabilization = c("none", "MAD", "L2")) {
  
  Family(
    ngradient = function(y, f, w = 1) {
      ngr = (y - f)
      ngr = stabilize_ngradient(ngr, w, stabilization)
      return(ngr)
    },
    loss = function(y, f) {
      (y - f)^2 
    },
    offset = weighted.mean,
    check_y = function(y) {
      if (!is.numeric(y) || !is.null(dim(y)))
        stop("response is not a numeric vector but ", sQuote("family = SEM"))
      y
    },
    name = "Stabilized L2 loss",
    fW = function(f) rep(1, length(f)),
    response = function(f) f
  )
}

stabilize_ngradient = function(ngr, w = 1, stabilization) {
  
  if (stabilization == "none") {
    ngr = ngr
  }
  
  if (stabilization == "MAD") {
    div = weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                          w = w, na.rm = TRUE)
    div = ifelse(div < 0.0001, 0.0001, div)
    ngr = ngr / div
  }
  
  if (stabilization == "L2") {
    div = sqrt(weighted.mean(ngr^2, w = w,  na.rm = TRUE))
    div = ifelse(div < 1e-04, 1e-04, div)
    div = ifelse(div > 1e+04, 1e+04, div)
    ngr = ngr / div
  }
  
  ngr
}

cvspat = function(n, r, Y = NULL, B, type = c("random", "spatial", "blocks", "lto", "llo", "llto"), map = NULL) {
  
  if(type == "random") {
    
    folds_list = vector("list", r)
    
    for (rep in seq_len(r)) {
      
      # Randomly shuffle indices of units
      idx = sample.int(n)
      
      # Assign fold numbers (roughly equal sized)
      fold_id = cut(
        seq_len(n),
        breaks = B,
        labels = FALSE
      )
      
      # Reorder fold assignment back to original ordering
      unit_fold = integer(n)
      unit_fold[idx] = fold_id
      
      # Create CV indicator matrix for this repetition
      folds_rep = matrix(1, nrow = n, ncol = B)
      
      for (fold in seq_len(B)) {
        test_ids = which(unit_fold == fold)
        folds_rep[test_ids, fold] = 0
      }
      
      # Store
      folds_list[[rep]] = folds_rep
    }
    
    # Stack repetitions vertically 
    folds = do.call(rbind, folds_list)
    
    attr(folds, "type") = paste(B, "-fold random cross-validation", sep = "")
    
  }
  
  if(type == "spatial") {
    
    # Calculate centroid coordinates for all geometries
    coords = suppressWarnings(st_coordinates(st_centroid(map)))
    
    # Add centroid print()# Add centroid coordinates safely
    map = map %>%
      mutate(x = coords[, 1],
             y = coords[, 2]) %>%
      st_set_geometry(NULL)
    
    # Run spatial partitioning on spatial units
    resamp = partition_kmeans(
      data = map,
      coords = c("x", "y"),
      nfold = B,
      repetition = r
    )
    
    folds_list = vector("list", r)
    
    for (rep in seq_len(r)) {
      folds_rep = matrix(1, nrow = n, ncol = B)
      for (fold in seq_len(B)) {
        test_ids = resamp[[rep]][[fold]]$test
        folds_rep[test_ids, fold] = 0
      }
      folds_list[[rep]] = folds_rep
    }
    
    # Combine all repetitions by row-binding
    folds = do.call(rbind, folds_list)
    
    attr(folds, "type") = paste(B, "-fold ", "k-means spatial clustering", sep = "")
    
  }
  
  if(type == "blocks") {
    
    folds_list = vector("list", r)
    
    for (rep in seq_len(r)) {
      folds_rep = matrix(1, nrow = n, ncol = B)
      blocks = cv_spatial(x = map,
                          k = B,
                          selection = "random",
                          hexagon = TRUE,
                          iteration = 100,
                          progress = FALSE,
                          report = FALSE,
                          plot = FALSE)
      folds_list[[rep]] = unname(blocks$biomod_table * 1)
    }
    
    # Combine all repetitions by row-binding
    folds = do.call(rbind, folds_list)
    
    attr(folds, "type") = paste(B, "-fold ", "rectangular blocking", sep = "")
  }
  
  if(type == "lto") {
    
    if(B > r) {
      stop("The number of folds is higher than the number of temporal units")
    }
    
    folds = sapply(1:B, function(b) {
      fold_vec = rep(1, r * n)
      fold_vec[((b-1)*n + 1):(b*n)] = 0
      return(fold_vec)
    })
    
    attr(folds, "type") = paste(B,"-fold ", "leave-time-out cross-validation",  sep = "")
    
  }
  
  if(type == "llo") {
    locus = rep(map$NAME, times = r)   
    tempus = rep(1:r, each = n)  
    
    LT = data.frame(Y,locus, tempus)
    
    coords = suppressWarnings(st_coordinates(st_centroid(map)))
    map = map %>%
      mutate(x = coords[, 1],
             y = coords[, 2]) %>%
      st_set_geometry(NULL)
    
    LT = left_join(LT, map, by = c("locus" = "NAME"))
    
    LT = LT[,c("locus", "tempus", "Y", "x", "y")]
    
    # Remove units if present
    LT[] = lapply(LT, function(col) {
      if (inherits(col, "units")) as.numeric(col) else col
    })
    
    
    # Define task
    task_spt = as_task_regr_st(
      LT, target = "Y",
      coordinate_names = c("x", "y"), 
      coords_as_features = FALSE
    )
    
    # Set roles
    task_spt$set_col_roles("locus", roles = "space")
    
    # Define resampling
    rsmp_cstf_loc = rsmp("sptcv_cstf", folds = B)
    rsmp_cstf_loc$instantiate(task_spt)
    
    # Initialize folds matrix with 1s (train by default)
    folds = matrix(1, nrow = n * r, ncol = B)
    
    # Iterate over resampling splits
    for (i in seq_len(B)) {
      test_ids = rsmp_cstf_loc$test_set(i = i)
      folds[test_ids, i] = 0
    }
    
    attr(folds, "type") = paste(B, "-fold ", "leave-location-out cross-validation", sep = "")
    
  }
  
  if(type == "llto") {
    locus = rep(map$NAME, times = r)   
    tempus = rep(1:r, each = n)  
    
    LT = data.frame(Y,locus, tempus)
    
    coords = suppressWarnings(st_coordinates(st_centroid(map)))
    map = map %>%
      mutate(x = coords[, 1],
             y = coords[, 2]) %>%
      st_set_geometry(NULL)
    
    LT = left_join(LT, map, by = c("locus" = "NAME"))
    
    LT = LT[,c("locus", "tempus", "Y", "x", "y")]
    
    # Remove units if present
    LT[] = lapply(LT, function(col) {
      if (inherits(col, "units")) as.numeric(col) else col
    })
    
    
    # Define task
    task_spt = as_task_regr_st(
      LT, target = "Y",
      coordinate_names = c("x", "y"), 
      coords_as_features = FALSE
    )
    
    # Set roles
    task_spt$set_col_roles("locus", roles = "space")
    task_spt$set_col_roles("tempus", roles = "time")
    
    # Define resampling
    rsmp_cstf_loc = rsmp("sptcv_cstf", folds = B)
    rsmp_cstf_loc$instantiate(task_spt)
    
    # Initialize folds matrix with 1s (train by default)
    folds = matrix(1, nrow = n * r, ncol = B)
    
    # Iterate over resampling splits
    for (i in seq_len(B)) {
      test_ids = rsmp_cstf_loc$test_set(i = i)
      folds[test_ids, i] = 0
    }
    
    attr(folds, "type") = paste(B, "-fold ", "leave-location-leave-time-out cross-validation", sep = "")
    
  }
  
  return(folds)
  
}

DeselectBoost = function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  
  tau = ifelse(is.null(tau), 0.01, tau)
  
  if(is.null(data) && class(object$model.frame()) == 'list'){return(stop("Please enter the data."))
  } else if(!is.null(data)){
    data = data
  } else{
    data = object$model.frame()
  }
  
  nameVar = names(coef(object, which = ''))[-1]
  which.response = which(sapply(1:(dim(data)[2]), function(x){identical(as.numeric(data[,x]), as.numeric(object$response))}))
  name.response = colnames(data)[which.response]
  
  mstop = object$mstop()
  RiskRed = object$risk()
  totalRiskRed = RiskRed[1] - RiskRed[mstop+1] 
  diffRiskRed = sapply(seq(1:(mstop)), function(k){RiskRed[k]-RiskRed[k+1]})
  
  select = selected(object) - 1
  diffRiskRed = diffRiskRed[selected(object) - 1 != 0]
  
  select = select[select != 0]
  Var = plyr::count(select)[[1]]
  Risk.Var = lapply(1:length(Var),function(j){sum(diffRiskRed[which(plyr::count(select)[[1]][j] == select)])})
  
  n.parameter = c(names(object$coef()))
  if('(Intercept)' %in% n.parameter) n.parameter = n.parameter[-which(n.parameter == '(Intercept)')]
  if('alpha' %in% n.parameter) n.parameter = n.parameter[-which(n.parameter == 'alpha')]
  
  Risk.order = data.frame(Var,n.parameter, as.numeric(Risk.Var))
  Risk.order = Risk.order[order(Risk.order$as.numeric.Risk.Var.),]
  Risk.order$CumRisk = cumsum(Risk.order$as.numeric.Risk.Var.)
  colnames(Risk.order) = c( 'Var', 'VarName', 'Risk', 'CumRisk')
  
  perc = ifelse(is.null(tau), 0.01, tau) 
  percRiskRed = totalRiskRed * perc
  if(method[1] == 'attributable'){RiskRedOver = Risk.order[which(Risk.order$Risk > percRiskRed),]
  } else{RiskRedOver = Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(plyr::empty(RiskRedOver)){form2 = as.formula(paste(name.response, "~ 1"))
  } else {form2 = as.formula(paste(name.response, " ~ 0 + alpha + ", paste(RiskRedOver$VarName, collapse= "+")))
  if(!is.null(environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]])){
    dfbase = environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]]
  }}
  
  model_after = glmboost(form2, data = data, weights = model.weights(object), offset = 0, center = FALSE,
                         family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  
  out = model_after
  out$tau  = tau
  out$deselectmethod = method[1] 
  class(out) = c(class(out))
  
  return(out)
  
} 

DeselectBoostFE = function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  
  tau = ifelse(is.null(tau), 0.01, tau)
  
  if(is.null(data) && class(object$model.frame()) == 'list'){return(stop("Please enter the data."))
  } else if(!is.null(data)){
    data = data
  } else{
    data = object$model.frame()
  }
  
  nameVar = names(coef(object, which = ''))
  which.response = which(sapply(1:(dim(data)[2]), function(x){identical(as.numeric(data[,x]), as.numeric(object$response))}))
  name.response = colnames(data)[which.response]
  
  mstop = object$mstop()
  RiskRed = object$risk()
  totalRiskRed = RiskRed[1] - RiskRed[mstop+1] 
  diffRiskRed = sapply(seq(1:(mstop)), function(k){RiskRed[k]-RiskRed[k+1]})
  
  select = selected(object) 
  diffRiskRed = diffRiskRed[selected(object) != 0]
  
  select = select[select != 0]
  Var = plyr::count(select)[[1]]
  Risk.Var = lapply(1:length(Var),function(j){sum(diffRiskRed[which(plyr::count(select)[[1]][j] == select)])})
  
  n.parameter = c(names(object$coef()))
  if('(Intercept)' %in% n.parameter) n.parameter = n.parameter[-which(n.parameter == '(Intercept)')]
  if('alpha' %in% n.parameter) n.parameter = n.parameter[-which(n.parameter == 'alpha')]
  
  Risk.order = data.frame(Var,n.parameter, as.numeric(Risk.Var))
  Risk.order = Risk.order[order(Risk.order$as.numeric.Risk.Var.),]
  Risk.order$CumRisk = cumsum(Risk.order$as.numeric.Risk.Var.)
  colnames(Risk.order) = c( 'Var', 'VarName', 'Risk', 'CumRisk')
  
  perc = ifelse(is.null(tau), 0.01, tau) 
  percRiskRed = totalRiskRed * perc
  if(method[1] == 'attributable'){RiskRedOver = Risk.order[which(Risk.order$Risk > percRiskRed),]
  } else{RiskRedOver = Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(plyr::empty(RiskRedOver)){form2 = as.formula(paste(name.response, "~ 1"))
  } else {form2 = as.formula(paste(name.response, " ~ 0 + ", paste(RiskRedOver$VarName, collapse= "+")))
  if(!is.null(environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]])){
    dfbase = environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]]
  }}
  
  model_after = glmboost(form2, data = data, weights = model.weights(object), offset = 0, center = FALSE,
                         family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  
  out = model_after
  out$tau  = tau
  out$deselectmethod = method[1] 
  class(out) = c(class(out))
  
  return(out)
  
} 