library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(tictoc)
library(actuar)
library(dirmult)

ALHR_AIPW <- function(data, beta0 = 0, threshold = 1e-8, max_iter = 100,T_AZ = 'fail1', C_AZ = 'cens1', min.S = 0.007, remove.last = T, run_boot_se = F, boot_runs = 500, bayesian_boot = F, tau = 1){
  # inputs:
  # data: nx(p+) matrix, first 3 columns (X, Delta, A) are fixed, the rest are baseline covariates Z
  #   X - nx1, the censored event time for each patient
  #   Delta - nx1, the event indicator for each patient, 1 = Died, 0 = Censored.
  #   A - nx1, 0/1 treatment variable
  #   Z - nxp, the baseline covariates matrix
  # beta0 - initial value for beta
  # threshold - convergence threshold for the difference in beta
  # max_iter - maximum iterations of Newton-Raphson allowed.
  # T_AZ - working model for the T|A,Z, available ones include ('Cox', 'Spline', 'RSF', 'llogis')
  # C_AZ - working model for the C|A,Z, available ones include ('Cox', 'Spline', 'RSF', 'llogis')
  n = nrow(data)
  X = data[,1]
  n_res = length(unique(X))
  Delta = data[,2]
  Delta_c = (1 - Delta)*(X < tau)
  A = data[,3]
  Z = as.matrix(data[,-(1:3)])
  p = ncol(Z)
  X.sort = sort(unique(X))
  event_rank = pmin(n_res,rank(X, ties.method = 'first'))
  Y = t(sapply(event_rank, function(x) c(rep(1,x), rep(0,n_res-x))))   # n x n_res  (sample x time)
  dNt = t(sapply(event_rank, function(x) tabulate(x, nbins = n_res))) * Delta    # n x n_res  (sample x time)
  model_se = boot_se = 0
  
  if (T_AZ == 'fail1'){
    X0 = X[A==0]
    X1 = X[A==1]
    Delta0 = Delta[A==0]
    Delta1 = Delta[A==1]
    Z0 = Z[A==0,]
    Z1 = Z[A==1,]
    X.sort0 = sort(unique(X0))
    X.sort1 = sort(unique(X1))
    S_t = data.frame(matrix(0, nrow = n_res, ncol = n))
    model0 = flexsurvreg(Surv(X0, Delta0)~ Z0, dist ='llogis')
    model1 = flexsurvreg(Surv(X1, Delta1)~ Z1, dist ='llogis')
    pred0 = predict(model0, type = 'survival', times = X.sort0)$.pred
    pred1 = predict(model1, type = 'survival', times = X.sort1)$.pred
    S_t[which(X.sort %in% X.sort0), A==0] = sapply(pred0, function(x) x$.pred)
    S_t[which(X.sort %in% X.sort1), A==1] = sapply(pred1, function(x) x$.pred)
    S_t[1,] = ifelse(S_t[1,] == 0, 1, S_t[1,])
    S_t[S_t==0] = NA
    S_t = t(as.matrix(tidyr::fill(S_t, names(S_t))))
  }
  if (T_AZ == 'Cox'){
    model = coxph(Surv(X, Delta) ~ A + Z, timefix = F)
    Lambda = exp(cbind(A,Z)%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
    S_t = exp(-Lambda)
    S_t[,1] = ifelse(S_t[,1] == 0, 1, S_t[,1])
    S_t[S_t==0] = NA
    S_t = t(as.matrix(tidyr::fill(data.frame(t(S_t)), names(data.frame(t(S_t))))))
  }
  if (T_AZ == 'Spline'){
    model = hare(X, Delta, cbind(A,Z))
    S_t = 1 - sapply(sort(unique(X)), function(x) phare(q = x, cov = cbind(A,Z), fit = model))
    S_t[is.nan(S_t)] = NA
    S_t[,1] = ifelse(S_t[,1] == 0 | is.na(S_t[,1]), 1, S_t[,1])
    if(any(is.na(S_t))) S_t = t(as.matrix(tidyr::fill(data.frame(t(S_t)), names(data.frame(t(S_t))))))
  }
  if (T_AZ == 'RSF'){
    dat = data.frame(cbind(X,Delta,A,Z))
    colnames(dat) = c('X','Delta','A',paste0('Z',1:p))
    model = rfsrc(Surv(X, Delta) ~ ., data = dat, ntree = 1000, ntime = X.sort, mtry = 2, splitrule = 'bs.gradient')
    S_t = matrix(0, nrow = n, ncol = n_res)
    S_t[, which(X.sort %in% model$time.interest)] = model$survival.oob
    S_t[,1] = ifelse(S_t[,1] == 0, 1, S_t[,1])
    S_t[S_t==0] = NA
    if(any(is.na(S_t))) S_t = t(as.matrix(tidyr::fill(data.frame(t(S_t)), names(data.frame(t(S_t))))))
  }
  # if (C_AZ == 'cens1'){
  #   model = coxph(Surv(X, Delta_c) ~ A + Z, timefix = F)
  #   Lambda_c = exp(cbind(A,Z)%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
  #   S_c = exp(-Lambda_c)
  #   S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
  #   S_c[S_c==0] = NA
  #   S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  # }
  if (C_AZ == 'cens1'){
    model = coxph(Surv(X, Delta_c) ~ Z, timefix = F)
    Lambda_c = exp(Z%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
    S_c = exp(-Lambda_c)
    S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  if (C_AZ == 'Cox'){
    model = coxph(Surv(X, Delta_c) ~ A + Z, timefix = F)
    Lambda_c = exp(cbind(A,Z)%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
    S_c = exp(-Lambda_c)
    S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  if (C_AZ == 'Spline'){
    model = hare(X, Delta_c, cbind(A,Z))
    S_c = 1 - sapply(sort(unique(X)), function(x) phare(q = x, cov = cbind(A,Z), fit = model))
    S_c[is.nan(S_c)] = NA
    S_c[,1] = ifelse(S_c[,1] == 0 | is.na(S_c[,1]), 1, S_c[,1])
    if(any(is.na(S_c))) S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  if (C_AZ == 'RSF'){
    dat = data.frame(cbind(X,Delta_c,A,Z))
    colnames(dat) = c('X','Delta_c','A',paste0('Z',1:p))
    model = rfsrc(Surv(X, Delta_c) ~ ., data = dat, ntree = 1000, ntime = X.sort, mtry = 2, splitrule = 'bs.gradient')
    S_c = matrix(0, nrow = n, ncol = n_res)
    S_c[, which(X.sort %in% model$time.interest)] = model$survival.oob
    S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
    S_c[S_c==0] = NA
    if(any(is.na(S_c))) S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  S_c = as.matrix(S_c)
  S_t = as.matrix(S_t)
  check = T
  calc_min.S = function(x) tanh(x)/2+0.5
  min.S.base = atanh((min.S - 0.5)*2)
  # while B2_0 is negative, increase min.S
  while (check){
    min.S.base = min.S.base + 0.25
    min.S = calc_min.S(min.S.base)
    if (min(S_t) < min.S) S_t = 1 - (1 - S_t)*(1 - min.S)/(1-min(S_t))
    if (min(S_c) < min.S) S_c = 1 - (1 - S_c)*(1 - min.S)/(1-min(S_c))
    Lambda_c = -log(S_c) # n x n_res  (sample x time)
    lambda_c = t(apply(Lambda_c, 1, function(x) c(x[1], diff(x))))    # n x n_res  (sample x time)
    # calculate J
    dNc = t(sapply(event_rank, function(x) tabulate(x, nbins = n_res))) * (1-Delta) # n x n_res  (sample x time)
    dMc = dNc - Y * lambda_c   # n x n_res  (sample x time)
    J = t(apply(dMc/S_t/S_c, 1, cumsum))   # n x n_res  (sample x time)
    # J = J - rowMeans(J)
    # J = apply(dMc/S_t/S_c, 1, cumsum)
    # J = t(J - rowMeans(J))
    dF_t = t(apply(1-S_t, 1, function(x) c(x[1], diff(x))))   # n x n_res  (sample x time)
    dB1 = (dNt/S_c + J*dF_t)  # n x n_res  (sample x time)
    if (remove.last) dB1 = dB1[,-n_res]
    
    # function for calculating B2_l
    calc_B2 = function(beta, l = 0)  colMeans(A^l*exp(beta*A) * (Y/S_c + J*S_t)) # n_res x 1   (time x 1)
    if (remove.last) calc_B2 = function(beta, l = 0)  colMeans(A^l*exp(beta*A) * (Y/S_c + J*S_t))[-n_res]
    calc_A_bar = function(beta)  (calc_B2(beta,1)/calc_B2(beta,0)) # n_res x 1   (time x 1)
    calc_U = function(A_bar, weight=1)  sum((dB1*outer(A, A_bar, '-'))*weight)/n   #  1x1
    calc_dU = function(A_bar, weight=1)  - sum((t(dB1*weight)*(A_bar - A_bar^2)))/n   #  1x1
    calc_PL = function(beta, weight=1) sum((dB1*outer(A*beta, log(abs(calc_B2(beta,0))), '-')*weight))/n   #  1x1
    calc_beta = function(beta0, weight=1){
      beta = beta0
      diff = 999
      iter_count = 0
      while(diff > threshold){
        beta_old = beta
        A_bar = calc_A_bar(beta)
        U = calc_U(A_bar, weight = weight)
        dU = calc_dU(A_bar, weight = weight)
        beta = beta - U/dU
        diff = abs(beta - beta_old)
        iter_count = iter_count + 1
        if (iter_count > max_iter) stop(paste('Did not converge after',max_iter,'steps.'))
      }  # 16.48s per 100 runs
      return(beta)
    }
    
    B2.min = min(min(calc_B2(-10)), min(calc_B2(10)))
    U_sign = calc_U(calc_A_bar(-10))*calc_U(calc_A_bar(10))
    not.converged = F
    tryCatch({beta = calc_beta(beta0)},
             error = function(e) not.converged <<- T)
    check = (B2.min < 1e-4 | is.nan(B2.min) | U_sign > 0 | not.converged)
    if (check == F) min.S.base = min.S.base + 0.25
    min.S = calc_min.S(min.S.base)
  }
  
  # points = seq(-10,10,length.out = 100)
  # plot(points, sapply(points, function(x) calc_PL(x, weight=weight)))
  # plot(points, sapply(points, function(x) calc_U(calc_A_bar(x), weight=1)))
  # plot(points, sapply(points, function(x) calc_dU(calc_A_bar(x))))
  # plot(points, sapply(points, function(x) min(calc_B2(x))))
  
  # model se
  A_bar = calc_A_bar(beta)
  K = mean(rowSums((dB1   -   t(t(exp(beta*A) * (Y/S_c + J*S_t))[-n_res,]/calc_B2(beta)*colMeans(dB1)))*outer(A, A_bar, '-'))^2)
  nu = mean(colSums(t(dB1)*(A_bar - A_bar^2)))
  model_se = sqrt(K/nu^2/n)
  # boot se
  if (run_boot_se) {
    if (bayesian_boot)  weights = rdirichlet(n = boot_runs, alpha = rep(1,n))
    else weights = t(rmultinom(n = boot_runs, size = n, prob = rep(1,n)))
    boot_beta = rep(0, boot_runs)
    for (i in 1:boot_runs){
      boot_beta[i] = calc_beta(beta, weight = weights[i,])
    }
    boot_se = sd(boot_beta)
  }
  return(list(beta = beta, model_se = model_se, boot_se = boot_se, min.S = min.S))
}

# beta0 = beta_true; threshold = 1e-6; max_iter = 100; C_AZ = 'exp-AZ'; min.S = 1e-8; run_model_se = F; tau = tau
ALHR_IPW <- function(data, beta0 = 0, threshold = 1e-6, max_iter = 100, C_AZ = 'cens1', remove.last = T, min.S = 0.05, run_boot_se = F, boot_runs = 500, bayesian_boot = F, tau = 1){
  # inputs:
  # data: nx(p+3) matrix, first 3 columns (X, Delta, A) are fixed, the rest are baseline covariates Z
  #   X - nx1, the censored event time for each patient
  #   Delta - nx1, the event indicator for each patient, 1 = Died, 0 = Censored.
  #   A - nx1, 0/1 treatment variable
  #   Z - nxp, the baseline covariates matrix
  # beta0 - initial value for beta
  # threshold - convergence threshold for the difference in beta
  # max_iter - maximum iterations of Newton-Raphson allowed.
  # C_AZ - working model for C|AZ,  'exp-AZ' or 'exp-A' or 'exp'
  n = nrow(data)
  X = data[,1]
  n_res = length(unique(X))
  Delta = data[,2]
  Delta_c = (1 - Delta)*(X<tau)
  A = data[,3]
  Z = as.matrix(data[,-(1:3)])
  X.sort = sort(unique(X))
  event_rank = pmin(n_res,rank(X, ties.method = 'first'))
  boot_se = model_se = 0
  Y = t(sapply(event_rank, function(x) c(rep(1,x), rep(0,n_res-x))))   # n x n_res  (sample x time)
  dNt = t(sapply(event_rank, function(x) tabulate(x, nbins = n_res))) * Delta    # n x n_res  (sample x time) 
  # estimate S_c using Cox for each A
  # if (C_AZ == 'cens1'){
  #   model = coxph(Surv(X, Delta_c) ~ A + Z, timefix = F)
  #   Lambda_c = exp(cbind(A,Z)%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
  #   S_c = exp(-Lambda_c)
  #   S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
  #   S_c[S_c==0] = NA
  #   S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  # }
  if (C_AZ == 'cens1'){
    model = coxph(Surv(X, Delta_c) ~ Z, timefix = F)
    Lambda_c = exp(Z%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
    S_c = exp(-Lambda_c)
    S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  if (C_AZ == 'Cox'){
    model = coxph(Surv(X, Delta_c) ~ A + Z, timefix = F)
    Lambda_c = exp(cbind(A,Z)%*%model$coefficients) %*% t(basehaz(model, centered = F)$hazard)
    S_c = exp(-Lambda_c)
    S_c[,1] = ifelse(S_c[,1] == 0, 1, S_c[,1])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(data.frame(t(S_c)), names(data.frame(t(S_c))))))
  }
  if (C_AZ == 'W2'){
    S_c = data.frame(matrix(0, nrow = n_res, ncol = n))
    model0 = survfit(Surv(X[A==0], Delta_c[A==0])~1, timefix = F)
    model1 = survfit(Surv(X[A==1], Delta_c[A==1])~1, timefix = F)
    X.sort = sort(unique(X))
    S_c[which(X.sort %in% model0$time), A==0] = cbind(exp(-model0$cumhaz))[,rep(1,sum(A==0))]
    S_c[which(X.sort %in% model1$time), A==1] = cbind(exp(-model1$cumhaz))[,rep(1,sum(A==1))]
    S_c[1,] = ifelse(S_c[1,] == 0, 1, S_c[1,])
    S_c[S_c==0] = NA
    S_c = t(as.matrix(tidyr::fill(S_c, names(S_c))))
  }
  if (C_AZ == 'W1'){
    model = survfit(Surv(X, Delta) ~ 1, timefix = F)
    S_t = model$surv
    S_t[S_t == 0] = min(S_t[S_t>0])
    S_c = matrix(replicate(n, colSums(Y)/S_t), nrow = n, byrow = T)/n
  }
  S_c = as.matrix(S_c)
  if (min(S_c) < min.S) S_c = 1 - (1 - S_c)*(1 - min.S)/(1-min(S_c))
  dB1 = dNt/S_c; if (remove.last) dB1 = dB1[,-n_res] # n x n_res  (sample x time) 
  # function for calculating B2_l
  calc_B2 = function(beta, l = 0)  colMeans(A^l*exp(beta*A) *  Y/S_c)  # n_res x 1   (time x 1) 
  if (remove.last)  calc_B2 = function(beta, l = 0)  colMeans(A^l*exp(beta*A) *  Y/S_c)[-n_res]
  calc_A_bar = function(beta)  calc_B2(beta,1)/calc_B2(beta,0) # n_res x 1   (time x 1)   
  calc_U = function(A_bar, weight=1)  sum(dB1*outer(A, A_bar, '-')*weight)/n   #  1x1         
  calc_dU = function(A_bar, weight=1)  - sum(t(dB1*weight)*(A_bar - A_bar^2))/n   #  1x1      
  calc_PL = function(beta)  sum((dB1*outer(A*beta, log(calc_B2(beta,0)), '-')))/n   #  1x1 
  calc_beta = function(beta0, weight=1){
    beta = beta0
    diff = 999
    iter_count = 0
    while(diff > threshold){
      beta_old = beta
      A_bar = calc_A_bar(beta)
      U = calc_U(A_bar, weight = weight)
      dU = calc_dU(A_bar, weight = weight)
      beta = beta - U/dU
      diff = abs(beta - beta_old)
      iter_count = iter_count + 1
      if (iter_count > max_iter) stop(paste('Did not converge after',max_iter,'steps.'))
    }  # 16.48s per 100 runs
    return(beta)
  }
  # points = seq(-10,10,length.out = 100)
  # plot(points, sapply(points, calc_PL))
  # plot(points, sapply(points, function(x) calc_U(calc_A_bar(x))))
  # plot(points, sapply(points, function(x) calc_dU(calc_A_bar(x))))
  
  # calculate beta
  beta = calc_beta(beta0)
  # model se using Gillen's robust estimator sqrt(|U.star|_2^2 / dU^2)
  A_bar = calc_A_bar(beta)
  B2.0 = calc_B2(beta)
  B2.1 = calc_B2(beta, l=1)
  dU = calc_dU(A_bar)
  T.Y_S_c = t(Y/S_c)
  if (remove.last) T.Y_S_c = T.Y_S_c[-n_res,]
  U.star = (rowSums(dB1*outer(A, A_bar, '-')) - A*exp(beta*A)*colSums(T.Y_S_c/B2.0) + exp(beta*A)*colSums(T.Y_S_c*B2.1/B2.0^2))/n
  model_se = sqrt(sum(U.star^2)/n/dU^2 / n)
  # boot se
  if (run_boot_se) {
    if (bayesian_boot)  weights = rdirichlet(n = boot_runs, alpha = rep(1,n))
    else weights = t(rmultinom(n = boot_runs, size = n, prob = rep(1,n)))
    boot_beta = rep(0, boot_runs)
    for (i in 1:boot_runs){
      boot_beta[i] = calc_beta(beta, weight = weights[i,])
    }
    boot_se = sd(boot_beta)
  }
  return(list(beta = beta, model_se = model_se, boot_se = boot_se, min.S = min.S))
}



simulate <- function(T_AZ_sim = 'exp-AZ', C_AZ_sim = 'exp-AZ', n = 1000, p = 2, A_prob = 0.5){
  # T_AZ_sim takes 3 models:  'llogis', 'exp-AZ', 'piecewise-exp' and 'llogis-lnorm'
  # C_AZ_sim takes 4 models:  'exp', 'exp-A', 'exp-AZ', 'piecewise-exp' and 'llogis'
  n_original = n
  n = ceiling(n*1.01)
  A = rbinom(n, 1, A_prob)
  Z = matrix(rnorm(n*p, 1, 1), nrow = n)  # Z1, Z2 ~ N(0,1)
  if (T_AZ_sim == 'exp-AZ') t = -log(runif(n))/exp(0.6*A+as.vector(Z%*%c(-0.3,0.5)))
  if (T_AZ_sim == 'llogis') t = exp(-0.4*A + as.vector(Z%*%c(-0.4,-0.4)) + 0.3*rlogis(n))
  if (T_AZ_sim == 'piecewise-exp'){
    t = sapply(1:n, function(x) exp_cdfsim(n=1, endtime=9999, 
                                           theta=c(exp(A[x]), exp(0.2*Z[x,1]), exp(0.2*Z[x,2])) ,tau=c(0.2, 0.4))$time)
  }
  if (T_AZ_sim == 'exp-lnorm'){
    t = rep(0,n)
    t[A==0] = -log(runif(sum(A==0)))/exp(0.7+as.vector(Z[A==0,]%*%c(-0.3,0.2)))
    t[A==1] = exp(-0.5 + as.vector(Z[A==1,]%*%c(-0.3, -0.3)) + 0.2*rnorm(sum(A==1)))
    #data = sim.survdata(X = cbind(A[1:10]), T = n, type = 'tvbeta', beta = cbind(runif(10)), censor = F)
    #t0 = sim.survdata(T = 20, X = cbind(A[1:n]), type = 'tvbeta', beta = data.frame(beta1 = (1:20)/20), censor = F)$data
    # ww = exp(0 + as.vector(matrix(rnorm(20000, 1, 1), nrow = 10000)%*%c(-0.3, -0.6)) + log(rexp(10000)))
    # hist(ww[ww<0.9], breaks= 100)
    # ee = exp(-0.4 + as.vector(matrix(rnorm(20000, 1, 1), nrow = 10000)%*%c(-0.3, -0.3)) + 0.2*rnorm(10000))
    # hist(ee[ee<0.9],breaks = 100)
  }
  if (C_AZ_sim == 'exp') c = -log(runif(n))/exp(0.3)
  if (C_AZ_sim == 'exp-A') c = -log(runif(n))/exp(0.6*A)
  if (C_AZ_sim == 'exp-AZ') c = -log(runif(n))/exp(0.6*A - as.vector(Z%*%c(0.2,-0.2)))
  if (C_AZ_sim == 'piecewise-exp'){
    c = sapply(1:n, function(x) exp_cdfsim(n=1, endtime=9999, 
                                           theta=c(exp(-0.2*A[x]), exp(0.8*A[x] + 0.4*Z[x,1] + 0.4*Z[x,2])) ,tau=0.35)$time)
  }
  if (C_AZ_sim == 'exp-lnorm') {
    # c = exp(-0.5*A + as.vector(Z%*%c(-0.3, -0.3)) + 0.2*rlogis(n))
    c = rep(0,n)
    c[A==0] = -log(runif(sum(A==0)))/exp(0.5+as.vector(Z[A==0,]%*%c(-0.3,0.1)))
    c[A==1] = exp(-0.4 + as.vector(Z[A==1,]%*%c(-0.3, -0.3)) + 0.2*rnorm(sum(A==1)))
  }
  X = pmin(t,c)
  Delta = 1*(t<c)
  duplicated_index = unique(c(which(duplicated(X)), which(duplicated(t)), which(duplicated(c))))
  if (length(duplicated_index) != 0) {
    cat(paste(length(duplicated_index) ,'duplicates detected.\n'))
    return(cbind(t, c, X, Delta, A, Z)[-duplicated_index,][1:n_original,])
  }else return(cbind(t, c, X, Delta, A, Z)[1:n_original,])
}



one_setting_simulation = function(runs = 500, n_sim = 1000, p = 2, tau = 0.7, T_AZ_sim = 'exp-lnorm', C_AZ_sim = 'llogis', min_S = 5e-2){
  set.seed(2021)
  runs_max = round(runs*2)
  all_data_ori = simulate(n = n_sim*runs_max, T_AZ_sim = T_AZ_sim, C_AZ_sim = C_AZ_sim, p = p, A_prob = 0.5)
  
  all_data = all_data_ori[,-(1:2)]
  all_data[,2] = all_data[,2]*(all_data[,1] < tau)
  all_data[,1] = ifelse(all_data[,1] < tau, all_data[,1], tau)
  estimators = c('aipw-llogis-cox', 'aipw-spline-cox', 'aipw-rsf-cox', 'aipw-spline-spline', 'aipw-rsf-rsf', 'ipw-az', 'ipw-a', 'ipw-1', 'pl')
  beta_estimates = model_se = boot_se = 
    matrix(0, nrow = runs_max, ncol = length(estimators), dimnames = list(NULL, estimators))
  
  i=0
  each_core = function(){
    
  }
  converged_runs = rep(FALSE, runs_max)
  while(sum(converged_runs) < runs & i < runs_max){
    i = i + 1
    errored <<- F
    tryCatch({
      data = all_data[(1 + (i-1)*n_sim):(i*n_sim),]
      beta_estimates[i,1] = ALHR_AIPW(data = data, T_AZ = 'llogis', C_AZ = 'Cox', min.S=min_S)
      beta_estimates[i,2] = ALHR_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'Cox', min.S=min_S)
      beta_estimates[i,3] = ALHR_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'Cox', min.S=min_S)
      beta_estimates[i,4] = ALHR_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'Spline', min.S=min_S)
      beta_estimates[i,5] = ALHR_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'RSF', min.S=min_S)  
      beta_estimates[i,6] = ALHR_IPW(data = data, C_AZ = 'exp-AZ', min.S=min_S)
      beta_estimates[i,7] = ALHR_IPW(data = data, C_AZ = 'exp-A', min.S=min_S)
      beta_estimates[i,8] = ALHR_IPW(data = data, C_AZ = 'exp', min.S=min_S)
      beta_estimates[i,9] = as.numeric(coxph(Surv(data[,1], data[,2])~data[,3], timefix = F)$coefficients)
      if (i %% 100 == 0) cat(paste0(i, ' runs done.\n'))
    },
    error = function(e){
      cat(paste0('error occurred in ', i,'-th run when T_AZ_sim = ', T_AZ_sim, ', C_AZ_sim = ',C_AZ_sim,'.\n'))
      error_data <<- data
      error_i <<- i
      errored <<- T
    }
    )
    converged_runs[i] = !errored
  }
  
  # beta_true_estimates = sapply(1:runs_max, function(x) as.numeric(coxph(Surv(all_data_ori[(1 + (x-1)*n_sim):(x*n_sim),1]*(all_data_ori[(1 + (x-1)*n_sim):(x*n_sim),1] < tau) + tau*(all_data_ori[(1 + (x-1)*n_sim):(x*n_sim),1] >= tau), 1*(all_data_ori[(1 + (x-1)*n_sim):(x*n_sim),1] < tau))~all_data_ori[(1 + (x-1)*n_sim):(x*n_sim),5], timefix = F)$coefficients))
  beta_true = as.numeric(coxph(Surv(all_data_ori[,1]*(all_data_ori[,1] < tau) + tau*(all_data_ori[,1] >= tau), 1*(all_data_ori[,1] < tau))~all_data_ori[,5], timefix = F)$coefficients)
  # mean(beta_true_estimates)
  if (i == runs_max & sum(converged_runs) < runs) cat(paste0('Only ', sum(converged_runs),' out of ',runs, 'required runs converged.\n'))
  
  return(list(beta_true = beta_true, beta_estimates = beta_estimates[converged_runs,], model_se = model_se[converged_runs,], boot_se = boot_se[converged_runs,], converged_runs = converged_runs))
}



output = function(results){
  bias = round(c(beta_true = results$beta_true, colMeans(results$beta_estimates) - results$beta_true),4)
  sd = round(c(0, apply(results$beta_estimates, 2, sd)),4)
  rmse = round(c(0, sqrt(colMeans((results$beta_estimates - results$beta_true)^2))),4)
  return(rbind(bias, sd,rmse))
}