library(survival)
library(parallel)
library(randomForestSRC)
library(polspline)
library(flexsurv)
library(tictoc)
library(actuar)
library(dirmult)
library(cpsurvsim)
library(coxed)



TATE_AIPW <- function(data, beta0 = 0, threshold = 1e-8, max_iter = 100,T_AZ = 'fail1', C_AZ = 'cens1', remove.last = T, min.S = 0.007, run_boot_se = T, boot_runs = 500, bayesian_boot = F, tau = 1){
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
    calc_PL = function(beta) sum((dB1*outer(A*beta, log(abs(calc_B2(beta,0))), '-')))/n   #  1x1
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
  # plot(points, sapply(points, calc_PL))
  # plot(points, sapply(points, function(x) calc_U(calc_A_bar(x))))
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


TATE_IPW <- function(data, beta0 = 0, threshold = 1e-6, max_iter = 100, C_AZ = 'cens1', remove.last = T, min.S = 0.05, run_boot_se = T, boot_runs = 500, bayesian_boot = F, tau = 1){
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
  dB1 = dNt/S_c; if (remove.last) dB1 = (dB1[,-n_res]) # n x n_res  (sample x time) 
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


simulate <- function(T_AZ_sim = 'fail1', C_AZ_sim = 'cens1', n = 1000, p = 2, A_prob = 0.5){
  # T_AZ_sim takes 3 models:  'llogis', 'exp-AZ', 'piecewise-exp' and 'llogis-lnorm'
  # C_AZ_sim takes 4 models:  'exp', 'exp-A', 'exp-AZ', 'piecewise-exp' and 'llogis'
  n_original = n
  n = ceiling(n*1.01)
  A = rbinom(n, 1, A_prob)
  Z = matrix(runif(n*p, -1, 1), nrow = n)  # Z1, Z2 ~ N(0,1)
  
  # want mean(T>tau)/mean(C>tau) > 0.05 for any A,Z
  # if (T_AZ_sim == 'fail1'){
  #   t = rep(0,n)
  #   t[A==0] = exp(-0+ as.vector(Z[A==0,]%*%c(0.15,0.15)) + 0.25*rnorm(sum(A==0)))
  #   t[A==1] = exp(-1.2 + as.vector(Z[A==1,]%*%c(-1,-1)) + 2.5*rnorm(sum(A==1)))
  # } # min inside exp is -0.4 when 0.25norm; -0.8 when 0.5rnorm; -1.6 when 1rnorm; 2.4 when 1.5rnorm; 3.2 when 2rnorm; 4.0 when 2.5rnorm; 4.8 when 3rnorm
  # # if (T_AZ_sim == 'fail2'){
  #   t = rep(0,n)
  #   t[A==0] = exp(-0.0+ as.vector(Z[A==0,]%*%c(0.1,0.2)) + 0.125*rlogis(sum(A==0)))
  #   t[A==1] = exp(-0.9 + as.vector(Z[A==1,]%*%c(-1, -1)) + 1*rlogis(sum(A==1)))
  # } # min inside exp is -0.35 when 0.125*rnorm;-0.7 when 0.25*rlogis; -1.1 when 0.5*rlogis; -2.9 when 1*rlogis
  if (T_AZ_sim == 'fail1'){
    t = rep(0,n)
    t[A==0] = exp(-0+ as.vector(Z[A==0,]%*%c(0.1,0.15)) + 0.125*rlogis(sum(A==0)))
    t[A==1] = exp(-1 + as.vector(Z[A==1,]%*%c(-0.9, -1)) + 1*rlogis(sum(A==1)))
  } # min inside exp is -0.35 when 0.125*rnorm;-0.7 when 0.25*rlogis; -1.1 when 0.5*rlogis; -2.9 when 1*rlogis
  if (T_AZ_sim == 'fail2'){
    t = rep(0,n)
    t[A==0] = exp(-0.5 - 0.3 + as.vector(Z[A==0,]%*%c(0.2,0.3)) + 0.5*rlgamma(sum(A==0), 2, 5))
    t[A==1] = exp(-2 - 1.7 + as.vector(Z[A==1,]%*%c(-0.6, -0.7)) + 2*rlgamma(sum(A==1), 2, 5))
  } # -2.5 when 1*rlg(2,5); -1.7 when 0.3*rlg(5,5); -1 when 0.4*rlg(2,5); -0.5 when 0.2*rlg(2,5)
  # -5 when 2*rlg(2,5), -1.3 when 0.5*rlg(2,5)
  # if (T_AZ_sim == 'fail2'){
  #   t = rep(0,n)
  #   t[A==0] = exp(-0.3+ as.vector(Z[A==0,]%*%c(0.3,0.2)) + 0.2*rlnorm(sum(A==0)))
  #   t[A==1] = exp(-1.8 + as.vector(Z[A==1,]%*%c(-0.7, 1.3)) + 0.8*rlnorm(sum(A==1)))
  # }  # min is -0.95 when 0.2rlnorm; -3.8 when 0.8rlnorm
  
  # if (T_AZ_sim == 'fail2'){
  #   t = rep(0,n)
  #   t[A==0] = exp(-1.2 + as.vector(abs(Z[A==0,])%*%c(0.3,0.3)) + 0.3*rlgamma(sum(A==0), 5, 5))
  #   t[A==1] = exp(-1.4 + as.vector(abs(Z[A==1,])%*%c(-0.5, -0.5)) + 1*rlgamma(sum(A==1), 2, 5))
  # } # -2.5 when 1*rlg(2,5); -1.7 when 0.3*rlg(5,5)
  # # 
  # # mean(exp(-1.7 + 0.3*rlgamma(sum(A==0), 5, 5))>1)
  # # mean(exp(-1 + as.vector(abs(Z[A==1,])%*%c(-0.6, -0.8)) + 1*rlgamma(sum(A==1), 2, 5))>1)

  if (C_AZ_sim == 'cens1') c = -log(runif(n))/exp(-0.5 + as.vector(Z%*%c(-0.8, -0.6))) # max inside exp is 0.9            
  if (C_AZ_sim == 'cens2'){
    c = rep(0,n)
    c[Z[,1]>0] = exp(-0.5 - 1*A[Z[,1]>0] + as.vector(sqrt(abs(Z[Z[,1]>0,]))%*%c(0.4, 0.4)) + 0.3*rlnorm(sum(Z[,1]>0)))
    c[Z[,1]<=0] = exp(2.4 - 1*A[Z[,1]<=0] + as.vector(sqrt(abs(Z[Z[,1]<=0,]))%*%c(0, -1.2)) - 1*rlnorm(sum(Z[,1]<=0)))
  }   # min inside exp is -1.5 when 0.3*rlnorm  ; min is 0.2 when -1*rlnorm 
  
  # if (T_AZ_sim == 'div_hazard'){
  #   a = 0.8
  #   b = 2.5
  #   uni = runif(n)
  #   t = -1/a*log(1-uni)
  #   t[A==1] = ifelse(uni[A==1] > 0.5, 1/b*( 0.5*(b-a) - log(1-uni[A==1])),t[A==1])
  # }
  X = pmin(t,c)
  Delta = 1*(t<c)
  duplicated_index = unique(c(which(duplicated(X)), which(duplicated(t)), which(duplicated(c))))
  if (length(duplicated_index) != 0) {
    cat(paste(length(duplicated_index) ,'duplicates detected.\n'))
    return(cbind(t, c, X, Delta, A, Z)[-duplicated_index,][1:n_original,])
  }else return(cbind(t, c, X, Delta, A, Z)[1:n_original,])
}


one_setting_simulation = function(runs = 500, n_sim = 1000, p = 2, tau = 1, T_AZ_sim = 'fail1', C_AZ_sim = 'cens1', run_boot_se = c(T,T,T,T,T,F,F,F,F,F,F), boot_runs = 500, bayesian_boot = F){
  set.seed(2021)
  all_data_ori = simulate(n = 2*n_sim*runs, T_AZ_sim = T_AZ_sim, C_AZ_sim = C_AZ_sim, p = p, A_prob = 0.5)
  beta_true = as.numeric(coxph(Surv(all_data_ori[,1]*(all_data_ori[,1] < tau) + tau*(all_data_ori[,1] >= tau), 1*(all_data_ori[,1] < tau))~all_data_ori[,5], timefix = F)$coefficients)
  
  all_data = all_data_ori[,-(1:2)]
  all_data[,2] = all_data[,2]*(all_data[,1] < tau)
  all_data[,1] = ifelse(all_data[,1] < tau, all_data[,1], tau)
  estimators = c('AIPW-fail1-cens1', 'AIPW-spline-cens1', 'AIPW-fail1-spline', 'AIPW-spline-spline', 'AIPW-rsf-cens1', 'AIPW-fail1-rsf', 'AIPW-rsf-rsf', 'IPW-cens1', 'W2', 'W1', 'PL')
  beta_estimates = model_se = boot_se = min_S_list = model_CP = boot_CP =
    matrix(0, nrow = runs, ncol = length(estimators), dimnames = list(NULL, estimators))
  
  # Test each data set for whether spline converges on both T and C.
  if (IPW_only == F){
    test_spline_converge = function(i){
      data = all_data[(1 + (i-1)*n_sim):(i*n_sim),]
      out1 = capture.output(hare(data[,1], data[,2], data[,3:5]))
      out2 = capture.output(hare(data[,1], (1 - data[,2])*(data[,1]<tau), data[,3:5]))
      return(ifelse("Convergence problems.... stopping addition" %in% c(out1, out2),F , T))
    }
    test_results = mclapply(1:(2*runs), test_spline_converge, mc.cores = cores)
    all_data = all_data[sort(sapply(which(unlist(test_results)), function(x) (1 + (x-1)*n_sim):(x*n_sim))),][1:(n_sim*runs),]
  }else{
    all_data = all_data[1:(n_sim*runs),]
  }
  # estimate beta and model se
  each_core = function(i){
    data = all_data[(1 + (i-1)*n_sim):(i*n_sim),]
    beta_est = model_se = boot_se = min_S_list = boot.ci.low = boot.ci.high = coverred = rep(0, length(estimators))
    # beta0 = 0;  T_AZ = 'fail1'; C_AZ = 'cens1';threshold = 1e-10; run_boot_se = F; boot_runs = boot_runs; max_iter = 100; remove.last = T; min.S = 0.007;boot_runs = 500; bayesian_boot = F;tau = tau
    if (IPW_only == T){
      res8 = TATE_IPW(data = data, C_AZ = 'cens1' , beta0 = 0, run_boot_se = run_boot_se[8], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res9 = TATE_IPW(data = data, C_AZ = 'W2' , beta0 = 0, run_boot_se = run_boot_se[9], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res10 = TATE_IPW(data = data, C_AZ = 'W1' , beta0 = 0, run_boot_se = run_boot_se[10], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      cox.model = coxph(Surv(data[,1], data[,2])~data[,3], timefix = F)
      beta_est[11] = as.numeric(cox.model$coefficient)
      model_se[11] = as.numeric(summary(cox.model)$coefficients[3])
      for (i in 8:(length(estimators) - 1)){
        eval(parse(text = paste0('beta_est[',i,'] = res',i,'$beta')))
        eval(parse(text = paste0('model_se[',i,'] = res',i,'$model_se')))
        eval(parse(text = paste0('boot_se[',i,'] = res',i,'$boot_se')))
        eval(parse(text = paste0('min_S_list[',i,'] = res',i,'$min.S')))
      }
    }else{
      res1 = TATE_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[1], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res2 = TATE_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[2], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res3 = TATE_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'Spline', beta0 = 0, run_boot_se = run_boot_se[3], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res4 = TATE_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'Spline', beta0 = 0, run_boot_se = run_boot_se[4], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res5 = TATE_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[5], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res6 = TATE_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'RSF', beta0 = 0, run_boot_se = run_boot_se[6], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res7 = TATE_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'RSF', beta0 = 0, run_boot_se = run_boot_se[7], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res8 = TATE_IPW(data = data, C_AZ = 'cens1' , beta0 = 0, run_boot_se = run_boot_se[8], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res9 = TATE_IPW(data = data, C_AZ = 'W2' , beta0 = 0, run_boot_se = run_boot_se[9], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res10 = TATE_IPW(data = data, C_AZ = 'W1' , beta0 = 0, run_boot_se = run_boot_se[10], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      cox.model = coxph(Surv(data[,1], data[,2])~data[,3], timefix = F)
      beta_est[11] = as.numeric(cox.model$coefficient)
      model_se[11] = as.numeric(summary(cox.model)$coefficients[3])
      for (i in 1:(length(estimators) - 1)){
        eval(parse(text = paste0('beta_est[',i,'] = res',i,'$beta')))
        eval(parse(text = paste0('model_se[',i,'] = res',i,'$model_se')))
        eval(parse(text = paste0('boot_se[',i,'] = res',i,'$boot_se')))
        eval(parse(text = paste0('min_S_list[',i,'] = res',i,'$min.S')))
      }
    }
    return(list(beta_estimates = beta_est, model_se = model_se, boot_se = boot_se, min_S_list = min_S_list))
  }
  parallel_results = mclapply(1:runs, each_core, mc.cores = cores)
  for (out in c('beta_estimates','model_se', 'boot_se', 'min_S_list')){
    eval(parse(text=paste0(out, ' = t(sapply(parallel_results, function(x) x$', out, '))')))
    eval(parse(text=paste0('colnames(',out,')=estimators')))
  }
  
  # calculate CP using both model_se and boot_se
  boot.ci.low = beta_estimates - boot_se*qnorm(0.975)
  boot.ci.high = beta_estimates + boot_se*qnorm(0.975)
  model.ci.low = beta_estimates - model_se*qnorm(0.975)
  model.ci.high = beta_estimates + model_se*qnorm(0.975)
  boot.coverred = (boot.ci.low < beta_true)*(boot.ci.high > beta_true)
  model.coverred = (model.ci.low < beta_true)*(model.ci.high > beta_true)
  
  return(list(beta_true = beta_true, beta_estimates = beta_estimates, model_se = model_se, boot_se = boot_se, min_S_list = min_S_list, boot.coverred = boot.coverred, model.coverred = model.coverred))
}


C_table = function(T_AZ_list= c('fail1', 'fail2'), C_AZ_list = c('cens1','cens2'), tau = 1){
  c.table = matrix(0, ncol = length(C_AZ_list), nrow = length(T_AZ_list), dimnames = list(T_AZ_list, C_AZ_list))
  for (i in 1:length(T_AZ_list)){
    for (j in 1:length(C_AZ_list)){
      set.seed(2021)
      data0 = simulate(n=600000, T_AZ_sim = T_AZ_list[i], C_AZ_sim = C_AZ_list[j])
      # tau = round(quantile(data0[,3], 0.9),1)
      A0 = data0[,5]==0
      A1 = data0[,5] == 1
      if (j==1) {
        beta_true = as.numeric(coxph(Surv(data0[,1]*(data0[,1] < tau) + tau*(data0[,1] >= tau), 1*(data0[,1] < tau))~data0[,5], timefix = FALSE)$coefficients)
        print(beta_true)
      }
      hist( (data0[,3]*(data0[,3]<tau) + tau*(data0[,3]>=tau))[data0[,5]==0], breaks= 100)
      hist( (data0[,3]*(data0[,3]<tau) + tau*(data0[,3]>=tau))[data0[,5]==1], breaks= 100)
      c = paste0(round(1-mean(data0[,4]*(data0[,1]<tau)), 3)*100, '%')
      c0 = paste0(round(1-mean(data0[data0[,5] == 0,4]*(data0[data0[,5] == 0,1]<tau)), 3)*100, '%')
      c1 = paste0(round(1-mean(data0[data0[,5] == 1,4]*(data0[data0[,5] == 1,1]<tau)), 3)*100, '%')
      t.true = paste0(round(mean(data0[,4]*(data0[,1]<tau)), 3)*100, '%')
      t.true0 = paste0(round(mean(data0[A0,4]*(data0[A0,1]<tau)), 3)*100, '%')
      t.true1 = paste0(round(mean(data0[A1,4]*(data0[A1,1]<tau)), 3)*100, '%')
      c.true0 = paste0(round(mean((data0[A0,2] < data0[A0,1])*(data0[A0,2]<tau)), 3)*100, '%')
      c.true1 = paste0(round(mean((data0[A1,2] < data0[A1,1])*(data0[A1,2]<tau)), 3)*100, '%')
      c.true = paste0(round(mean((data0[,2] < data0[,1])*(data0[,2]<tau)), 3)*100, '%')
      c.table[i,j] = paste0(c,'/',c0,'/',c1,'|',t.true0,'/',t.true1, '/',c.true0, '/',c.true1)
    }
  }
  return(c.table)
}


