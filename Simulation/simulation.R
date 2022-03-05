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
      res8 = ALHR_IPW(data = data, C_AZ = 'cens1' , beta0 = 0, run_boot_se = run_boot_se[8], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res9 = ALHR_IPW(data = data, C_AZ = 'W2' , beta0 = 0, run_boot_se = run_boot_se[9], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res10 = ALHR_IPW(data = data, C_AZ = 'W1' , beta0 = 0, run_boot_se = run_boot_se[10], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
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
      res1 = ALHR_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[1], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res2 = ALHR_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[2], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res3 = ALHR_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'Spline', beta0 = 0, run_boot_se = run_boot_se[3], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res4 = ALHR_AIPW(data = data, T_AZ = 'Spline', C_AZ = 'Spline', beta0 = 0, run_boot_se = run_boot_se[4], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res5 = ALHR_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'cens1', beta0 = 0, run_boot_se = run_boot_se[5], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res6 = ALHR_AIPW(data = data, T_AZ = 'fail1', C_AZ = 'RSF', beta0 = 0, run_boot_se = run_boot_se[6], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res7 = ALHR_AIPW(data = data, T_AZ = 'RSF', C_AZ = 'RSF', beta0 = 0, run_boot_se = run_boot_se[7], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res8 = ALHR_IPW(data = data, C_AZ = 'cens1' , beta0 = 0, run_boot_se = run_boot_se[8], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res9 = ALHR_IPW(data = data, C_AZ = 'W2' , beta0 = 0, run_boot_se = run_boot_se[9], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
      res10 = ALHR_IPW(data = data, C_AZ = 'W1' , beta0 = 0, run_boot_se = run_boot_se[10], boot_runs = boot_runs, bayesian_boot = bayesian_boot, tau = tau)
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


