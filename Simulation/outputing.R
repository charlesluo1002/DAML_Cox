source('simulation.R')
output = function(results){
  bias = formatC(colMeans(results$beta_estimates) - results$beta_true, format = 'f', digits= 3)
  sd = formatC(apply(results$beta_estimates, 2, sd), format = 'f', digits= 3)
  bias.sd = paste0(bias, ' (',sd, ')')
  model_se = formatC(apply(results$model_se, 2, mean), format = 'f', digits= 3)
  boot_se = formatC(apply(results$boot_se, 2, mean), format = 'f', digits= 3)
  rmse = formatC(sqrt(colMeans((results$beta_estimates - results$beta_true)^2)), format = 'f', digits= 3)
  boot_cp = formatC(colMeans(results$boot.coverred), format = 'f', digits= 2)
  model_cp = formatC(colMeans(results$model.coverred), format = 'f', digits= 2)
  se = paste0(model_se,'/',boot_se)
  cp = c(0,paste0(model_cp, '/',boot_cp))
  return(list(beta_true = formatC(results$beta_true, format = 'f', digits= 3), out = rbind(bias, sd, model_se, boot_se, model_cp, boot_cp)))
}

latex_out = function(results1, results2, results3, results4, n0, p0, runs0){
  q2=output(results1)
  q3=output(results2)
  q4=output(results3)
  q5=output(results4) 
  # latex outputing  for q2 and q3
  res1 = rbind(t(q2$out), t(q3$out))
  res1 = apply(res1,c(1,2), function(x) ifelse(substr(x,1,1) == '-', paste0(x, '~'), paste0('~',x, '~')))
  for (i in 1:nrow(res1)){
    n.est = ncol(q2$out)
    if (i == 1) {
      cat('\n\n')
      cat(r'(\begin{table}[htbp] \begin{center} \begin{tabular}{llllllll} \hline )')
      cat(r'(\multicolumn{1}{c}{$T/C$ distribution} & \multicolumn{1}{c}{Estimators} & \multicolumn{1}{c}{Bias} & \multicolumn{1}{c}{SD} & \multicolumn{1}{c}{Model SE} & \multicolumn{1}{c}{Boot SE} & \multicolumn{1}{c}{Model CP} & \multicolumn{1}{c}{Boot CP} \\  )')
      cat(r'(\hline   \multirow{11}{*}{fail1-cens1} )')
    }
    if (i == 1 + n.est){
      cat(r'(\hline   \multirow{11}{*}{fail1-cens2} )')
    }
    string = paste0(r'(\multicolumn{1}{c}{)', res1[i,], '}',  collapse = ' & ')
    cat(r'( & )', paste0(colnames(q2$out)[ifelse(i>n.est,i-n.est,i)], r'( & )',string,' \\\\ ', ifelse(i %in% c(1, 7, 1+n.est, 7+n.est), '\\cline{2-8}', ''), ' \n'))
    if (i == nrow(res1)){
      cat(paste0(r'(\hline \end{tabular} \end{center} \captionsetup{width=.9\linewidth}  \caption{Simulation based on )', runs0, ' data sets. $(n,p) = (', n0, ',',p0, r'()$. Failure time $T$ follows distribution from fail1 scenario with $\beta^* = )', q2$beta_true, r'($.} \end{table})'))
    }
  }
  cat('\n\n\n\n')
  # latex outputing fr q4 and q5
  res1 = rbind(t(q4$out), t(q5$out))
  res1 = apply(res1,c(1,2), function(x) ifelse(substr(x,1,1) == '-', paste0(x, '~'), paste0('~',x, '~')))
  for (i in 1:nrow(res1)){
    n.est = ncol(q4$out)
    if (i == 1) {
      cat(r'(\begin{table}[htbp] \begin{center} \begin{tabular}{llllllll} \hline )')
      cat(r'(\multicolumn{1}{c}{$T/C$ distribution} & \multicolumn{1}{c}{Estimators} & \multicolumn{1}{c}{Bias} & \multicolumn{1}{c}{SD} & \multicolumn{1}{c}{Model SE} & \multicolumn{1}{c}{Boot SE} & \multicolumn{1}{c}{Model CP} & \multicolumn{1}{c}{Boot CP} \\  )')
      cat(r'(\hline   \multirow{11}{*}{fail2-cens1} )')
    }
    if (i == 1 + n.est){
      cat(r'(\hline   \multirow{11}{*}{fail2-cens2} )')
    }
    string = paste0(r'(\multicolumn{1}{c}{)', res1[i,], '}',  collapse = ' & ')
    cat(r'( & )', paste0(colnames(q4$out)[ifelse(i>n.est,i-n.est,i)], r'( & )',string,' \\\\ ', ifelse(i %in% c(1, 7, 1+n.est, 7+n.est), '\\cline{2-8}', ''), ' \n'))
    if (i == nrow(res1)){
      cat(paste0(r'(\hline \end{tabular} \end{center} \captionsetup{width=.9\linewidth}  \caption{Simulation based on )', runs0, ' data sets. $(n,p) = (', n0, ',',p0, r'()$. Failure time $T$ follows distribution from fail2 scenario with $\beta^* = )', q4$beta_true, r'($.} \end{table})'))
      cat('\n')
    }
  }
}


############
path = ''


IPW_only = F
runs0 = 1000
n0 = 500
p0 = 2
tau0 = 1
run_boot_se0 = rep(F,11)  #c(F,F,F,F,F,T,T,T,T,T,F)
boot_runs0 = 200
bayesian_boot = T
cores = ifelse(startsWith(getwd(), 'C'), 1, 16)
tic()
# T_AZ_sim='fail1';C_AZ_sim='cens1';runs = runs0; n_sim = n0; p = p0; tau = tau0; boot_runs = boot_runs0; run_boot_se = run_boot_se0
results1 = one_setting_simulation(runs = runs0, n_sim = n0, p = p0, tau = tau0, T_AZ_sim='fail1',C_AZ_sim='cens1',boot_runs = boot_runs0, run_boot_se = run_boot_se0, bayesian_boot = bayesian_boot)
results2 = one_setting_simulation(runs = runs0, n_sim = n0, p = p0, tau = tau0, T_AZ_sim='fail1', C_AZ_sim='cens2',boot_runs = boot_runs0, run_boot_se = run_boot_se0, bayesian_boot = bayesian_boot)
results3 = one_setting_simulation(runs = runs0, n_sim = n0, p = p0, tau = tau0, T_AZ_sim='fail2', C_AZ_sim='cens1',boot_runs = boot_runs0, run_boot_se = run_boot_se0, bayesian_boot = bayesian_boot)
results4 = one_setting_simulation(runs = runs0, n_sim = n0, p = p0, tau = tau0, T_AZ_sim='fail2', C_AZ_sim='cens2',boot_runs = boot_runs0, run_boot_se = run_boot_se0, bayesian_boot = bayesian_boot)
toc()
# save(results1, results2, results3, results4, file = paste0(path,  gsub('-','_', Sys.Date()),'_simul_results_final1.RData'))

latex_out(results1, results2, results3, results4, n0 = n0, p0 = p0, runs0 = runs0)




