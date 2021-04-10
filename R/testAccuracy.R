##################################################
# TEST ACCURACY AFTER RUNNING MCMC

# input:
#   gammaOut - gamma posterior samples, size: iterations x (2pt+2pe)
#   burnin - number of burn-ins
#   index - indicates which gammas we want to look at
#   trueGamma - the true of selection

##################################################

accuracy = function(
  gammaOut,
  burnin,
  index,
  trueGamma
){
  g_out = gammaOut[-(1:burnin),index]
  mppi = colMeans(g_out)
  trueGamma = as.numeric(trueGamma[index]!=0)
  Gindex1 = which(trueGamma==1)
  Gindex0 = which(trueGamma==0)
  pos = which(mppi>0.5)
  neg = which(mppi<0.5)
  FP = sum(pos %in% Gindex0)
  FN = sum(neg %in% Gindex1)
  TP = sum(pos %in% Gindex1)
  TN = sum(neg %in% Gindex0)
  # fpr
  fpr = FP/(FP+TN)
  # fnr
  fnr = FN/(FN+TP)
  
  tpr = 1 - fnr
  tnr = 1 - fpr
  
  # Matthews correlation coefficient
  mcc = (TP*TN - FP*FN)/sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)  )
  
  return(list(TPR = tpr, TNR = tnr, FPR = fpr, FNR = fnr, MCC = mcc))
}