sim_data = function(
  ind = 150,
  window = 30,
  lambda_0 = 0.5,
  mu_0 = 0.5,
  emisP0 = 0.95,
  emisP1 = 0.95,
  pt = 30,
  pe = 20,
  bval = c(0.7,1.0,1.5),
  cor = 0.0,
  biased = TRUE
  ){
  
  library(reshape)
  library(MASS)
  library(abind)
  
  # ind = 150
  # window = 30
  # emisP0 = 0.95 # approx. probability of telling the truth
  # emisP1 = 0.95
  # pt = 30 # pt is the number of transition coefficients
  # pe = 20 # pe is the number of emission coefficients
  # bval = c(0.7,1.0,1.5) # the absolute value of the coefficients are randomly selected from these 3 values
  
  # set true coefficients
  
  correl = 0


  intercept_0 = 1/(1/emisP0-1)
  intercept_1 = 1/(1/emisP1-1)
  
  beta_lambda = beta_mu = rep(0,pt); beta_0 = beta_1 = rep(0,pe)
  
  beta_lambda[c(1,7,10)] = sample(bval,3,replace = TRUE); beta_lambda[c(2,5,8)] = -sample(bval,3,replace = TRUE)
  beta_mu[c(1,13,15)] = sample(bval,3,replace = TRUE); beta_mu[c(2,9,11)] = -sample(bval,3,replace = TRUE)
  
  beta_0[c(1,7)] = sample(bval,2,replace = TRUE); beta_0[c(2,5)] = -sample(bval,2,replace = TRUE) 
  beta_1[c(1,5)] = sample(bval,2,replace = TRUE); beta_1[c(2,6)] = -sample(bval,2,replace = TRUE)
  
  betatrue = c(beta_lambda,beta_mu,beta_0,beta_1)
  # betatrue = c(log(lambda_0),beta_lambda,log(mu_0),beta_mu,log(intercept_0),beta_0,log(intercept_1),beta_1)
  
  # initialize for later use
  
  allSubjTab = data.frame()
  trueState = NULL
  
  for(id in 1:ind){
    
    # initialize temporary table to store transition information
    tempTrans = data.frame()
    
    # variables to put in transition data frame
    # binary variables
    nBi = 4 # number of binary variables
    biVar = NULL
    for (i in 1:nBi){
      biVar = c(biVar,rbinom(1,1,.5))
    }
    
    # continuous variables
    nCon = length(beta_lambda) - nBi # number of continuous variables
    
    # Make X - Correlation between covariates is cor
    sigma2 <- diag( nCon )
    
    for( i in 1:nCon ){
      for( j in 1:nCon ){
        if( i != j ){
          sigma2[ i , j ] = cor^abs(i - j)
        }
      }
    }
    
    conVar = mvrnorm(n = 1, mu = rep(0,nCon), Sigma = sigma2)
    curr = 1
    prev = NA  # fill this in retrospectively
    t = 0
    delta = 1 # always the case because of when the measurememt times are set (every unit)
    
    # initialize the transition 
    newTrans = c(id, biVar, conVar, curr, prev, delta, t, 1, NA) # vector of new transition for subject i
    tempTrans = rbind(tempTrans, newTrans) # append new transition
    
    # name the transition
    tempName = "ID"
    for (i in 2:(ncol(tempTrans)-6)){
      tempName = c(tempName, paste0("tVar",i-1))
    }
    tempName = c(tempName, "curr", "prev", "delta", "time", "firstobs", "length")
    names(tempTrans) = tempName
    
    trans_count = 0 # prevent too many transitions
    
    # build the transition
    # while(tempTrans$time[nrow(tempTrans)] < window && trans_count < 200){
    while(tempTrans$time[nrow(tempTrans)] < window){
      
      trans_count = trans_count+1
      
      lastRow = nrow(tempTrans) # how many rows does tempTrans currently have?
      
      if(tempTrans$curr[lastRow] == 0){
        tempTrans$prev[lastRow] = 1 # Make previous transition state 1
        # Set rate for transition to 1
        tranIndex = c(2:(ncol(tempTrans)-6))  # 1 is id, last 6 are cur, pre, del, tim, fir, len
        rate = lambda_0 * exp( as.matrix(tempTrans[lastRow, tranIndex]) %*% beta_lambda ) 
      }else if(tempTrans$curr[lastRow] == 1){
        tempTrans$prev[lastRow] = 0 # Make previous transition state 0
        # Set rate for transition to 0
        tranIndex = c(2:(ncol(tempTrans)-6))  # 1 is id, last 6 are cur, pre, del, tim, fir, len
        rate = mu_0 * exp( as.matrix(tempTrans[lastRow, tranIndex]) %*% beta_mu )
      }
      
      tempTrans$length[lastRow] = rexp(1, rate = rate) # record time of transition
      if(lastRow == 10){ biVar[1] = rbinom(1,1,.5) }
      if(lastRow == 10){ biVar[3] = rbinom(1,1,.5) }
      # 2 of the Binary covariates can change value if obs transitions more than 10 times, 
      conVar[c(1,3)] = as.numeric(tempTrans[lastRow,c(6,8)] + mvrnorm(n = 1, rep(0,2), .001*diag(2)))
      # 2 continous covariates that follow random walk
      
      # Update states and time 
      curr = tempTrans$prev[lastRow]   
      prev = NA  # fill this in retrospectively
      t = sum(tempTrans$length)
      delta = 1 # always the case because of when the measurememt times are set (every unit)
      
      # Append the new transition to the temp data frame
      newTrans = c(id, biVar, conVar,
                   curr, prev, delta, t, 0, NA) # vector of new transition for subject i
      tempTrans = rbind(tempTrans, newTrans) # append new transition
    }
    
    # Parse through temp files to collect OBSERVED states with covariates
    # Use the random list of observation times
    obs = data.frame()
    # random_obs = sort(runif(window,0,window))   
    random_obs = sort(runif(window,0,floor(max(tempTrans$time))))
    tmp = cbind(tempTrans[1,],0)
    colnames(tmp)[ncol(tmp)]="i"
    obs = rbind(obs, tmp)
    for(i in random_obs){ # obs time
      for(j in 1:nrow(tempTrans)){ # real transition time
        if( (tempTrans$time[j] <= i) & (i <= tempTrans$time[j+1]) ){
          obs <- rbind(obs, cbind(tempTrans[j,],i))
        }
      }
    }
    names(obs)[ncol(obs)] = "obs_time"
    
    # Adjust next state given observations (not transitions) and adjust delta
    for(i in 2:nrow(obs)){
      obs$prev[i] <- obs$curr[(i-1)]
      obs$delta[i] <- obs$obs_time[i] - obs$obs_time[i-1]
    }  
    obs$prev[1] = NA
    
    # add true state to list  
    trueState = c(trueState,obs$curr)  
    
    # build emission covariate matrix
    nEmis = length(beta_0)-2
    
    # Make X - Correlation between covariates is cor
    sigma2p <- diag( nEmis )
    
    for( i in 1:nEmis ){
      for( j in 1:nEmis ){
        if( i != j ){
          sigma2p[ i , j ] = cor^abs(i - j)
        }
      }
    }
    
    cont2 = mvrnorm(n=1, rep(0,nEmis), sigma2p) 
    # borrowed 2 covariates from transition
    
    emisX = NULL
    for(i in 1:nrow(obs)){
      cont2[2] = cont2[2]+rnorm(1,0,0.1)
      cont2[4] = cont2[4]+rnorm(1,0,0.1)
      emisX = rbind(emisX,cont2)
    }
    emisX = cbind(obs$tVar1, obs$tVar3, emisX) # first two covariates for emission is from transition 
    
    tempName = c("tVar1","tVar3")
    for (i in 3:length(beta_0)){
      tempName = c(tempName, paste0("eVar",i))
    }
    
    colnames(emisX) = tempName
    rownames(emisX) = c(1:nrow(obs))
    
    # build emission
    ytmp = NULL # initialize observation
    for (i in 1:nrow(obs)){
      if (obs$curr[i] == 0){
        p = intercept_0*exp( as.numeric(emisX[i,]) %*% beta_0 ) / (  1 + intercept_0*exp( as.numeric(emisX[i,]) %*% beta_0 ) )
        p = 1-p
        ytmp = c(ytmp,rbinom(1,1,p)) # 1-p chance observe 1 <==> probability of telling truth is p
      }else if(obs$curr[i] == 1){
        p = intercept_1*exp( as.numeric(emisX[i,]) %*% beta_1 )/  (  1 + intercept_1*exp( as.numeric(emisX[i,]) %*% beta_1 ) )
        ytmp = c(ytmp,rbinom(1,1,p)) # again telling truth probability is p
      }
      
    }
    
    # attach emission covariates and observation to obs data frame
    obs = cbind(obs, emisX, ytmp)
    allSubjTab = rbind(allSubjTab,obs)
    
  }
  
  # clean up the data frame
  # adjust delta accordingly
  for (i in 1:nrow(allSubjTab)){
    if (is.na(allSubjTab$prev[i])){
      allSubjTab$delta[i] = NA
    }
  }
  
  x = allSubjTab
  obstime = x[,(pt+8)]
  id = x$ID
  delta = x$delta
  y = x$ytmp
  colnames(x)[pt+9] = "eVar1"
  colnames(x)[pt+10] = "eVar2"
  ##############################################################################
  ##############################################################################
  
  datmat = cbind(id,x[,2:(pt+1)],x[,(pt+9):(pt+9+pe-1)],y,obstime)
  # datmat = as.matrix(datmat)  
  trueBeta = rbind(c("baseline01",log(0.5)),
                   cbind(paste0("t01_",colnames(datmat)[2:(pt+1)]),betatrue[1:pt]),
                   c("baseline10",log(0.5)),
                   cbind(paste0("t10_",colnames(datmat)[2:(pt+1)]),betatrue[(pt+1):(2*pt)]),
                   c("baseline00",1/(1/emisP0-1)),
                   cbind(paste0("e00_",colnames(datmat)[(pt+2):(pt+pe+1)]),betatrue[(2*pt+1):(2*pt+pe)]),
                   c("baseline11",1/(1/emisP1-1)),
                   cbind(paste0("e11_",colnames(datmat)[(pt+2):(pt+pe+1)]),betatrue[(2*pt+pe+1):(2*pt+2*pe)]))
  ##############################################################################
  ##############################################################################
  
  if(!biased){
    datmat[,pt+pe+2] = trueState
    # datmat = datmat[,-c((pt+2):(pt+pe+1))]
    # trueBeta = trueBeta[-c((2*pt+3):(2*pt+2*pe+4)),]
  }
  
  out = list(DATA = datmat, BETA = trueBeta, STATE = trueState, DELTA = delta)
  return(out)
}
