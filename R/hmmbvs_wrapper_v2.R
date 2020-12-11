HMMbvs_R = function(data, 
                    tcova=NULL,
                    tforce=NULL,
                    ecova=NULL,
                    eforce=NULL,
                    standardize=NULL,
                    model="HMM", 
                    init="baseline", 
                    initvalue = NULL,
                    iter=5000, 
                    v1=5, 
                    v2=1, 
                    a=1, 
                    b=9, 
                    thin = 10, 
                    thin_hidden=10){
  
  library(msm)
  
  data = data[order(data$id,data$obstime),]
  data$delta = time2delta(data$obstime)
  
  stdcol = which(colnames(data) %in% standardize)
  data[stdcol] = scale(data[stdcol])
  
  # Preprocess the input and the dataset
  if(length(tcova)!=2){
    tcova = list("t01"=tcova,
                 "t10"=tcova)
  }
  tran01 = tcova$t01
  tran10 = tcova$t10
  fin01 = tforce$tf01
  fin10 = tforce$tf10
  xt = union(tran01,tran10)
  datat = as.data.frame(data)
  datat = datat[xt]
  datat = cbind(1,datat)
  tran01 = c("int",tran01)
  tran10 = c("int",tran10)
  xt = union(tran01,tran10)
  
  if(length(ecova)!=2){
    ecova = list("e00"=ecova,
                 "e11"=ecova)
  }
  emis00 = ecova$e00
  emis11 = ecova$e11
  fin00 = eforce$ef00
  fin11 = eforce$ef11
  xe = union(emis00,emis11)
  datae = as.data.frame(data)
  datae = datae[xe]
  datae = cbind(1,datae)
  emis00 = c("int",emis00)
  emis11 = c("int",emis11)
  xe = union(emis00,emis11)
  
  id = data$id
  y = data$y
  delta = data$delta
  datmat = cbind(id,datat,datae,y,delta)
  
  fout01 = setdiff(xt,tran01)
  fout10 = setdiff(xt,tran10)
  fout00 = setdiff(xe,emis00)
  fout11 = setdiff(xe,emis11)
  
  pt = length(xt)
  pe = length(xe)
  klist0 = seq(0,2*(pt+pe)-1)
  
  fixin = c(2+c(0,pt,2*pt,2*pt+pe),
            which(colnames(datmat) %in% fin01),
            pt + which(colnames(datmat) %in% fin10),
            pt + which(colnames(datmat) %in% fin00),
            pt+pe + which(colnames(datmat) %in% fin11))-1
  
  fixout = c(which(colnames(datmat) %in% fout01),
             pt + which(colnames(datmat) %in% fout10),
             pt + which(colnames(datmat) %in% fout00),
             pt+pe + which(colnames(datmat) %in% fout11))-1
  
  klist = klist0[-which(klist0 %in% sort(c(klist0[c(fixin,fixout)])))]
  
  
  if(init=="baseline"){
    
    InitGamma = InitBeta = rep(0,2*pt+2*pe)
    InitGamma[fixin] = 1
    InitBeta[c(0,pt,2*pt,2*pt+pe)+1]=1
    
  } else if (init=="warmstart"){
    
    obsTime = delta2time(data)
    
    cav = cbind(data,obsTime)
    cav = as.data.frame(cav)
    st = statetable.msm(y, subject=id, data=cav)
    twoway4.q <- rbind(c(0, st[3]/(st[1]+st[3])), c(st[2]/(st[2]+st[4]), 0)) # qmatrix
    rownames(twoway4.q) <- colnames(twoway4.q) <- c("Non-smoke", "Smoke")
    
    eq01 = genEq(tran01,hasint = TRUE)
    eq10 = genEq(tran10,hasint = TRUE)
    eq00 = genEq(emis00,hasint = TRUE)
    eq11 = genEq(emis11,hasint = TRUE)
    
    misc.msm <- msm(y ~ obsTime, subject = id, data = cav, qmatrix = twoway4.q,
                    covariates = list("1-2" = eq01, "2-1" = eq10),
                    hcovariates = list(eq00,eq11),
                    hmodel = list(hmmBinom(size=1, prob=0.2),hmmBinom(size=1, prob=0.8)),control = list(maxit = 20000,reltol = 1e-16))
    
    nest = length(misc.msm$estimates)
    start11 = nest-(length(emis11)-2)
    end00 = start11-1
    start00 = end00-(length(emis00)-2)
    
    index1 = start00:end00
    index2 = start11:nest
    int1 = (misc.msm$estimates[start00-3])
    int2 = (misc.msm$estimates[start00-1])
    state0 = -c(int1,misc.msm$estimates[index1])
    state1 = c(int2,misc.msm$estimates[index2])
    
    index00 = c(1+pt+pt,pt+which(colnames(datmat) %in% emis00)-1)
    index11 = c(1+pt+pt+pe,pt+pe+which(colnames(datmat) %in% emis11)-1)
    
    InitGamma = InitBeta = rep(0,2*pt+2*pe)
    InitGamma[fixin] = 1
    InitBeta[index00] = state0
    InitBeta[index11] = state1
    
  } else if (init=="manual"){
    
    InitGamma = InitBeta = rep(0,2*pt+2*pe)
    InitGamma[fixin] = 1
    
    index01 = which(colnames(datmat) %in% tran01)-1
    index10 = pt+which(colnames(datmat) %in% tran10)-1
    index00 = pt+which(colnames(datmat) %in% emis00)-1
    index11 = pt+pe+which(colnames(datmat) %in% emis11)-1
    
    InitBeta[
      c(1,
        index01,
        pt+1,
        index10,
        2*pt+1,
        index00,
        2*pt+pe+1,
        index11)
      ]=c(initvalue$init01,initvalue$init10,initvalue$init00,initvalue$init11)
    
  }
  
  
  # get obs length for each individual
  len = NULL
  for (i in 1:(nrow(datmat)-1)){
    if (datmat$id[i]!=datmat$id[i+1]){
      len = c(len,i)
    }
  }
  
  startIndex = c(0,len) # index for firstobs for each subject
  endIndex = c(len,nrow(datmat)) # index for last obs for each subject
  subjLen = endIndex-startIndex # number of observation for each subject
  
  ishmm = (model=="HMM")
  if(!ishmm){thin_hidden=iter/thin}

  output = HMMbvs(as.matrix(datmat), InitBeta, InitGamma, nobs=subjLen, sp=startIndex, klist=klist, hmm=ishmm, iterations=iter, v1=v1, v2=v2, pt=pt, pe=pe, a=a, b=b, thin=thin, thin_hidden=thin_hidden)
  
  rm.index = c(which(colnames(datmat) %in% fout01)-1,
               pt+which(colnames(datmat) %in% fout10)-1,
               pt+which(colnames(datmat) %in% fout00)-1,
               pt+pe+which(colnames(datmat) %in% fout11)-1)
  
  bout = output[[1]]
  gout = output[[2]]
  hout = output[[3]]
  ll = output[[4]]
  
  cova_name = c("baseline01",
                paste0("t01_",colnames(datmat)[3:(3+pt-2)]),
                "baseline10",
                paste0("t10_",colnames(datmat)[3:(3+pt-2)]),
                "baseline00",
                paste0("e00_",colnames(datmat)[(3+pt):(3+pt+pe-2)]),
                "baseline11",
                paste0("e11_",colnames(datmat)[(3+pt):(3+pt+pe-2)]))
  
  if(length(rm.index)>0){
    bout = bout[,-rm.index]
    gout = gout[,-rm.index]
    cova_name = cova_name[-rm.index]
  }
  
  if(!ishmm){
    cova_name=cova_name[1:(2*pt)]
    bout=bout[,1:(2*pt)]
    gout=gout[,1:(2*pt)]
  }
  
  hmmout = list()
  hmmout$bout = bout
  hmmout$gout = gout
  if(ishmm){hmmout$hout = hout}
  hmmout$ll = ll
  hmmout$name = cova_name
  hmmout$data = datmat
  
  return(hmmout)
}