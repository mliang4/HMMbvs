ppc = function(output,file,type="sum_one", burnin=0, postsample = NULL){
  suppressMessages(library(bayesplot))
  suppressMessages(library(msm))
  suppressMessages(library(ggplot2)  )
  sigmoid = function(z){
    return(exp(z)/(1+exp(z)))
  }
  
  bout = output$bout
  hout = output$hout
  hind = seq(from=1,by = nrow(bout)/nrow(hout),length.out = nrow(hout))
  hout = hout[which(hind>burnin),]
  hind = hind[hind>burnin]
  beta_post = bout[hind,]
  
  data=output$data
  name=output$name
  name=substr(name,5,10)
  bind1 = 1:(which(name=="line10")-1)
  bind2 = (which(name=="line10")):(which(name=="line00")-1)
  bind3 = (which(name=="line00")):(which(name=="line11")-1)
  bind4 = (which(name=="line11")):length(name)
  xind3 = c(2,which(colnames(data)%in%name[bind3]))  # column of intercept and emis00 covariates
  xind4 = c(2,which(colnames(data)%in%name[bind4]))  # column of intercept and emis11 covariates
  
  uidvec = unique(data$id)
  ymat = rep(0,length(data$y))
  
  if(is.null(postsample)){
    message("Sampling from posterior predictive distribution, It might take a while...")
    
    for(sample in 1:nrow(beta_post)){
      
      z_temp0 = as.matrix(data[,xind3]) %*% beta_post[sample,bind3]
      z_temp1 = as.matrix(data[,xind4]) %*% beta_post[sample,bind4]
      p0 = 1 - sigmoid(z_temp0)
      p1 = sigmoid(z_temp1)
      y_temp = rbinom(length(p0),1,as.numeric(p0))
      y_temp = cbind(y_temp, rbinom(length(p1),1,as.numeric(p1)))
      
      yvec = NULL
      for(idx in 1:ncol(hout)){
        yvec = c(yvec,y_temp[idx,1+hout[sample,idx]])
      }
      
      if(sample %% 10 ==0){
        message(paste0(round(sample/nrow(beta_post)*100),"%"))
      }
      # 
      # yvec = NULL
      # hvec = hout[sample,]
      # for (i in 1:length(uidvec)){
      #   
      #   uid = uidvec[i]
      #   data0 = data[which(data$id==uid),]
      #   
      #   hvec_single = hvec[which(data$id==uid)]
      #   yvec_single = NULL
      #   for (t in 1:length(data0$y)){
      #     if (hvec_single[t] == 0){
      #       z = sum(data0[t,xind3] * beta_post[sample,bind3])
      #       p = 1-sigmoid(z) # emission probability to 1
      #       yvec_single = c(yvec_single,rbinom(1,1,p))
      #     }else{
      #       z = sum(data0[t,xind4] * beta_post[sample,bind4])
      #       p = sigmoid(z) # emission probability to 1
      #       yvec_single = c(yvec_single,rbinom(1,1,p))
      #     }
      #   }
      #   yvec = c(yvec,yvec_single)
      # }
      ymat = cbind(ymat,yvec)
    
    } # end for sample in 1:nrow(beta_post)
    
    # save(yrep,file="ppc_prequit.RData")
    yrep = ymat[,-1]
    yrep = t(yrep)
  }else{
    yrep = postsample
  }
  
  y = data$y
  
  if(!is.null(file)){
    pdf(file)
  }
  
  if("sum_one" %in% type){
    # total number of smoking state predicted
    py = sum(y)
    pyrep = rowSums(x = yrep)
    # pyrep_med = quantile(pyrep,0.5)
    pval = 2*min(sum(pyrep>py),sum(pyrep<py))/length(pyrep)
    sum_one <- function(x) sum(x == 1)
    g=ppc_stat(y, yrep, stat = "sum_one")+
      legend_none() +
      xaxis_title(size = 13, family = "sans") +
      labs(x="Count", title="Total number of state '1' predicted", subtitle=paste0("Two-sided p-value = ", round(pval,3)))
    suppressMessages(print(g))

  }
  yy = NULL
  yyrep = matrix(ncol = length(unique(data$id)), nrow = dim(yrep)[1])
  for (uid in uidvec){
    uidx = which(data$id==uid)
    yy = c(yy,mean(y[uidx]==1))
  }
  for (j in 1:dim(yrep)[1]){
    ytmp = NULL
    for (uid in uidvec){
      uidx = which(data$id==uid)
      ytmp = c(ytmp,mean(yrep[j,uidx]==1))
    }
    yyrep[j,] = ytmp
  }

  if("avg" %in% type){
    py = mean(yy)
    pyrep = apply(yyrep,1,mean)
    pval = 2*min(sum(pyrep>py),sum(pyrep<py))/length(pyrep)
    # mean of % of smoking states for individual
    g=ppc_stat(yy, yyrep, stat = "mean")+
      legend_none() +
      xaxis_title(size = 13, family = "sans") +
      labs(x="Value", title="The mean of averaged % of state '1' per subject", subtitle=paste0("Two-sided p-value = ", round(pval,3)))
    suppressMessages(print(g))
  }
  
  if("var" %in% type){
    bern_var = function(x){
      p = mean(x)
      return(p*(1-p))
    }
    py = bern_var(yy)
    pyrep = apply(yyrep,1,bern_var)
    pval = 2*min(sum(pyrep>py),sum(pyrep<py))/length(pyrep)
    g=ppc_stat(yy, yyrep, stat = "bern_var")+
      legend_none() +
      xaxis_title(size = 13, family = "sans") +
      labs(x="Value", title="The variance of averaged % of state '1' per subject", subtitle=paste0("Two-sided p-value = ", round(pval,3)))
    suppressMessages(print(g))

  }
  
  table = statetable.msm(y, subject=id, data=data)
  ns = table[1,2]
  ss = table[2,2]
  nn = table[1,1]
  sn = table[2,1]
  yrep_msm = cbind(t(yrep), data$id)
  colnames(yrep_msm)[dim(yrep_msm)[2]] = "id"
  nsvec = snvec = nnvec = ssvec = NULL
  
  message("Getting transition info, might take a while...")
  for (i in 1:dim(yrep)[1]){
    table = statetable.msm(yrep[i,], subject=id, data = yrep_msm)
    ns = table[1,2]
    nn = table[1,1]
    ss = table[2,2]
    sn = table[2,1]
    nsvec = c(nsvec,ns)
    nnvec = c(nnvec,nn)
    snvec = c(snvec,sn)
    ssvec = c(ssvec,ss)
  }
  
  ppc_trans = function(py,pyrep,state1,state2){
    pval = 2*sum(pyrep>py)/length(pyrep)
    if (pval>1){
      pval = 2*min(sum(pyrep>py),sum(pyrep<py))/length(pyrep)
    }
    g=ppc_stat(as.vector(py), as.matrix(pyrep))+
      legend_none() +
      xaxis_title(size = 13, family = "sans") +
      labs(x="Count", title=paste0("Total number of transitions from ",state1," to ",state2," predicted"), subtitle=paste0("Two-sided p-value = ", round(pval,3)))
    suppressMessages(print(g))
  }
  
  
  if("tran01" %in% type){
    ppc_trans(ns,nsvec,"0", "1")
  }
  
  if("tran10" %in% type){
    ppc_trans(sn,snvec,"1", "0")
  }
  
  if("tran00" %in% type){
    ppc_trans(nn,nnvec,"0", "0")
  }
  
  if("tran11" %in% type){
    ppc_trans(ss,ssvec,"1", "1")
  }

  if(!is.null(file)){
    dev.off()
  }
  
  if(is.null(postsample)){
    return(yrep)
  }
  
}