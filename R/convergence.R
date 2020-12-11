convergence = function(output,file){
  
  library(coda)
  
  ll = output$ll
  gout = output$gout
  postbeta = output$bout
  betaname = output$name
  
  g = NULL
  xlim = NULL
  for(i in 1:dim(gout)[1]){
    xlim = c(xlim,i)
    y = sum(gout[i,])
    g = c(g,y)
  }
  
  if(!is.null(file)){
    pdf(file)
  }
  
  plot(ll,type = "l",xlab = "MCMC Iteration", ylab = "Log-likelihood", main = "Log-likelihood")
  plot(xlim,g,type = "l", main = paste0("Included coefficients"), xlab = "MCMC Iteration", ylab = "Number of Covariates Included")
  
  for(i in 1:ncol(postbeta)){
    plot(postbeta[,i], type="l", main=betaname[i], xlab = "MCMC Iteration", ylab = "Sampled Coefficient")
  }
  
  if(!is.null(file)){
    dev.off()
  }
}





