delta2time = function(hmmDatmat){
  
  # convert time difference (delta) to time (starting from 0):
  # hmmDatmat: must be a dataframe with columns id, delta
  
  delta = as.numeric(hmmDatmat$delta)
  deltalist = c(which(!duplicated(hmmDatmat$id)), nrow(hmmDatmat)+1) # firstobs
  timevec = rep(0,nrow(hmmDatmat))
  for (i in 1:(length(deltalist)-1)){
    from = deltalist[i]+1
    to = deltalist[i+1]-1
    timevec[from:to] = cumsum(hmmDatmat$delta[from:to])
  }
  return(timevec)
}


