time2delta = function(time){
  
  # convert time (sorted in each id) to time difference (delta):
  
  # time needs to be sorted
  
  deltavec = NA
  deltavec = c(deltavec,time[2:length(time)]-time[1:(length(time)-1)])
  deltavec[which(deltavec<0)] = NA
  return(deltavec)
}
