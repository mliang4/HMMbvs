findtruegamma = function(data=df,output=hmmout){
  trueBetaSel = data$BETA[which(data$BETA[,1] %in% output$name),2]
  return(as.numeric(trueBetaSel!=0))
}