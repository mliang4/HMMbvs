genEq = function(covavec, hasint = TRUE){
  
  # convert a vector of covariate names to a R formula
  # can remove the "int" name in the name vector by setting 'hasint=TRUE'
  
  if(hasint){
    covavec = covavec[covavec != "int"]
  }
  
  output = NULL
  
  if(length(covavec)>0){
    sstr = "~"
    for (item in covavec){
      if (item == covavec[1]){
        sstr = paste0(sstr,item)
      } else {
        sstr = paste0(sstr,"+",item)
      }
    }
    output = as.formula(sstr)
  }

  return(output)
}

