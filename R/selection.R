selection = function(output, file, burnin, cred=0.95, threshold=0.5, trueGamma=NULL, plotting=FALSE){
  
  library(ggplot2)
  
  gout = output$gout
  gl = gout[(burnin+1):nrow(gout),]
  gl_na = gl
  gl_na[gl_na==0]=NA #change gamma to NA if they are 0
  
  bout = output$bout
  bl = bout[(burnin+1):nrow(bout),]
  
  blgl=bl*gl_na
  
  
  ddff = data.frame(x=c(1:ncol(gout)),z=colMeans(gl))
  name = output$name
  ddff$y = name
  
  bind1 = 1:(which(name=="baseline10")-1)
  bind2 = (which(name=="baseline10")):length(name)
  bind3 = bind4 = NULL
  if(length(which(name=="baseline00"))>0){
    bind2 = (which(name=="baseline10")):(which(name=="baseline00")-1)
    bind3 = (which(name=="baseline00")):(which(name=="baseline11")-1)
    bind4 = (which(name=="baseline11")):length(name)
  }

  
  bindlist=list(bind1,bind2,bind3,bind4)
  mppivec=c("MPPI01","MPPI10","MPPI00","MPPI11")
  
  colnames(ddff) = c("id","MPPI","Terms")
  
  if(!is.null(trueGamma)){
    ddff$truegamma=trueGamma
  }
  
  # mppi plot
  if (plotting){
    
    if(!is.null(file)){
      pdf(file, width=8, height=8)
    }
    
    for(ind in 1:4){
      
      if(!is.null(bindlist[[ind]])){
        p=ggplot(aes(x=Terms,y=MPPI), data=ddff[bindlist[[ind]],])+
          geom_bar(stat="identity", width=0.1)+
          geom_hline(yintercept=threshold, linetype="dashed", color = "black",size=1)+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
          ggtitle(mppivec[[ind]])
        if(!is.null(trueGamma)){
          p=p+geom_point(data=ddff[bindlist[[ind]],], aes(x = Terms, y= truegamma),colour = "red", size = 5)
        }
        print(p)
      }
      


    }
    
    if(!is.null(file)){
      dev.off()
    }

  }
  
  # table of all terms
  ftable = NULL
  lobd = (1-cred)/2
  upbd = 1-lobd
  for (sel in 1:ncol(gl)){
    row = c(name[sel], round(colMeans(gl)[sel],3), round(colMeans(blgl,na.rm = T)[sel],3), round(quantile(blgl[,sel],c(lobd,upbd),na.rm = T),3))
    ftable = rbind(ftable,row)
  }
  colnames(ftable) = c("Terms", "MPPI", "Est", paste0("CI_",lobd*100), paste0("CI_",upbd*100))
  rownames(ftable) = ftable[,1]
  ftable=ftable[,-1]
  ftable=as.data.frame(ftable)
  ftable=data.frame(apply(ftable,2,as.numeric),row.names = rownames(ftable))
  
  # table of selected terms
  table = NULL
  lobd = (1-cred)/2
  upbd = 1-lobd
  for (sel in which(colMeans(gl)>threshold)){
    row = c(name[sel], round(colMeans(gl)[sel],3), round(colMeans(blgl,na.rm = T)[sel],3), round(quantile(blgl[,sel],c(lobd,upbd),na.rm = T),3))
    table = rbind(table,row)
  }
  colnames(table) = c("Terms", "MPPI", "Est", paste0("CI_",lobd*100), paste0("CI_",upbd*100))
  rownames(table) = table[,1]
  table=table[,-1]
  table=as.data.frame(table)
  table=data.frame(apply(table,2,as.numeric),row.names = rownames(table))
  
  mppi = colMeans(gl)
  
  return(list(full_table = ftable, table = table,mppi = mppi))
  
}