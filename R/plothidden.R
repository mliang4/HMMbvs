plot_hidden=function(

output=hmmout,
data,
burnin=0, # burnin for coefficients, not for hidden states
subject=1,
true_hidden=NULL

){
  
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(dplyr)
  
  bout = output$bout
  hout = output$hout
  hind = seq(from=1,by = nrow(bout)/nrow(hout),length.out = nrow(hout))
  hout = hout[which(hind>burnin),]
  hind = hind[hind>burnin]
  
  nsub=length(unique(data$id))
  truehidden=true_hidden

  if(subject=="all" && is.null(truehidden)){
    hacc=NULL
    simhidden=as.numeric(colMeans(hout)>0.5)
    for(sub in 1:nsub){
      idx=which(data$id==sub)
      yidx=data$y[idx]
      hidx=simhidden[idx]
      hacc=c(hacc,mean(yidx==hidx))
    }
    
    hdf=data.frame(index=c(1:nsub),id=c(1:nsub)[order(hacc,decreasing = T)],accu=sort(hacc,decreasing = T),check.names = F)
    
    p=ggplot(aes(x=index,y=accu),data=hdf)+
      geom_bar(stat="identity", width=0.2)+
      theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))+
      labs(y="Proportion", x = "Subjects (relabeled according to the proportion)")+
      scale_x_continuous(limits = c(0.9,nsub+.1),breaks = round(seq(from = 1,to = nsub,length.out = round(nsub/10))))+
      ylim(-0.05, 1.05)
    print(p)
  }else if(subject=="all" && !is.null(truehidden)){
    hacc=NULL
    for(sub in 1:nsub){
      idx=which(data$id==sub)
      yidx=data$y[idx]
      hidx=truehidden[idx]
      df=data.frame(timepoint=c(1:length(idx)),sim_h=hidx,obs=yidx)
      hacc=c(hacc,mean(yidx==hidx))
    }
    
    hdf=data.frame(index=c(1:nsub),id=c(1:nsub)[order(hacc,decreasing = T)],accu=sort(hacc,decreasing = T),check.names = F)
    
    p=ggplot(aes(x=index,y=accu),data=hdf)+
      geom_bar(stat="identity", width=0.2)+
      theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))+
      labs(y="Proportion", x = "Subjects (relabeled according to the proportion)")+
      scale_x_continuous(limits = c(0.9,nsub+.1),breaks = round(seq(from = 1,to = nsub,length.out = round(nsub/10))))+
      ylim(-0.05, 1.05)
    print(p)
  }else{
    sub=subject
    idx=which(data$id==sub)
    yidx=data$y[idx]
    hidx=colMeans(hout)[idx]
    
    if(!is.null(truehidden)){
      htrueidx=truehidden[idx]
    }
    
    
    df=data.frame(timepoint=c(1:length(idx)),
                  sim_h=as.character(as.numeric(hidx>0.5)),
                  sim_state=hidx, 
                  obs=as.character(yidx))
    
    if(!is.null(truehidden)){
      df=cbind(df,true_h=as.character(htrueidx))
    }
    
                  
    
    tp2=rep(c(1:length(idx)),2)
    var2=as.factor(
      c(
        rep(1,length(idx)),rep(0,length(idx))
      )
    )
    est2=c(hidx,1-hidx)
    df2=data.frame(timepoint=tp2, variable=var2,sim_state=est2)
    
    mycols <- c("0"="darkmagenta", "1"="turquoise")
    
    g0 <- (ggplot(df, aes(x = timepoint, y = obs, fill = obs, col = obs)) + 
             geom_bar(stat = "identity", alpha = I(0.7)) + 
             scale_fill_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1")) +
             scale_color_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1")) +
             theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
             labs(y = "Reported state")) %>% ggplotGrob
    
    if(!is.null(truehidden)){
    g1 <- (ggplot(df, aes(x = timepoint, y = true_h, fill = true_h, col = true_h)) + 
             geom_bar(stat = "identity", alpha = I(0.7)) + 
             scale_fill_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1"),drop=FALSE) +
             scale_color_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1"),drop=FALSE) +
             theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
             labs(y = "True State")) %>% ggplotGrob
    }
    
    g2 <- (ggplot(df, aes(x = timepoint, y = sim_h, fill = sim_h, col = sim_h)) + 
             geom_bar(stat = "identity", alpha = I(0.7)) +
             scale_fill_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1")) +
             scale_color_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1")) +
             theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
             labs(y = "Estimated State")) %>% ggplotGrob
    g3 <- (ggplot(df2, aes(x = timepoint, y = sim_state, col = variable)) + geom_line() +
             scale_color_manual(values = mycols, name = "State", labels = c("0"="0", "1"="1")) +
             theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
             labs(y = "Posterior Prob.")) %>%
      ggplotGrob()
    
    if(!is.null(truehidden)){
    g0$widths <- g1$widths
    grid.arrange(g0, g1, g2, g3, widths = 1, nrow = 4,top = paste0("Plot for subject ",subject))
    }else{
      grid.arrange(g0, g2, g3, widths = 1, nrow = 3,top = paste0("Plot for subject ",subject))
    }
  }
  
}



