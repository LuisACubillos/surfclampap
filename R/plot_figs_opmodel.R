### Funcion grafica series de tiempo de mediana y con límite inferior y superior
plot_ts <- function(data=syr,data_ev=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500){
  dfp <- data
  Fstr <- rep(Str,nsim)
  dfp <- as.data.frame(dfp)
  colnames(dfp) <- yrs
  dfp <- cbind(dfp,Str=Fstr)
  dfe <- as.data.frame(data_ev)
  colnames(dfe) <- yr
  df <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfp[dfp$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    tmp_df <- t(apply(tmp_df[,1:(dim(tmp_df)[2]-1)],2,quantile,c(0.025,0.5,0.975),na.rm=TRUE))
    df <- rbind(df,tmp_df)
  }
  colnames(df) <- c("Li","Mediana","Ls")
  iyr <- c(yr,yrs)
  #hcr <- c(rep(Str[1],length(iyr)),rep(Str[2],length(iyr)),rep(Str[3],length(iyr)),rep(Str[4],length(iyr)),rep(Str[5],length(iyr)),rep(Str[6],length(iyr)))
  iyr2 <- rep(iyr,nstr)
  hcr <- NULL
  for(i in 1:nstr){
    tmp_hcr <- rep(Str[i],length(iyr))
    hcr <- c(hcr,tmp_hcr)
  }
  df <- data.frame(df,Year=iyr2,Str=hcr)

  p <- ggplot(data=df,aes(x=Year,y=Mediana,group=Str))+
    geom_ribbon(aes(ymin=Li,ymax=Ls),fill="grey")+
    geom_rect(xmin = yr[1],
              xmax = yr[length(yr)],
              ymin = ymin, ymax = max(df$Ls,na.rm = TRUE),
              fill = "lightblue", alpha = 0.01)+
    geom_line()+facet_wrap(~Str,ncol=nstr)+
    scale_y_continuous(name=ylabel,breaks=seq(ymin,max(df$Ls,na.rm = TRUE),by=ytick),limits = c(ymin,max(df$Ls,na.rm=T)))+mi.tema()
  p
}

#### Grafica boxplot
boxplot_ts <- function(data=syr,data_ev=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500,breaks_x =c(2011,2020,2030),yref=Srms){
  dfp <- data
  Fstr <- rep(Str,nsim)
  dfp <- as.data.frame(dfp)
  colnames(dfp) <- yrs
  dfp <- cbind(dfp,Str=Fstr)
  dfe <- as.data.frame(data_ev)
  colnames(dfe) <- yr
  df <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfp[dfp$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    df <- rbind(df,tmp_df)
  }
  df2 <- melt(df)
  colnames(df2) <- c("u","Year","Var")
  
  p <- ggplot(aes(factor(Year),Var),data=df2)+geom_boxplot(outlier.shape = NA)+facet_wrap(vars(u),nrow=1)+
    xlab("Años")+ ylab(ylabel) + geom_hline(yintercept=yref,linetype="dashed")+ theme(axis.title.y = element_text(size = rel(0.1), angle = 90))+ theme(axis.title.x = element_text(size = rel(0.1), angle = 00))+scale_x_discrete(breaks = breaks_x)+mi.tema()
  p
}
######## Gafica una muestra de realizaciones
plot_ts_fideos <- function(data=syr,data_ev=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500,breaks_x =c(2015,2020,2025,2030,2035),yref=Srms,nview=6){
  dfp <- data
  Fstr <- rep(Str,nsim)
  dfp <- as.data.frame(dfp)
  colnames(dfp) <- yrs
  dfp <- cbind(dfp,Str=Fstr)
  dfe <- as.data.frame(data_ev)
  colnames(dfe) <- yr
  df <- NULL
  for(i in 1:nstr){
    tmp_dfp <- dfp[dfp$Str==Str[i],,] 
    tmp_df <- cbind(dfe,tmp_dfp)
    df <- rbind(df,tmp_df)
  }
  selec_sim <- NULL
  for(i in 1:nstr){
    ind1 <- nsim*i-nsim+1
    ind2 <- nsim*i
    tmp_sel <- floor(runif(1,min=ind1,max=ind2))
    selec_sim <- c(selec_sim,tmp_sel)
  }
  #selec_sim <- floor(runif(nview,min=1,max=nsim))
  #selec_sim <- c(3,505,1015,)
  df2 <- df[selec_sim,]
  df2 <- melt(df2)
  colnames(df2) <- c("u","Year","Var")
  
  p <- ggplot(data=df2)+
    geom_line(aes(x=Year,y=Var,group=u,linetype=u))+
    xlab("Años")+
    ylab(ylabel)+
    geom_hline(yintercept=yref,linetype="dashed")+
    theme(axis.title.y = element_text(size = rel(0.1), angle = 90))+
    theme(axis.title.x = element_text(size = rel(0.1), angle = 00))+scale_x_discrete(breaks = breaks_x)+mi.tema()
  p
}

plot_prob_status <- function(data=nrpr,depl=0.2,targt=0.4,Str=Str,nstr=nstr,nsim=nsim,breaks_x=c(2015,2020,2025,2030,2035)){
  df_tmp <- data
  Fstr <- rep(Str,nsim)
  dfp_tmp <- as.data.frame(df_tmp)
  colnames(dfp_tmp) <- yrs
  dfp <- cbind(dfp_tmp,Str=Fstr)
  dfd <- NULL
  dft <- NULL
  nstr=6
  for(i in 1:nstr){
    tmp_dfd <- dfp[dfp$Str==Str[i],,]
    tmp_dfd <- tmp_dfd[,1:length(yrs)]
    tmp_dfd <- t(apply(tmp_dfd, 2,function(x) sum(ifelse(x<depl,1,0))/nsim))
    dfd <- rbind(dfd,tmp_dfd)
    tmp_dft <- dfp[dfp$Str==Str[i],,]
    tmp_dft <- tmp_dft[,1:length(yrs)]
    tmp_dft <- t(apply(tmp_dft, 2,function(x) sum(ifelse(x>targt,1,0))/nsim))
    dft <- rbind(dft,tmp_dft)
  }
  dat_prob <- data.frame(dfd)
  colnames(dat_prob) <- yrs
  dat_prob$u <- Str
  dat_prob_depletion <- melt(dat_prob)
  p1 <- ggplot(dat_prob_depletion)+xlab("Year")+
    geom_line(aes(x=variable,y=value,group=u,linetype=u))+
    mi.tema()+theme(legend.position = c(0.8,0.8))+
    scale_x_discrete(breaks = breaks_x)+
    scale_y_continuous(name="Probability",breaks=seq(0,1,by=0.2),limits = c(0,1))
  
  dat_prob <- data.frame(dft)
  colnames(dat_prob) <- yrs
  dat_prob$u <- Str
  dat_prob_target <- melt(dat_prob)
  p2 <- ggplot(dat_prob_target)+xlab("Year")+
    geom_line(aes(x=variable,y=value,group=u,linetype=u))+
    mi.tema()+theme(legend.position = "none")+
    scale_x_discrete(breaks = breaks_x)+
    scale_y_continuous(name="Probability",breaks=seq(0,1,by=0.2),limits = c(0,1))
  
  p3 <- ggarrange(p1,p2,labels = c("A","B"),nrow=1,ncol = 2)
}

