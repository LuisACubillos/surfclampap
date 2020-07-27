plot_ts <- function(data=syr,data_ev=NULL,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500,yref=Target){
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
    geom_hline(yintercept=yref,linetype="dashed")+
    scale_y_continuous(name=ylabel,breaks=seq(ymin,max(df$Ls,na.rm = TRUE),by=ytick),limits = c(ymin,max(df$Ls,na.rm=T)))+mi.tema()
  p
}
