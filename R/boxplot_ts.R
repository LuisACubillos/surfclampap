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
