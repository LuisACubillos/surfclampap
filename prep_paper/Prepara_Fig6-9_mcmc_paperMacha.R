#MSE_macha: Prepara la figura proyecciones
code_dir <- "./R/"
file_path <- list.files(code_dir,full.names = T)
for(f_path in file_path){source(f_path,encoding = "UTF-8")}

library(ggplot2)
library(ggpubr)
library(reshape2)
library(lattice)
theme_set(theme_bw())

#Prepara lectura de datos
md <- read.admb("models/opmacha")
nstr = 6
nsim= 500
nr=nstr*nsim #numero de filas
Str<-c("u=0","u=10","u=15","u=20","u=25","u=30")
yrs <- seq(2017+1,2017+20,1)
yr  <- seq(2011,2017,1)
nyr <- length(yrs)
low <- 0.05
upr <- 0.95
So=md$S0
Srms <- 0.4*So

#Lee los datos:
#Biomasa desovante
syr=as.matrix(read.table(paste0("data/","05Bdes_mcmc.txt"),nrow=nr,fill=T))
seval=as.matrix(read.table(paste0("data/","10SD_historia.txt"),nrow=nr,fill=T))
#Reclutamiento
ryr=as.matrix(read.table(paste0("data/","13Recluta_mcmc.txt"),nrow=nr,fill=T))
reval=as.matrix(read.table(paste0("data/","12Reclutas_historia.txt"),nrow=nr,fill=T))
#Mortalidad por pesca
fyr=as.matrix(read.table(paste0("data/","09Fmort_mcmc.txt"),nrow=nr,fill=T))
feval=as.matrix(read.table(paste0("data/","11Fmort_historia.txt"),nrow=nr,fill=T))
#Capturas
cyr=as.matrix(read.table(paste0("data/","09Fmort_mcmc.txt"),nrow=nr,fill=T))

#Graficos
############## FIGURA 6 RECLUTAMIENTO
# Figura 6: Reclutamientos
Rprom <- mean(md$Reclutas)
f6 <- plot_ts_fideos(data=ryr,data_ev=reval,ylabel="Recruitment (millions)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=5,yref=Rprom,breaks_x =c(2015,2020,2025,2030,2035),nview=10,seed = 4645)
f6
ggsave(paste0("figs/","Fig6_recruits_realization.jpg"),plot=f6,dpi=300,width = 20,height = 16,units = "cm")


############## FIGURA 7 
## Fig. 7 Reclutamiento - Biomasa - Mort. Pesca
f1 <- plot_ts(data=ryr,data_ev=reval,ylabel="Recruitment (millions)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=50,yref=NULL)
f2 <- plot_ts(data=syr,data_ev=seval,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500,yref=Srms)
f3 <- plot_ts(data=fyr,data_ev=feval,ylabel="Fishing mortality",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.2,yref=NULL)

Fig7 <- ggarrange(f1,f2,f3,labels = c("A","B","C"),nrow=3,ncol = 1)
Fig7
ggsave(paste0("figs/","Fig7_RecSSBFishMort.jpg"),plot=Fig7,dpi=300,width = 30,height = 27,units = "cm")


############## FIGURA 8 AGOTAMIENTO
depletion <- seval/So
Target <- 0.4
dyr <- as.matrix(read.table(paste0("data/","06RPR_mcmc.txt"),nrow=nr,fill=T))
f8 <- plot_ts(data=dyr,data_ev=depletion,ylabel="Depletion",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.4,yref = Target)
f8
ggsave(paste0("figs/","Fig8_Depletion.jpg"),plot=f8,dpi=300,width = 30,height = 9,units = "cm")



############## FIGURA 8 AGOTAMIENTO (ALTERNATIVA CON BOXPLOT 
Target <- 0.4
f9 <- boxplot_ts(data=dyr,data_ev=depletion,ylabel="Depletion",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.2,yref=Target,breaks_x =c(2011,2020,2030))
#ggsave(paste0("figs/","Fig8_Depletion_v2.jpg"),plot=f9,dpi=300,width = 30,height = 9,units = "cm")


############# FIGURA 9 Probabilidada de estaus ######
fig9 <- plot_prob_status(data=dyr,depl=0.2,targt=0.4,Str=Str,nstr=nstr,nsim=nsim,breaks_x=c(2015,2020,2025,2030,2035))
fig9
ggsave(paste0("figs/","Fig9_Probability.jpg"),plot=fig9,dpi=300,width = 20,height = 15,units = "cm")


#################   NO USADOS  ##################################
boxplot_ts(data=syr,data_ev=seval,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=200,yref=Srms,breaks_x =c(2011,2020,2030))

Rprom <- mean(md$Reclutas)
p5 <- boxplot_ts(data=ryr,data_ev=reval,ylabel="Recruitment (n)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=5,yref=Rprom,breaks_x =c(2011,2020,2030))