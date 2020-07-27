#source(paste0('R/','my.theme.ggplot.R')) #personalizacion ggplot2
#source(paste0('R/','read.admb.R'))
#source(paste0('R/','plot_figs_opmodel.R'))
code_dir <- "./R/"
file_path <- list.files(code_dir,full.names = T)
for(f_path in file_path){source(f_path,encoding = "UTF-8")}

library(ggplot2)
library(ggpubr)
library(reshape2)
library(lattice)
theme_set(theme_bw())

md <- read.admb("models/opmacha")
# 01Btota_mcmc.txt: Biomasa total
# 02Bvuln_mcmc.txt: Biomasa Explotable
# 03Bcru_mcmc.txt : Biomasa evaluada directa proyectada
# 04BcruVuln_mcmc.txt: Biomasa evaluada explotable
# 05Bdes_mcmc.txt
# 06RPR_mcmc.txt
# 07CTP_mcmc.txt
# 08YTP_mcmc.txt
# 09Fmort_mcmc.txt

#Lee los datos
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
########################3
#Fitted to length composition
pobs_cut <- md$pobs_cru
yrs_cut <- yr
pobs <- as.data.frame(pobs_cut)
talla <- md$Tallas
colnames(pobs)<- talla
pobs2 <- melt(pobs)
pobs2$yrs <- rep(yrs_cut,length(talla))
colnames(pobs2) <- c("talla","frq","yrs")
ppred_cut <- md$ppred_cru
ppred <- as.data.frame(ppred_cut)
talla <- md$Tallas
colnames(ppred)<- talla
ppred2 <- melt(ppred)
ppred2$yrs <- rep(yrs_cut,length(talla))
colnames(ppred2) <- c("talla","frq","yrs")

pobs2$Data <- rep("Observed",length(pobs2$talla))
head(pobs2)
ppred2$Data <- rep("Predicted",length(ppred2$talla))
head(ppred2)
comptalla <- rbind(pobs2,ppred2)
head(comptalla)
p3 <- ggplot(data=pobs2,aes(x=talla,y=frq))+
  geom_col(data=pobs2,aes(x=talla,y=frq))+
  facet_wrap(~yrs,ncol=1)+ xlab("Length (mm)")+
  scale_x_discrete(breaks=c(5, 25,50,75,95))+mi.tema()
#p3 <- p3 + geom_col(data=ppred2,aes(x=talla,y=frq,fill="Fitted"))#

p4 <- ggplot(data=comptalla)+
  geom_col(aes(x=talla,y=frq,fill=Data),position = "dodge")+mi.tema()+scale_x_discrete(breaks=c(5, 25,50,75,95))+ylab("Proportion")+xlab("Length (mm)")+
  facet_wrap(~yrs,ncol=1)+scale_colour_grey()

ggsave(paste0("figs/","Fig04_Ajuste_Tallas.png"),p4,height=25,width=12,units="cm",dpi=300)
#################################
Total <- md$Bt
Spawning <- md$SD
Vulnerable <- md$Bv
Catch <- md$Cobs
df <- data.frame(Total,Spawning,Vulnerable,Catch)
#df <- melt(df[,c("bt","bd","bv","yt")],id=yr)
df <- melt(df)
df$yr <- rep(yr,4)
p2 <- ggplot(df)+
  geom_line(aes(x=yr,y=value,linetype=variable))+
  ylab("Biomass (t)")+xlab("Year")+mi.tema()+theme(legend.position = c(0.8,0.9))
p2

# Reclutamiento
Recruits <- md$Reclutas
Rec_devs <- md$log_desv_Rt
Fishing <- md$fyr
df3 <- data.frame(Recruits,Rec_devs,Fishing,yr)

p3b <- ggplot(data=df3)+
  geom_col(aes(x=yr,y=Recruits))+mi.tema()+
  geom_hline(yintercept = md$R0)+
  ylab("Recruits (millions)")+xlab("Year")
p3b

p3c <- ggplot(data=df3)+
  geom_line(aes(x=yr,y=Fishing))+
  geom_point(aes(x=yr,y=Fishing))+
  mi.tema()+
  ylab("Fishing mortality")+xlab("Year")
p3c

Fig5 <- ggarrange(p2,p3b,p3c,labels = c("A","B","C"),nrow=1,ncol = 3)
Fig5
ggsave(paste0("figs/","Fig5_PopIndicadores.jpg"),plot=Fig5,dpi=300,width = 30,height = 18,units = "cm")



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
# Figura 6: Reclutamientos
Rprom <- mean(md$Reclutas)
f6 <- plot_ts_fideos(data=ryr,data_ev=reval,ylabel="Recruitment (n)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=5,yref=Rprom,breaks_x =c(2015,2020,2025,2030,2035),nview=10)
f6
ggsave(paste0("figs/","Fig6_recruits_realization.jpg"),plot=f6,dpi=300,width = 20,height = 16,units = "cm")

## Fig. 7 Reclutamiento - Biomasa - Mort. Pesca
f1 <- plot_ts(data=ryr,data_ev=reval,ylabel="Recruitment (n)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=50)
f2 <- plot_ts(data=syr,data_ev=seval,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=500)
f3 <- plot_ts(data=fyr,data_ev=feval,ylabel="Fishing mortality",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.2)

Fig7 <- ggarrange(f1,f2,f3,labels = c("A","B","C"),nrow=3,ncol = 1)
Fig7
ggsave(paste0("figs/","Fig7_RecSSBFishMort.jpg"),plot=Fig7,dpi=300,width = 30,height = 27,units = "cm")

depletion <- seval/So
dyr <- as.matrix(read.table(paste0("data/","06RPR_mcmc.txt"),nrow=nr,fill=T))
f8 <- plot_ts(data=dyr,data_ev=depletion,ylabel="Depletion",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.2)
f8
ggsave(paste0("figs/","Fig8_Depletion.jpg"),plot=f8,dpi=300,width = 30,height = 9,units = "cm")

Target <- 0.4
f9 <- boxplot_ts(data=dyr,data_ev=depletion,ylabel="Depletion",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=0.2,yref=Target,breaks_x =c(2011,2020,2030))
f9
#ggsave(paste0("figs/","Fig8_Depletion_v2.jpg"),plot=f9,dpi=300,width = 30,height = 9,units = "cm")

fig9 <- plot_prob_status(data=dyr,depl=0.2,targt=0.4,Str=Str,nstr=nstr,nsim=nsim,breaks_x=c(2015,2020,2025,2030,2035))
fig9
ggsave(paste0("figs/","Fig9_Probability.jpg"),plot=fig9,dpi=300,width = 20,height = 15,units = "cm")


###################################################
boxplot_ts(data=syr,data_ev=seval,ylabel="Spawning biomass (t)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=200,yref=Srms,breaks_x =c(2011,2020,2030))

Rprom <- mean(md$Reclutas)
p5 <- boxplot_ts(data=ryr,data_ev=reval,ylabel="Recruitment (n)",nstr=nstr,nsim=nsim,Str=Str,ymin=0,ytick=5,yref=Rprom,breaks_x =c(2011,2020,2030))

