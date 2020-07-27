#MSE_macha: Prepara la figura de indicadores poblacionales
code_dir <- "./R/"
file_path <- list.files(code_dir,full.names = T)
for(f_path in file_path){source(f_path,encoding = "UTF-8")}

library(ggplot2)
library(ggpubr)
library(reshape2)
library(lattice)
theme_set(theme_bw())

# Lee el reporte del modelo operativo
md <- read.admb("models/opmacha")

#Prepara los datos a graficar Panel A
Total <- md$Bt    # biomasa total
Spawning <- md$SD # biomasa desovante
Vulnerable <- md$Bv #biomasa vulnerable
Catch <- md$Cobs  #Capturas
df <- data.frame(Total,Spawning,Vulnerable,Catch)
df <- melt(df)
df$yr <- rep(yr,4) #agreaga años por 4 indicadores a graficar

# Reclutamiento a graficar Panel B
Recruits <- md$Reclutas
Rec_devs <- md$log_desv_Rt
Fishing <- md$fyr
df3 <- data.frame(Recruits,Rec_devs,Fishing,yr)

# Figura Panel A: Biomasas
p3a <- ggplot(df)+
  geom_line(aes(x=yr,y=value,linetype=variable))+
  ylab("Biomass (t)")+xlab("Year")+mi.tema()+theme(legend.position = c(0.8,0.9))

# Figura Panel B: Reclutamiento
p3b <- ggplot(data=df3)+
  geom_col(aes(x=yr,y=Recruits))+mi.tema()+
  geom_hline(yintercept = md$R0)+
  ylab("Recruits (millions)")+xlab("Year")

#Figura Panel C: Mortalidad por pesca
p3c <- ggplot(data=df3)+
  geom_line(aes(x=yr,y=Fishing))+
  geom_point(aes(x=yr,y=Fishing))+
  mi.tema()+
  ylab("Fishing mortality")+xlab("Year")

Fig5 <- ggarrange(p3a,p3b,p3c,labels = c("A","B","C"),nrow=1,ncol = 3)

ggsave(paste0("figs/","Fig5_PopIndicadores.jpg"),plot=Fig5,dpi=300,width = 30,height = 18,units = "cm")