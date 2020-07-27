#MSE_macha: Prepara la Composicion de Tallas Observada y Estimada

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
names(md)
yr <- seq(2011,2017)
#Prepara los datos para graficar
pobs_cut <- md$pobs_cru #proporcion de tallas observadas
yrs_cut <- yr           # años
pobs <- as.data.frame(pobs_cut) 
talla <- md$Tallas      # Tallas
colnames(pobs)<- talla
pobs2 <- melt(pobs)     # melting
pobs2$yrs <- rep(yrs_cut,length(talla))
colnames(pobs2) <- c("talla","frq","yrs")
#Propoecion de tallas ajustadas
ppred_cut <- md$ppred_cru  #proporcion de tallas predichas
ppred <- as.data.frame(ppred_cut)
colnames(ppred)<- talla
ppred2 <- melt(ppred) #melting data
ppred2$yrs <- rep(yrs_cut,length(talla))
colnames(ppred2) <- c("talla","frq","yrs")

# Combina los datos observado y estimados
pobs2$Data <- rep("Observed",length(pobs2$talla))
ppred2$Data <- rep("Predicted",length(ppred2$talla))
comptalla <- rbind(pobs2,ppred2) 

#Composicion de tallas observadas
p3 <- ggplot(data=pobs2,aes(x=talla,y=frq))+
  geom_col(data=pobs2,aes(x=talla,y=frq))+
  facet_wrap(~yrs,ncol=1)+ xlab("Length (mm)")+
  scale_x_discrete(breaks=c(5, 25,50,75,95))+mi.tema()

p3 <- ggplot(data=comptalla)+
  geom_col(aes(x=talla,y=frq,fill=Data),position = "dodge")+mi.tema()+scale_x_discrete(breaks=c(5, 25,50,75,95))+ylab("Proportion")+xlab("Length (mm)")+
  facet_wrap(~yrs,ncol=1)+scale_colour_grey()

ggsave(paste0("figs/","Fig04_Ajuste_Tallas.png"),p4,height=25,width=12,units="cm",dpi=300)
