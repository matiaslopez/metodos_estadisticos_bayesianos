setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Teorico05/")

options(digits=3)
mc.cores=parallel::detectCores()
packages.needed=c("ggplot2","reshape2", "runjags","future","coda", "cowplot","arm", "bayesplot", "parallel")
lapply(packages.needed,FUN=require,character.only=T)

#########################################
# Occupancy models
DF=read.csv(file="Carbonero Montano Teo05.csv", header=T)
str(DF)

# Estandardiza las vars explicativas numericas
DF$elev=as.vector(scale(DF$elev, center=T, scale=T))
DF$forest=as.vector(scale(DF$forest, center=T, scale=T))

# Otros calculos necesarios
DF[,c("y.1","y.2", "y.3")] # las 3 visitas a los 237 sitios
visitas=rowSums(!is.na(DF[,c("y.1","y.2", "y.3")])) # Nb de visitas por sitio
detec=rowSums(DF[,c("y.1","y.2","y.3")],na.rm=T) # Nb de detecciones por sitio


# Creamos una "variable latente binaria" occ (o efecto aleatorio binario): 
# cuando ocup=0, hay incertidumbre sobre su valor real y por ende NA, 
# cuando detect>0, hay ocupacion
ocup =  ifelse(detec > 0, 1, NA)

# los datos en JAGS deben ser entrados como una lista
jagsData <- list(y = detec, visitas = visitas, nSites = length(visitas),
                 ocup =  ifelse(detec > 0, 1, NA),
                 forest = DF$forest, ele = DF$elev, ele2 = DF$elev*DF$elev)

# El modelo a ajustar: 
m1="model{
  for(i in 1:nSites) {
    logit(pres[i]) <- b0 + bFor * forest[i] + 
                     bElev * ele[i] + bElev2 * ele2[i]
    ocup[i] ~ dbern(pres[i])
    Y[i] ~ dbin(p * ocup[i],visitas[i])  
  }
  
  # Previas
  b0 ~ dnorm(0, 0.5) # intercepto
  bFor ~ dnorm(0, 0.5)    # pendiente de forest
  bElev ~ dnorm(0, 0.5)   # pendiente de elevation
  bElev2 ~ dnorm(0, 0.5)  # pendiente de elevation^2
  p ~ dbeta(3, 3) # prob de observacion
  
  # variable calculada: 
  N <- sum(ocup) #abundancia global a partir de ocupaciones locales 
}"

# Visualización de la distribucion beta usada como previa de p.
a=data.frame(x=seq(from=0, to=1, length.out = 100))
a$Beta1=dbeta(x=a$x,shape1=3, shape2=1) # mean=3/5 phi=4
a$Beta2=dbeta(x=a$x,shape1=1, shape2=3) # mean=1/4 phi=4
a$Beta3=dbeta(x=a$x,shape1=0.5, shape2=0.5) # mean=0.5 phi=3
a$Beta4=dbeta(x=a$x,shape1=1, shape2=1) # mean=0.5 phi=2
a$Beta5=dbeta(x=a$x,shape1=3, shape2=3) # mean=0.5 phi=6

# reformateando el dataframe a en "long format"
a=melt(a, id="x", value.name = "Probability")
a$line=rep(c("solid", "dashed","dotted","longdash", "dotdash"), each=100)
ggplot(data=a, aes(x=x, y=Probability, linetype=variable))+ 
  geom_line(aes(linetype=line), linewidth=1.1)+
  labs(y="Probability (Y=y)", x= "Y")+
  theme_bw()+
  theme (axis.title=element_text(size=18), 
         axis.text=element_text(size=16),
         legend.position = "null")+
  annotate("text", x=.125, y=3, label= expression(paste(mu)["Y"]*"= 0.33, "*paste(Phi)*" =3"), size=6)+
  geom_segment(aes(x = 0.123, y = 2.9, xend = 0.1, yend = 2.5),size=0.5,
               arrow = arrow(length = unit(0.2,"cm")))+
  annotate("text", x=.825, y=3, label= expression(paste(mu)["Y"]*"= 0.67,"*paste(Phi)*" =3"), size=6)+
  geom_segment(aes(x = 0.85, y = 2.9, xend = 0.875, yend = 2.5),size=0.5,
               arrow = arrow(length = unit(0.2,"cm")))+
  annotate("text", x=.25, y=2.5, label= expression(paste(mu)["Y"]*"= 0.5,"*paste(Phi)*" =1"), size=6)+
  geom_segment(aes(x = 0.25, y = 2.4, xend = 0.08, yend = 1.2),size=0.5,
               arrow = arrow(length = unit(0.2,"cm")))+
  annotate("text", x=.4, y=3, label= expression(paste(mu)["Y"]*"= 0.5,"*paste(Phi)*" =2"), size=6)+
  geom_segment(aes(x = 0.4, y = 2.9, xend = 0.5, yend = 1.1),size=0.5,
               arrow = arrow(length = unit(0.2,"cm")))+
  annotate("text", x=.6, y=2.5, label= expression(paste(mu)["Y"]*"= 0.5,"*paste(Phi)*" =6"), size=6)+
  geom_segment(aes(x = 0.55, y = 2.4, xend = 0.55, yend = 2),size=0.5,
               arrow = arrow(length = unit(0.2,"cm"))) # Fig 12.1

# Ajustando el modelo
out.m1=run.jags(data=jagsData,
                monitor=c("p", "b0", "bFor", "bElev", "bElev2","N", "deviance"), 
                model=m1, n.chains=3, thin=5, sample=5000, burnin = 100, method="rjparallel") 
# Resultados del modelo
summary(out.m1)

# Convierte el output de JAGS en objeto tipo MCMC
m1.mcmc <- as.mcmc.list(out.m1)

# Convergencia de las cadenas
mcmc_dens_overlay(m1.mcmc)+
  geom_density(linewidth=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 16)) 

# graficos de trazas de las cadenas
mcmc_trace(m1.mcmc)+
  theme_bw()+
  ylab("valor del parámetro")+
  scale_x_continuous(breaks=seq(from=0, to=5000, by=1000))+
  theme(legend.position="none",
        axis.text=element_text(size = 12),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12))

# Autocorrelacion de los parámetros en las cadenas
mcmc_acf(m1.mcmc)+
  geom_line(size=0.8)+
  theme_classic()+
  scale_x_continuous(limits=c(2,10))+
  scale_y_continuous(limits=c(-0.1,0.1))+
  theme(axis.text = element_text(size=14),
        axis.title =  element_text(size=18),
        strip.text = element_text(size=14),
        strip.background=element_rect(colour = "black", fill = "white"))

# distribuciones posteriores de los params que afectan la prob ocupacion
mcmc_areas(m1.mcmc, prob_outer = 0.95, regex_pars = c("b"))+
  geom_density(linewidth=1, alpha=0.9)+ 
  scale_x_continuous(breaks=seq(from=-2, to=5, by=0.5))+
  theme_bw()+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))  

# distribuciones posteriores de la prob de observación
mcmc_areas(m1.mcmc, prob_outer = 0.95, regex_pars = c("^p"))+
  geom_density(linewidth=1.2, alpha=0.9)+ 
  scale_x_continuous(breaks=seq(from=0, to=1, by=0.05))+
  theme_bw()+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))  

# Obtiene la Distr predictiva posterior del modelo.
DPP.m1=run.jags(data=jagsData,monitor=c("Y"),model=m1, n.chains=3, thin=5, 
                sample=5000, burnin = 100, method="rjparallel") 

# transformando la Distr Pred Posterior en un dataframe:
DPP=as.data.frame(unlist(DPP.m1$mcmc [[1]]))

# Distr predictiva posterior y la 1era visita 
ppc_dens_overlay(y=DF$y.1,yrep=as.matrix(DPP), trim = F, size = 0.5,alpha = 1)+
  xlab("1era visita")+ 
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")

# Suma de los valores en las 3 visitas, hay NAs!
DF$suma=DF$y.1+DF$y.2+DF$y.3
# filas para las cuales DF$suma contiene un NA
which(is.na(DF$suma))

# Distr predictiva posterior SIN las filas para las que DF$suma tiene NA
DPP1=DPP[, which(!is.na(DF$suma))]

ppc_dens_overlay(y=DF$suma[which(!is.na(DF$suma))],
                 yrep=as.matrix(DPP[, which(!is.na(DF$suma))]), 
                 trim = F, size = 0.5,alpha = 1)+
  xlab("Detecciones en 3 visitas")+ 
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")

# Media de detecciones en 3 visitas
ppc_stat(y=DF$suma[which(!is.na(DF$suma))],
          yrep=as.matrix(DPP[, which(!is.na(DF$suma))]),stat = "mean", binwidth = 0.05)+  
  xlab("Media de detecciones en 3 visitas")+
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")

# Uso de los pp-checks a partir de la dist predictiva posterior
cero=function (x){sum(x==0)/length(x)} # proporcion de ceros en un vector
unos=function (x){sum(x==1)/length(x)} # proporcion de unos en un vector

# proporcion de ceros en 3 visitas
ppc1=ppc_stat(y=DF$suma[which(!is.na(DF$suma))],
         yrep=as.matrix(DPP[, which(!is.na(DF$suma))]),stat = cero)+
  xlab("Prob (Y=0) en las 3 visitas")+ 
  theme_bw()+
  scale_x_continuous(breaks=seq(from=0, to=1, by=0.1))+
  theme(axis.text = element_text(size=14,hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=16),
        legend.position = "none")

# proporcion de unos en 3 visitas
ppc2=ppc_stat(y=DF$suma[which(!is.na(DF$suma))],
          yrep=as.matrix(DPP[, which(!is.na(DF$suma))]),stat = unos)+
   xlab("Prob (Y=1) en 3 visitas")+ 
  theme_bw()+
  scale_x_continuous(breaks=seq(from=0, to=1, by=0.1))+
  theme(axis.text = element_text(size=14,hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=16),
        legend.position = "none")
plot_grid(ppc1, ppc2, ncol=2,labels = LETTERS[1:2],align="hv",label_x=0.85, label_y=0.95) 

########################################################################################
#             Binomial -mixture model
# 

DF1=read.csv(file="ReinitaHornera Teo05.csv", header=T)
DF1$sitio=as.factor(DF1$sitio)
summary(DF1)

# Tabulaciones descriptivas
ftable(table(DF1$y.1, DF1$y.2, DF1$y.3))

# Modelo SIN variables explicativas
Bimix="model{
          for(i in 1:nSites) {
               # modelo de abundancia real (no observable)
               Conteo[i] ~ dpois(lambda)
               # modelo de observación
               for(j in 1:nOcc) {
                  Y[i,j] ~ dbin(p, Conteo[i])
               }
            }
  # Previas
  lambda ~ dunif(0, 10) # media de dist Poisson
  p ~ dbeta(1,1) # probabilidad de deteccion
  # variable derivada
  Ntotal <- sum(Conteo) # abundancia global
}"

# los datos
JAGSdata1<-list(nSites=70, nOcc=4, 
                Y=as.matrix(DF1[,c("y.1","y.2","y.3","y.4")]))
str(JAGSdata1)
# el primer modelo
m2.out=run.jags(data=JAGSdata1,model=Bimix, monitor=c("p","lambda","Ntotal"),
                n.chains=3, thin=5, sample=5000, burnin = 100, method="rjparallel")
summary(m2.out)

# Convierte el output de JAGS en objeto tipo MCMC
m2.mcmc <- as.mcmc.list(m2.out)

# Convergencia de las cadenas
mcmc_dens_overlay(m2.mcmc)+
  geom_density(linewidth=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 16)) 

# No verifiqué los graficos de trazas y de autocorrelacion de los estimados.
# ni hice una mínima interpretación de estos resultados

# Modelo CON variables explicativas
BiMix2="model{
  for(i in 1:nSites) {
    # modelo de abundancia real (no observable)
    Conteo[i] ~ dpois(lambda)
    #probabilidad de deteccion
    p[i]<- 1/(1+exp(-1*(b0 + b.ufc*ufc[i] + b.trba*trba[i])))
    for(j in 1:nOcc) {
      Y[i,j] ~ dbin(p[i], Conteo[i])
    }
  }
  # Distribuciones previas
  lambda ~ dunif(0, 10) # media de dist Poisson
  b0 ~dnorm(0,0.5) # previa intercepto
  b.ufc~dnorm(0,0.5)  # previa pendiente parcial
  b.trba~dnorm(0,0.5) # previa pendiente parcial
  # variable derivada
  Ntotal <- sum(Conteo)# abundancia global
}"

# los datos para el 2do Binomial mixture model
JAGSdata2 <-list(nSites=nrow(DF1), nOcc=4,
             Y=as.matrix(DF1[,c("y.1","y.2", "y.3", "y.4")]),
             ufc = DF1$ufc, trba = DF1$trba)

# El modelo ajustado
m3.out=run.jags(data=JAGSdata2,model=BiMix2, 
                monitor=c("b0","b.trba","lambda","b.ufc","Ntotal","p"),
                n.chains=3, thin=5, sample=5000, burnin = 100, method="rjparallel")
summary(m3.out)

# Convierte el output de JAGS en objeto tipo MCMC
m3.mcmc <- as.mcmc.list(m3.out)

# Convergencia de las cadenas para los parámetros
mcmc_dens_overlay(m3.mcmc, regex_pars = c("^b", "Ntotal", "lambda"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 16)) 

# Convergencia de las cadenas para algunas prob de deteccion [i]
mcmc_dens_overlay(m3.mcmc, pars = c("p[1]","p[2]","p[3]","p[27]","p[42]","p[66]"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 16)) 

# graficos de trazas: OK
mcmc_trace(m3.mcmc, regex_pars = c("^b", "Ntotal"))+
  theme_bw()+
  ylab("valor del parámetro")+
  scale_x_continuous(breaks=seq(from=0, to=5000, by=1000))+
  theme(legend.position="none",
        axis.text=element_text(size = 12),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12))

# Autocorrelacion of sampled parameter values per chain
mcmc_acf(m3.mcmc, regex_pars = c("^b", "Ntotal"))+
  geom_line(size=0.8)+
  theme_bw()+
  scale_x_continuous(limits=c(1,20))+
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18),
        strip.text = element_text(size=14))

# Aumentando el thinning para disminuir la autocorrelación observada en el grafico anterior
m3.out1=run.jags(data=JAGSdata2,model=BiMix2, 
                monitor=c("b0","b.trba","lambda","b.ufc","Ntotal","p"),
                n.chains=3, thin=30, sample=5000, burnin = 100, method="rjparallel")
summary(m3.out1)
m31.mcmc <- as.mcmc.list(m3.out1)

# Autocorrelacion of sampled parameter values per chain
mcmc_acf(m31.mcmc, regex_pars = c("^b", "Ntotal"))+
  geom_line(size=0.8)+
  theme_bw()+
  scale_x_continuous(limits=c(1,20))+
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18),
        strip.text = element_text(size=14))

#distribuciones posteriores
mcmc_areas(m31.mcmc, regex_pars = c("^b"), prob=0.79)+
  geom_density(lwd=1.2, alpha=0.9)+ 
  scale_x_continuous(breaks=seq(from=-3.5, to=1, by=0.5))+
  theme_bw()+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))  

# distribuciones posteriores
mcmc_areas(m31.mcmc, prob=0.71, regex_pars = c("Ntotal"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))  

# Obtiene la distr predictva posterior de cada los conteos en cada sitio y ocasion
m4.out=run.jags(data=JAGSdata2,model=BiMix2, monitor=c("Y"),
                n.chains=3, thin=5, sample=5000, burnin = 100, method="rjparallel")

# Convierte el output del modelo en un objeto mcmc.list
m4.mcmc <- as.mcmc.list(m4.out)
# Se unen las 3 cadenas -- ya que el modelo convergió -- y se convierte en data frame
m4.dist.pred.post=as.data.frame(combine.mcmc(m4.mcmc))
str(m4.dist.pred.post)
names(m4.dist.pred.post)
