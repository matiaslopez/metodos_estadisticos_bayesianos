options (digits=3)


setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Practico01/")

# Paquetes requeridos en este practico
paquetes<-c("ggplot2","reshape2","cowplot", "reshape2", "ggExtra", "brms")
lapply(paquetes,FUN=require,character.only=T)

#############################################################

# Explicacion de max verosimilitud para conteos 
DF=read.csv("Practico_01_counts.csv", header=T, stringsAsFactors=T) # importa los datos
str(DF)
table(DF$insectos)

ggplot(data=DF, aes(x=insectos))+
  geom_histogram(fill="red")+
  theme_bw()+
  labs(x="# insectos", y="Frecuencia")+
  scale_x_continuous(breaks=0:8)+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16))

# Construccion progresiva de la funcion de verosimilitud
dpois(DF$insectos, lambda=1)  # dens Poisson para mu=1
log(dpois(DF$insectos, lambda=1)) # log dens Poisson para mu=1
sum(log(dpois(DF$insectos, lambda=1))) # suma de log dens Poisson para mu=1

# Funcion de max verosimilitud 
LPoisson=function(lambda) {sum(log(dpois(DF$insectos,lambda)))}

# Grafico de la Funcion de max verosimilitud 
DF1=data.frame(media=seq(0.1,10,0.01))
DF1$Ltot=sapply(DF1$media,LPoisson)

ggplot(data=DF1, aes(x=media, y=Ltot))+
  geom_line(col="red", size=1)+
  theme_bw()+
  labs(x=expression(mu[Y]), y="Log Ltot")+
  scale_x_continuous(breaks=seq(from=0.1,to=10, by=1.1))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16))
mean(DF$insectos)

#########################################

# Funciones de Distr previas plausibles:
DF2=data.frame(media=seq(from=0.5,to=10, by=0.01))
DF2$"medlog=1, sd=0.5"=dlnorm(DF2$media,meanlog=1,sdlog=0.5)
DF2$"medlog=2, sd=0.5"=dlnorm(DF2$media,meanlog=2,sdlog=0.5)
DF2$"medlog=1, sd=1.0"=dlnorm(DF2$media,meanlog=1,sdlog=1)
DF2$"medlog=2, sd=1.0"=dlnorm(DF2$media,meanlog=2,sdlog=1)

DF2=melt(DF2,id="media") # transforma el dataframe en "long format"
# Graficos de las Funciones de Distr previas plausibles:
ggplot(data=DF2, aes(x=media, y=value, col=variable))+
  geom_line(lwd=1.3)+
  theme_bw()+
  labs(x=expression(mu[Y]), y="Dens. Probabilidad")+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.position = c(0.6,0.75),
        legend.title = element_blank(),
        legend.text=element_text(size=14)) 

# Distribucion previa seleccionada "a ciegas"
ggplot(data=DF2[DF2$variable=="medlog=2, sd=1.0",], aes(x=media, y=value))+
  geom_line(lwd=1.3)+
  labs(x=expression(mu[Y]), y="Dens. Probabilidad")+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.position =  "none") 

# Un algoritmo MCMC "básico" con 5000 iteraciones
MCMC1= data.frame(matrix(nrow=5000,ncol=1)); names(MCMC1)="mean"
MCMC1$mean[1]=3.5
for(i in 2:5000){
  current=MCMC1$mean[i-1]
  prop=rnorm(n=1,mean=current, sd=0.5)
  post.prop=sapply(prop, LPoisson)+ dlnorm(prop,meanlog=0.5,sdlog=0.5, log=T)
  post.curr=sapply(current, LPoisson)+ dlnorm(current,meanlog=0.5,sdlog=0.5, log=T)
  R=min(1,exp(post.prop)/exp(post.curr))
  MCMC1$mean[i]=ifelse(R >runif(n=1, min=0, max=1), prop, current)
}  
summary(MCMC1$mean)
quantile(MCMC1$mean, probs=c(0.025,0.975))

# Distribucion posterior obtenida con el algoritmo anterior:
ggplot(data=MCMC1, aes(x=mean))+
  geom_density(size=1)+ 
  theme_bw()+ 
  scale_x_continuous(breaks=seq(from=2, to=4, by=0.25))+
  labs(x=(expression(mu[Y])), y="Dens. probabilidad")+ 
  geom_point(aes(x=2.619, y=0), size=4)+ 
  geom_point(aes(x=3.559, y=0), size=4)+ 
  geom_segment(aes(x = 2.619, y = 0, xend = 3.559, yend = 0))+
  theme(axis.title =element_text(size=18),
        axis.text = element_text(size=16))

# Distribucion empleada para obtener el proximo valor candidato de la cadena de Markov
ggplot(data.frame(x = c(1, 6)), aes(x)) +  
  stat_function(fun = dnorm, n = 1000, args = list(mean=3.5, sd=0.5), size=1, col="black") +
  labs(x=expression(mu[Y]), y="Dens. Probabilidad")+
  scale_x_continuous(breaks=seq(from=0, to=10, by=0.5))+
  geom_vline(xintercept=3.5, lty=2, size=1.2, col="dark blue")+
  theme_bw()+
  geom_area(stat = 'function',
            fun = dnorm, args= list(mean=3.5, sd=0.5),
            fill = 'blue',xlim = c(qnorm(0.025, mean=3.5, sd=0.5),qnorm(0.975, mean=3.5, sd=0.5)),
            alpha = 0.4)+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16)) 
