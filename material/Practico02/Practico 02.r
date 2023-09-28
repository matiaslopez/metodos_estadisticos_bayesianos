options(digits=3)
#To obtain the figures you must have first imported the data and ran the statistical models. 
packages.needed=c("ggplot2", "bayesplot","fitdistrplus","qqplotr", "ggeffects", "GGally" , "broom", 
                  "qqplotr", "brms",  "bayesplot", "doBy","multcomp", "cowplot")
lapply(packages.needed,FUN=require,character.only=T)
#####################################################################################################################################################  
setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Practico02/")

DF2=read.csv(file="Pr 02 sleep mammals.csv", header=T,stringsAsFactors=T)
str(DF2)

#  Analisis exploratorio de datos
summary(DF2[,c("vore", "sleep_total")])
DF2=na.omit(DF2) # elimina los datos faltantes

# Estatisticos descriptivos
desc.stats=function(x){c(mean=mean(x), median=median(x),sd=sd(x), n=length(x))} 
summaryBy(sleep_total~vore, data=DF2, FUN=desc.stats)

# Relacion entre las vars de respuesta y explicativas 
ggplot(data=DF2, aes(x=vore, y=sleep_total))+
  geom_boxplot(size=0.5, alpha=0.9)+
  theme_bw()+
  labs(x="Dieta", y="Horas de sueño")+
  stat_summary(fun = "mean", colour = "dark red", size = 4, pch=15, geom = "point")+
  geom_jitter(size=2, width=0.1)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16)) 

#  Distribucion de Prob para la verosimilitud
norm=fitdist(DF2$sleep_total, "norm")
lognorm=fitdist(DF2$sleep_total,"lnorm")
cdf.sleep=cdfcomp(list(norm,lognorm), main="", legendtext =c("Normal", "Lognormal"), fitcol = c("black", "grey"), 
                  plotstyle ="ggplot")+
  geom_line(size=1.2)+
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16),
        legend.position = c(0.7,0.25),
        legend.text= element_text (size=14))
qq.sleep=qqcomp(list(norm,lognorm), main="",fitcol = c("black","grey"), 
                plotstyle ="ggplot")+
  geom_line(size=1.2)+theme_bw()+
  theme(axis.title= element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position ="none")
plot_grid(cdf.sleep, qq.sleep, ncol=2) 

unique(model.matrix(sleep_total~vore, data=DF2)) # Matriz de diseño

get_prior(formula=sleep_total~vore, data=DF2,family=gaussian) #las distr previas a definir 
# las distr previas propuestas
prior.m3 = c(set_prior("lognormal(log(9), 0.5)", class = "Intercept"),
             set_prior("normal(0, 1.8)", class = "b"),
             set_prior("lognormal(log(6), 0.5)", class = "sigma"))
# Graficos de las distr previas propuestas 
mean.ref.m3=ggplot(data.frame(x = c(0, 24)), aes (x))+ 
  stat_function(fun=dlnorm,n=1000,args=list(meanlog = log(9), sdlog =0.5)) +
  labs(x="Mean (ref)", y="Prob. density")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
diff.means.m3=ggplot(data.frame(x = c(-8, 8)), aes (x))+   
  stat_function(fun = dnorm,n=1000, args = list(mean=0, sd=1.8)) + 
  labs(x="Diff. betw. means")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),              
        axis.title.y = element_blank())
sd.Y.m3=ggplot(data.frame(x = c(0, 20)), aes (x))+  
  stat_function(fun=dlnorm, n=1000, args=list(meanlog=log(6),sdlog =0.5)) +
  labs(x=expression(sigma))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),              
        axis.title.y = element_blank())
plot_grid(mean.ref.m3, diff.means.m3, sd.Y.m3, ncol=3,
          labels = LETTERS[1:3],align="hv",label_x=0.9, label_y=0.95)

# El modelo a ajustar
options(mc.cores=parallel::detectCores())
m3.brms=brm(formula=sleep_total~vore, data=DF2, family=gaussian, prior = prior.m3, 
            warmup = 1000, future=TRUE, chains=3, iter=2000, thin=3)
summary(m3.brms)

### Evaluacion de la convergencia del modelo 
# Dist posterior de cada parametro para cada cadena de m3.brms 
mcmc_dens_overlay(m3.brms, regex_pars = c("^b", "sigma"))+
  geom_density(lwd=1.2, alpha=0.1)+
  theme_bw()+
  ylab("Dens. Probabilidad")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position="none")
 
# Trace plots: degree of mixing of the three chains
mcmc_trace(m3.brms, regex_pars = c("^b", "sigma"), size=0.3)+
  theme_bw()+
  ylab("Valor del Parametro")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),              
        strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position="none") 

# Autocorrelacion de los parametros estimados por cadena: 
mcmc_acf(m3.brms, regex_pars = c("^b", "sigma"))+
  geom_line(size=1)+
  scale_x_continuous(limits=c(2,20))+
  scale_y_continuous(limits=c(-0.25,0.25))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 16)) 

# Una mejor visualizacion del summary del modelo
mcmc_plot(m3.brms, regex_pars = c("^b", "sigma"), type="intervals", prob_outer = 0.95)+ 
  scale_y_discrete(labels=c("Mean carn.","diff herb.","diff insect.",
                            "diff omni.", "sigma"))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),              
        axis.text.y = element_text (angle=45,hjust=0.4)) 
mcmc_areas_ridges(m3.brms, regex_pars = c("^b", "sigma"), prob_outer = 0.95)+ 
  scale_y_discrete(labels=c("Mean carn.","diff herb.","diff insect.",
                            "diff omni.", "sigma"))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        axis.text.y = element_text (angle=45,hjust=0.4)) 
  theme_bw()+
  scale_x_continuous(breaks=seq(from=-10, to=16, by=2))+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        axis.text.y = element_text (angle=45,hjust=0.4)) 

# Extrae muestras de las distr posteriores marginales de cada parametro
post.m3brms=as_draws_df(m3.brms)
names(post.m3brms)

#Grafico de las distr previas y de las distr posteriores de cada parametro
post.pr.mean.ref.m3=mean.ref.m3+  
  geom_density(data=post.m3brms,aes(x=b_Intercept), linetype=3)
post.pr.diff.m3=diff.means.m3+
  geom_density(data=post.m3brms,aes(x=b_voreherbi), linetype=2)+
  geom_density(data=post.m3brms,aes(x=b_voreinsecti), linetype=3)+
  geom_density(data=post.m3brms,aes(x=b_voreomni), linetype=4)
post.pr.sigma.m3=sd.Y.m3+
  geom_density(data=post.m3brms,aes(x=sigma), linetype=2)
plot_grid(post.pr.mean.ref.m3, post.pr.diff.m3, post.pr.sigma.m3, ncol=3,
          labels = LETTERS[1:3],align="hv",label_x=0.9, label_y=0.95) 

# medias de cada grupo a partir de las dist posteriores
mediasm3=data.frame(carni=post.m3brms$b_Intercept,
                    herbi=post.m3brms$b_Intercept+post.m3brms$b_voreherbi,
                    insecti=post.m3brms$b_Intercept+post.m3brms$b_voreinsecti,
                    omni=post.m3brms$b_Intercept+post.m3brms$b_voreomni)
ggplot(data=mediasm3, aes(x=carni))+
  geom_density()+
  geom_density(aes(x=herbi), col="red")+
  geom_density(aes(x=insecti), col="blue")+
  geom_density(aes(x=omni), col="dark green")+
  theme_bw()+
  labs(x="horas de sueño", y="Dens. probabilidad")+
  scale_x_continuous(breaks=seq(from=-10, to=16, by=2))+
  annotate("text",x = 7, y = 0.3, label = "Carnívoros", size=5, col="black" )+
  annotate("text",x = 7, y = 0.28, label = "Herbívoros", size=5, col="red" )+
  annotate("text",x = 7, y = 0.26, label = "Insectívoros", size=5, col="blue" )+
  annotate("text",x = 7, y = 0.24, label = "Omnívoros", size=5, col="dark green" )+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))

bayes_R2(m3.brms) # R2 de Gelman et al (2018)

# Residual analysis of Bayesian One way ANOVA: 
# Obteniendo los residuales de Pearson 
res.m3.brms=data.frame(residuals(m3.brms, type="pearson", ndraws=1000, summary=T))
head(res.m3.brms,2)
# Obteniendo los valores predichos (a partir de 1000 muestras)de la var respuesta  
fit.m3.brms=data.frame(fitted(m3.brms, scale="linear", summary=T, ndraws=1000))
head(fit.m3.brms,2)
DF2$fit=fit.m3.brms$Estimate # medias de valores spredichos para cada dato 
DF2$resid=res.m3.brms$Estimate # medias de los residuales de Pearson para cada dato

# Residuals vs fitted 
res.fit.m3.brms=ggplot(data=DF2, aes(x=fit,y=resid))+ 
  geom_jitter(col="black", width=0.15, size=3)+
  geom_hline(yintercept =0,linetype = 2, size=1.5)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuals vs vore
res.vore.m3.brms=ggplot(data=DF2, aes(x=vore,y=resid))+ 
  geom_jitter(col="black", width=0.15, size=3)+
  geom_boxplot(alpha=0.1)+
  geom_hline(yintercept =0,linetype = 2, size=1.5)+
  theme_bw()+
  labs(x="Diet",y="Residuals")+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        axis.title.y=element_blank())
# QQplot
qq.m3.brms=ggplot(data=DF2, mapping=aes(sample = resid)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(res.fit.m3.brms, res.vore.m3.brms, qq.m3.brms, ncol=3,
          labels = LETTERS[1:3],align="hv",label_x=0.9, label_y=0.95) 

# Plot de efectos condicionales predichos por m3.brms 
plot(conditional_effects(m3.brms, effects="vore",prob=0.89))[[1]]+ 
  labs(y="Sleep total", x="Diet")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16)) 

###################

# Test de hipotesis usando los factores de Bayes.
# Media car - Media herb>1.4
hypothesis(m3.brms, "voreherbi>1.4", class="b")
plot(hypothesis(m3.brms, "voreherbi>1.4", class="b"))[[1]]+
  theme_bw()+
  labs(x="Carn-herb>1.4", y="Dens Probabilidad")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.position = "none",
        strip.background = element_blank(), 
        strip.text = element_text(size = 16)) 


# Tambien se pueden hacer varios tests simultaneamente
h = c("carn>herb"="Intercept> Intercept+voreherbi",
      "carn>insect"="Intercept> Intercept+voreinsecti",
      "carn>prom.herb.insect"="Intercept> (Intercept+voreherbi+Intercept+voreinsecti)/2")
hypothesis(m3.brms, h, class="b")

plot(hypothesis(m3.brms, h, class="b"))[[1]] +
  theme_bw()+
  labs(y="Dens Prob")+
  theme(axis.title=element_text(size=18), 
       axis.text=element_text(size=16),
       legend.position = "none",
       strip.background = element_blank(), 
       strip.text = element_text(size = 16)) 
