options(digits=3)
packages.needed=c("ggplot2", "bayesplot","fitdistrplus", "ggbreak", "gridExtra", "qqplotr", "ggeffects", "GGally" , "broom", "qqplotr", "brms",
            "bayesplot","posterior", "doBy", "cowplot")
lapply(packages.needed,FUN=require,character.only=T)
#################################################################################################################################################################

options(mc.cores=parallel::detectCores())
#################################################################################################
#              REGRESION MULTIPLE 
DF1= read.csv("Teo 02 regr multiple.csv", header=T)
str(DF1)
summary(DF1)

# a) Distribucion de prob de la var respuesta
norm.m2=fitdist(DF1$SO2,"norm")
lognorm.m2=fitdist(DF1$SO2,"lnorm")
cdf.m2=cdfcomp(list(norm.m2,lognorm.m2), xlogscale=T,  main="",
               legendtext = c("Normal", "Lognormal"),fitcol = c("black", "grey"), plotstyle ="ggplot")+
  geom_line(size=1.2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position = c(0.70,0.25),
        legend.text=element_text(size=16))
qq.m2=qqcomp(list(norm.m2,lognorm.m2), main="",fitcol = c("black", "grey"),plotstyle ="ggplot")+
  geom_line(size=1.2)+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position ="none")
grid.arrange(cdf.m2, qq.m2, ncol=2) # 

DF1$logSO2=log(DF1$SO2) # transformacion logaritmica de la var respuesta 

#b) Relacion entre var respuesta y  vars explicativas 
ggpairs(DF1[,c("Temp","ManufEnter","Population1970","AvgWindSpeed","AvgPrecip","AvgRainyDays","SO2")], 
        lower=list(continuous = wrap("smooth_loess", alpha = 0.3,size=1)),
        upper = list(continuous = wrap("cor", size=6)))+
  theme_bw()+ 
  theme(strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=14),
        axis.text = element_text(size=12)) # 

# Estandardiza las vars explicativas
DF1s=scale(DF1[,c("Temp","ManufEnter","Population1970","AvgWindSpeed","AvgPrecip","AvgRainyDays")], center=T, scale=T)
summary(DF1s) # verificando
DF1s=as.data.frame(cbind(logSO2=DF1$logSO2, DF1s)) # añadiendo la var de respuesta al DF recien creado

# Distr previas que hay que definir para el modelo a ajustar
get_prior(formula=logSO2~Temp+ManufEnter+Population1970+AvgWindSpeed+AvgPrecip+AvgRainyDays, 
          data=DF1s,family=gaussian) #muestra las distr previas que debe ser especificadas

# Distr previas para el modelo a ajustar
prior.m2 = c(set_prior("normal(0, 2)", class = "Intercept"),
             set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0,3)", class = "sigma"))

# Graficos de las distr previas 
int.m2=ggplot(data.frame(x = c(-4, 4)), aes(x)) +  
  stat_function(fun = dnorm, n = 1000, args = list(mean=0, sd=2)) +
  labs(x="Intercepto", y="Dens. probabilidad")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
slope.m2=ggplot(data.frame(x = c(-2, 2)), aes(x)) +  
  stat_function(fun = dnorm, n = 1000, args = list(mean=0, sd=1)) + 
  labs(x="Pendiente parcial")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        axis.title.y = element_blank())
sd.m2=ggplot(data.frame(x = c(0, 5)), aes(x)) +  
  stat_function(fun = dcauchy, n = 1000, args = list(location=0, scale=3)) +
  labs(x=expression(sigma))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        axis.title.y = element_blank())
plot_grid(int.m2, slope.m2, sd.m2, ncol=3,labels = LETTERS[1:3],
          align="hv",label_x=0.90, label_y=0.95) 

# the statistical model of Bayesian multiple regression
m2.brms=brm(formula=logSO2~Temp+ManufEnter+Population1970+AvgWindSpeed+AvgPrecip+AvgRainyDays, 
            data=DF1s,family=gaussian, prior = prior.m2, warmup = 900, chains=4, iter=2000, thin=3,
            future=TRUE)
summary(m2.brms)

### Evaluacion de la convergencia del modelo
#  distribuciones marginales posteriores para cada cadena y parámetro
par.names.m2=variables(m2.brms)[1:length(variables(m2.brms))-2]
mcmc_dens_overlay(m2.brms,pars=par.names.m2)+ 
  facet_text(on =T)+
  theme_bw()+
  theme(axis.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        strip.background = element_blank()) # fig 4.20

# trace plots: convergencia de las cadenas  
mcmc_trace(m2.brms,size=0.3, pars=par.names.m2)+
  facet_text(on = F)+ 
  theme_bw()+
  theme(axis.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=14),
        strip.background = element_blank()) # not shown 

# autocorrelacion de  parametros 
mcmc_acf(m2.brms, pars=par.names.m2)+
  geom_line(size=1)+
  scale_x_continuous(limits=c(2,10))+
  scale_y_continuous(limits=c(-0.2,0.2))+
  theme_bw()+
  theme(axis.title =  element_text(size=18),
        axis.text = element_text(size=14), 
        strip.text = element_text(size=14,hjust=0.5),
        strip.background = element_blank())

# distributiones posteriores de los parameteros 
dist.post.m2=as_draws_df(m2.brms,pars=pars.m2.brms, add_chain = F)
names(dist.post.m2) 

# Estatisticos decsriptivos de los parametros por cadena 
summaryBy(b_Intercept+b_Temp+b_ManufEnter+b_Population1970+b_AvgWindSpeed+b_AvgPrecip+
            b_AvgRainyDays+sigma~.chain, FUN=c(mean), data=dist.post.m2)

# distribuciones posteriores de los parameteros uniendo las cadenas
dist.post.m2=as_draws_df(m2.brms,pars=pars.m2.brms, add_chain = T)
names(dist.post.m2) 

# Dist posteriores de los parameteros (uniendo las cadenas) 
mcmc_areas(m2.brms,pars=par.names.m2,point_est="mean", prob = 0.9)+
  scale_x_break(c(0.75,2.75))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=15)) 

# R2 de Gelman et al (2018)
str(posterior_predict(m2.brms)) # distr predictiva posterior
R2=data.frame(R2=bayes_R2(m2.brms, summary=F))  #  valores del R2 Bayesiano
ggplot(data=R2, aes(x=R2))+
  geom_density(size=1)+
  labs(x=expression(paste("R"^2)), y="Dens. probabilidad")+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16)) 
bayes_R2(m2.brms) # estadisticos descriptivos del R2

# Analisis de residuales 
# Obteniendo los residuales de Pearson 
res.m2.brms=as.data.frame(residuals(m2.brms, type="pearson", ndraws=1000, summary=T))
head(res.m2.brms,2)
DF1s$resid=res.m2.brms$Estimate # means of Pearson residuals for each data point 
DF1s$resid.Q2.5=res.m2.brms$Q2.5 # credible Intervals of residuals
DF1s$resid.Q97.5=res.m2.brms$Q97.5 # CI residuales
# 1000 valores predichos de la var de respuesta  
fit.m2.brms=as.data.frame(fitted(m2.brms, scale="linear", summary=T, ndraws=1000))
head(fit.m2.brms,2)
DF1s$fit=fit.m2.brms$Estimate # means of fitted values for each data point 
DF1s$fit.Q2.5=fit.m2.brms$Q2.5 # credible Intervals of fitted values
DF1s$fit.Q97.5=fit.m2.brms$Q97.5 # credible Intervals of fitted values

# Residuals vs fitted
res.fit.m2.brms=ggplot(data=DF1s, aes(x=fit,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
res.Temp.m2.brms=ggplot(data=DF1s, aes(x=Temp,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+
  labs(x="Temp",y="Residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(),
        axis.text = element_text(size=16))
res.ManufEnter.m2.brms=ggplot(data=DF1s, aes(x=ManufEnter,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="ManufEnter",y="Residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(),
        axis.text = element_text(size=16))
res.Population1970.m2.brms=ggplot(data=DF1s, aes(x=Population1970,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="Population1970",y="Residuals")+
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16))
res.AvgPrecip.m2.brms=ggplot(data=DF1s, aes(x=AvgPrecip,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="AvgPrecip",y="Residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))
res.AvgWindSpeed.m2.brms=ggplot(data=DF1s, aes(x=AvgWindSpeed,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="AvgWindSpeed",y="Residuals")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_blank(), 
        axis.text = element_text(size=16))
res.AvgRainyDays.m2.brms=ggplot(data=DF1s, aes(x=AvgRainyDays,y=resid))+ 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.1)+
  theme_bw()+ 
  labs(x="AvgRainyDays",y="Residuals")+
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16))
qq.m2.brms=ggplot(data=DF1s, mapping=aes(sample = resid)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(res.fit.m2.brms, res.Temp.m2.brms,res.ManufEnter.m2.brms,
          res.Population1970.m2.brms,res.AvgPrecip.m2.brms,res.AvgWindSpeed.m2.brms,
          res.AvgRainyDays.m2.brms, qq.m2.brms,ncol=3,labels = LETTERS[1:8],
          align="hv",label_x=0.90, label_y=0.95) 

# Graficos de los efectos conditionales 
m2.brms.cond.eff=conditional_effects(m2.brms) # lista de las relaciones marginales predichas 
names(m2.brms.cond.eff) # nombres de los efectos conditionales 
# plots of conditional effects
m2.brms.Temp=ggplot(data=m2.brms.cond.eff$Temp, aes(x=Temp, y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
  labs(y="logSO2",x="Temperature")+ 
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
m2.brms.ManufEnter=ggplot(data=m2.brms.cond.eff$ManufEnter, aes(x=ManufEnter,y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
   labs(y="logSO2",x="# Manufact. Enterp.")+ 
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_null (), 
        axis.text=element_text(size=16))
m2.brms.Population1970=ggplot(data=m2.brms.cond.eff$Population1970,aes(x=Population1970, y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
   labs(y="logSO2",x="Population1970")+ 
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_null (), 
        axis.text=element_text(size=16))
m2.brms.AvgWindSpeed=ggplot(data=m2.brms.cond.eff$AvgWindSpeed,aes(x=AvgWindSpeed, y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
  labs(y="logSO2",x="AvgWindSpeed")+ 
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
m2.brms.AvgPrecip=ggplot(data=m2.brms.cond.eff$AvgPrecip,aes(x=AvgPrecip, y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
  labs(y="logSO2",x="Avg Annual Precipitation")+ 
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_null (), 
        axis.text=element_text(size=16))
m2.brms.AvgRainyDays=ggplot(data=m2.brms.cond.eff$AvgRainyDays,aes(x=AvgRainyDays, y=estimate__))+ 
  geom_line(size=1)+
  theme_bw()+
  geom_ribbon(aes(ymin=lower__,ymax=upper__),fill = "grey90", alpha=0.5)+
  labs(y="logSO2",x="Avg # Rainy Days")+ 
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_null (), 
        axis.text=element_text(size=16))
plot_grid(m2.brms.Temp, m2.brms.ManufEnter, m2.brms.Population1970, 
          m2.brms.AvgWindSpeed, m2.brms.AvgPrecip,m2.brms.AvgRainyDays,  
          ncol=3,labels = LETTERS[1:6],align="hv",label_x=0.90, label_y=0.95) 

#Dist previas con distinto grado de informacion
ggplot(data.frame(x = c(-7, 7)), aes(x)) +  
  stat_function(fun = dunif, n = 1000, args = list(min=-7, max=7), size=1, col="red") +  
  stat_function(fun = dnorm, n = 1000, args = list(mean=0, sd=100), size=1, col="dark blue") +
  stat_function(fun = dnorm, n = 1000, args = list(mean=0, sd=10), size=1, col="dark green") +
  stat_function(fun = dnorm, n = 1000, args = list(mean=0, sd=1), size=1, col="dark magenta") +
  stat_function(fun = dnorm, n = 1000, args = list(mean=1.4, sd=0.2), size=1, col="brown") +
  labs(x="Parámetro", y="Dens. probabilidad (sqrt)")+
  annotate("text",x = -4, y = 0.25, label = "Uniforme(-7,7)", size=6, col="red" )+
  annotate("text",x = 4, y = 0.55, label = "Normal(0,10)", size=6, col="dark green" )+
  annotate("text",x = 5.9, y = 0.25, label = "Normal(0,100)", size=6, col="dark blue" )+
  annotate("text",x = -4, y = 0.5, label = "Normal(0,10)", size=6, col="dark magenta" )+
  annotate("text",x = -3.5, y = 0.85, label = "Normal(0.4,0.2)", size=6, col="brown" )+
  geom_segment(aes(x = -4, y = 0.2, xend = -4, yend = 0.08),col="red", arrow = arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 6, y = 0.2, xend = 6, yend = 0.01),col="dark blue", arrow = arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 4, y = 0.5, xend = 4, yend = 0.04),col="dark green", arrow = arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = -2.5, y = 0.42, xend = -1.2, yend = 0.25),col="dark magenta", arrow = arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = -1.5, y = 0.8, xend = 1, yend = 0.7),col="brown", arrow = arrow(length = unit(0.1,"cm")))+
  scale_x_continuous(breaks = seq(from=-7, to=7, by=1))+
  scale_y_sqrt()+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
