mc.cores=parallel::detectCores()
options(digits=3)
packages.needed=c("ggplot2","fitdistrplus","brms","future","GGally",  "ggeffects","performance", "cowplot",
                  "ggdist", "arm", "DHARMa","qqplotr","reshape2", "bayesplot", "parallel", "doBy","ggbreak", "ggridges", "ggdist","pROC")
lapply(packages.needed,FUN=require,character.only=T)

##########################################################################################################################
DF=read.csv(file="actdiest Pr05.csv", header=T, stringsAsFactors = T) # importa archivo
str(DF)
# convierte variables a factores
DF$Fieldnumber=as.factor(DF$Fieldnumber)
DF$Plot=as.factor(DF$Plot)
DF$rootdiamscore=as.factor(DF$rootdiamscore)
str(DF)
names(DF)
####################################################################################################################################

####################################################################################################################################
#   0. Caracterización básica de las variables del DF. 
####################################################################################################################################
summary(DF)
# Cuenta el número de ceros en la var de respuesta
table(DF$actdiest==0) # 11.5% de ceros
#531    69

table(DF$spsab)
# cuenta el numero de especies para las que hay 1,2,3,4,...26 datos en cada dataframe
table(table(DF$spsab)) 

table(DF$spsab)[table(DF$spsab)>=9] # especies para las que hay al menos 10 datos

# nombres de las 25 especies para las que hay al menos 10 datos: va a ser usado en gráficos luego
sp.10.diest=names(table(DF$spsab)[table(DF$spsab)>=9])

#############################################################################################
#   1. Distr de probdablidades de la variable de respuesta
#############################################################################################
# Grafico de densidades de la var de respuesta
ggplot(data=DF, aes(x=actdiest))+
  geom_density(size=1.4, col="red")+
  theme_bw()+
  labs(x="Act. diesterasa", y="Densidad")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# La var de respuesta actdiest son reales >= ceros. Se podría considerar una
# distribución de mezcla como alternativa a Tweedie (tiene el mismo número de parámetros).
# Una posibilidad seria Zero augmented Gama, que combina una distr binomial (para los ceros) 
# y una distr Gama para los valores reales estrictamente positivos.
# 
# Para diesterasa: usando los valores de actdiest estrictamente positivos
dens.diest=ggplot(data=DF[DF$actdiest>0,], aes(x=actdiest))+
  geom_density(linewidth=1)+
  theme_bw()+
  labs(x="Act.diestesterasa", y="Densidad")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16))
diest.gama=fitdist(data=DF[DF$actdiest>0,]$actdiest, distr = "gamma")
diest.lognorm=fitdist(DF[DF$actdiest>0,]$actdiest, distr = "lnorm")
cdf.diest=cdfcomp(list(diest.lognorm,diest.gama), xlogscale=T, ylogscale=F, main="",
                  legendtext =c("Lognormal", "Gama"),fitcol = c("blue", "red"), 
                  plotstyle ="ggplot")+
  geom_line(linewidth=1, col="black")+
  geom_point(size=1)+
  xlab("Act. diesterasa")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position = c(0.3,0.75))
qq.diest=qqcomp(list(diest.lognorm,diest.gama), main="", fitcol = c("blue", "red"), plotstyle ="ggplot")+
  geom_line(linewidth=1)+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position = "none")
plot_grid(dens.diest, cdf.diest,qq.diest, nrow=1)
# CONCLUSION: la dist Gama ajusta mejor para los valores positivos de diestest. 
# Ninguna distr ajusta muy bien a valores extremos de cola superior. Se usara ZA Gama 

#############################################################################################
#   3. Analisis exploratorio de datos
#############################################################################################
# Covariacion de las vars explicativas
ggpairs(DF[,c("pH", "Ntotal","ResinP", "Clay", "Silt", "C")],
        lower = list(continuous = "smooth"))+
  theme(strip.background =element_rect(fill="white"),
        strip.text=element_text(size=12))

# Centra y estandardiza las variables explicativas a incluir para comparar sus efectos relativos 
DF$pH.s=as.vector(scale(DF$pH, center=T, scale=T))
DF$Clay.s=as.vector(scale(DF$Clay, center=T, scale=T))
DF$Resin.P.s=as.vector(scale(DF$ResinP, center=T, scale=T))
DF$C.s=as.vector(scale(DF$C, center=T, scale=T))
DF$Silt.s=as.vector(scale(DF$Silt, center=T, scale=T))

summary_by(actdiest~rootdiamscore, data=DF, FUN=c(mean,sd,length))
# hay diferencias de medias de actdiest segun rootdiamscore

# Relación entre cada var respuesta y las 5 variables explicativas a incluir
diest1=ggplot(data=DF, aes(x=ResinP, y=actdiest, color=rootdiamscore))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  ylim(0,60)+
  geom_smooth(se=F)+ 
  labs(x="Resin P",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position=c(0.75,0.75)) 
diest2=ggplot(data=DF, aes(x=pH, y=actdiest, color=rootdiamscore))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  geom_smooth(se=F)+ 
  ylim(0,60)+
  labs(x="pH",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none") 
diest3=ggplot(data=DF, aes(x=Clay, y=actdiest, color=rootdiamscore))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  ylim(0,60)+
  geom_smooth(se=F)+ 
  labs(x="Clay",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none") 
diest4=ggplot(data=DF, aes(x=Silt, y=actdiest, color=rootdiamscore))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  ylim(0,60)+
  geom_smooth(se=F)+ 
  labs(x="Silt",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")  
diest5=ggplot(data=DF, aes(x=C, y=actdiest, color=rootdiamscore))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  ylim(0,60)+
  geom_smooth(se=F)+ 
  labs(x="Carbon",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")  
diest6=ggplot(data=DF, aes(x=rootdiamscore, y=actdiest))+ 
  geom_jitter(size=1.2, width=0.1)+
  theme_bw()+
  ylim(0,60)+
  geom_boxplot(alpha=0.1)+ 
  labs(x="rootdiamscore",y="actdiest")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16)) 
plot_grid(diest1, diest2, diest3, diest4, diest5, diest6,nrow=2)
          
# Visualizando la relacion entre actdiest y ResinP para las sp con más de 10 datos 
ggplot(data=DF[DF$spsab %in% sp.10.diest,], aes(x=ResinP, y=actdiest))+ 
  geom_point(shape=1)+
  theme_bw()+
  geom_smooth(method="lm", se=F)+ 
  labs(x="Resin P",y="actdiest")+
  facet_wrap (~ spsab, scale="free_y") + 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=14),
        strip.background = element_blank(),
        strip.text=element_text(size=16)) 
# hay una debil heterogeneidad de pendientes, lo que justifica usar pendientes aleatorias 
# en el modelo de actdiest

#############################################################################################
#   4. Modelo estadístico a ajustar
#############################################################################################

formula.m1 = bf(actdiest~rootdiamscore+ Resin.P.s+Clay.s+pH.s+Silt.s+C.s +
                        (1|Fieldnumber)+(1+Resin.P.s|spsab)+(1|fam),
                hu~Resin.P.s+Clay.s+pH.s+Silt.s+C.s+(1|fam))
get_prior(formula=formula.m1,data=DF,family='hurdle_gamma')


#############################################################################################
#   5.Ajustando el modelo
#############################################################################################
m1=brm(formula=formula.m1,data=DF,family='hurdle_gamma',warmup = 1000, chains=3, 
       iter=5000, thin=2, control = list(adapt_delta = 0.95))
summary(m1)

# Pendientes de las vars explicativas en la parte log(muY):
fixef(m1)[c(5:9)]
# denotan el efecto de aumentar  una SD de cada var explicativa sobre log(muY) 

# Desviaciones estandard de las vars explicativas numéricas
sqrt(diag(var(DF[,c("pH","ResinP","Clay","Silt","C")])))

# Efecto relativo (%) de cambiar 1 SD sobre muY:
100*(exp(fixef(m1)[c(5:9)])-1)

# Pendientes de la parte Pr(Y=0)
fixef(m1)[10:14]

# máximo efecto de cambio en 1 Desviaciones estandard de las vars explicativas sobre Pr(Y=0)
fixef(m1)[10:14]/4

# Convergencia de las cadenas
# Dist posterior de cada parametro para cada cadena de m3.brms 
mcmc_dens_overlay(m1, regex_pars = c("^b"))+
  geom_density(lwd=1.2, alpha=0.1)+
  theme_bw()+
  ylab("Dens. Probabilidad")+
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position="none")

# Dist posterior de cada parametro para cada cadena de m3.brms 
mcmc_dens_overlay(m1, regex_pars = c("^s"))+
  geom_density(lwd=1.2, alpha=0.1)+
  theme_bw()+
  ylab("Dens. Probabilidad")+
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 14),
        legend.position="none")

# Trace plots: degree of mixing of the three chains
mcmc_trace(m1, regex_pars = c("^b"), size=0.3)+
  theme_bw()+
  ylab("Valor del Parametro")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=14),              
        strip.background = element_blank(), 
        strip.text = element_text(size = 14),
        legend.position="none") 

# Autocorrelacion de los parametros estimados por cadena: 
mcmc_acf(m1, regex_pars = c("^b"))+
  geom_line(size=1)+
  scale_x_continuous(limits=c(1,20))+
  scale_y_continuous(limits=c(-0.25,0.25))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) 

#  Distribuciones posteriores de (algunos) parametros
variables(m1)[c(1,3:9)] # parametros relacionados con log(mediaY)
mcmc_areas(m1, regex_pars =variables(m1)[c(1,3:9)])+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))
 
# parametros relacionados con logit(Y)
mcmc_areas(m1, regex_pars = c("^b_hu"))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16))

# Algunas distribuciones posteriores del modelo m1 
post.m1=as_draws_df(m1,variable = "^b_", regex = T)

bayes_R2(m1)
# Graficos de distr posteriores de los efectos poblacionales
mcmc_plot(m1, regex_pars = "^b", type="intervals", prob_outer = 0.95)+ 
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16)) 

bayes_R2(m1) # % de varianza explicada

#############################################################################################
#   6. Evaluar la calida de ajuste del modelo estadístico
#############################################################################################
# (a) distribucion predictiva posterior y pp-checks
dist.pred.post.m1=predict(m1, ndraws=1e3, summary=F) 

ppc.density.m1=ppc_dens_overlay(y=DF$actdiest,yrep=dist.pred.post.m1, trim = F, size = 0.5,alpha = 1)+
  xlab("Act. diesterasa")+ 
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
ppc.mean.m1=ppc_stat(y=DF$actdiest,yrep=dist.pred.post.m1,stat = "mean", binwidth = 0.05)+  
  xlab("Media Act. diesterasa")+
  theme_bw()+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")

# Uso de los pp-checks a partir de la dist predictiva posterior
cero=function (x){sum(x==0)/length(x)} # proporcion de ceros en un vector
ppc_ceros.m1=ppc_stat(y=DF$actdiest,yrep=dist.pred.post.m1,stat = cero)+
  xlab("Prob (Act. diest=0)")+ 
  theme_bw()+
  theme(axis.text = element_text(size=14,hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=16),
        legend.position = "none")
plot_grid(ppc.density.m1, ppc.mean.m1,ppc_ceros.m1, ncol=3,
          labels = LETTERS[1:3],align="hv",label_x=0.85, label_y=0.95) 

# (b) ANALISIS de residuos con RQR
# simulación de los RQR
qres.m1=createDHARMa(simulatedResponse = t(dist.pred.post.m1), 
                     observedResponse = DF$actdiest, 
                     fittedPredictedResponse = apply(dist.pred.post.m1, 2, median), 
                     integerResponse = T) 
# Convierte los RQR a dist Normal
res.m1=data.frame(res=qnorm(residuals(qres.m1))) 
# Pone todo en un DF para hacer los graficos
res.m1=cbind(res.m1, 
             DF[,c("Clay.s","Resin.P.s","C.s","Silt.s","pH.s","rootdiamscore")],
             fitted=fitted(m1)[,1], # average of fitted values
             pareto=loo(m1, pointwise=T)$diagnostics$pareto_k) #LOO_CV Pareto k

# Los graficos
qq.m1=ggplot(data=res.m1, mapping=aes(sample = res)) + 
  stat_qq_point(aes(col=rootdiamscore))+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position = "none")
# Residuals vs fitted
res.fit.m1= ggplot(data=res.m1, aes(x=fitted,y=res))+ 
  geom_jitter(aes(color=rootdiamscore), size=1)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position = "none")
# Residuals vs explanatory variables
res.Resin.m1=ggplot(data=res.m1, aes(x=Resin.P.s,y=res, col=rootdiamscore))+ 
  geom_jitter(size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="Resin.P (std)",y="Residuals")+
  theme(axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
res.Clay.m1=ggplot(data=res.m1, aes(x=Clay.s,y=res, color=rootdiamscore))+ 
  geom_jitter(size=1)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+
  labs(x="Clay (std)",y="Residuals")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
res.Silt.m1=ggplot(data=res.m1, aes(x=Silt.s,y=res, color=rootdiamscore))+ 
  geom_jitter(size=1)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+
  labs(x="Silt (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
res.pH.m1=ggplot(data=res.m1, aes(x=pH.s,y=res, color=rootdiamscore))+ 
  geom_jitter(size=1)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+
  labs(x="pH (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
res.C.m1=ggplot(data=res.m1, aes(x=C.s,y=res, color=rootdiamscore))+ 
  geom_jitter(size=1)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+
  labs(x="C (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
res.Root.m1=ggplot(data=res.m1, aes(x=rootdiamscore,y=res, color=rootdiamscore))+ 
  geom_jitter(size=1, width=0.05, alpha=0.5)+
  geom_boxplot(alpha=0.2)+
  geom_hline(yintercept =0,linetype = 2, linewidth=1.2)+
  theme_bw()+
  labs(x="Rootdiamscore",y="Residuals")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
Pareto.m1=ggplot(data=res.m1, aes(x=1:nrow(res.m1),y=pareto))+ 
  geom_jitter(aes(color=rootdiamscore), size=1, height = 0.05)+ 
  theme_bw()+
  labs(y="Pareto's k")+
  geom_hline(yintercept =0.5,linetype = 2)+ 
  geom_hline(yintercept =0.7,linetype = 2)+
  geom_hline(yintercept =1,linetype = 2)+
  theme(axis.title.x=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18),
        legend.position = "none")
plot_grid(res.fit.m1,res.Resin.m1,res.C.m1, 
          res.Clay.m1, res.Silt.m1,res.pH.m1,
          res.Root.m1, qq.m1,Pareto.m1,
          ncol=3,labels = LETTERS[1:9],
          align="hv",label_x=0.85, label_y=0.95) 

# busca valores de residuales > +2 y < -2:
DF[res.m1$res>2 | res.m1$res< (-2),1:12]

summary(DF[res.m1$res>2 | res.m1$res< (-2),1:12]) # estadisticos descriptivos de estos 17 valores
summary(DF[,1:12])

# Efectos de grupo en m1
str(ranef(m1)) # estructura de la lista - hay que extraer y separar los 3 componentes
ranef.fam=as.data.frame(ranef(m1)$fam) # efectos de grupo para familia
names(ranef.fam)

ranef.sp=as.data.frame(ranef(m1)$spsab) # efectos de grupo para especie
names(ranef.sp)

ranef.field=as.data.frame(ranef(m1)$Fieldnumber) # efectos de grupo para Fieldnumber
names(ranef.field)

# plots de visualización de los efectos de grupo para familia
Int.fam1=ggplot(data=ranef.fam, aes(x=reorder(rownames(ranef.fam),Estimate.Intercept), y=Estimate.Intercept))+
 geom_point()+
 geom_linerange(aes(ymin=Q2.5.Intercept, ymax= Q97.5.Intercept))+
 theme_bw()+
 ylab("Intercepto log(media)")+
 coord_flip()+
 theme(axis.title.y=element_blank(), 
       axis.title.x=element_text(size=16), 
       axis.text.y = element_text(size=6),
       axis.text.x = element_text(size=14))

Int.fam2=ggplot(data=ranef.fam, aes(x=reorder(rownames(ranef.fam),Estimate.hu_Intercept), y=Estimate.hu_Intercept))+
  geom_point()+
  geom_linerange(aes(ymin=Q2.5.hu_Intercept, ymax= Q97.5.hu_Intercept))+
  theme_bw()+
  scale_y_continuous(breaks=seq(from=-7, to=7, by=1))+
  ylab("Intercepto logit(p)")+
  coord_flip()+
  theme(axis.title.y=element_blank(), 
        axis.title.x=element_text(size=16), 
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=14))

# plots de visualización de los efectos de grupo para especie
Int.sp=ggplot(data=ranef.sp, aes(x=reorder(rownames(ranef.sp),Estimate.Intercept), y=Estimate.Intercept))+
  geom_point()+
  geom_linerange(aes(ymin=Q2.5.Intercept, ymax= Q97.5.Intercept))+
  theme_bw()+
  labs(y="Intercepto", x="Especies")+
  scale_y_continuous(breaks=seq(from=-2, to=2, by=0.5))+
  coord_flip()+
  theme(axis.title=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank())
Pend.sp=ggplot(data=ranef.sp, aes(x=reorder(rownames(ranef.sp),Estimate.Intercept), y=Estimate.Resin.P.s))+
  geom_point()+
  geom_linerange(aes(ymin=Q2.5.Resin.P.s, ymax= Q97.5.Resin.P.s))+
  theme_bw()+
  labs(y="Pendiente", x="Especies")+
  scale_y_continuous(breaks=seq(from=-1, to=1, by=0.25))+
  coord_flip()+
  theme(axis.title=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank())
Int.int.pend=ggplot(data=ranef.sp, aes(x=Estimate.Intercept, y=Estimate.Resin.P.s))+
  geom_point(size=2)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_bw()+
  labs(y="Pend aleatoria", x="Intercepto aleatorio")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14))

# plots de visualización de los efectos de grupo para Fieldnumber
Int.field=ggplot(data=ranef.field, aes(x=reorder(rownames(ranef.field),Estimate.Intercept), y=Estimate.Intercept))+
  geom_point()+
  geom_linerange(aes(ymin=Q2.5.Intercept, ymax= Q97.5.Intercept))+
  theme_bw()+
  labs(y="Intercepto  Field", x="Fieldnumber")+
  coord_flip()+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14))
plot_grid(Int.fam1,Int.fam2,Int.field,
          Int.sp,Pend.sp, Int.int.pend, 
          ncol=3,labels = LETTERS[1:6],
          align="hv",label_x=0.9, label_y=0.2) 


# Curvas condicionales predichas por el modelo
# Condicionales  para cada var explicativa/interaccion:
m1.cond.eff=conditional_effects(m1) 
names(m1.cond.eff) # nombre de los efectos conditionales 

plot(m1.cond.eff, plot = F, points=T)[[1]]+
  theme_bw()+
  labs(y="Act. diesterasa",x="rootdiamscore")+ 
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position = "none")
plot(m1.cond.eff, plot = F,points=T)[[2]]+
  theme_bw()+
  labs(y="Act. diesterasa",x="ResinP (std)")+ 
  ylim(0,60)+
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank())
plot(m1.cond.eff, plot = F,points=T)[[3]]+
  theme_bw()+
  labs(y="Act. diesterasa",x="Clay (std)")+ 
  ylim(0,60)+
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank())
plot(m1.cond.eff, plot = F,points=T)[[4]]+
  theme_bw()+
  labs(y="Act. diesterasa",x="pH (std)")+ 
  ylim(0,60)+
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank())
plot(m1.cond.eff, plot = F,points=T)[[5]]+
  theme_bw()+
  labs(y="Act. diesterasa",x="C (std)")+ 
  ylim(0,60)+
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank())










# visualizacion de los efectos aleatorios 
ranef.Fieldnumber.Plot=ggplot(data=ranefFieldnumber.Plot, aes(x=row.names(ranefFieldnumber.Plot), y=int))+
  geom_point(size=3, col="red")+
  theme_bw()+
  labs(x="Site.Plot", y="Intercept")+
  geom_hline(yintercept = 0,linetype = 2, size=1.2)+
  coord_flip()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
ranef.Spsab.int=ggplot(data=ranefSpsab, aes(x=reorder(row.names(ranefSpsab), int), y=int))+
  geom_point(size=2, col="red")+
  theme_bw()+
  labs(x="Spp name", y="Intercepts")+
  geom_hline(yintercept = 0,linetype = 2, size=1.2)+
  coord_flip()+
  theme(axis.title=element_text(size=18), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=8))
ranef.Spsab.slope=ggplot(data=ranefSpsab, aes(x=reorder(row.names(ranefSpsab), int), y=slope))+
  geom_point(size=2, col="red")+
  theme_bw()+
  labs(x="Spp name", y="Slopes")+
  geom_hline(yintercept = 0,linetype = 2, size=1.2)+
  coord_flip()+
  theme(axis.title=element_text(size=18), 
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=8))
ranef.Fam=ggplot(data=ranefFam, aes(x=reorder(row.names(ranefFam), int), y=int))+
  geom_point(size=3, col="red")+
  theme_bw()+
  labs(x="Family", y="Intercept")+
  geom_hline(yintercept = 0,linetype = 2, size=1.2)+
  coord_flip()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(ranef.Spsab.int, ranef.Spsab.slope, ranef.Fieldnumber.Plot, ranef.Fam,nrow=2)

# Data frame con los interceptos y pendientes por especie (efectos fijos + aleatorios)
int.sl.spp=coef(m1)$cond$spsab [,c('(Intercept)',"Resin.P.s" )]

# Curvas predichas por spp: act diest vs resinP(std)
pred.spp=data.frame(Resip=seq(from=-3, to=3, by=0.1)) # set de valores de ResinP (std)

# Este loop genera valores predichos por el modelo de media (actdiest) de cada una de las 105 spp 
# solo usando el valor de ResinP (std), i.e. suponiendo las otras variables =0 (en sus medias)
# Estas predicciones NO incorporan la parte ZI (los ceros!) del modelo ZAGamma 
for (i in 1:length(rownames(int.sl.spp))) 
  pred.spp[,i+1]=exp(int.sl.spp[i,1]+(int.sl.spp[i,2]*pred.spp$Resip))
# pone nombre spp a cada columna de valores predichos
names(pred.spp)[2:(length(rownames(int.sl.spp))+1)]= rownames(int.sl.spp) 
# el exp es para invertir la funcion de enlace log empleada en el modelo ZAGamma

# Reformatea el dataframe de predicciones a "long format", excluyendo (por ahora) ResinP
pred.spp=melt(data=pred.spp[,-1],variable.name="spsab",value.name="act.diest")
# Reincopora variable ResinP en el dataframe reformateado: los 61 valores * 105 spp= 6405 filas
pred.spp$Resip=rep(seq(from=-3, to=3, by=0.1),times=length(rownames(int.sl.spp)))
dim(pred.spp)
#[1] 6405    3

# dataframe de predicciones de act.diest para la "especie promedio" usando los efectos fijos
Spp.prom=data.frame(Resip=seq(from=-3, to=3, by=0.1)) # set de valores de ResinP (std)
Spp.prom$act.diest=exp(fixef(m1)$cond[1]+ (fixef(m1)$cond[4]*Spp.prom$Resip))

Curvas.SPP=ggplot(data=pred.spp, aes(x=Resip, y=act.diest, col=spsab))+
  geom_line(linewidth=1, alpha=0.3)+
  geom_line(data=Spp.prom, aes(x=Resip, y=act.diest), linewidth=1, linetype="dashed",color="black")+
  theme_bw()+
  labs (x="ResinP (std)", y="Avg(act diestesterasa)")+
  theme(legend.position="none", 
        axis.title = element_text(size=18),
        axis.text = element_text(size=16))

# Curvas condicionales predichas por el modelo (ambos componentes: con y ZI + los efectos aleatorios)
pred.m1=ggpredict(model=m1,terms=c("Resin.P.s[all]"), type="zi_random")
Curv.cond.m1=plot(pred.m1, add.data = T)+
  theme_bw()+
  labs(x="Resin.P (estand.)", y="diestesterase activity")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        plot.title=element_blank())
Curv2.cond.m1=plot(pred.m1, add.data = T)+
  theme_bw()+
  scale_y_log10()+ # mejora la visualizacion
  labs(x="Resin.P (estand.)", y="diestesterase activity")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        plot.title=element_blank())
