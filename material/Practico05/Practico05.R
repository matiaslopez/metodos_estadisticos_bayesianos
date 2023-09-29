setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Practico05/")

mc.cores=parallel::detectCores()
options(digits=3)
packages.needed=c("ggplot2","fitdistrplus","brms","future","GGally",  "ggeffects","performance", "cowplot",
                  "ggdist", "arm", "DHARMa","qqplotr","reshape2", "bayesplot", "parallel", "doBy","ggbreak", "ggridges", "ggdist","pROC")
lapply(packages.needed,FUN=require,character.only=T)
##############################

library(dplyr)
# importar datos
DF=read.csv("actdiest Pr05.csv", header = T,stringsAsFactors=T) %>%
  mutate(Plot=factor(Plot),
         Fieldnumber=factor(Fieldnumber),
         rootdiamscore=factor(rootdiamscore))
str(DF)
summary(DF)

DF.s=data.frame(
  Fieldnumber=DF$Fieldnumber,
                Plot=DF$Plot,
                actdiest=DF$actdiest,
                fam=DF$fam,
                spsab=DF$spsab,
                scale(DF[ c("pH","Clay", "Silt", "C", "Ntotal", "ResinP")], scale=T, center=T), 
                rootdiamscore=DF$rootdiamscore
                )
str(DF.s)

ggplot(DF, aes(x=Plot, fill=fam)) + geom_bar(color="black") + facet_wrap(.~Fieldnumber)


# Distribucion de prob de la var de respuesta: actdiest
ggplot(DF %>% filter(actdiest>0), aes(x=actdiest)) + geom_density()
ggplot(DF, aes(x=actdiest)) + geom_density()

sin_ceros = DF.s[DF.s$actdiest>0,]$actdiest

norm_fit=fitdist(sin_ceros,"lnorm")
exp_fit=fitdist(sin_ceros,"exp")
gamma_fit=fitdist(sin_ceros,"gamma")


CDF.ZI=cdfcomp(list(norm_fit,exp_fit,gamma_fit),
               addlegend=T,main="",legendtext=c("LogNormal","Exponencial", "Gamma"),
               plotstyle = "ggplot")+
  xlab("Value of actdiest")+
  geom_line(linewidth=0.8)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position = c(0.75,0.25),
        legend.text=element_text(size=14))
QQ.ZI=qqcomp(list(norm_fit,exp_fit, gamma_fit),addlegend=F,main="",legendtext=c("LogNormal","Exponencial", "Gamma"),
             plotstyle = "ggplot")+
  theme_bw()+
  geom_jitter(size=2, height=0.2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        title=element_blank(),
        legend.position = c(0.75,0.25),
        legend.text=element_text(size=14))
plot_grid(CDF.ZI, QQ.ZI, ncol=2)


# Miro los datos a ver si incluyo todos o no
ggpairs(data=DF.s %>% select(-fam, -spsab, -Fieldnumber, -Plot), 
        lower=list(continuous = wrap("smooth", method = "lm")))+
  theme_bw()+ 
  theme(strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=12)) 

summary(DF.s[c("pH", "Silt", "C", "ResinP")])
summaryBy(actdiest~fam+rootdiamscore, data=DF.s, FUN=c(mean, sd))

f1 = bf(actdiest~pH + Silt + C + ResinP + rootdiamscore +
          (1| Fieldnumber) +
          (1+ ResinP | fam)
)

f2 = bf(actdiest~pH + Silt + C + ResinP + rootdiamscore +
          (1| Fieldnumber) +
          (1+ ResinP | spsab),
        hu~ResinP + Clay + pH + C + (1|fam),
        shape~1
)

print(f2)

get_prior(formula=f2,data=DF.s,family=hurdle_gamma()) 

prior.m11 = c(set_prior("normal(0,2)", class = "b")) 
# prior.m12 = c(set_prior("normal(1.5,0.5)",class = "b")) # más informativos
get_prior(formula=f2,data=DF.s,family=hurdle_gamma(), prior=prior.m11) 


# Modelos estadísticos
m11=brm(f2,
        data=DF.s,family=hurdle_gamma(),cores=mc.cores, 
        prior = prior.m11, warmup = 1000,
        chains=3, iter=2000, thin=2,control = list(adapt_delta = 0.95))

summary(m11)

prior_summary(m11)

fixef(m11)

# Cuanto cambia la media por cada unidad que me muevo de la unidad
100*(exp(fixef(m11)[c(5:9)])-1)

a








# Estructura de los datos
table(DF$Nest)
with (DF, table(SexParent,FoodTreatment))
with (DF, ftable(Nest, SexParent,FoodTreatment))

# Usar offset?
ggplot(data=DF, aes(x=BroodSize, y=SiblingNegotiation)) +
  geom_boxplot(aes(x=as.factor(BroodSize)))+
  geom_jitter(alpha=0.3, size=2,position = position_jitter(width = .2))+
  labs(x="BroodSize", y="SiblingNegotiation")+ 
  theme_bw()+
  stat_summary(fun=mean, geom="point", shape=19, size=4,color="black")+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16))
  
summaryBy(SiblingNegotiation~BroodSize, data=DF,FUN=mean)

# Graficos exploratorios 
ggplot(DF, aes(x=interaction(SexParent, FoodTreatment), y=SiblingNegotiation/BroodSize)) + 
  geom_boxplot(alpha=0.1)+
  theme_bw()+
  labs(x="sexo y trat", y="Tasa")+ 
  geom_jitter(alpha=0.3, size=2,width = 0.2)+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16))

summaryBy(SiblingNegotiation~SexParent+FoodTreatment, data=DF, FUN=c(mean, sd))
DF$Tasa = DF$SiblingNegotiation/DF$BroodSize
summaryBy(Tasa~SexParent+FoodTreatment, data=DF, FUN=c(mean, sd))


# variacion entre zonas: Nest
ggplot(DF, aes(x=Nest, y=SiblingNegotiation/BroodSize)) + 
  geom_boxplot(alpha=0.1)+
  theme_bw()+
  stat_summary(fun=mean, geom="point", shape=19, size=2,color="black")+
  ylab("Tasa")+xlab("Nest")+ 
  theme(axis.text.x = element_text(size = 12,angle =25,hjust = 1),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=16))


# Modelos estadisticos
# Que priors se requieren especificar y cuales son sus valores por defecto?
get_prior(formula=SiblingNegotiation~SexParent+FoodTreatment+offset(log(BroodSize))+
          (1|Nest),data=DF,family=negbinomial()) 
prior.m11 = c(set_prior("normal(0,5)", class = "b")) # debilmente informativos
prior.m12 = c(set_prior("normal(1.5,0.5)",class = "b")) # más informativos

# Qué significa el prior para los coeficientes b?
ggplot(data = data.frame(x = seq(from = -2, to = 2, by = .01)),aes(x = x)) + 
  geom_ribbon(aes(ymin = 0, ymax = dnorm(x, 1.5, 0.1)),fill = "blue", alpha=0.5)+
  geom_ribbon(aes(ymin = 0, ymax = dnorm(x, 0, 1)),fill = "red", alpha=0.5)+
  labs(x="coeficiente en escala g(E(Y))", y="Dens Prob")+
  theme_bw()+
  scale_x_continuous(breaks = -2:2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
        
# Modelos estadísticos
m11=brm(SiblingNegotiation~SexParent+FoodTreatment+offset(log(BroodSize))+(1|Nest),
        data=DF,family=negbinomial(),cores=mc.cores, prior = prior.m11,warmup = 1000,
        chains=3, iter=2000, thin=2,control = list(adapt_delta = 0.95))
m12=brm(SiblingNegotiation~SexParent+FoodTreatment+offset(log(BroodSize))+(1|Nest),
        data=DF,family=negbinomial(),cores=mc.cores, prior = prior.m12,warmup = 1000,
        chains=3, iter=2000, thin=2,control = list(adapt_delta = 0.95))

# Resultados
summary(m11)
summary(m12)

# Sigamos con m11
# Dist posteriores del modelo m11
post.m11=as_draws_df(m11)
head(post.m11,2)

# # Convergencia de las cadenas para efectos poblacionales
mcmc_dens_overlay(m11,regex_pars = c("^b", "sd", "shape"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 14)) 

# # Convergencia de las cadenas para efectos de grupo
mcmc_dens_overlay(m11,regex_pars = c("r_Nest"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 12),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

#Grafico de trazas
mcmc_trace(m11, regex_pars = c("^b", "sd", "shape"))+
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

# Autocorrelacion de los parametros estimados por cadena: 
mcmc_acf(m11, regex_pars = c("^b", "sd", "shape"))+
  geom_line(size=1)+
  scale_x_continuous(limits=c(1,20))+
  scale_y_continuous(limits=c(-0.25,0.25))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) 

# medias e IntCred 95% de los parametros poblacionales
mcmc_areas(m11, prob=0.62, prob_outer = 0.95, regex_pars = c("^b", "sd", "shape"))+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  scale_y_discrete(labels=c("intercepto","dif.Food","dif.Sex","sd grupos", "shape"))+
  scale_x_continuous(breaks=seq(from=-1.5, to=1, by=0.25))+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))

interc=variables(m11)[grep("Intercept]", variables(m11))]# interceptos por sitio
mcmc_areas(m11, prob=0.62, prob_outer = 0.95, pars = interc)+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  geom_vline(xintercept = 0, col="red", size=1)+
  scale_y_discrete(labels=rev(unique(DF$Nest)))+
  scale_x_continuous(breaks=seq(from=-1, to=1, by=0.25))+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))

# Efectos de grupo : no graficado pero usado más abajo
ef.gr.m11=as.data.frame(ranef(m11))
names(ef.gr.m11)=c("int", "int.SE", "int.Q2.5","int.Q97.5")
head(ef.gr.m11)
# Interceptos
ggplot(ef.gr.m11, aes(x=rownames(ef.gr.m11), y=int))+geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = int.Q2.5, ymax = int.Q97.5), width = 0.5)+
  theme_bw()+
  ylab("intercepto")+xlab("Nest")+ 
  coord_flip()+
  geom_hline(yintercept =0,linetype = 2, size=1, col="red")
theme(axis.text = element_text(size= 10),
      axis.text.x = element_text(angle =90,hjust = 1))+
  
# % de varianza explicada por el modelo ajustado
bayes_R2(m11)

## Validación de m11 por analisis de residuales
post.pred.m11=predict(m11, ndraws=1e3, summary=F)
qres.m11=createDHARMa(simulatedResponse = t(post.pred.m11), 
                     observedResponse = DF$SiblingNegotiation, 
                     fittedPredictedResponse = apply(post.pred.m11,2, median),integerResponse=T)
res.m11=data.frame(res=qnorm(residuals(qres.m11))) 

res.m11=cbind(res.m11, 
             DF[,c("FoodTreatment","SexParent")],
             fitted=fitted(m11, ndraws=1000)[,1], # average of fitted values
             pareto=loo(m11, pointwise=T)$diagnostics$pareto_k) #LOO_CV Pareto k

# Residuales vs SexParent 
res2=ggplot(data=res.m11, aes(x=SexParent,y=res))+ 
  geom_boxplot(alpha=0.2)+ 
  geom_jitter(alpha=0.3,size=2,position = position_jitter(width = .2))+
  theme_bw()+
  labs(y="Residuales", x="SexParent")+
  coord_flip()+
  theme(axis.text = element_text(size= 16),
        axis.title = element_text(size= 18))
  
# Residuales vs FoodTreatment
res3=ggplot(data=res.m11, aes(x=FoodTreatment,y=res))+ 
  geom_boxplot(alpha=0.2)+ 
  geom_jitter(alpha=0.3, size=2,width = 0.2)+
  theme_bw()+
  coord_flip()+
  theme(axis.text = element_text(size= 16),
        axis.title = element_text(size= 18))

# QQ plot Residuales 
res4=ggplot(data=res.m11, mapping=aes(sample = res)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")

# Leave-one-out cross validation
res5=ggplot(data=res.m11, aes(x=1:nrow(res.m11),y=pareto))+ 
  geom_point(size=1)+ 
  theme_bw()+
  labs(y="Pareto's k")+
  geom_hline(yintercept =0.5,linetype = 2)+ 
  geom_hline(yintercept =0.7,linetype = 2)+
  geom_hline(yintercept =1,linetype = 2)+
  theme(axis.title.x=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18),
        legend.position="none")

# QQ plot interceptos de efectos grupales 
qq.int.m11=ggplot(data=ef.gr.m11, mapping=aes(sample=int))+ 
  stat_qq_point(size=2)+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(y="mean of group effects")+
  theme_bw()+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

plot_grid(res2,res3,res5,NULL,res4,  qq.int.m11, ncol=3,
             labels = c(LETTERS[1:3],"",LETTERS[4:5]),
             align="hv",label_x=0.9, label_y=0.95)


# Distr predictiva posterior de m11
m11.pred.post=posterior_predict(m11, ndraws = 500) 

# Densidad de valores de la var de respuesta
ppc.density=ppc_dens_overlay(y=DF$SiblingNegotiation,yrep=m11.pred.post)+
  scale_x_continuous(limits = c(0,40))+
  xlab("SiblingNegotiation")+ 
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
ppc.media=ppc_stat(y=DF$SiblingNegotiation, yrep=m11.pred.post, stat = mean,binwidth = 0.2)+
  xlab("Media")+
  theme(legend.position = "none",
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
ppc.sd=ppc_stat(y=DF$SiblingNegotiation, yrep=m11.pred.post, stat = sd,binwidth = 0.2)+
  xlab("Sd")+
  theme(legend.position = "none",
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
plot_grid(ppc.density, ppc.media, ppc.sd, ncol=3,
          labels = LETTERS[1:3],align="hv",label_x=0.85, label_y=0.95) 

# pp checks por grupo
ppc_stat_grouped(y=DF$SiblingNegotiation, yrep=m11.pred.post, stat="mean", 
                              group=DF$Nest,binwidth = 1)+
  theme(legend.position = "none")
ppc_stat_grouped(y=DF$SiblingNegotiation, yrep=m11.pred.post, stat = sd,
                 group=DF$Nest,binwidth = 1)+
  theme(legend.position = "none")

# Efectos Conditionales 
m11.cond.eff=conditional_effects(m11) 
names(m11.cond.eff) # names of the conditional effects fitted
cond.m11.sexo=plot(m11.cond.eff, plot = F, points=T,
                  point_args = list(width =0.05, col="blue"))[[1]]+
  theme_bw()+
  labs(x="Sexo",y="Sibling Negotation")+ 
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18))
cond.m11.trat=plot(m11.cond.eff, plot = F, points=T,
                   point_args = list(width =0.05, col="blue"))[[2]]+
  theme_bw()+
  labs(x="Tratamiento",y="Sibling Negotation")+ 
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18))
plot_grid(cond.m11.trat, cond.m11.sexo,ncol=2,
          labels = LETTERS[1:2],align="hv",label_x=0.2, label_y=0.95)
