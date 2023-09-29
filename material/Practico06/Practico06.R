setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Practico06/")

mc.cores=parallel::detectCores()
options(digits=3)
packages.needed=c("ggplot2","fitdistrplus","brms","future","GGally",  "ggeffects","performance", "cowplot",
                  "ggdist", "arm", "DHARMa","qqplotr","reshape2", "bayesplot", "parallel", "doBy","ggbreak", "ggridges", "ggdist","pROC")
lapply(packages.needed,FUN=require,character.only=T)
##############################

DF2=read.csv("bangladesh Pr06.csv", header=T, stringsAsFactors = T)
str(DF2)
DF2$distrito=as.factor(DF2$distrito)

# Analisis exploratorio de datos
with(DF2, ftable(metodo,educ))
with(DF2, ftable(religion,educ,metodo))

# Crea dataframe para hacer ggplot a partir de tabulaciones
a=as.data.frame(with(DF2, ftable(religion,educ,metodo)))
ggplot(data=a)+ 
  geom_bar(mapping = aes(x=educ,y=Freq, fill=metodo),position="fill", stat="identity")+ 
  facet_wrap(~religion)+
  theme_bw()+ 
  labs(x="Nivel de educación",y="Proporción")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18),
        axis.text=element_text(size = 14), 
        axis.title = element_text(size = 18))

# Crea dataframe para hacer ggplot a partir de tabulaciones
a2=as.data.frame(with(DF2, ftable(vive, religion, metodo)))
ggplot(data=a2)+ 
  geom_bar(mapping = aes(x=metodo,y=Freq, fill=vive),position="fill", stat="identity")+ 
  facet_wrap(~religion)+
  theme_bw()+ 
  labs(x="Método",y="Proporción")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18),
        axis.text=element_text(size = 14),
        axis.title = element_text(size = 18))
ggplot(data=DF2,aes(x=as.factor(metodo),y=desc.viva))+ 
  geom_boxplot()+
  theme_bw()+ 
  facet_wrap(vive~religion)+
  labs(x="Método",y="# hijos")+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18),
        axis.text=element_text(size = 14),
        axis.title = element_text(size = 18))

summaryBy(desc.viva~metodo+vive+religion,data=DF2,FUN=mean)

ggplot(DF2, aes(x=metodo, y=desc.viva)) +
  geom_jitter(alpha=0.3, size=2, aes(color=vive), position_jitter=0.2)

# MODELO ESTADISTICO
# previas SIN "efectos de grupos" 
get_prior(formula=metodo~educ*religion+vive+religion*desc.viva, data=DF2,
          family=categorical())

# previas CON efectos de grupos 
get_prior(formula=metodo~educ*religion+vive+religion*desc.viva+(1|distrito),data=DF2,
          family=categorical())

prior.m3 <- c(#set_prior("normal(0,3)",class = "b"),
              set_prior("normal(0,2)", class = "Intercept", dpar="mu2"))

m3=brm(formula=metodo~educ*religion+vive+religion*desc.viva+(1|distrito),
       family=categorical(),data=DF2,prior = prior.m3, future=T,
       warmup = 1000,chains=3, iter=2000,thin=2)

# resultados del modelo m3
summary(m3)
levels(DF2$educ)
levels(DF2$religion)
levels(DF2$vive)

# Nombre de los parametros del modelo
variables(m3) 
# Nombre de parametros con efectos poblacionales
mu2=variables(m3)[grep("b_mu2", variables(m3))]
mu3=variables(m3)[grep("b_mu3", variables(m3))]
mu4=variables(m3)[grep("b_mu4", variables(m3))]
# convergencia de las cadenas a una distr posterior estacionaria
mcmc_dens_overlay(m3,pars = mu2)+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 12),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 
# autocorrelacion de los estimados muestrados
mcmc_acf(m3, pars = mu3)+
  geom_line(size=1)+
  scale_x_continuous(limits=c(1,20))+
  scale_y_continuous(limits=c(-0.25,0.25))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=14),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) 

# distribuciones posteriores de algunos parametros poblacionales (no mostrada)
# medias e IntCred 95% de los parametros poblacionales
mcmc_areas(m3, prob=0.62, prob_outer = 0.95, pars=mu4)+
  geom_density(lwd=1.2, alpha=0.9)+ 
  theme_bw()+
  theme(axis.text.y=element_text(size = 14, hjust=0.4),
        axis.text.x=element_text(size = 16))


### ANALISIS DE RESIDUALES
# Distr predictiva posterior de modelo m3
post.pred.m3=predict(m3, ndraws=1e3, summary=F)
qres.m3=createDHARMa(simulatedResponse = t(post.pred.m3), 
                     observedResponse = DF2$metodo, 
                     fittedPredictedResponse = apply(post.pred.m3,2, median),integerResponse=T)
res.m3=data.frame(res=qnorm(residuals(qres.m3))) 

# estructura de los valores predichos
str(as.data.frame(fitted(m3, ndraws=1))) # average of fitted values

# creates a dataframe with avg values of residuals, of fitted, values, Pareto'k of LOO-CV, and explanatory variables
res.m3=cbind(res.m3, 
             DF2[,c("educ","religion", "vive", "desc.viva")],
             fitted1=as.data.frame(fitted(m3, ndraws=1000))[,1],
             fitted2=as.data.frame(fitted(m3, ndraws=1000))[,5],
             fitted3=as.data.frame(fitted(m3, ndraws=1000))[,9],
             fitted4=as.data.frame(fitted(m3, ndraws=1000))[,13],
             pareto=loo(m3, pointwise=T)$diagnostics$pareto_k)

res.fit1.m3= ggplot(data=res.m3, aes(x=fitted1,y=res))+ 
  geom_point(size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted1",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
res.fit2.m3= ggplot(data=res.m3, aes(x=fitted2,y=res))+ 
  geom_point(size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted2",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")
res.fit3.m3= ggplot(data=res.m3, aes(x=fitted3,y=res))+ 
  geom_point(size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted3",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")
res.fit4.m3= ggplot(data=res.m3, aes(x=fitted4,y=res))+ 
  geom_point(size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted4",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")
res.fit.musul= ggplot(data=res.m3, aes(x=religion,y=res))+ 
  geom_boxplot()+
  geom_jitter(width=0.05)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
res.fit.educ= ggplot(data=res.m3, aes(x=educ,y=res))+ 
  geom_boxplot()+
  geom_jitter(width=0.05)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
res.fit.vive= ggplot(data=res.m3, aes(x=vive,y=res))+ 
  geom_boxplot()+
  geom_jitter(width=0.05)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
res.fit.desc= ggplot(data=res.m3, aes(x=desc.viva,y=res))+ 
  geom_jitter(width=0.05)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="desc.viva",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position="none")
plot_grid(res.fit1.m3, res.fit2.m3,res.fit3.m3, res.fit4.m3,
          res.fit.musul, res.fit.educ,res.fit.vive,res.fit.desc, ncol=4,
          labels = LETTERS[1:8],align="hv",label_x=0.8, label_y=0.95)

# qqplot para los residuales y los efectos de grupo 
qq.m3=ggplot(data=res.m3, mapping=aes(sample = res)) + 
  stat_qq_point()+
  theme_bw()+
  annotate("text",x=0, y=-2.5, label="datos", size=6)+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# Efectos de grupo : hay un efecto de grupo para 3 valores de la var respuesta vectorial
ef.gr.m3=as.data.frame(ranef(m3))
ef.gr.m3=data.frame(gr2=ef.gr.m3[,1],
                    gr3=ef.gr.m3[,5],
                    gr4=ef.gr.m3[,9])
qq.gr2=ggplot(data=ef.gr.m3, mapping=aes(sample = gr2)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  annotate("text",x=0, y=-0.25, label="gr2", size=6)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
qq.gr3=ggplot(data=ef.gr.m3, mapping=aes(sample = gr3)) + 
  stat_qq_point()+
  theme_bw()+
  annotate("text",x=0, y=-0.5, label="gr3", size=6)+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
qq.gr4=ggplot(data=ef.gr.m3, mapping=aes(sample = gr4)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  annotate("text",x=0, y=-1, label="gr4", size=6)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
pareto.m3=ggplot(data=res.m3, aes(x=1:nrow(res.m3),y=pareto))+ 
  geom_point(size=1)+ 
  theme_bw()+
  labs(y="Pareto's k")+
  geom_hline(yintercept =0.5,linetype = 2)+ 
  geom_hline(yintercept =0.7,linetype = 2)+
  geom_hline(yintercept =1,linetype = 2)+
  theme(axis.title.x=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
plot_grid(qq.m3,qq.gr2, qq.gr3,qq.gr4, pareto.m3,
          nrow=2,labels = LETTERS[1:5],align="hv",label_x=0.8, label_y=0.95)

# pp-checks a partir de distr predictivas posteriores
# Densidad de valores de la var de respuesta
ppc.density=ppc_dens_overlay(y=DF2$metodo,yrep=post.pred.m3)+
  xlab("Método")+ 
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")

prop_uno=function(x) mean(x == 1) # función que calcula la prop de 1
prop_dos=function(x) mean(x == 2) # función que calcula la prop de 2
prop_tres=function(x) mean(x == 3) # función que calcula la prop de 3
prop_cuatro=function(x) mean(x == 4) # función que calcula la prop de 4

prepos1=ppc_stat(y=DF2$metodo, yrep=post.pred.m3, stat = "prop_uno")+
  annotate("text",x=0.09, y=60, label="Pr(Y=1)", size=6)+
  theme(legend.position="none",
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
prepos2=ppc_stat(y=DF2$metodo, yrep=post.pred.m3, stat = "prop_dos")+
  annotate("text",x=0.18, y=60, label="Pr(Y=2)", size=6)+
  theme(legend.position="none",
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
prepos3=ppc_stat(y=DF2$metodo, yrep=post.pred.m3, stat = "prop_tres")+
  annotate("text",x=0.08, y=60, label="Pr(Y=3)", size=6)+
  theme(legend.position="none",
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
prepos4=ppc_stat(y=DF2$metodo, yrep=post.pred.m3, stat = "prop_cuatro")+
  annotate("text",x=0.575, y=60, label="Pr(Y=4)", size=6)+
  theme(legend.position="none",
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
plot_grid(prepos1,prepos2, prepos3, prepos4, ncol=2)

ppc_stat_grouped(y=DF2$metodo, yrep=post.pred.m3, stat="prop_cuatro", group=DF2$distrito,binwidth = 0.05)+
  xlab("Probabilidad método 4")+ 
  scale_x_continuous(breaks=seq(from=0, to=1, by=0.25))+
  theme(axis.text = element_text(size=12),
        axis.title=element_text(size=18),
        legend.position = "none")


# Distribuciones marginales:efectos poblacionales  
result.m3=marginal_effects(m3, categorical=T)
length(result.m3)

resm3.1= plot(result.m3, plot=F)[[1]]+ 
  theme_bw()+ 
  labs(y="Probabilidad", x="Educación")+
  theme(legend.position = "none",
        axis.text.y= element_text(size = 16),
        axis.text.x= element_text(size = 16, vjust=0.3), 
        axis.title = element_text(size = 18)) 
  # Educacion
resm3.2=plot(result.m3, plot=F)[[2]]+ 
  theme_bw()+ 
  theme(legend.position=c(0.9,0.5),
        axis.title.y = element_blank(),
        axis.text=element_text(size = 16), 
        axis.title = element_text(size = 18))# Religion
resm3.3=plot(result.m3, plot=F)[[3]]+ 
  theme_bw()+
  labs(y="Probabilidad")+
  theme(legend.position = "none",
        axis.text=element_text(size = 16), 
        axis.title = element_text(size = 18))# Vive
resm3.4=plot(result.m3, plot=F)[[4]]+ 
  theme_bw()+ 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text=element_text(size = 16), 
        axis.title = element_text(size = 18))# Descendencia viva
plot_grid(resm3.1, resm3.2,resm3.3, resm3.4, ncol=2,
          labels = LETTERS[1:4],align="hv",label_x=0.9, label_y=0.95)
