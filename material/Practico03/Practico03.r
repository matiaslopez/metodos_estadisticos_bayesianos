options(digits=3)
paquetes<-c("ggplot2","fitdistrplus","brms","future","car", "GGally","gamlss.dist",
           "ggeffects","ggbreak", "cowplot", "doBy", "DHARMa","qqplotr", "bayesplot","arm")
lapply(paquetes,FUN=require,character.only=T)
options(mc.cores=parallel::detectCores())
setwd("/home/mlopez/git/metodos_estadisticos_bayesianos/material/Practico03/")

#############################

# Curvas de la distribucion beta para diferentes combinaciones de parametros 
ggplot(data=data.frame(x=seq(from=0.01, to=0.99, by=0.01)), aes(x))+
  stat_function(fun = dbeta, n = 1000, args = list(shape1=3, shape2=1), size=1, col="red")+
  stat_function(fun = dbeta, n = 1000, args = list(shape1=1, shape2=3), size=1, col="blue")+
  stat_function(fun = dbeta, n = 1000, args = list(shape1=0.5, shape2=0.5), size=1, col="dark green")+
  stat_function(fun = dbeta, n = 1000, args = list(shape1=1, shape2=1), size=1, col="black")+
  stat_function(fun = dbeta, n = 1000, args = list(shape1=3, shape2=4), size=1, col="magenta")+
  labs(y="Dens Probabilidad", x= "Y")+
  theme_bw()+
  theme (axis.title=element_text(size=18), 
         axis.text=element_text(size=16))+
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

################################################################################################
# los datos
DF=read.csv ("Plant cover Pr03.csv", header=T, stringsAsFactors=T)
summary(DF)


# Estandarzando las variables explicativas numericas 
DF.s=data.frame(scale(DF[ c("Gr.pressure","patch.size","inter.patch.dist")], scale=T, center=T), 
                Totalcover=DF$Totalcover, Prod=DF$Prod)

# Graficos exploratorios
ggpairs(data=DF,mapping = aes(col = Prod), 
        lower=list(continuous = wrap("smooth", method = "lm")))+
  theme_bw()+ 
  theme(strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=12)) 

# Modelos estadisticos
#  a) sin vars explicativas en el parametro de dispersion  phi:
formula.m1=bf(Totalcover~Gr.pressure*Prod + patch.size*Prod+inter.patch.dist)+
                   lf(phi~1)
get_prior(formula.m1, data=DF.s,family='Beta')

#  b) with explanatory variables for the dispersion parameter phi:
formula.m2=bf(Totalcover~Gr.pressure*Prod + patch.size*Prod+inter.patch.dist)+
                lf(phi~Gr.pressure*Prod + patch.size*Prod+inter.patch.dist)
get_prior(formula.m2, data=DF.s,family='Beta')

# Defining priors for the more complex model:
prior.m1 = c(set_prior("normal(0,0.5)",class = "b"),
             set_prior("normal(logit(0.5),(logit(0.8)-logit(0.2))/4", class = "Intercept"))

# Modelos estadisticos a ajustar:
m1=brm(formula=formula.m1, family=Beta (link = "logit", link_phi = "log"), 
        warmup = 1000,data=DF.s,prior=prior.m1,chains=3, iter=3000, future=T,
        control = list(adapt_delta =0.9999, max_treedepth=15))
m2=brm(formula=formula.m2, family=Beta(link = "logit", link_phi = "log"), 
            warmup = 1000,data=DF.s, prior=prior.m1,chains=3, iter=3000, future=T,
            control = list(adapt_delta =0.9999, max_treedepth=15))

# añadir el criterio LOO para hacer la selección de modelos:
m1=add_criterion(m1, "loo")
m2=add_criterion(m2, "loo")
# compara modelos:
loo_compare(m1,m2, criterion= "loo")

# examina el resultado global del modelo seleccionado
summary(m1)

# Dist posteriores del modelo m1
post.m1=as_draws_df(m1,variable = c("^b_"), regex = T)
head(post.m1,2)

# Convergencia de las cadenas de m1
mcmc_dens_overlay(m1, regex_pars = c("^b"))+
  geom_density(lwd=1.2, alpha=0.1)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

# Autocorrelacion de los valores muestrados de los parametros por cadena
mcmc_acf(m1, regex_pars = "b")+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(1,10), breaks=seq(1,9,2))+
  scale_y_continuous(limits=c(-0.2,0.6))+
  theme(axis.text.y=element_text(size = 16),
        axis.text.x=element_text(size = 14, angle=90),
        axis.title = element_text(size = 18))
        
# Plot de medias e Int Credibilidad de los parámetros
mcmc_plot(m1, type="intervals", regex_pars = "b", prob_outer = 0.95)+
  theme_bw()+
  scale_x_break(c(0.4, 6.8),scales = 0.5)+ 
  geom_vline(xintercept=0, lty=2, size=1.2)+
  theme(axis.text.y=element_text(size = 16),
        axis.text.x=element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y=element_blank())

bayes_R2(m1) # % of varianza explicada por el modelo


############ ANALiSIS DE RESIDUOS  ##################
# distribucion predictiva posterior de m1
post.dist.m1=posterior_predict(m1, ndraws=1000) 
head(post.dist.m1,2) # 53 columnas orque hay 53 datos; 1000 filas

# Randomized quantile residuals 
qres.m1=createDHARMa(simulatedResponse = t(post.dist.m1), 
                          observedResponse = DF.s$Totalcover, 
                          fittedPredictedResponse = apply(post.dist.m1, 2, median), integerResponse = T) # calculates uniform residuals 
res.m1=data.frame(res=qnorm(residuals(qres.m1))) # convierte RQR a Normal

# Crea un DF con mediana de residuales, avg fitted. Pareto'k of LOO-CV, and vars explicativas 
res.m1=cbind(res.m1, 
             DF.s[,c("Gr.pressure","patch.size","inter.patch.dist","Prod")],
             fitted=fitted(m1, ndraws=1000)[,1], # average of fitted values
             pareto=loo(m1, pointwise=T)$diagnostics$pareto_k) #LOO_CV Pareto k

# Los graficos del analisis de residuos
qq.m1=ggplot(data=res.m1, mapping=aes(sample=res)) + 
  stat_qq_point(aes(col=Prod), size=2)+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position=c(0.2,0.7))
# Residuales vs fitted
res.fit.m1= ggplot(data=res.m1, aes(x=fitted,y=res))+ 
  geom_point(aes(col=Prod), size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position = "none")
# Residuals vs explanatory variables
res.Gr.pressure.m1=ggplot(data=res.m1, aes(x=Gr.pressure,y=res))+ 
  geom_point(aes(col=Prod), size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="Grazing pressure (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
res.patch.size.m1=ggplot(data=res.m1, aes(x=patch.size,y=res))+ 
  geom_point(aes(col=Prod), size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="Patch size (std)",y="Residuals")+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
res.inter.patch.dist.m1=ggplot(data=res.m1, aes(x=inter.patch.dist,y=res))+ 
  geom_point(aes(col=Prod), size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="Inter patch distance (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18),
        legend.position = "none")
Pareto.m1=ggplot(data=res.m1, aes(x=1:nrow(res.m1),y=pareto))+ 
  geom_point(aes(col=Prod),size=2)+ 
  theme_bw()+
  labs(y="Pareto's k")+
  scale_x_continuous(breaks=seq(from=1, to=53, by=5))+
geom_hline(yintercept =0.5,linetype = 2)+ 
  geom_hline(yintercept =0.7,linetype = 2)+
  geom_hline(yintercept =1,linetype = 2)+
  theme(axis.title.x=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18),
        legend.position = "none")
plot_grid(res.fit.m1,res.Gr.pressure.m1,qq.m1, 
          res.patch.size.m1,res.inter.patch.dist.m1,Pareto.m1,
          ncol=3,labels = LETTERS[1:8],
          align="hv",label_x=0.85, label_y=0.95) 

# pp checks for density, mean per level of Prod and max Total plat cover
ppc.density=ppc_dens_overlay(y=DF.s$Totalcover,yrep=post.dist.m1, trim = F, size = 0.5, alpha = 1)+
              xlab("Total Plant cover")+ 
              theme(axis.text = element_text(size=16),
                    axis.title=element_text(size=18),
                    legend.position = "none")
ppc.mean=ppc_stat(y=DF.s$Totalcover,yrep=post.dist.m1,stat = "mean", binwidth = 0.001)+  
          xlab("Average Total Plant cover")+ 
          theme(axis.text = element_text(size=16),
                axis.title=element_text(size=18),
                legend.position = "none")
ppc.interv=ppc_intervals(y=DF.s$Totalcover,yrep=post.dist.m1)+  
  xlab("data point")+ 
  scale_x_continuous(breaks=seq(from=1, to=nrow(DF.s), by=5))+
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
ppc_max_Prod=ppc_stat_grouped(y=DF.s$Totalcover,yrep=post.dist.m1, 
                  group = DF.s$Prod,stat = "max",binwidth=0.005)+
            xlab("Maximum Total Plant cover")+ 
            theme(axis.text = element_text(size=14,angle=90, hjust=1),
                  axis.title=element_text(size=18),
                  strip.text=element_text(size=16),
                  legend.position = "none")
plot_grid(ppc.density, ppc.interv, ppc.mean,ppc_max_Prod,  ncol=2,
          labels = LETTERS[1:4],align="hv",label_x=0.95, label_y=0.95) 

# Condicionales  para cada var explicativa/interaccion:
m1.cond.eff=conditional_effects(m1) 
names(m1.cond.eff) # names of the conditional effects fitted

# plots of Conditional effects
m1.Grazing.Prod=plot(m1.cond.eff, plot = F, points=T)[[5]]+
  theme_bw()+
  labs(y="Total Plant cover",x="Grazing pressure (std)")+ 
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position = "none")
m1.Patchsize.Prod=plot(m1.cond.eff, plot = F, points=T)[[6]]+
  theme_bw()+
  labs(y="Total Plant cover",x="Patch.size (std)")+ 
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank(),
        legend.position = c(0.20,0.85))
m1.inter.patch.dist=plot(m1.cond.eff, plot = F, points=T)[[4]]+
  theme_bw()+
  labs(y="Total Plant cover",x="Inter patch distance (std)")+ 
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =element_blank(),
        legend.position = "none")
plot_grid(m1.Grazing.Prod, m1.Patchsize.Prod,m1.inter.patch.dist,
          ncol=3,labels = LETTERS[1:3],align="hv",label_x=0.85, label_y=0.95) 

########################################################################################

####   Conteos con exceso de ceros  #########
# Importa datos
DF1=read.csv("warblers Pr03.csv",header=T,stringsAsFactors=T)
summary(DF1)
DF1$Year=as.factor(DF1$Year)

table(DF1$NumFledged) # Distr de frecuencias de los valores de var respuesta
table(DF1$NumFledged==0)/length(DF1$NumFledged) # proporcion de ceros

# Dist Probabilidad de la var de respuesta
poiss=fitdist(DF1$NumFledged,"pois")
negbin=fitdist(DF1$NumFledged,"nbinom")
ZAP = fitdist(DF1$NumFledged, "ZAP",discrete=T,  
              start = list(mu=mean(DF1$NumFledged), 
                           sigma=mean(DF1$NumFledged == 0)))
ZANBI = fitdist(DF1$NumFledged, "ZANBI", 
               start = list(mu=mean(DF1$NumFledged),
                            sigma=(mean(DF1$NumFledged)^2) /
                                  (var(DF1$NumFledged)-mean(DF1$NumFledged)),
                            nu=mean(DF1$NumFledged == 0)),
               method = "mge",optim.method = "L-BFGS-B", lower = c(0.01, Inf), upper = c(0, Inf))
CDF.ZI=cdfcomp(list(poiss,negbin,ZAP,ZANBI),
               addlegend=T,main="",legendtext=c("Poisson","NegBin","ZAPoisson","ZANegBin"),
               plotstyle = "ggplot")+
  xlab("Number of Fledglings")+
  geom_line(size=0.8)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position = c(0.75,0.25),
        legend.text=element_text(size=14))
QQ.ZI=qqcomp(list(poiss,negbin,ZAP,ZANBI),addlegend=F,main="",legendtext=c("Poisson","NegBin","ZAPoisson","ZANegBin"),
             plotstyle = "ggplot")+
  theme_bw()+
  geom_jitter(size=2, height=0.2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        title=element_blank(),
        legend.position = c(0.75,0.25),
        legend.text=element_text(size=14))
plot_grid(CDF.ZI, QQ.ZI, ncol=2)

# Analisis exploratorio de datos 
ggpairs(DF1[,c("BreedingDensity","Precip","NumFledged")],
        lower = list(continuous=wrap("points", position=position_jitter(height=0.2, width=0.2)))) +
  theme_bw()+ 
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size=14))
# Estandardizando las vars explicativas
DF1s=data.frame(scale(DF1[,c("BreedingDensity","Precip")],center=T,scale=T),NumFledged=DF1$NumFledged)

# Distr previas para los modelos
formula.m3=bf(NumFledged~BreedingDensity+Precip)+
              lf(hu~1)
formula.m4=bf(NumFledged~BreedingDensity+Precip)+
              lf(hu~BreedingDensity+Precip)
get_prior(formula=formula.m3,data=DF1s,family='hurdle_poisson')
get_prior(formula=formula.m4,data=DF1s,family='hurdle_poisson')

# grafico de la distr previa por defecto de Pr(Y=0)
ggplot(data.frame(x=c(-5, 5)), aes(x))+
  stat_function(fun=dlogis,n=1000, args=list(location=0, scale=1), col="black", size=1.1)+
  stat_function(fun=dnorm,n=1000, args=list(mean=0, sd=1), col="blue", size=1.1)+
  theme_bw()+
  scale_x_continuous(breaks=-5:5)+
  labs (x="LOGIT(Prob (Y=0))", y="Dens. Prob.")+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16))# Prob (Y=0)          

prior.m3 = c(set_prior("normal(0,0.5)",class = "b",coef = "Precip"))
prior.m4 = c(set_prior("normal(0,0.5)",class = "b",coef = "Precip"),
             set_prior("normal(0,0.45)",class = "b",coef = "Precip",dpar="hu" ))

# Modelos estadísticos a ajustar:
m3=brm(formula=formula.m3,data=DF1s,family='hurdle_poisson',prior=prior.m3,
       warmup = 1000, chains=3, iter=2000, thin=3)
m4=brm(formula=formula.m4,data=DF1s,family='hurdle_poisson',prior=prior.m4,
       warmup = 1000, chains=3, iter=2000, thin=3)

# Comparacion de los modelos ajustados
m3=add_criterion(m3, criterio="loo")
m4=add_criterion(m4, criterio="loo")
loo_compare(m3,m4)

# Resultados del modelo seleccionado
summary(m4)

# Convergencia de las cadenas
# Dist posterior de cada parametro para cada cadena de m3.brms 
mcmc_dens_overlay(m4, regex_pars = c("^b"))+
  geom_density(lwd=1.2, alpha=0.1)+
  theme_bw()+
  ylab("Dens. Probabilidad")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position="none")

# Trace plots: degree of mixing of the three chains
mcmc_trace(m4, regex_pars = c("^b"), size=0.3)+
  theme_bw()+
  ylab("Valor del Parametro")+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),              
        strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        legend.position="none") 

# Autocorrelacion de los parametros estimados por cadena: 
mcmc_acf(m4, regex_pars = c("^b"))+
  geom_line(size=1)+
  scale_x_continuous(limits=c(1,20))+
  scale_y_continuous(limits=c(-0.25,0.25))+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),             
        strip.background = element_blank(), 
        strip.text = element_text(size = 12)) 


# distribuciones posteriores del modelo m4 
post.m4=as_draws_df(m4,variable = "^b_", regex = T)

# Graficos de distr posteriores
mcmc_plot(m4, regex_pars = "^b", type="intervals", prob_outer = 0.95)+ 
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16)) 

bayes_R2(m4) # % de varianza explicada


############ ANALISIS de residuos de ZA-Poisson (m4)  ##################
# 1) distribucion predictiva posterior
dist.pred.post.m4=predict(m4, ndraws=1e3, summary=F) 
# 2) simlación de los RQR
qres.m4=createDHARMa(simulatedResponse = t(dist.pred.post.m4), 
                     observedResponse = DF1s$NumFledged, 
                     fittedPredictedResponse = apply(dist.pred.post.m4, 2, median), integerResponse = T) # calculates uniform residuals 
#3) Conviert los RQR a dist Normal
res.m4=data.frame(res=qnorm(residuals(qres.m4))) 
#4) POne todo en un DF para haer los graficos
res.m4=cbind(res.m4, 
             DF1s[,c("BreedingDensity","Precip")],
             fitted=fitted(m4, ndraws=1000)[,1], # average of fitted values
             pareto=loo(m4, pointwise=T)$diagnostics$pareto_k) #LOO_CV Pareto k
# Los graficos
qq.m4=ggplot(data=res.m4, mapping=aes(sample = res)) + 
  stat_qq_point()+
  theme_bw()+
  stat_qq_line()+
  stat_qq_band(alpha=0.3)+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuals vs fitted
res.fit.m4= ggplot(data=res.m4, aes(x=fitted,y=res))+ 
  geom_jitter(col="black", size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))
# Residuals vs explanatory variables
res.Breed.m4=ggplot(data=res.m4, aes(x=BreedingDensity,y=res))+ 
  geom_jitter(col="black", size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="BreedingDensity (std)",y="Residuals")+
  theme(axis.text = element_text(size=16),
        axis.title.x=element_text(size=18))
res.Precip.m4=ggplot(data=res.m4, aes(x=Precip,y=res))+ 
  geom_jitter(col="black", size=1)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_bw()+
  labs(x="Precip (std)",y="Residuals")+
  theme(axis.title.y=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.x=element_text(size=18))
Pareto.m4=ggplot(data=res.m4, aes(x=1:nrow(res.m4),y=pareto))+ 
  geom_jitter(size=1, height = 0.05, col="black")+ 
  theme_bw()+
  labs(y="Pareto's k")+
  geom_hline(yintercept =0.5,linetype = 2)+ 
  geom_hline(yintercept =0.7,linetype = 2)+
  geom_hline(yintercept =1,linetype = 2)+
  theme(axis.title.x=element_blank(), 
        axis.text = element_text(size=16),
        axis.title.y=element_text(size=18))
plot_grid(res.fit.m4,res.Precip.m4,res.Breed.m4, 
          qq.m4,NULL, Pareto.m4,
          ncol=2,labels = c(LETTERS[1:4],"",LETTERS[5]),
          align="hv",label_x=0.85, label_y=0.95) 

# Uso de los pp-checks a partir de la dist predictiva posterior
cero=function (x){sum(x==0)/length(x)} # proporcion de ceros en un vector
ppc.density.m4=ppc_dens_overlay(y=DF1s$NumFledged,yrep=dist.pred.post.m4, trim = F, size = 0.5, alpha = 1)+
  xlab("NumFledged")+ 
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
ppc.mean.m4=ppc_stat(y=DF1s$NumFledged,yrep=dist.pred.post.m4,stat = "mean", binwidth = 0.05)+  
  xlab("media NumFledged")+ 
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
ppc_ceros.m4=ppc_stat(y=DF1s$NumFledged,yrep=dist.pred.post.m4,stat = cero)+
  xlab("Prob (NumFledged=0)")+ 
  theme(axis.text = element_text(size=14,angle=90, hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=16),
        legend.position = "none")
ppc_media_sd.m4=ppc_stat_2d(y=DF1s$NumFledged,yrep=dist.pred.post.m4, stat = c("mean", "sd"))+
  labs(x="Media (NumFledged)",y="sd (NumFledged)")+ 
  theme(axis.text = element_text(size=16),
        axis.title=element_text(size=18),
        legend.position = "none")
plot_grid(ppc.density.m4, ppc.mean.m4,ppc_ceros.m4,ppc_media_sd.m4,  ncol=2,
          labels = LETTERS[1:4],align="hv",label_x=0.35, label_y=0.95) 

# Curvas condicionales predichas por el modelo

# Condicionales  para cada var explicativa/interaccion:
m4.cond.eff=conditional_effects(m4) 
names(m4.cond.eff) # nombre de los efectos conditionales 

m4.Breed=plot(m4.cond.eff, plot = F)[[1]]+
  theme_bw()+
  labs(y="NumFledged",x="BreedingDensity (std)")+ 
  theme(axis.text = element_text(size=16),
        axis.title =  element_text(size=18),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position = "none")
m4.Precip=plot(m4.cond.eff, plot = F)[[2]]+
  theme_bw()+
  labs(y="NumFledged",x="Precip (std)")+ 
  theme(axis.text = element_text(size=16),
        axis.title.x =  element_text(size=18),
        axis.title.y =  element_blank(),
        legend.position = c(0.20,0.85))
plot_grid(m4.Breed,m4.Precip,ncol=2) 
