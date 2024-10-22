#############################333
# Modified code by M. Schuckers
#   based, in part, on code by CJ Turtoro
# this code generages simulated data for
# player aging and fits models proposed in 
# Schuckers, Lopez and Macdonald
#
#


library(tidyverse);
#library(zoo);
library(mgcv);
#library(ggpubr);
library(modelr);
library(splines)
#library(gam);
#library(caret)
library(dplyr)
library(fitdistrplus)
library(truncnorm)
library(stringr)
library(BART)
library(lme4)
#library(optimx)

knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
basefolder=""
year_length_days=365


# Range of Decay is from x to 1
# y is the "true" baseline aging curve

quadform = function(x){
  -(1/7)*(x-26)^2
}

ageform2=function(x,a0=0,a=-1/7,b=0,agemax=26,age2=30)
{
 a*(x-agemax)^2+-1*b*as.numeric(x>age2)*(x-age2)^2
}
ageform3=function(x,a0=0,a=-1/7,b=0,agemax=26,age2=30)
{
  a*(x-agemax)^2+-1*b*as.numeric(x>age2)*(x-age2)^2+a0*(x-age2)*as.numeric(x>age2)
}
ageform4=function(x,a0=0,a=-0.1,b=0.8,agemax=25,age2=30)
{
  a*(x-agemax)^2+b*(x-age2)*as.numeric(x>age2)
}

ageform5=function(x,a=-1/9,b=-0.006,c=0.0045,agemax=25)
{
  a*(x-agemax)^2+
    b*(x-agemax)^2*as.numeric(x>agemax)+
    c*(x-agemax)^3*as.numeric(x>agemax)        
}

######################
# full_data_gen is function to generate data 
#     for players 
#     number of players will be 26*n_players
minage0=18
maxage0=40
#y1=ageform4(minage:maxage,a=-0.12,b=0.8,agemax=25, age2=30)
y1=ageform5(x=minage0:maxage0,agemax=25)

y2=c(-5:0,-1:-10)
numb.letters=10

full_data_gen = function(y=y1,n_players=np,minage=18,
                         sd.players=3,sd.noise=1,sd.quad=0.02,
                         maxage=40,overall_mean=2,missing_param=0){
  nyears=length(y)
  agemax=25
  player_data<- data.frame(
      p1 = c(rep(letters[1:numb.letters],each=nyears*n_players)),
      p2 = c(rep(rep(1:n_players,each=nyears),numb.letters)), ## Create Players
      age = rep(minage:maxage,n_players*numb.letters), # Ages 
      traj = rep(y,numb.letters*n_players), # Baseline Aging
      r_goodness = c(rep(rnorm(n_players*numb.letters,overall_mean,sd.players),
                     each=nyears))) %>% # Random "Quality"
      mutate(r_quad=c(rep(rnorm(n_players*numb.letters,0,sd.quad),
                      each=nyears))*(age-agemax)^2) %>%
  #r_rate = rep(runif(n_players*numb.letters,x,1),each=nyears), # Random "Decay Rate"
  mutate(r_spray = rep(rnorm((n_players*numb.letters)*nyears,
                             0,sd.noise*25/age)) # Random Yearly Noise
  ) %>% 
  mutate(gar = r_goodness+traj+r_quad+r_spray, # Rolls Dice
         #gar = r_goodness+r_spray+traj+r_quad,
         player = paste0(p1,p2))%>%
  dplyr::select(-p1) %>%
  dplyr::select(-p2) %>% # Name Player
  group_by(player) %>%
    
    # "Dropout Mechanism" -- Based on 3-year GAR 
    # (Avg player career = 6.6 seasons)
    mutate(lag1 = lag(gar),lag2 = lag(gar,2)) %>%
    rowwise() %>% 
    mutate(threeyr = mean(c(lag1,lag2,gar),na.rm = T)) %>%
    group_by(player) %>%
    mutate(fcode = ifelse(lag(player)!=player,0,1),
           include = ifelse(threeyr<4,0,1),
           d = gar - lag(gar)) 
  
  
  trugar = data.frame(age=minage:maxage,true_gar = overall_mean+y)
  player_data<-player_data %>%
  left_join(trugar,by="age") 
  
  return(player_data)
}

############################33
# inverse logit function

ilogit<-function(x){
  return(exp(x)/(1+exp(x)))
}


# function to make some of the observations in the df
# missing values, i.e not in the league
# option parameter in the data_trim_miss
# is as follows:
#   0  : no missingness, data is complete
#   1  : missingness is deterministic by player value
#   2  : missingness is probabilistic by player value (gar)
#   3  : missingness is probabilistic by lag function of CJ Torturo
#
#   Note that missing_param has different functions depending on
#       the option value above


data_add_miss = function(player_data,option=0,missing_param=1,minage,maxage)
{
  #option 0, no missingness
  if (option==0) {
    player_data<-player_data %>% mutate(gar_miss=gar)
  }
  #option 1, missingness deterministic below a threshold
  if (option==1)
  {
    player_data <-player_data %>%
      mutate(gar_miss=if_else(gar<missing_param,NA_real_,gar))
  }
  # option 2 logistic prob 
  if(option==2)
  {
    kappa=0.01
    player_data$cumsum=player_data %>% 
      group_by(player) %>% 
      group_map(~cumsum(.$gar)) %>% 
      unlist
    player_data<- player_data %>%
      mutate(miss_prob=1-ilogit(gar-missing_param+0.03*(age-20))) %>%
      mutate(miss_true=rbinom(n=length(gar),size=1,prob=miss_prob)) %>%
      mutate(gar_miss=if_else(as.logical(miss_true),NA_real_,gar))
  }
  
  if(option==3)
  {
    player_data<-player_data %>%
      
    group_by(player) %>%
      # "Dropout Mechanism" -- Based on 3-year GAR 
      # (Avg player career = 6.6 seasons)
      mutate(lag1 = lag(gar),lag2 = lag(gar,2)) %>%
      rowwise() %>% 
      mutate(threeyr = mean(c(lag1,lag2,gar),na.rm = T)) %>%
      group_by(player) %>%
      mutate(fcode = ifelse(lag(player)!=player,0,1),
             gar_miss = ifelse(threeyr<missing_param,NA_real_,gar),
             d = gar - lag(gar)) 
    
  }
  if(option==4)
  {
    player_data$cumsum=player_data %>% 
      group_by(player) %>% 
      group_map(~cumsum(.$gar)) %>% 
      unlist
    perc_in_byage=read.csv("~/PlayerAging/aging_forward_pctbyage.csv") %>%
      dplyr::select(c("pct","age")) %>%
      filter(age %in% minage:maxage)
    player_data$gar_miss=NA
    ages<-minage:maxage
    #kappa=missing_param
    for(i in 1:length(minage:maxage)){
       which.age=which(player_data$age==ages[i])
       #print(which.age)
       #beta.p=rbeta(1,kappa*perc_in_byage$percent[i],kappa*(1-perc_in_byage$percent[i]))
       #print(beta.p)
       n.samp=round(length(which.age)*perc_in_byage$pct[i],0)
       #print(n.samp)
       which.in=sample(which.age,n.samp,replace=FALSE,prob=exp(player_data$cumsum[which.age])/missing_param)
       #print(which.in)
       player_data$gar_miss[which.in]=player_data$gar[which.in]
    }
    
  }
  return(player_data)
}

## need to drop players who do not appear in league, i.e 
## players with all values of include==0

###############################################
# this function eliminates any player who do not appear at 
# all i.e. they do not have any years of observed data
# 
#  this is for model fitting and prediction
#       need to have at least one value of a player to make that work.
sift_players=function(player_data)
{
  players_in<-player_data %>%
    filter(!is.na(gar_miss)) %>%
    group_by(player) %>%
    summarize(count=n())%>%
    filter(count>0) %>%
    dplyr::select(-count)
  
  player_data<-player_data %>%
    filter(player %in% players_in$player)
  
return(player_data)
}




  

  #print(mean_max)

####################################33
#  Original Schuckers method for fitting model for aging
#  - create base model
#  - create lower threshold
#  - for missing values (outside league) impute below threshold
#       following the base model
#  
##############################################################3
missingness_method=function(player_data,perc_use=0.75){
  
  #######################
  # initial model with just observed data
  # 
  model0<-mgcv::gam(gar_miss~s(age,k=6)+player,data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  cat("m_m0")
  
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 
  minage<-min(player_data$age)
  maxage<-max(player_data$age)
  ages<- minage: maxage
  cat(ages)
  min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
  q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}
  z.quant=tapply(player_data$gar,player_data$age,q.2)
  cat("c-1")
  z_df<-data.frame(z.quant,ageyrs=as.numeric(names(z.quant)))
  # smooth the lower percentile curve
  cat("c-2")
  quant_df=data.frame(ageyrs=ages) %>%
    full_join(z_df, by = "ageyrs")
  
  min.reg=loess(z.quant~ageyrs,data=quant_df)
  min.fit=predict(min.reg,ages)
  cat(min.fit)
  cat("\n")
  # since we are making an upper bound if cannot calculate the upper bound then 
  # take the max of the predicted ones
  min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)
  #write_csv(min.fit,"~/PlayerAging/miss_opt4delta/bounds_0_600_0.8.csv")
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  player_data$ubound<-min.fit[player_data$age-minage+1]
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  

  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
  for(i in 1:n.gen)
      {
      #  if(i %% 50 ==0) cat(".")
        z_gen[i]<-etruncnorm( a=-Inf, 
                             b=player_data$ubound[which.gen[i]], 
                             mean = mean0+player_data$eta[which.gen[i]], 
                             sd = sd0)
        
      }
      
  player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  
  modelx<-mgcv::gam(gar_impute~as.factor(player)+s(age,k=6),data=player_data)
  player_data$eta2<-predict(modelx,newdata=player_data)
  player_data$gar_impute2=player_data$gar
  if(n.gen>0){
            for(i in 1:n.gen)
            
            {
              #  if(i %% 50 ==0) cat(".")
              z_gen[i]<-etruncnorm( a=-Inf, 
                                    b=player_data$ubound[which.gen[i]], 
                                    mean = mean0+player_data$eta2[which.gen[i]], 
                                    sd = sd0)
              
            }
            
            player_data$gar_impute2[which.gen]<-z_gen
            }
  k2=6 # of knots
  model.out<-mgcv::gam(gar_impute2~player-1+s(age,k=k2),data=player_data)
  tempname=str_replace_all(names(model.out$coefficients[1]),"player","")
  #player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  
  #pred.out<-predict(model.out,type="terms", exclude="player",newdata=player_1)
  cat("c0")
  pred.out<-data.frame(pred=predict(model.out,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  #print(str(player_data))
  #write.csv(player_data,"~/PlayerAging/miss_opt4delta/miss1_trunc_0_600_0.8_single.csv")
  cat("c1")
  return(pred.out)
  
}  




  
####################################
# function for estimation of underlying aging curve that 
#     uses percentiles of the observed data at a given year
#     and the total number of players to estimate
#     percentiles of a given year

miss_quantile_method=function(player_data, perc_use=0.2,perc_estimand=0.5)
{
  minage=min(player_data$age)
  maxage=max(player_data$age)
  library(truncnorm)
  # calculate the total number of players in the dataset
  n.players<-length(table(player_data$player))+4
  player_data <- player_data %>%
    mutate(miss_true=as.numeric(is.na(gar_miss)))
  
  #  calculate the total number of players observed in a given year
  nperyear<- tapply(1-player_data$miss_true,player_data$age,sum)+2
  
  # get the observed ages 
  agevalues<-as.numeric(names(nperyear))

  # calculate perc_use percentile
  quant.na=function(x){return(quantile(x,perc_use,na.rm="TRUE"))}
  
  # calculate perc_use percentile by year
  quant_by_year<-tapply(player_data$gar_miss,player_data$age,quant.na)
  
  # calcualted adjusted quantile for all players by year
  quant_by_year_all <- 1-(1-perc_use)*nperyear/n.players
 mintemp<-cbind(1e-3,quant_by_year_all)
 maxtemp<-cbind(apply(mintemp,1,max),1-1e-3)
 quant_out<-apply(maxtemp,1,min)
 
   
  
  #fit simple model to get estimated sd of residuals
  model0<-mgcv::gam(gar_miss~bs(age,knots=seq((minage-2),(maxage+2),length=7))+player,
                    data=player_data)
  sd0<-sd(model0$residuals)
  
  # adjust sd by amount of truncation in the distribution
  #      of residuals
  sd1<-sd0*1/sqrt(vtruncnorm(qnorm(quant_out),Inf,0,1))
  
  # get adjusted value for quantile we want to estimate
  ### XXX ###
  out.estimated<- quant_by_year+
    (qnorm(perc_estimand)-qnorm(quant_out))*sd1
  age.estimated=as.numeric(names(out.estimated))
    
  knots0=seq((minage-1),(maxage+1),length.out=5)
  mean.smooth=mgcv::gam(out.estimated~bs(age.estimated),
                        weights=nperyear+0.1)
  #mean.smooth=loess(mean.post~ages,weights=nperyear)
  #print(predict(mean.smooth))
  #print(ages)
  qe=data.frame(age=minage:maxage,
                quant.predicted=predict(mean.smooth,newdata=data.frame(age.estimated=minage:maxage)))
  #print(qe)
  #mean.qe$q2=mean.qe$q2-max(mean.qe$q2)
  return(qe)
}  

####################################33
#  Second  method for fitting model for aging
#  - create base model via miss_quantile_method
#  - for missing values (outside league) impute 
#       following the miss_quantile_method + player effects
#  
##############################################################3
missingness_method2=function(player_data,perc_use=0.2){
  
  #######################
  # initial model with just observed data
  # 
  model0<-mgcv::gam(gar_miss~s(age,k=6)+player,data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)

  minage=min(player_data$age)
  maxage=max(player_data$age)
# get estimates of average aging curve from miss_quantile_method    
quant.out<-  miss_quantile_method(player_data,perc_use=perc_use,perc_estimand=0.5)

# join average aging curve 
player_data<-player_data %>%
  right_join(quant.out,by=c("age"))
#print(table(player_data$player))
player_data2<-player_data %>%
  mutate(minus_aging=gar_miss-quant.predicted) 
#print(player_data2[player_data2$player=="b19",])
#print(dim(player_data))
#cat("d1")
# fit curve to data with aging removed to get player factors
player_effects_model<-mgcv::gam(minus_aging~player,data=player_data2)
#cat("d1.5")
#print(predict(player_effects_model)[1:50])
#print(length(player_data$quant.predicted))
player_data$nu=predict(player_effects_model,newdata=player_data)+player_data$quant.predicted
#cat("d2")


ages<- minage: maxage

min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}

z.quant=tapply(player_data$gar,player_data$age,q.2)

# smooth the lower percentile curve
quant_df=data.frame(ageyrs=ages,z.quant)
min.reg=loess(z.quant~ageyrs,data=quant_df)
min.fit=predict(min.reg,ages)

# since we are making an upper bound if cannot calculate the upper bound then 
# take the max of the predicted ones
min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)

# create variable for imputed_gar values
#
player_data$gar_impute_nu<-player_data$gar
player_data$ubound<-min.fit[player_data$age-minage+1]
which.gen<-which(is.na(player_data$gar_miss))
n.gen<-length(which.gen)

library(truncnorm)
if(n.gen>0){
  z_gen<-rep(NA,n.gen)
  
        for(i in 1:n.gen)
        {
          #  if(i %% 50 ==0) cat(".")
          z_gen[i]<-etruncnorm( a=-Inf, 
                                b=player_data$ubound[which.gen[i]], 
                                mean = mean0+player_data$nu[which.gen[i]], 
                                sd = sd0)
          
        }
        
        player_data$gar_impute_nu[which.gen]<-z_gen
      }
##############################
# fit a new curve to get the mean average aging based upon 
#    these imputed values
#
k2=6
model.out<-mgcv::gam(gar_impute_nu~player+s(age,k=k2),data=player_data)
#tempname=str_replace_all(names(model.out$coefficients[2]),"player","")
#player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
cat("d3")
#pred.out<-predict(model.out,type="terms", exclude="player",newdata=player_1)
pred.out<-data.frame(pred=predict(model.out,newdata=player_data),age=player_data$age) %>%
  group_by(age) %>%
  summarise(out=mean(pred))


return(pred.out)
}


missingness_method3=function(player_data,perc_use=0.2){
  
  #######################
  # initial model with just observed data
  # 
  model0<-mgcv::gam(gar_miss~s(age,k=6)+as.factor(player),data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  
  minage=min(player_data$age)
  maxage=max(player_data$age)
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 

  ages<- minage: maxage
  min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
  q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}
  z.quant=tapply(player_data$gar,player_data$age,q.2)
  
  z_df<-data.frame(z.quant,ageyrs=as.numeric(names(z.quant)))
  # smooth the lower percentile curve
  quant_df=data.frame(ageyrs=ages) %>%
    full_join(z_df, by = "ageyrs")
  
  min.reg=loess(z.quant~ageyrs,data=quant_df)
  min.fit=predict(min.reg,ages)
  
  # since we are making an upper bound if cannot calculate the upper bound then 
  # take the max of the predicted ones
  min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)
  
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  player_data$ubound<-min.fit[player_data$age-minage+1]+
    rnorm(length(player_data$age),mean=0,sd=sd0/3)
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  
  
  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
    for(i in 1:n.gen)
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  
  modelx<-mgcv::gam(gar_impute~s(age,k=6)+player,data=player_data)
  player_data$ubound<-min.fit[player_data$age-minage+1]+
    rnorm(length(player_data$age),mean=0,sd=sd0/3)
  player_data$eta2<-predict(modelx,newdata=player_data)
  player_data$gar_impute2=player_data$gar
  
  if(n.gen>0){
    for(i in 1:n.gen)
      
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta2[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute2[which.gen]<-z_gen
  }
  k2=6 # of knots - 3
  model.out<-mgcv::gam(gar_impute2~player+s(age,k=k2),data=player_data)
  tempname=str_replace_all(names(model.out$coefficients[2]),"player","")
  player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  
  #pred.out<-predict(model.out,type="terms", exclude="player",newdata=player_1)
  # type="terms" changed to "iterms" on 1/24/21
  pred.out<-data.frame(pred=predict(model.out,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  return(pred.out)
  
}  


####################################
#  This function runs BART and returns an estimated age curve
#
#

eval_BART <- function(truecurve = y1, input_data=player_data_sifted){
  # Observed data
  player_data_observed <- input_data %>% filter(include == 1)
 
  minage=min(input_data$age)
  maxage=max(input_data$age)
  # treatment is their age
  # x is their character ID
  trt <- player_data_observed$age    
  x <- as.factor(player_data_observed$player)
  
  ## outcome variable
  y <- player_data_observed$gar
  
  ## combining treatment variable with x
  xt = cbind(trt, x)
  
  #### #### #### #### #### 
  #### Impute across all players
  #### #### #### #### #### 
  
  trt_full <- input_data$age    
  #table(trt_full) ## one at each age for each of the players where we observe at least one year
  x_full <- as.factor(input_data$player)
  xt_full <- cbind(trt_full, x_full)
  
  
  ## run BART (wbart for continuous outcomes)
   n.imps = 1000 # number of imputations
  bart_mod = wbart(x.train = xt, y.train = y, 
                   k=2, ntree=100, ndpost=n.imps, 
                   nskip=500, printevery=100L)
  
  ## posterior predictions
  bart_pred = pwbart(xt_full, bart_mod$treedraws)
  
  #dim(bart_pred)  ## nsim rows x number of playerage X player predictions
  
  # Average prediction at each age for each person
  age_est_all = NULL
  for (m in 1:n.imps) {
    # potential outcomes for each age
    # computes average age curve 
    xt_fake <- as.data.frame(xt_full) 
    xt_fake$gar_hat <- bart_pred[m,] ## the unobserved at each m
    xt_fake$gar <- input_data$gar ## not sure this is needed
    xt_fake$is_obs <- input_data$include ## not sure this is needed
    
    age_ave <- suppressMessages(xt_fake %>% 
      group_by(trt_full) %>% 
      summarise(ave_age = mean(gar_hat)) %>% 
      mutate(sim_num = m))
    
    age_est_all <- bind_rows(age_est_all, age_ave)
  }
  
  bart_est <- age_est_all %>% 
    group_by(trt_full) %>% 
    summarise(age_est_overall = mean(ave_age))
  return(bart_est)  
}


missingness_lmer=function(player_data,perc_use=0.2){
  
  #######################
  # initial model with just observed data
  # 
  player_data$age2=player_data$age^2
  model0<-mgcv::gam(gar_miss~s(age,k=6)+player,data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  
  
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 
  minage<-min(player_data$age)
  maxage<-max(player_data$age)
  ages<- minage: maxage
  min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
  q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}
  z.quant=tapply(player_data$gar,player_data$age,q.2)
  
  z_df<-data.frame(z.quant,ageyrs=as.numeric(names(z.quant)))
  # smooth the lower percentile curve
  quant_df=data.frame(ageyrs=ages) %>%
    full_join(z_df, by = "ageyrs")
  
  min.reg=loess(z.quant~ageyrs,data=quant_df)
  min.fit=predict(min.reg,ages)
  
  # since we are making an upper bound if cannot calculate the upper bound then 
  # take the max of the predicted ones
  min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)
  
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  player_data$ubound<-min.fit[player_data$age-minage+1]
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  
  
  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
    for(i in 1:n.gen)
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  
  modelx<-mgcv::gam(gar_impute~as.factor(player)+s(age,k=6),data=player_data)
  player_data$eta2<-predict(modelx,newdata=player_data)
  player_data$gar_impute2=player_data$gar
  if(n.gen>0){
    for(i in 1:n.gen)
      
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta2[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute2[which.gen]<-z_gen
  }
  
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  m.imp.poly5  = lmer(gar_impute2  ~  1+ns(age,df=6) + (1+age|player)+(1+age2| player), 
                      REML = FALSE,data=player_data)
  pred.out<-data.frame(pred=predict(m.imp.poly5,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  
  return(pred.out)
}

missingness_lmer2=function(player_data,perc_use=0.2){
  
  #######################
  # initial model with just observed data
  # 
  model0<-mgcv::gam(gar_miss~s(age,k=6)+player,data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  
  
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 
  minage<-min(player_data$age)
  maxage<-max(player_data$age)
  ages<- minage: maxage
  min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
  q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}
  z.quant=tapply(player_data$gar,player_data$age,q.2)
  
  z_df<-data.frame(z.quant,ageyrs=as.numeric(names(z.quant)))
  # smooth the lower percentile curve
  quant_df=data.frame(ageyrs=ages) %>%
    full_join(z_df, by = "ageyrs")
  
  min.reg=loess(z.quant~ageyrs,data=quant_df)
  min.fit=predict(min.reg,ages)
  
  # since we are making an upper bound if cannot calculate the upper bound then 
  # take the max of the predicted ones
  min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)
  
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  player_data$ubound<-min.fit[player_data$age-minage+1]
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  
  
  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
    for(i in 1:n.gen)
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  
  modelx<-mgcv::gam(gar_impute~as.factor(player)+s(age,k=6),data=player_data)
  player_data$eta2<-predict(modelx,newdata=player_data)
  player_data$gar_impute2=player_data$gar
  if(n.gen>0){
    for(i in 1:n.gen)
      
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta2[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute2[which.gen]<-z_gen
  }
  
  m.imp.ns2    = lmer(gar_impute2  ~ 1 + ns(age,df=6) + (1+ns(age, df=6)      ||player), data=player_data)
 
  pred.out<-data.frame(pred=predict(m.imp.ns2,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  return(pred.out)
  
}


missingness_notrunc=function(player_data){
  
  #######################
  # initial model with just observed data
  # 
  model0<-mgcv::gam(gar_miss~s(age,k=6)+as.factor(player),data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  
  minage=min(player_data$age)
  maxage=max(player_data$age)
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 
  ages<- minage: maxage
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  #player_data$ubound<-min.fit[player_data$age-minage+1]
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  
  
  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
    for(i in 1:n.gen)
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-rnorm(1,mean = mean0+player_data$eta[which.gen[i]], 
                       sd = sd0)
      
    }
    
    player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
 k2=6
   model.out<-mgcv::gam(gar_impute~player+s(age,k=k2),data=player_data)
  tempname=str_replace_all(names(model.out$coefficients[2]),"player","")
  player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  pred.out<-data.frame(pred=predict(model.out,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  #write_csv(player_data,"miss_no_trunc_single.csv")
  #pred.out<-predict(model.out,type="terms", exclude="player",newdata=player_1)
  print(pred.out)
  #write_csv(data.frame(player_data),"miss_notrunc_0_600_0.8_single.csv")
  
  #write_csv(player_data,"~/PlayerAging/miss_opt4delta/miss1_trunc_0_600_0.8_single.csv")
  return(pred.out)
  
}

missingness_quad=function(player_data,perc_use=0.75){
  
  #######################
  # initial model with just observed data
  # 
  model0<-lm(gar_miss~age+I(age^2)+player,data=player_data)
  mean0<-mean(model0$residuals)
  sd0<-sd(model0$residuals)
  
  
  ######
  # Create predicted values for all observations
  #         use these that  
  player_data$eta<-predict(model0,newdata=player_data)
  
  # Next set of code builds a lower (permeable) threshold for 
  #     data being observed, something like replacement level
  #     but 
  minage<-min(player_data$age)
  maxage<-max(player_data$age)
  ages<- minage: maxage
  min.na.rm.true=function(x){return(min(x,na.rm=TRUE))}
  q.2=function(x){return(quantile(x,perc_use,na.rm=TRUE))}
  z.quant=tapply(player_data$gar,player_data$age,q.2)
  
  z_df<-data.frame(z.quant,ageyrs=as.numeric(names(z.quant)))
  # smooth the lower percentile curve
  quant_df=data.frame(ageyrs=ages) %>%
    full_join(z_df, by = "ageyrs")
  
  min.reg=loess(z.quant~ageyrs,data=quant_df)
  min.fit=predict(min.reg,ages)
  
  # since we are making an upper bound if cannot calculate the upper bound then 
  # take the max of the predicted ones
  min.fit[is.na(min.fit)]=max(min.fit,na.rm=TRUE)
  
  # create variable for imputed_gar values
  #
  player_data$gar_impute<-player_data$gar
  player_data$ubound<-min.fit[player_data$age-minage+1]
  which.gen<-which(is.na(player_data$gar_miss))
  n.gen<-length(which.gen)
  
  
  library(truncnorm)
  if(n.gen>0)
  {
    z_gen<-rep(NA,n.gen)  
    for(i in 1:n.gen)
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute[which.gen]<-z_gen
  }
  ##############################
  # impute a second time using the imputed values
  # to get a different curve for the mean per age
  # 
  
  modelx<-lm(gar_impute~age+I(age^2)+player,data=player_data)
  player_data$eta2<-predict(modelx,newdata=player_data)
  player_data$gar_impute2=player_data$gar
  if(n.gen>0){
    for(i in 1:n.gen)
      
    {
      #  if(i %% 50 ==0) cat(".")
      z_gen[i]<-etruncnorm( a=-Inf, 
                            b=player_data$ubound[which.gen[i]], 
                            mean = mean0+player_data$eta2[which.gen[i]], 
                            sd = sd0)
      
    }
    
    player_data$gar_impute2[which.gen]<-z_gen
  }
  k2=6 # of knots
  model.out<-lm(gar_impute2~player+age+I(age^2),data=player_data)
  tempname=str_replace_all(names(model.out$coefficients[1]),"player","")
  player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  
  #pred.out<-predict(model.out,type="terms", exclude="player",newdata=player_1)
  pred.out<-data.frame(pred=predict(model.out,newdata=player_data),age=player_data$age) %>%
    group_by(age) %>%
    summarise(out=mean(pred))
  return(pred.out)
  
}  


###################################333
#   this function evaluates the average player aging function and 
#      and reports rmse from original function as well as the 
#      percent of times that 


eval_model_meancurve=function(truecurve=y1,player_data_sifted,player_data)
{
  minage=min(player_data_sifted$age)
  maxage=max(player_data_sifted$age)
  #print(truecurve)
  ##########################
  # goal of this function is to evaluate the different methods
  #    for player aging
  #    delta method, spline method, loess method
  #     and 2 proposed by Schuckers et al 2020
  #     missingness_method and missingness_method2
  #     straight quanile estimation
  #     delta_plus
  #numb_methods=9
 
  
  #######################3
  # empty vectors for returning evaluation results
  
  ## ML -- why is this not include == 1?
  ## 
  player_sub <- player_data %>%
    filter(!is.na(gar_miss))
  
  ##############################
  # code for calculating delta methods by CJ Turtoro
  # This is method #1
  
  gar_means<-player_sub %>%
    group_by(age) %>%
    summarise(mean_gar = mean(gar_miss))
  
  delta =  player_data %>%
   mutate(d=ifelse(is.na(lag(gar_miss)),NA,d)) %>%
      filter(!is.na(gar_miss)) %>%
    group_by(age) %>%
    summarise(change = mean(d,na.rm=TRUE)) %>%
    mutate(change = ifelse(is.na(change),0,change),
           chained = cumsum(change),
           delta_gar = chained - max(chained)) %>%
    dplyr::select(age,delta_gar)
  
  
   #print(delta$delta_gar)
   #print(age_max)
  out_df=data.frame(method="delta",t(delta$delta_gar))

      #out_max01[1]=sum(as.numeric(which.max(delta$delta_gar)==age_max))
    #out_rmse[1]=sqrt(sum((delta$delta_gar-truecurve)^2))
    #out_name[1]="delta"
  
    
  ##################################3
  #
  # create a full data frame for player 'a1' so make predictions from
  # 
  #pnames=  names(table(player_sub$player))
  #player_a1<-data.frame(age=minage:maxage,player=rep(pnames[1],length(y1)))

  cat("a")
  ###################################  
  #get predicted average aging curve for spline model
  # This is method 2
    k2=6 # of knots in the spline
  spline_mod = mgcv::gam(gar_miss ~ player  + s(age,k=k2),data=player_sub)
  n2=length(spline_mod$coefficients)-(k2-1)
  mean.player.spline=mean(spline_mod$coefficients[1:n2])
  tempname=str_replace_all(names(spline_mod$coefficients[2]),"player","")
  player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  
  spline_pred<-data.frame(pred=predict(spline_mod,newdata=player_data_sifted),age=player_data_sifted$age) %>% 
    group_by(age) %>%
    summarise(out=mean(pred))
  
  #out_max01[2]=sum(as.numeric(which.max(spline_pred)==age_max))
  #out_rmse[2]=sqrt(sum((spline_pred-truecurve)^2))
  #out_name[2]="spline_nomiss"
  out_df<-out_df %>%
    bind_rows(data.frame(method="spline_nomiss",t(spline_pred)))
 cat("b")
 ###################################  
 #get predicted average aging curve for spline model
 # This is method 2
 k2=6 # of knots in the spline
 spline_mod2 = mgcv::gam(gar_miss ~ s(age,k=k2),data=player_sub)
 #n2=length(spline_mod$coefficients)-(k2-1)
 #mean.player.spline=mean(spline_mod$coefficients[1:n2])
 #tempname=str_replace_all(names(spline_mod$coefficients[2]),"player","")
 #player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
 
 spline_pred2<-data.frame(pred=predict(spline_mod2,newdata=player_data_sifted),age=player_data_sifted$age) %>% 
   group_by(age) %>%
   summarise(out=mean(pred))
 
 #out_max01[2]=sum(as.numeric(which.max(spline_pred)==age_max))
 #out_rmse[2]=sqrt(sum((spline_pred-truecurve)^2))
 #out_name[2]="spline_nomiss"
 out_df<-out_df %>%
   bind_rows(data.frame(method="spline_nomiss_noplayer",t(spline_pred2)))
 cat("c")
   ######################################
  # get predicted average aging curve for loess model
  # This is method 3
  #loess_mod = mgcv::gam(gar_miss ~ player -1 + lo(age,span = 0.3,degree = 2),data=player_sub)
  #n2=length(loess_mod$coefficients)-1
  #mean.player.loess=mean(loess_mod$coefficients[1:n2])
  #tempname=str_replace_all(names(loess_mod$coefficients[2]),"player","")
  #player_1<-data.frame(age=minage:maxage,player=rep(tempname,length(minage:maxage)))
  #loess_pred<-predict(loess_mod,type="terms", exclude="player",newdata=player_1)
  #out_max01[3]=sum(as.numeric(which.max(loess_pred)==age_max))
  #out_rmse[3]=sqrt(sum((loess_pred-truecurve)^2))
  #out_name[3]="loess_nomiss"
  #out_df<-out_df %>%
  #  bind_rows(data.frame(method="loess_nomiss",t(loess_pred)))
#cat("c")
    #######################################
  #
  #  This is method 4
  # missingness method 1
  
  miss1_pred<-missingness_method(player_data=player_data_sifted,perc_use=0.75)
  #print(pred.miss1)
  #out_max01[4]=sum(as.numeric(which.max(pred.miss1)==age_max))
  #out_rmse[4]=sqrt(sum((pred.miss1-truecurve)^2))
  #out_name[4]="missing1"
 cat("c3")
  out_df<-out_df %>%
    bind_rows(data.frame(method="missing1",t(miss1_pred)))
  cat("d")
  #######################################
  #
  #  This is method 5
  # missingness method 2
  miss2_pred<-missingness_method2(player_data_sifted,perc_use=0.75)
  #out_max01[5]=sum(as.numeric(which.max(pred.miss2)==age_max))
  #out_rmse[5]=sqrt(sum((pred.miss2-truecurve)^2))
  #out_name[5]="missing2"
  out_df<-out_df %>%
    bind_rows(data.frame(method="missing2",t(miss2_pred)))
  cat("e")
  #######################################
  #
  #  This is method 5
  # missingness method 2
  #miss3_pred<-missingness_method3(player_data_sifted)
  #out_max01[5]=sum(as.numeric(which.max(pred.miss2)==age_max))
  #out_rmse[5]=sqrt(sum((pred.miss2-truecurve)^2))
  #out_name[5]="missing2"
  #out_df<-out_df %>%
  #  bind_rows(data.frame(method="missing3",t(miss3_pred)))
  #cat("f")
  #######################################3
  #
  #  This is method 6
  mean.qe=miss_quantile_method(player_data,perc_use=0.75)
  mean.qe_pred=mean.qe$quant.predicted
  #out_max01[6]=sum(as.numeric(which.max(mean.qe$quant.predicted)==age_max))
  #out_rmse[6]=sqrt(sum((mean.qe$quant.predicted-truecurve)^2))
  #out_name[6]="miss_quantile"  
  out_df<-out_df %>%
    bind_rows(data.frame(method="miss_quantile",t(mean.qe_pred)))
  cat("g")
  
  
  ########################################33
  #
  # This is method 7
  delta$delta2=max(gar_means$mean_gar)+delta$delta_gar
  #out_max01[7]=sum(as.numeric(which.max(delta$delta_gar)==age_max))
  #out_rmse[7]=sqrt(sum((delta$delta2-truecurve)^2))
  #out_name[7]="delta plus"
  out_df<-out_df %>%
    bind_rows(data.frame(method="delta_plus",t(delta$delta2)))
  #return(data.frame(out_max01,out_rmse,out_name))
  cat("*\n")
  #cat(dim(player_data_sifted))
  ########################################33
  #
  # This is method 8
  ### Mike to write a BART function that runs here

#  tt<-tapply(!is.na(player_data_sifted$gar_miss),player_data_sifted$age,sum)
#tt<-tt[tt>0]
#
#if (length(tt)==length(minage:maxage)){
#    bart_est = eval_BART(input_data=player_data_sifted)
#  bart_est$method <- "bart"
#  out_df<-out_df %>%
#    bind_rows(data.frame(method="bart",t(bart_est$age_est_overall)))
#}
  cat("*\n")
  ################################################
  #
  # This is method 9
  ### BMac to write function with splines+ Random effects
  #m.imp.poly5  = lmer(gar_impute  ~ 1 + bs(age,bs="cr",k=5) + (1+age|player)+(1+age2| player), data=player_data_sifted)
  # 
  # 
   predmat<-missingness_lmer(player_data_sifted,perc_use=0.75)
   out_df<-out_df %>%
     bind_rows(data.frame(method="lmer_poly5",t(as.vector(predmat))))
  cat("**\n")
  
  predmat<-missingness_lmer2(player_data_sifted,perc_use=0.75)
  out_df<-out_df %>%
    bind_rows(data.frame(method="lmer_m.imp.ns2",t(as.vector(predmat))))
  cat("**\n")
  
  ################################################
  #
  # This is method 10?
  # 
  pred.notrunc<-missingness_notrunc(player_data_sifted)
  out_df<-out_df %>%
    bind_rows(data.frame(method="miss_notrunc",t(pred.notrunc)))
  #return(data.frame(out_max01,out_rmse,out_name))
  
 
  
  
  ###############################################3
  pred.quad<-missingness_quad(player_data_sifted,perc_use=0.75)
  out_df<-out_df %>%
    bind_rows(data.frame(method="miss_quad",t(pred.quad)))
 
  names(out_df)=c("method",paste0("age",minage:maxage))   
  out_df<-out_df %>% filter(age18!=18)
  return(out_df)
  }
 
