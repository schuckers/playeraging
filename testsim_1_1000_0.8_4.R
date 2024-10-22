source("NewAgingModelsTesting6.R")
player_aging_sim=function(
#overall mean
mu=1,
# number of players
nn=25,
# number of simulations
nreps=1,
# standard deviation of player values
sd.pl=4,
# missingness_option
miss_opt=2,
#missingness_parameter
miss_param=0.5,
sd.quad=0.005,
minage=18,
maxage=40,
seed0=2222021
)
{
  set.seed(seed0)
  outsim=NULL
for (i in 1:nreps)
{
player_data=full_data_gen(n_players=nn, sd.players=sd.pl, sd.quad=sd.quad,overall_mean=mu)
########################
# the function data_add_miss selects players to be missing from the full data
# NB: gar is the players value (goals above replacement)
#
# option =0 ==> no missing
# option =1 ==> missing is deterministic for values below missing_param
# option =2 ==> missing is stochastic based upon logistic prob = 1-ilogit(gar-missing_param))
# option =3 ==> missing is function of cumulative 3yr sum < missing_param

player_data=data_add_miss(player_data,option=miss_opt,missing_param=miss_param,minage,maxage)
# make data set with unobserved/missing obs eliminated
player_data_sifted=sift_players(player_data)
if(i==1) {
  write.csv(player_data_sifted,paste0("outsingle_",mu0,"_",nn1*10,
                            "_",sdp0,"_",miss0,"q.csv"))
}
modelmeans=tapply(player_data$true_gar,
                  player_data$age,
                  mean)
modelmeans=data.frame(modelmeans,age=as.numeric(names(modelmeans)))
allmeans=tapply(player_data_sifted$gar,
                player_data_sifted$age,
                mean)
allmeans=data.frame(allmeans,age=as.numeric(names(allmeans)))

obsmeans=tapply(player_data_sifted$gar[!is.na(player_data_sifted$gar_miss)],
                player_data_sifted$age[!is.na(player_data_sifted$gar_miss)],
                mean)
obsmeans=data.frame(obsmeans,age=as.numeric(names(obsmeans)))

pct_obs<-table(player_data_sifted$age[!is.na(player_data_sifted$gar_miss)])/(nn*10)
pct_obs=data.frame(pct_obs=as.vector(pct_obs),age=as.numeric(names(pct_obs)))
  
cat.means<-right_join(modelmeans,allmeans,by="age") %>%
  left_join(obsmeans,by="age")%>%
  left_join(pct_obs,by="age") %>%
  dplyr::select(-age)

datameans<-data.frame(method=c("true_gen","true_data","true_obs","pct_obs"),
                      t(cat.means))
colnames(datameans)=c("method",paste0("age",minage:maxage))

out_eval<-eval_model_meancurve(truecurve=y1+mu,player_data_sifted,player_data)
out_eval<-out_eval %>%
  bind_rows(datameans)%>%
  mutate(sim.number=i,
         mu=mu, 
         numb_players=10*nn,
         sd_players=sd.pl,
         missing_option=miss_opt,
         missing_param=miss_param,
         seed=seed0)
  outsim<-bind_rows(outsim,out_eval)
}

#write.csv(x=player_data_sifted,"outsim_0_1000_2_2_1.5b.csv")
return(outsim)
}

nn1=100; sdp0=0.8; mp0=1; mu0=1; miss0=4;
outtest1<-player_aging_sim(  mu=mu0,nn=nn1,   nreps=200,
 sd.pl=sdp0,   miss_opt=miss0,miss_param=mp0,sd.quad=0.005, minage=18,  maxage=40)


#outtest<-outtest1 
write.csv(outtest1,paste0("outtest_",mu0,"_",nn1*10,
                         "_",sdp0,"_",miss0,"q.csv"))