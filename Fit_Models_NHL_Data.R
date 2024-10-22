source("PlayerAging/NewAgingModelsTesting.R")

player_data<-read_csv("forwards_zppg_missing_long_2022.csv")%>%
  rename(age=ageyrs) %>%
  rename(gar_miss=z_ppg) %>%
  mutate(gar=gar_miss) %>%
  mutate(player=as.factor(playerid)) %>%
  group_by(playerid) %>%
  #filter(playerid<20000)%>%
  # "Dropout Mechanism" -- Based on 3-year GAR 
  # (Avg player career = 6.6 seasons)
  mutate(lag1 = lag(gar),lag2 = lag(gar,2)) %>%
  rowwise() %>% 
  mutate(threeyr = mean(c(lag1,lag2,gar),na.rm = T)) %>%
  group_by(player) %>%
  mutate(fcode = ifelse(lag(playerid)!=playerid,0,1),
         include = ifelse(threeyr<4,0,1),
         d = gar - lag(gar)) 

player_data_sifted<-sift_players(player_data)
player_data<-player_data_sifted
########################
# the function data_add_miss selects players to be missing from the full data
# NB: gar is the players value (goals above replacement)
#
# option =0 ==> no missing
# option =1 ==> missing is deterministic for values below missing_param
# option =2 ==> missing is stochastic based upon logistic prob = 1-ilogit(gar-missing_param))
# option =3 ==> missing is function of cumulative 3yr sum < missing_param

minage=min(player_data$age)
maxage=max(player_data$age)

numb_players=length(table(player_data$playerid))
pct_obs<-table(player_data_sifted$age[!is.na(player_data_sifted$gar_miss)])/(numb_players)
pct_obs=data.frame(pct_obs=as.vector(pct_obs),age=as.numeric(names(pct_obs)))

datameans<-data.frame(method=c("pct_obs"),
                      t(pct_obs))
colnames(datameans)=c("method",paste0("age",minage:maxage))

out_eval<-eval_model_meancurve(truecurve=0,player_data_sifted,player_data)
out_eval<-out_eval %>%
  bind_rows(datameans) 

names(out_eval)=c("method",paste0("age",minage:maxage))   
out_eval<-out_eval %>% filter(age18!=18)

write.csv(out_eval,paste0("Forwards_models_out_071022.csv"))