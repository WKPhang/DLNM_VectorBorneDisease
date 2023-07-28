rain_avg <- mydata$rain_avg # 2 3 -37.6 -33.7
rain_sd <- mydata$rain_sd # 1 4 -39.3 -36.3
rain_kurt <- mydata$rain_kurt # 1 1 -34.1 -31.9
temp_avg <- mydata$temp_avg # 1 1 -36.9 -34.6 # 1658.5
temp_kurt	 <- mydata$temp_kurt # 3 1 -36.8 -33.5
t_max_avg	<- mydata$t_max_avg # 5 1 # 1645.3
t_min_avg	<- mydata$t_min_avg # 1 1 # 1680.7
t_min_sd <- mydata$t_min_sd # 1 1
rh_avg	<- mydata$rh_avg # 1 1
rh_sd <- mydata$rh_sd # 5 1
csSOI <- mydata$csSOI # 1 1
csIOD <- mydata$csIOD # 1 1
?crossbasis

rain_avg.cb <- crossbasis(mydata$rain_avg[2:417], lag=16, 
                         argvar=list(fun="ns", df = 5),
                         arglag=list(fun="poly", degree = 2))

rain_sd.cb <- crossbasis(mydata$rain_sd[2:417], lag=16, 
                          argvar=list(fun="ns", df= 2),
                          arglag=list(fun="poly", degree = 1))

rain_kurt.cb <- crossbasis(mydata$rain_kurt[2:417], lag=16, 
                          argvar=list(fun="ns", df= 1),
                          arglag=list(fun="poly", degree = 2))

temp_avg.cb <- crossbasis(mydata$temp_avg[2:417], lag=16, 
                          argvar=list(fun="ns", df= 3),
                          arglag=list(fun="poly", degree = 1))

rh_avg.cb <- crossbasis(mydata$rh_avg[2:417], lag=16, 
                             argvar=list(fun="ns", df= 2),
                             arglag=list(fun="poly", degree =1))

rh_sd.cb <- crossbasis(mydata$rh_sd[2:417], lag=16, 
                             argvar=list(fun="ns", df= 5),
                             arglag=list(fun="poly", degree = 1))

csSOI.cb <- crossbasis(mydata$csSOI[2:417], lag=16, 
                          argvar=list(fun="ns", df= 2),
                          arglag=list(fun="poly", degree = 1))

csIOD.cb <- crossbasis(mydata$csIOD[2:417], lag=16, 
                          argvar=list(fun="ns", df = 5),
                          arglag=list(fun="poly", degree = 1))

#Based on BIC

rain_avg.cb <- crossbasis(mydata$rain_avg[2:417], lag=16, 
                          argvar=list(fun="ns", df = 1),
                          arglag=list(fun="poly", degree = 1))

rain_sd.cb <- crossbasis(mydata$rain_sd[2:417], lag=16, 
                         argvar=list(fun="ns", df= 1),
                         arglag=list(fun="poly", degree = 1))

rain_kurt.cb <- crossbasis(mydata$rain_kurt[2:417], lag=16, 
                           argvar=list(fun="ns", df= 1),
                           arglag=list(fun="poly", degree = 1))

temp_avg.cb <- crossbasis(mydata$temp_avg[2:417], lag=16, 
                          argvar=list(fun="ns", df= 1),
                          arglag=list(fun="poly", degree = 1))

rh_avg.cb <- crossbasis(mydata$rh_avg[2:417], lag=16, 
                        argvar=list(fun="ns", df= 1),
                        arglag=list(fun="poly", degree =1))

rh_sd.cb <- crossbasis(mydata$rh_sd[2:417], lag=16, 
                       argvar=list(fun="ns", df= 1),
                       arglag=list(fun="poly", degree = 1))

csSOI.cb <- crossbasis(mydata$csSOI[2:417], lag=16, 
                       argvar=list(fun="ns", df= 2),
                       arglag=list(fun="poly", degree = 1))

csIOD.cb <- crossbasis(mydata$csIOD[2:417], lag=16, 
                       argvar=list(fun="ns", df = 1),
                       arglag=list(fun="poly", degree = 1))




#Test all


cores=detectCores()
cores
cl <- makeCluster(cores)
cl
registerDoParallel(cl)
Parm_check <- foreach (ai = 1:5, .combine=rbind) %do% {
  foreach( aj= 1:5, .combine = rbind) %do% {
    foreach(bi = 1:5, .combine = rbind) %do% {
      foreach(bj = 1:5, .combine = rbind) %do% {
        foreach(ci = 1:5, .combine = rbind) %do% {
          foreach(cj = 1:5, .combine = rbind) %do% {
            rain_avg.cb <- crossbasis(mydata$rain_avg[2:417], lag=15, 
                                      argvar=list(fun="ns", df = ai),
                                      arglag=list(fun="poly", degree = aj))
            temp_avg.cb <- crossbasis(mydata$temp_avg[2:417], lag=15, 
                                      argvar=list(fun="ns", df= bi),
                                      arglag=list(fun="poly", degree = bj))
            rh_avg.cb <- crossbasis(mydata$rh_avg[2:417], lag=15, 
                                    argvar=list(fun="ns", df= ci),
                                    arglag=list(fun="poly", degree =cj))
      model_1var <- glm.nb(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) + 
                          + rain_avg.cb + rh_avg.cb + temp_avg.cb, data=mydata)
      
      Sys.sleep(0.00000000001)
      xyz <- c(ai, aj, bi, bj, ci, cj)
      print(xyz)
      # update GUI console
      flush.console()
      Parm_check <- c(ai, aj, bi, bj, ci, cj,
                      AIC(model_1var), AICc(model_1var), BIC(model_1var), 
                      Predictive_Rsqr(model_1var))
}}}}}}
stopCluster(cl)

Parm_check_df <- as.data.frame(Parm_check)
Parm_check[which.min(Parm_check[,8]),]

rain_avg.cb <- crossbasis(mydata$rain_avg[2:417], lag=15, 
                          argvar=list(fun="ns", df = 3),
                          arglag=list(fun="poly", degree = 1))
temp_avg.cb <- crossbasis(mydata$temp_avg[2:417], lag=15, 
                          argvar=list(fun="ns", df= 1),
                          arglag=list(fun="poly", degree = 1))
rh_avg.cb <- crossbasis(mydata$rh_avg[2:417], lag=15, 
                        argvar=list(fun="ns", df= 1),
                        arglag=list(fun="poly", degree =1))

best_mod <- glm.nb(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417]))+
                     rain_avg.cb + rh_avg.cb, data=mydata)
