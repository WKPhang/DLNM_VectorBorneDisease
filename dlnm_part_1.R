library(ggplot2)
library(lattice)
library(scales)
library(dlnm)
library(npreg)
library(MASS)
library(splines)
library(RColorBrewer)
library(reshape2)
library(MuMIn)
library(pscl) #zero inflated model
library(forecast) # need this for autoplot ts object

#### Actual work start here
getwd()
setwd("C:/YOUR/WORK/DIRECTORY")

#Data wrangling
mydata <- read.csv("dlnm_dataset.csv", header=T, stringsAsFactors = T)
mydata$Month <- as.factor(mydata$Month) #optional
head(mydata)
column_names <- colnames(mydata)

#Distribution of dependent varriable
hist(mydata$pk_case, breaks=20)
mean(mydata$pk_case)
var(mydata$pk_case)
100*sum(mydata$pk_case == 0)/length(mydata$pk_case)

#cubic spline interpolation for population
cspop <- ss(mydata$epi_week, mydata$pop, nknots=8, m=2)$y
plot(mydata$pop, type="l")
lines(cspop, type="l", col="blue")
mydata$pop <- cspop

#Create Incidence rate per 100 000 ppl
IR <- ((mydata$pk_case)/(mydata$pop))*1000000 
mydata$IR <- IR

#Smoothing the SOI and IOD
library(npreg)
mod.cs <-  ss(mydata$epi_week, mydata$SOI, all.knots=T, m=2)
mod.cs$fit
plot(mod.cs)

csSOI <- ss(mydata$epi_week, mydata$SOI, nknots=96, m=2)
csIOD <- ss(mydata$epi_week, mydata$IOD, nknots=96, m=2)

mydata$csSOI <- csSOI$y
mydata$csIOD <- csIOD$y

plot(mydata$IOD, type = "l")
lines(mydata$csIOD, type="l", col="red")

spline(mydata$Month)

#Comparison and mann-whitney test
colnames(mydata)

ggplot(mydata, aes(x =csSOI, y =  csIOD,
                   group=factor(group), colour=factor(group))) +
  geom_point() +
  scale_color_manual(values = c("grey", "black"))+
  theme_classic() +
  theme(text = element_text(size = 12, colour="black"), 
        axis.text = element_text(size = 15, colour="black")) +
  labs(x= "SOI", y = "DMI") +
  stat_ellipse(type="norm", lwd=1.2)

#Plot time-series graphs
tsactual <- ts(mydata$rain_sd, start=c(2012,1), frequency=52.25)
plot(tsactual, xaxt="n", yaxt="n")
axis(4)
axis(side=1, at=2012:2020)

#Non-parametric Spearman's cross-correlation
library(wmtsa)

mydata$pk_case[1:417]

#Cross-correlation analysis
## Correlation between climate variables
cor.test(mydata$rh_avg, mydata$csIOD,  method = "spearman")

## Negative lag means x (predictor) leads y (target) at lag of x
ccfvalues <- ccf(mydata$rain_avg,mydata$pk_case, lag.max = 24, ci=0.95)
plot(ccfvalues[-24:0], ylim=c(-0.3,0.4))
ccfvalues$series

#### Reversed the x and y variables
ccfvalues <- ccf(mydata$pk_case,mydata$csIOD, lag.max = 24, ci=0.95)
plot(ccfvalues[0:24], ylim=c(-0.3,0.4))



ccfvalues <- ccf(mydata$rain_avg ,mydata$pk_case, lag.max = 24, ci=0.95)
ccf(mydata$csSOI[18:417] ,mydata$resid[18:417], lag.max = 24, ci=0.95)

library(astsa)
lag2.plot (mydata$rain_avg, mydata$pk_case, 15)
?lag2.plot
## Run all ccf at once
str(mydata)

ccf_matrix <- matrix(0, 37, ncol(mydata[7:26]))
ccf_all <- as.data.frame(ccf_matrix)
k <- 1
for (i in mydata[,7:26]){
  ccfvalues <- ccf(mydata$logIR_1, i , lag.max = 36, ci=0.95)
  ccf_all[,k] <- ccfvalues$acf[-c(1:36)]
  k<-k+1
  print(k)
}

row_x <-matrix(0,37,1)
r <- 1
for (i in 0:36){
  row_x[r,] <- paste("lag_", i , sep ="")
  r <- r+1
  print(r)
}
row_x

colnames(ccf_all) <- colnames(mydata[7:26])
rownames(ccf_all) <- row_x
ccf_all

write.csv(ccf_all, file="ccf_all.csv", row.names=T, col.names=T)

##create heatmap from ccf_all
library(gplots)
col<- colorRampPalette(c("green", "yellow", "white", "blue", "purple"))(256)

heatmap.2(as.matrix(ccf_all), scale="none",
          col=col,
          dendrogram="none",
          trace="none",
          Rowv=F, 
          Colv=T,
          key=T,
          density.info = "none"
)

# Granger causality test
## HO: p > 0.05 , X does not granger cause Y
library(lmtest)
M <- 36
grgPvalue <- array(NA, M)

for ( j in 1:M){
  grgPvalue[j] <- grangertest(mydata$pk_case ~ mydata$rain_sd, order = j)$`Pr(>F)`[2]
}

## Plot the p-values against lag order to see the impact of lag order
plot(grgPvalue, type = "b", main = "P-values of the Granger test for 
     IR ~ IV.", xlab = "Order", ylab = "p-value")
abline(h = 0.1, col = "blue", lty = 2)
abline(h = 0.05, col = "red", lty = 2)

## Matrix of granger causality test

M <- 36
grgPvalue <- array(NA, M)
grgPvalue_matrix <- matrix(0, 36, ncol(mydata[7:26]))
grgPvalue_all <- as.data.frame(grgPvalue_matrix)
k <- 1

for (i in mydata[,7:26]){
  for (j in 1:M){
    grgPvalue[j] <- grangertest(mydata$logIR_1 ~ i , order = j)$`Pr(>F)`[2]
    grgPvalue_all[,k] <- grgPvalue
  }
  k <- k+1
  print(k)
}


row_x <-matrix(0,36,1)
r <- 1
for (i in 1:36){
  row_x[r,] <- paste("order_", i , sep ="")
  r <- r+1
  print(r)
}
row_x

colnames(grgPvalue_all) <- colnames(mydata[7:26])
rownames(grgPvalue_all) <- row_x
grgPvalue_all

write.csv(grgPvalue_all, file="grgPvalue_all.csv", row.names=T, col.names=T)

##create heatmap from grgPvalue_all
library(gplots)

heatmap.2(as.matrix(grgPvalue_all), scale="none",
          col=ifelse(grgPvalue < 0.05, "red","grey"),
          dendrogram="none",
          trace="none",
          Rowv=F, 
          Colv=T,
          key=T,
          density.info = "none"
)

#Construct null model

str(mydata$Month)
ggplot(mydata, aes(pk_case)) + geom_histogram(binwidth = 1)

summary(model_null)
model_null_1 <- glm.nb(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) + Month[2:417], data=mydata)
summary(model_null_1)
model_null_2 <- glm(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) + Month[2:417], family=poisson, data=mydata)
model_null_3 <- glm(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) + Month[2:417], family=quasipoisson, data=mydata)
model_null_4 <- zeroinfl(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) 
                         + Month[2:417] | pk_case[1:416] + offset(log(pop[2:417])) + Month[2:417] , data=mydata)

model_null_1.2 <- glm.nb(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])), data=mydata)

Predictive_Rsqr <- function(x) {
  numerator <- (sum((x$y - x$fitted.value)^2))
  denominator <- (nrow(x$model)*var(x$y))
  Predictive_Rsqr <- 1 - (numerator/denominator)
  return(Predictive_Rsqr)
}

summary(model_null)
AIC(model_null)
AICc(model_null)
BIC(model_null)
Predictive_Rsqr(model_null)

model_null <- model_null_1.2

plot(model_null$y, type="l")
lines(model_null$fitted.values, type="l", col="red")
lines(model_null$fitted.values, type="l", col="blue")

#MOdel fitting test
1 - pchisq(summary(model_null_1)$deviance, summary(model_null_1)$df.residual)

#Dispersion statistics
library(performance)
check_overdispersion(model_null_1)

library(pscl)
odTest(model_null_1)
?odTest
summary(model_null_2)
#Dispersion ratio
E2 <- resid(model_null_1, type = "pearson")
N  <- length(model_null_1$y)
p  <- attributes(logLik(model_null_1))$df
sum(E2^2) / (N - p)

#Run single parameter
weather_cb <- crossbasis(mydata$temp_avg[2:417], lag=15, 
                         argvar=list(fun="ns", df = 10),
                         arglag=list(fun="poly", degree = 1))
model_1var <- glm(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417])) + 
                    Month[2:417] + weather_cb, data=mydata)


summary(model_1var)
weather_cb
1151*2 
###### Test ########
model_1var_noAR <- glm.nb(pk_case[2:417] ~ offset(log(pop[2:417])) + 
                            weather_cb, data=mydata)

plot(model_1var$y, type="l", col="grey")
lines(model_1var$fitted.values, type="l", col="red")
lines(model_1var_noAR$fitted.values, type="l", col="blue")
lines(model_null$fitted.values[16:416], type="l", col="blue")
lines(model_null_1.2$fitted.values[15:415], type="l", col="green")

#Fix parameter automatically for all crossbasis
### Go to Part 2 ###

#Dominance analysis to identify important variables
library(dominanceanalysis)
mod_test <- glm.nb(pk_case[1:417] ~ offset(log(pop[1:417])) + Month[1:417] +
                     rain_avg.cb + temp_avg.cb + rh_avg.cb + csSOI.cb + csIOD.cb,
                   data=mydata)
da.glm.fit(mod_test)

dapres<-dominanceAnalysis(mod_test)
?da.negbin.fit
?dominanceAnalysis

#run DLNM
## remove stringasFactor for calling formula
## lag 0 at week 30
Modellist<- read.csv("exhaustive_var_test/exhaustive_var_combinelist.csv", header=T)
head(Modellist[1,])

#Prediction error index for dlnm
#MAPE can't be used if actual value is zero
MAE_dlnm <- function(x) {
  abs_val <- abs(x$y - x$fitted.value)
  MAE <- (1/nrow(x$model))*sum(abs_val)
  return(MAE)
}

RMSE_dlnm <- function(x) {
  sum_val <- sum((x$y - x$fitted.value)^2)
  RMSE <- sqrt(sum_val/nrow(x$model))
  return(RMSE)
}

Predictive_Rsqr <- function(x) {
  numerator <- (sum((x$y - x$fitted.value)^2))
  denominator <- (nrow(x$model)*var(x$y))
  Predictive_Rsqr <- 1 - (numerator/denominator)
  return(Predictive_Rsqr)
}

MDirAcc_lag1 <- function(x) {
  return( mean(sign(diff(x$y, lag=1))==sign(diff(x$fitted.value, lag=1))) )
}


cores=detectCores()
cores

cl <- makeCluster(cores)
cl
registerDoParallel(cl)
model_output <- foreach (i = 1:nrow(Modellist), .combine=rbind) %do% {
  fm <- Modellist[i,]
  dlnm_mod <- glm.nb(fm, data=mydata)
  AIC_out <- AIC(dlnm_mod)
  BIC_out <- BIC(dlnm_mod)
  AICc_out <- AICc(dlnm_mod)
  MASE_out <- MASE_dlnm(dlnm_mod)
  RMSE_out <- RMSE_dlnm(dlnm_mod)
  Predictive_Rsqr_out <- Predictive_Rsqr(dlnm_mod)
  MDA_out <- MDirAcc_lag1(dlnm_mod)
  LB_pvalue <- Box.test(ts(dlnm_mod$y - dlnm_mod$fitted.values),lag=16, type = "Ljung-Box")$p.value
  diff <- AICc(model_null) - AICc_out
  Sys.sleep(0.00000000001)
  print(i)
  # update GUI console
  flush.console()
  model_output <- c(fm, AIC_out, BIC_out, AICc_out, MASE_out, RMSE_out, 
                    Predictive_Rsqr_out, MDA_out,LB_pvalue , diff)
}
stopCluster(cl)
model_output_df <- as.data.frame(model_output)
model_output_df
col_name <- c("formula", "AIC", "BIC", "AICc", "MASE", "RMSE",
              "Predictive_R^2","MDA_out", "LB_pvalue", "null_diff")

colnames(model_output_df) <- col_name
head(model_output_df)
write.csv(model_output_df, file="exhaustive_model_output2.csv", row.names=T, col.names=T)
Modellist[1,]
#Cross-predict
##Select best mode
best_combination <- model_output_df[which.min(model_output_df$AICc),]
best_combination <- model_output_df[which.max(model_output_df$MDA_out),]
best_combination$formula
best_combination

best_mod <- glm.nb(best_combination$formula, data=mydata)
best_mod <- glm.nb(pk_case[2:417] ~ pk_case[1:416]+offset(log(pop[2:417]))+rain_avg.cb+rh_avg.cb+csSOI.cb, data=mydata)
best_mod <- glm.nb(pk_case[2:417] ~ pk_case[1:416] + offset(log(pop[2:417]))+rh_avg.cb, data=mydata)
best_mod <- glm.nb(pk_case[2:417] ~ offset(log(pop[2:417]))+rain_avg.cb + rh_avg.cb + csSOI.cb, data=mydata)
best_mod <- glm.nb(pk_case[2:417] ~ offset(log(pop[2:417]))+Month[2:417]+rain_avg.cb, data=mydata)
best_mod <- glm.nb(pk_case[2:417] ~ pk_case[1:416]+offset(log(pop[2:417]))+rain_avg.cb+rh_avg.cb+csSOI.cb+temp_avg.cb, data=mydata)
summary(best_mod)

AICc(best_mod)


###Goodness of fit test
with(best_mod, cbind(res.deviance = deviance, df = df.residual,
                     p = pchisq(deviance, df.residual, lower.tail=FALSE)))

#Ljung-Box test
ts_residuals <- ts(best_mod$y - best_mod$fitted.values)

library(portes)
LjungBox(model_null, lags=seq(5,30,5),order=0,season=1,squared.residuals=FALSE)

Box.test(ts_residuals,lag=16, type = "Ljung-Box")
p <- Box.test(ts_residuals, lag=16, type = "Ljung-Box", fitdf=1)
p$p.value

acf(ts_residuals)
pacf(ts_residuals)
plot(ts_residuals)

#produce residual vs. fitted plot
plot(best_mod$fitted.values, best_mod$y - best_mod$fitted.values)
abline(0,0)

qqnorm(best_mod$y - best_mod$fitted.values)
qqline(best_mod$y - best_mod$fitted.values)

library(car)
qqPlot(best_mod$y - best_mod$fitted.values)

shapiro.test(best_mod$y - best_mod$fitted.values)

#Plot
nrow(mydata)
mydata$pk_case
best_mod$y

#Time series plot
tsactual <- ts(best_mod$y, start=c(2012,18), frequency=52)
tspredict <- ts(best_mod$fitted.values, start=c(2012,18), frequency=52)
plot(tsactual, lwd=1, col="darkgrey")
lines(tspredict, col="red", lwd=1)

autoplot(tsactual, series="Observed",size=0.6)+
  autolayer(tspredict,  series="Predicted",size=0.7) +
  ylab("Number of cases") + xlab("Time")+
  ggtitle("") +
  scale_color_manual(labels = c("Observed", "Predicted"), values = c("gray26", "red"))+
  theme_bw()+
  theme(panel.grid.major = element_blank())+
  theme (panel.grid.minor = element_blank())

#Create confidence interval
library(tibble)
library(dplyr)
fam <- family(best_mod)
fam
str(fam)
ilink <- fam$linkinv
ilink
## add fit and se.fit on the **link** scale
ndata <- data.frame(best_mod$fitted.values)

ndata <- data.frame(mydata$pk_case[2:417], mydata$pop[2:417], mydata$Month[2:417])
colnames(ndata) <- c("pk_case", "pop", "Month")
ndata <- bind_cols(ndata, setNames(as_tibble(predict(best_mod, ndata, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))
ndata

ndata$index <- c(1:416)
ndata$fit[17:416] <- best_mod$fitted.values

plot(ndata$pk_case, type="l", col="grey")
lines(ndata$fit, type="l", col="red")
lines(ndata$right_upr, type = "l",lty=2, col="green")
lines(ndata$right_lwr, type = "l",lty=2, col="green")

#Insert residuals back into original dataset
nrow(best_mod$model)
mydata$residual[18:417] <- best_mod$y - best_mod$fitted.values


#Scatter plot
library("ggpubr")
scatter_data <- data.frame(fit = best_mod$fitted.values, observed= best_mod$y)

p <- ggscatter(x = "observed", y="fit", data= scatter_data, 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlim = c(0,17), ylim = c(0,17),
               yaxs="i", xaxs="i",
               xlab = "observed", ylab = "fit")
p + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

#Contour plot preparation
## Rainfall

percentiles <- round(quantile(mydata$IOD, c(0.05, 0.5, 0.95)), 2)
mean(mydata$csSOI)
percentiles

mean(mydata$rain_avg)

pred.prec <- crosspred(rain_avg.cb, best_mod, cen=7.62, by=0.1)
pred.temp <- crosspred(temp_avg.cb, best_mod, cen=25.31, by=0.1)
pred.rh <- crosspred(rh_avg.cb, best_mod, cen=83.52, by=0.1)
pred.csSOI <- crosspred(csSOI.cb, best_mod, cen=-0.29, by=0.1)
pred.csIOD <- crosspred(csIOD.cb, best_mod, cen=0.28, by=0.1)
pred.prec_sd <- crosspred(rain_sd.cb, best_mod, cen=4.3, by=0.1)
pred.prec_kurt <- crosspred(rain_kurt.cb, best_mod, cen=0.03, by=0.1)
pred.rh_sd <- crosspred(rh_sd.cb, best_mod, cen=2.54, by=0.1)

redall_PREC <- crossreduce(rain_avg.cb, best_mod, cen=83.52, by=0.1)

lab.prec <- "rain"

##Contour plot

plot(pred.prec, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Daily precipitation (mm)", ylab="Lag (weeks)"))

plot(pred.rh, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Relative humidity (%)", ylab="Lag (weeks)")) 

plot(pred.temp, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Mean temperature (degree Celcius)", ylab="Lag (weeks)")) 

plot(pred.csSOI, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="SOI", ylab="Lag (weeks)")) 

plot(pred.csIOD, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="DMI", ylab="Lag (weeks)"))

plot(pred.prec_sd, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Daily precipication SD", ylab="Lag (weeks)")) 

plot(pred.prec_kurt, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Daily precipication kurtosis", ylab="Lag (weeks)"))

plot(pred.rh_sd, "contour", key.title=title("RR"), cex.title=0.5, 
     plot.title=title("", xlab="Relative humidity SD", ylab="Lag (weeks)"))

#3-D plot
plot(pred.csSOI, zlab= "RR",
     plot.title=title("", xlab="Precipitation (mm)", ylab="Lag (weeks)"))


#Overall
plot(pred.prec, "overall", xlab=lab.prec, ylab="RR", ylim=c(0.0, 3.0), main="Overall Cumulative Association")
lines(redall_PREC, ci="lines", col=4, lty=3)
range(pred.prec$allRRfit)

##Sensitivity
plot.crosspred(pred.prec, "slices", var= c("20", "7", "1"),  
               lag=c(1, 10, 15), col="blue", ci.level=0.95,ylab="IR", 
               xlab="", cex.axis=1.25, lwd=1.5, lag.lab="Lag (weeks)", 
               var.lab="Precipitation (mm)", cex.main=0.8, cex.lab=1.25,  
               ci.arg=list(density=80, col=grey(0.7)) )


#Plot monthl=y seasonal
newdata1 <- data.frame(Month = factor(mydata$Month),
                       fit_val = mydata$rain_avg)

shapiro.test(mydata$rain_avg)

p <- ggplot(newdata1, aes(x=Month, y=fit_val)) + 
  geom_jitter(position=position_jitter(0.2), color="gray40")

p + stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
                 geom="errorbar", color="red", width =0.2) + 
  stat_summary(fun=mean, geom="point", color="red", cex=3) +
  theme_classic() +
  geom_hline(yintercept=0, lty="dashed", col="red")

p + stat_summary(fun.data=median_hilow,
                 geom="errorbar", color="red", width =0.2) + 
  stat_summary(fun=median, geom="point", color="red", cex=3) +
  theme_classic() +
  geom_hline(yintercept=0, lty="dashed", col="red")

cor(mydata$csSOI, mydata$rain_avg)
ccf(mydata$rain_avg, mydata$rh_avg)
?cor
