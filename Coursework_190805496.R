#====================================================================================

# Time Series for Business
# Coursework
# Amaya Syed - ID: 190805496

#====================================================================================

# install necessary packages
#install.packages("")

# loading necessary packages
library(rstudioapi)
library(naniar)
library(ggplot2)
library(imputeTS)
library(fpp)
library(fma)
library(tseries)
library(TSA) # fourrier transform
library(forecast)
library(PerformanceAnalytics)
library(dplyr)
library(vars)

# Setting work directory
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

#====================================================================================

# loading data
data_raw = read.csv("AirQuality.csv", header=TRUE)

# replacing the missing values coded as -200 with NA
data = data_raw %>% replace_with_na(replace = list(NOx = -200, NO2 = -200, Temp= -200, RH = -200, AH = -200))

# a) TRAINING/TESTING SETS

# "the training dataset should incude the first 296 observations (until the last observation of 2004). The test dataset should include the 31 daily observations for January 2005."

# splitting the data into training and testing sets
training = window(ts(data, frequency = 7), end = c(43, 2))
testing = window(ts(data, frequency = 7), start = c(43, 3), end = c(47, 5))

# Performing linear interpolation on all the training and testing ata
training_int = na_seadec(training, algorithm = "interpolation", find_frequency = TRUE)
testing_int = na_seadec(testing, algorithm = "interpolation", find_frequency = TRUE)

summary(training_int)
summary(testing_int)
# -----------------------------------------------------------------------------------

# I want to check the dominant frequency is in my dataset so I'll do a fourrier analysis. 

PGram = periodogram(training_int[,2])

PGramAsDataFrame = data.frame(freq=PGram$freq, spec=PGram$spec)
order = PGramAsDataFrame[order(-PGramAsDataFrame$spec),]
top2 = head(order, 2)

TimePeriod =1/top2[2,1]
TimePeriod # 6.976744 -> 7

# This quick analysis confirms that the dominant seasonality information in the dataset comes from the weekly pattern. 

# -----------------------------------------------------------------------------------

# Assigning each time series to individual variables for later use
ts_NOx = training_int[,2]
ts_NO2 = training_int[,3]
ts_Temp = training_int[,4]
ts_RH = training_int[,5]

test_NOx = testing_int[,2]
test_NO2 = testing_int[,3]
test_Temp = testing_int[,4]
test_RH = testing_int[,5]

# Plotting the main times series -
plot(training_int[,2:5], main="", xlab = "Time [weeks]")

# Checking the correlation between the series
chart.Correlation(training_int[,2:6], histogram=TRUE, pch=19)
# because AH does not correlate highly with NOx data we will remove it from our analysis.
chart.Correlation(training_int[,2:5], histogram=TRUE, pch=19)

# plotting the interpolated values in red on the original data in gray
autoplot(ts_NOx, series="I") +
  autolayer(training[,2], series="O") +
  scale_colour_manual(
    values=c(`I`="red",`O`="black"))

plot(ts_NOx)
plot(test_NOx)

# checking for outliers
tsoutliers(test_NOx) # supposedly none

# check for outliers visually
boxplot(ts_NOx, xlab = "NOx", ylab = "ppb")
# we still have outliers in the top values despite removing values over 1000. 
boxplot(test_NOx)

# autocorrolation in the data
acf(ts_NOx, main= "")
pacf(ts_NOx, main="")

# confirming the data is not stationary
adf.test(ts_NOx, alternative = "stationary") # p-value = 0.2519 -> non stationary

# checking out what the recommended box.cox transformation parameter is. 
lambda <- BoxCox.lambda(ts_NOx)
lambda # 0.08973438
# This is very close to zero, so we'll try a log transform

# decomposing the data
decomposed_data = decompose(ts_NOx)
plot(decomposed_data)

# effect of applying a log transformation and difference to the data
ln_NOx = log(ts_NOx)
diff_ln_NOx = diff(ln_NOx)
plot(ln_NOx)
plot(diff_ln_NOx, ylab= "diff(log(NO_x))") # the data now looks stationary

# checking if taking the log and differencing takes care of stationarity
adf.test(diff_ln_NOx, alternative = "stationary") # stationary
kpss.test(diff_ln_NOx) # stationary

# lambda = 0 and diff = 1 seems like a good choice for this data. 

acf(diff_ln_NOx, main = "") # lag in 1, 2, 3
pacf(diff_ln_NOx, main = "") # significan autocorrelation until lag 7 - 

# sARIMA MODELS

# checking what auto.arima thinks is the best sARIMA model. 
auto.arima(ts_NOx, trace = T, stepwise = F, approximation = F, lambda = 0) # Best model found: ARIMA(0,1,3)(2,0,0)[7] AIC=454.13   AICc=454.43   BIC=476.26

sarima1 = Arima(ts_NOx, order= c(0,1,3), seasonal = c(2, 0, 0), lambda= 0)
sarima1
sarima1_res = residuals(sarima1)

checkresiduals(sarima1)
pacf(sarima1_res)
qqnorm(sarima1_res)
qqline(sarima1_res, distribution=qnorm)

# can we find a lower AIC by tweaking the model ? 

sarima2 = Arima(ts_NOx, order= c(2,1,1), seasonal = c(1, 0, 3), lambda = 0)
sarima2
sarima2_res = residuals(sarima2) # AIC = 419.03 

sarima3 = Arima(ts_NOx, order= c(2,1,1), seasonal = c(1, 1, 1), lambda = 0)
sarima3
sarima3_res = residuals(sarima3) # AIC = 409.36


# Lowest AIC we could find was 409 corresponding to the sarima3 model. Lets check residuals, fit and forecast for all three models. 

# checking residuals
# sarima2
checkresiduals(sarima2)
pacf(sarima2_res)
qqnorm(sarima2_res)
qqline(sarima2_res, distribution=qnorm)

# sarima3
checkresiduals(sarima3)
pacf(sarima3_res)
qqnorm(sarima3_res)
qqline(sarima3_res, distribution=qnorm)

# fitting the model
fitted_values_sarima1 = fitted.values(sarima1)
plot(fitted_values_sarima1,  col="blue", lty=2)
lines(ts_NOx)

fitted_values_sarima2 = fitted.values(sarima2)
plot(fitted_values_sarima2,  col="blue", lty=2)
lines(ts_NOx)

fitted_values_sarima3 = fitted.values(sarima3)
plot(fitted_values_sarima3,  col="blue", lty=2)
lines(ts_NOx)

# RMSE for model fit
model_rmse1 = sqrt(sum((fitted_values_sarima1 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse1 # 9.089435

model_rmse2 = sqrt(sum((fitted_values_sarima2 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse2 # 8.227808

model_rmse3 = sqrt(sum((fitted_values_sarima3 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse3 # 8.120169

# forecasting 31 days in the future
forecast1 = forecast(sarima1, h=31)
plot(forecast1)
lines(test_NOx, col="red")

forecast2 = forecast(sarima2, h=31)
plot(forecast2)
lines(test_NOx, col="red")

forecast3 = forecast(sarima3, h=31)
plot(forecast3)
lines(test_NOx, col="red")

# RMSE for forecast fit
forecast_rmse1 = sqrt(sum((forecast1$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse1 # 41.10467

forecast_rmse2 = sqrt(sum((forecast2$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse2 # 36.41724

forecast_rmse3 = sqrt(sum((forecast3$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse3 # 36.46474

# model sarima2 is marginally better, but not significantly so. Overall the sARIMA models seem to provide adequate mean values which follows the fluctuactions in NOx values - except for the two before last values in the testing set!!! instrument erro? Weird climate?

# DYNAMIC REGRESSION MODELS

# NOx modeled as function of n explanatory variables, with the error term allowed to be autocorrelated instead of assumed to be white noise. 

reg_arima1 = auto.arima(ts_NOx, xreg= ts_NO2 + ts_Temp + ts_RH, lambda= 0)
reg_arima1 # Regression with ARIMA(1,1,4)(2,0,0)[7] errors - AIC = 227.48 

# we'll tweak that
reg_arima2 = Arima(ts_NOx, order= c(0, 1, 4), seasonal= c(0,1,1), xreg= ts_NO2 + ts_Temp + ts_RH, lambda= 0) # AIC = 192.88
reg_arima2

reg_arima3 = Arima(ts_NOx, order= c(0, 1, 3), seasonal= c(0,1,1), xreg= ts_NO2 + ts_RH, lambda= 0) # AIC = 196.71
reg_arima3

# checking residuals
# reg_arima1
checkresiduals(reg_arima1)
res_reg1 = residuals(reg_arima1)
pacf(res_reg1)
qqnorm(res_reg1)
qqline(res_reg1, distribution=qnorm)

# reg_arima2
checkresiduals(reg_arima2)
res_reg2 = residuals(reg_arima2)
pacf(res_reg2)
qqnorm(res_reg2)
qqline(res_reg2, distribution=qnorm)

# reg_arima3
checkresiduals(reg_arima3)
res_reg3 = residuals(reg_arima3)
pacf(res_reg3)
qqnorm(res_reg3)
qqline(res_reg3, distribution=qnorm)

# fitting the models
fitted_values_reg1 = fitted.values(reg_arima1)
plot(fitted_values_reg1,  col="blue", lty=2)
lines(ts_NOx)

fitted_values_reg2 = fitted.values(reg_arima2)
plot(fitted_values_reg2,  col="blue", lty=2)
lines(ts_NOx)

fitted_values_reg3 = fitted.values(reg_arima3)
plot(fitted_values_reg3,  col="blue", lty=2)
lines(ts_NOx)

# RMSE for model fit
model_rmse4 = sqrt(sum((fitted_values_reg1 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse4 # 6.226941

model_rmse5 = sqrt(sum((fitted_values_reg2 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse5 # 5.484693

model_rmse6 = sqrt(sum((fitted_values_reg3 - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse6 # 5.578701

# forecasting 31 days into the future
forecast4 = forecast(reg_arima1, xreg= rep(test_Temp + test_RH + test_NO2), h = 31)
plot(forecast4)
lines(test_NOx, col="red")

forecast5 = forecast(reg_arima2, xreg= rep(test_Temp + test_RH + test_NO2), h = 31)
plot(forecast5)
lines(test_NOx, col="red")

forecast6 = forecast(reg_arima3, xreg= rep(test_RH + test_NO2), h = 31)
plot(forecast6)
lines(test_NOx, col="red")

# RMSE for forecast
forecast_rmse4 = sqrt(sum((forecast4$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse4 # 41.26112

forecast_rmse5 = sqrt(sum((forecast5$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse5 # 42.91447

forecast_rmse6 = sqrt(sum((forecast6$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse6 # 46.41265

# EXPONENTIAL SMOOTHING MODELS
# Because of the nature of the data, I will apply Holt-Winter directly to the data - seasonal = multiplicative. 

NOx_hw <- hw(ts_NOx, h= 31, seasonal = "multiplicative")
plot(NOx_hw)
lines(NOx_hw$fitted, col='red',main="")
lines(test_NOx, col="red")

# checking residuals
plot(as.vector(NOx_hw$fitted),as.vector(NOx_hw$residuals), xlab="Fitted values",ylab="Residuals")
hist(NOx_hw$residuals,main="")
plot(NOx_hw$residuals,ylab="Residuals") # not great
acf(NOx_hw$residuals)

# model rmse
model_rmse7 = sqrt(sum((NOx_hw$fitted - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse7 # 7.748265

# forecast rmse
forecast_rmse7 = sqrt(sum((NOx_hw$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse7 # 39.12458

# NEURAL NETWORK MODEL
# Feed-forward neural networks with a single hidden layer and lagged inputs for forecasting univariate time series.

set.seed(32)
nn_model = nnetar(ts_NOx)
nn_forecast = forecast(nn_model, h = 31, PI =T)
plot(nn_forecast)
lines(test_NOx, col="red")

# checking residuals
checkresiduals(nn_forecast)
qqnorm(nn_forecast$residuals)
qqline(nn_forecast$residuals, distribution=qnorm)

# model rmse
model_rmse8 = sqrt(sum((nn_forecast$fitted - ts_NOx) ^ 2)) / length(ts_NOx)
model_rmse8 # NA

# forecast rmse
forecast_rmse8 = sqrt(sum((nn_forecast$mean - test_NOx) ^ 2)) / length(test_NOx)
forecast_rmse8 # 39.41189
