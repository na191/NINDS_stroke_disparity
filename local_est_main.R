rm(list=ls())
library(glmnet)
library(splines)
library(MASS)
library(dplyr)
library(survival)
library(statmod)
library(ggplot2)
source('function_local.R')
set.seed(123)
### The commented code is for simulating data only
# n_clinical_feature = 15
# num_bh_knot = 2
# d_phi = 3
# T0=300
# lambda=0.05
# Xcov=1*ar1_cor(n_clinical_feature,0.2)
# rateC=0.1
# missing_rate = c(rep(.1,5),rep(.3,5),rep(.5,5))
# betas = data.frame(X1=-.8,X2=.8,X3=.1,X4=-.1,
#                    M1=.8,M2=.8,M3=-.1,M4=.1,
#                    'fw1.X1'=.6,'fw1.X2'=-.6,'fw1.M1'=.6,'fw1.M2'=-.6,fw=1)
# fw_coef=1
# fw_intercept=1
# dist_family = 'gompertz'
# dat_tab=local_est_tab=sum_stat_tab=params_tab=list()
# 
# shape_noise = rnorm(2,0,0.1) # generate purterbation for coefficients and intercept of shape function
# local_fw_coef = shape_noise[1]+fw_coef
# local_fw_intercept = shape_noise[2]+fw_intercept
# rho=runif(1,0.3,0.6)
# n_patient = sample(c(3000,1000),1)
# beta_noise = rnorm(ncol(betas),0,0.1)
# local_betas = beta_noise+betas
# ### simulate data for each site given parameters
# dat_obj=generate_data(n_patient,n_clinical_feature,local_betas,T0,missing_rate,
#                       lambda,rho,rateC,Xcov,local_fw_coef,local_fw_intercept,d_phi,dist_family=dist_family)
# sample_data = dat_obj$dat.obs %>% select(id:X15,calendar_time)
# colnames(sample_data) = c('id','y','failed',
#                           'a','b','c','d','e',
#                           'f','g','h','i','j',
#                           'k','l','m','n','o','calendar_time')
# write.csv(sample_data,'sample_data.csv',row.names=F)


### The simulated data is stored in the sample_data.csv file
sample_data <- read.csv("sample_data.csv")

sample_data <- read.csv('~/111 Dropbox/Molei Liu/distributed_cox/code/sample_data.csv')

### user have to specify 
### x: clinical feature; 
### time: survival time; 
### event: censorship (1 for not censored); 
### calendar_time: calendar_time 
x = sample_data %>% select(a:o)
time = sample_data$y
event = sample_data$failed
calendar_time = sample_data$calendar_time


# see surv_local in function_local.R for details
sum_stat = surv_local(time=time, event=event, 
                      x=x, calendar_time=calendar_time, num_bh_knot=2, d_phi=3)

# User-specified mean and sd vectors 
# (NA in the j-th entry represents not standardizing xj):

mean_user <- c(rep(0.5, 13), NA, NA)
sd_user <- c(rep(2, 13), NA, NA)

# Use mean_user and sd_user for standardization

sum_stat = surv_local(time=time, event=event, 
                      x=x, calendar_time=calendar_time, num_bh_knot=2, d_phi=3,
                      mean = mean_user, sd = sd_user)

# To add
sum_stat$X <- NULL
sum_stat$cbhhat.obj <- NULL
sum_stat$alasso <- NULL

save(sum_stat, file = 'xx')


