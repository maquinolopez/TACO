   ###########################################
###          Bayesian Weighted Clymo          ###
###          By Marco A. Aquino-Lopez         ###
   ###########################################
library(Rtwalk)
rm(list=ls())
wkdir="~/GitHub/Bayesian_Carbon_Acc/Data/"
setwd(wkdir)
# MCMC parameters
th = 200
s_size = 1e+3
burnin = 2e+3
newrun = TRUE
# read data
Data_full <- read.csv('./Cpeat_final_0415.csv')
rownames <- names(Data_full)
sites <- unique(Data_full$site)  
# prepare data
site <- Data_full[which(Data_full$site == sites[46]),] # 46 work well - 7 does not work at all
Data <- site$cumulative_mass
Data <- diff(site$cumulative_mass) 
n_data <- length(Data)
Psi = 100 
dates <- site$zeroed_age
top_ages <- head(dates,n_data)
bot_ages <- tail(dates,n_data)
depths <- tail(site$depth,n_data) - (site$depth[2]-site$depth[1])/2
delay_p<- 5
boun_p <- 7
# plot age depth model and data over depth
par(mfrow=c(2,1))
plot(depths,top_ages,type = 'l',main='Age depth model')
plot(depths,Data,type = 'l',main='Mass accumulation per section')
which(Data==max(Data))
##
# Functions
##
w <- function(param){
  b_1 = log(Psi-1)/param[5] # param[5] is the delay
  b_0 = -b_1*param[6] # param[6] is the A/C boundery
  w <- 1/ (1 + exp(-b_0-b_1*depths))
}

simu <- function(param){
  # calculate for the cat
  z1 <-  (param[1] / param[2]) * (exp(-param[2] * top_ages) - exp(-param[2] * bot_ages))  
  # calculate from the acro
  z2 <-  (param[3] / param[4]) * (exp(-param[4] * top_ages) - exp(-param[4] * bot_ages))
  # calculate the logistic regression   
  w1 <- w(param)
  # calculate the expected carbon given the weight function
  z <- w1 * z1 + (1-w1) * z2
  return(z)
}

supp <- function (param){
  upper = c(rep(1, 4), # first 4 parameters should be lower than 1
            50, # delay parameter should be lower than 20
            tail(depths,1) # breaking point should be lower than the the highest depth
            )
  if (all(param>0) + all(param<upper) == 2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

ll <- function(param){
  z = simu(param)
  z = Data-z
  l = dnorm(z,mean = 0,sd = Data*.01,log = TRUE)
  mask = !(l == -Inf)
  l <- sum(l[mask]) 
  return(l)
}

prior_d = function(param){
  # priors for the decay model 
  d1 = sum(dnorm(head(param,4), rep(c(.01, .001),2),rep(c(.1,.1),2), log =TRUE))
  # prior for the weight function model
  d2 = sum(dnorm(tail(param,2), mean= c(delay_p,boun_p), sd = c(1,8), log =TRUE)) 
  return( d1 + d2)
}

obj <- function(param){
  return( -(prior_d(param) + ll(param)) )
}

ini = function(param){
  d1 = abs( rnorm(4, rep(c(.01, .001),2),rep(c(.01,.1),2)) )
  d2 = abs( rnorm(2, mean= c(delay_p,boun_p), sd = c(.5,4)) )
  return(c(d1, d2))
}
##
# Run twalk
##
if (newrun){
  twalk <- Runtwalk( dim=length(ini()),  Tr=th *(burnin+s_size),  Obj=obj, Supp=supp, x0=ini(), xp0=ini(),PlotObj = FALSE,PlotLogPost = FALSE) 
}else{
  state <- as.matrix(read.csv('state.csv'))
  twalk <- Runtwalk( dim=length(ini()),  Tr=th *(burnin+s_size),  Obj=obj, Supp=supp, x0=state[1,], xp0=state[2,],PlotObj = FALSE,PlotLogPost = FALSE) 
}
subsample = seq(1, length(info$output[,1]), th)
sample <- tail(info$output[subsample,],s_size)

##
# Plot results
##

# energy plot
par(mfrow=c(1,1))
plot(tail(info$Us[subsample],s_size),type = 'l')
# posterior distributions of parameters
par(mfrow=c(2,2))
plot(density(sample[,1]),main='Density a Cato')
plot(density(sample[,2]),main='Density b Cato')
plot(density(sample[,3]),main='Density a Acro')
plot(density(sample[,4]),main='Density b Acro')
# posterior distribution of boundery parameters
par(mfrow=c(3,1))
hist(sample[,5],main='Delay')
change_date <-  approx(site$depth,site$zeroed_age,sample[,6])$y
hist(change_date,main='A/C boundery age')
points(1960,0,pch=16,col=rgb(1,0,0,.9))
hist(sample[,6],main='A/C boundery depth ')
points(63,0,pch=16,col=rgb(1,0,0,.9))
# Cato limit function
par(mfrow=c(1,1))
plot(depths,w(sample[,1]),type = 'l',col=rgb(0,0,0,.1),xlab = 'Depth',
     ylab = 'Cato limit',ylim = c(0,1))
for (s in 1:dim(sample)[1]){
  lines(depths,w(sample[s,]),col=rgb(0,0,0,.1))
}
# data over depths 
par(mfrow=c(1,1))
plot(site$depth,site$cumulative_mass,type='l',col=rgb(1,0,0,.5))
polygon(c(site$depth,rev(site$depth) ), 
        c(c(0,cumsum(Data*(.9))), rev(c(0,cumsum(Data*(1.1))))),
        col=rgb(1,0,0,.5),border = NA)

for (s in 1:dim(sample)[1]){
  lines(head(site$depth,n_data),cumsum(simu(sample[s,])),col=rgb(0,0,0,.015))
  abline(v=sample[s,6],col=rgb(0,1,0,.01))
}
points(boun_p,0,col=rgb(1,0,0,.8),pch='|')
# data vs sumulated over ages
plot(top_ages,Data,type='l',col=rgb(1,0,0,.5))
polygon(c(top_ages,rev(top_ages) ), c(Data*(.9), rev(Data*(1.1))), col=rgb(1,0,0,.5),border = NA)
for (s in 1:dim(sample)[1]){
    lines(top_ages,simu(sample[s,]),col=rgb(0,0,0,.015))
    abline(v=approx(site$depth,site$zeroed_age,sample[s,6])$y,col=rgb(0,1,0,.01))
}
points(approx(site$depth,site$zeroed_age,boun_p)$y,0,col=rgb(1,0,0,.8),pch='|')


# This saves the state of the MCMC. 
#It can be use to continue the chain
state = matrix(c(tail(twalk$output,1),tail(twalk$outputp,1)),ncol=length(tail(twalk$output,1)),byrow = T)
write.csv(state,'state.csv',col.names = F,row.names = F)