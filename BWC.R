   ###########################################
###          Bayesian Weighted Clymo          ###
###          By Marco A. Aquino-Lopez         ###
   ###########################################
library(Rtwalk)
rm(list=ls())
wkdir="~/GitHub/Bayesian_Carbon_Acc/Data/"
setwd(wkdir)
# MCMC parameters
th = 100 # MCMC thinning 
s_size = 1e+3 # Final sample size
burnin = 2e+3 # burin-in (to be discarted)
newrun = FALSE
# read data
Data_full <- read.csv('./Cpeat_final_0415.csv')
rownames <- names(Data_full)
sites <- unique(Data_full$site)  

# prepare data
site <- sites[3] # 177 
print(site)
site <- Data_full[which(Data_full$site == site),] # 46 work well - 7 does not work at all
Data <- diff(site$cumulative_mass) 
# Data <- site$cumulative_mass[-1]
n_data <- length(Data)
Psi = 1000 # 1/psi is how quickly we get to a .5 proportion
dates <- site$zeroed_age
top_ages <- head(dates,n_data)
bot_ages <- tail(dates,n_data)
depths <- tail(site$depth,n_data) - (site$depth[2]-site$depth[1])/2
delay_p<- 2 # in depth
boun_p <- 8 # in depth
ADmodel <- function(x)(approx(depths,bot_ages - (top_ages-bot_ages )/2,x)$y)
# plot age depth model and data over depth
par(mfrow=c(2,1))
plot(depths,top_ages,type = 'l',main='Age depth model',ylab = 'Ages (yr)',xlab = '')
plot(depths,Data,type = 'l',main='Mass accumulation per section',xlab = 'Depth (cm)')

##
# Functions
##
w <- function(param){ # Pondering function (it is a logistic regression with a fix)
  b_1 = log(Psi-1)/param[5] # param[5] is the delay
  b_0 = -b_1*param[6] # param[6] is the A/C boundary
  w <- 1/ (1 + exp(-b_0-b_1*depths))
  return(w)
}

simu <- function(param){
  # param[1] and param[3] are peat addition  .001
  # param[2] and param[4] are decomposition  .0001
  # preivous formula cumulative_mass ~ (a / b) * (1 - exp(-b * zeroed_age)),
  # calculate for the cat
  # z1 <-  (param[1] / param[2]) * (exp(-param[2] * top_ages) - exp(-param[2] * bot_ages))
  # calculate from the acro
  # z2 <-  (param[3] / param[4]) * (exp(-param[4] * top_ages) - exp(-param[4] * bot_ages))
  # calculate for the cat
  ts <- bot_ages-ADmodel(param[6]-param[5])
  ts[ts<0] <- 0
  z1 <-  (param[1] / param[2]) * (1 - exp(-param[2] * ts))
  # # calculate from the acro
  z2 <-  (param[3] / param[4]) * (1 - exp(-param[4] * bot_ages))
  # calculate the logistic regression   
  w1 <- w(param)
  # calculate the expected carbon given the weight function
  z <- w1 * z1 + (1-w1) * z2
  return(z)
}

supp <- function (param){
  upper = c(rep(1, 4), # first 4 parameters should be lower than 1
            delay_p*4, # delay parameter should be lower than 2 * prior (boun_p)
            2*boun_p # AC boundary should be lower than the the highest depth
            )
  if (all(param>0) + all(param<upper) == 2){
    if (param[6]-param[5] > 0){
      if (all(param[c(1,3)] > param[c(2,4)] ))  {
        return(TRUE)   
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

ll <- function(param){
  z = simu(param)
  z = Data - c(z[1],diff(z))
  l = dnorm(z,mean = 0,sd = Data*.01,log = TRUE)
  mask = !(l == -Inf)
  l <- sum(l[mask]) 
  return(l)
}

prior_d = function(param){
  # param[1] and param[3] are peat addition  .001
  # param[2] and param[4] are decomposition  .0001
  # priors for the decay model 
  d1 = sum(dnorm(head(param,4), rep(c( 0.001 , 0.0001),2),rep(c(.005,.0005),2), log =TRUE))
  # prior for the weight function model, 
  # first is the delay (0) and then the boundary (boun)
  d2 = sum(dnorm(tail(param,2), mean= c(0,boun_p), sd = c(delay_p,4), log =TRUE)) 
  return( d1 + d2 )
}

obj <- function(param){
  return( -(prior_d(param) + ll(param)) )
}

ini = function(param){
  d1 = abs( rnorm(4, rep(c(.001, .0001),2),rep(c(.0005,.00005),2)) )
  d2 = abs( rnorm(2, mean= c(delay_p,boun_p), sd = c(.05,1)) )
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
# This saves the state of the MCMC. 
state = matrix(c(tail(twalk$output,1),tail(twalk$outputp,1)),ncol=length(tail(twalk$output,1)),byrow = T)
write.csv(state,'state.csv',row.names = F)
# subsampling
subsample = seq(1, length(twalk$output[,1]), th)
sample <- tail(twalk$output[subsample,],s_size)

##
# Plot results
##
# energy plot
par(mfrow=c(1,1))
plot(tail(twalk$Us[subsample],s_size),type = 'l',ylab = '- Log posterior',xlab='Iterations')

###########
# parameters plot
# posterior distributions of parameters
# param[1] and param[3] are peat addition  .001
# param[2] and param[4] are decomposition  .0001
par(mfrow=c(2,2))
plot(density(sample[,1]),main='Density of peat addition Cato')
plot(density(sample[,2]),main='Density of decomposition Cato')
plot(density(sample[,3]),main='Density of peat addition Acro')
plot(density(sample[,4]),main='Density of decomposition Acro')
###########
# posterior distribution of boundary parameters
par(mfrow=c(3,1))
plot(density(sample[,5]),main='Delay',xlab = 'cm')
change_date <-  approx(site$depth,site$zeroed_age,sample[,6])$y
plot(density(change_date),main='A/C boundery age',xlab = 'years before coring')
points(1960,0,pch=16,col=rgb(1,0,0,.9))
plot(density(sample[,6]),main='A/C boundery depth ',xlab = 'cm')
points(63,0,pch=16,col=rgb(1,0,0,.9))
###########
# Cato limit function
par(mfrow=c(1,1))
plot(depths,w(sample[,1]),type = 'l',col=rgb(0,0,0,.1),xlab = 'Depth',
     ylab = 'Cato limit',ylim = c(0,1))
for (s in 1:dim(sample)[1]){
  lines(depths,w(sample[s,]),col=rgb(0,0,0,.1))
}

###########
# data over prediction
# data over depths 
par(mfrow=c(1,1))
plot(tail(site$depth,length(site$site)-1),(Data),type='l',col=rgb(1,0,0,.5),
     ylim=c(min(Data*(.5)),max(Data*(1.4))))
polygon(c(site$depth[-1],rev(site$depth[-1]) ), 
        c(c((Data*(.97))), rev(c((Data*(1.03))))),
        col=rgb(1,0,0,.5))
for (s in 1:dim(sample)[1]){
  sim = simu(sample[s,])
  points(site$depth[-1],c(sim[1],diff(sim)),col=rgb(0,0,0,.015),pch='_')
  abline(v=sample[s,6],col=rgb(0,1,0,.01))
}
points(boun_p,0,col=rgb(1,0,0.5,.8),pch='|')

###########
# data vs simulated over ages
par(mfrow=c(1,1))
plot(bot_ages,Data,type='l',col=rgb(1,0,0,.5),ylim=c(min(Data*(.5)),max(Data*(1.4))))
polygon(c(bot_ages,rev(bot_ages) ), c(Data*(.97), rev(Data*(1.03))), col=rgb(1,0,0,.5),border = NA)
for (s in 1:dim(sample)[1]){
    sim = simu(sample[s,])
    points(top_ages,c(sim[1],diff(sim)),col=rgb(0,0,0,.015),pch='_')
    abline(v=approx(site$depth,site$zeroed_age,sample[s,6])$y,col=rgb(0,1,0,.01))
}
points(approx(site$depth,site$zeroed_age,boun_p)$y,0,col=rgb(1,0,0,.8),pch='|')

summary(sample[,1:4])

summary(sample[,5:6])
###########
# difference between simulation and data
par(mfrow=c(1,1))
sim = simu(colMeans(sample))
plot(top_ages,Data-c(sim[1],diff(sim)),col=rgb(0,0,0,.015),ylim=c(-max(Data*(.5)),max(Data*(1.4))),pch='_',
     xlab='age',ylab='difference between data and model')
for (s in 1:dim(sample)[1]){
  sim = simu(sample[s,])
  points(top_ages,Data-c(sim[1],diff(sim)),col=rgb(0,0,0,.015),pch='_')
  abline(v=approx(site$depth,site$zeroed_age,sample[s,6])$y,col=rgb(0,1,0,.01))
}
abline(h=0,col='red')


