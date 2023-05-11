###################################################
###  Testing Carbon Accumulations               ###
###  By Marco A. Aquino-Lopez                   ###
###################################################
library(Rtwalk)
rm(list=ls())
#for profiling delete before updating
Rprof(tmp <- tempfile())
wkdir="~/MEGA_A/PaleoStats 101/Carbon accumulation/code/"
source("Random-Bacon.R")
setwd(wkdir)
# MCMC parameters
th = 400
s_size = 1e+3
burnin = 1e+4
newrun = FALSE
# define Bacon paramebers
folder = "~/MEGA_A/PaleoStats 101/Carbon accumulation/code/Bacon_runs/"
core = 'Titus'
run = FALSE
#set up bacon
setup_RandBacon(core, folder ,run)  
# read data
Data_full <- read.csv('./Cpeat_final_0415.csv')
rownames <- names(Data_full)
sites <- unique(Data_full$site)  
# prepare data
site <- Data_full[which(Data_full$site == core),] # 46 work well - 7 does not work at all
Data <- site$cumulative_mass
Data <- diff(site$cumulative_mass) 
n_data <- length(Data)
depths <- tail(site$depth,n_data) - (site$depth[2]-site$depth[1])/2
Psi = 100 
delay_p<- 5
boun_p <- 7
delay_p_sd <- .1 
boun_p_sd <- 2

# plot age depth model and data over depth
par(mfrow=c(2,1))
plot(site$depth,site$zeroed_age,type = 'l',main='Age depth model')
plot(head(site$depth,n_data),Data,type = 'l',main='Mass accumulation per section')
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
  ages <- Random_bacon(site$depth) 
  ages <- ages - ages[1]
  top_ages <- head(ages , n_data)
  bot_ages <- tail(ages , n_data)
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
            tail(site$depth,.5*n_data)[1] # breaking point should be lower than the the highest depth
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
  d2 = sum(dnorm(tail(param,2), mean= c(delay_p,boun_p), sd = c(delay_p_sd,boun_p_sd), log =TRUE)) 
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
subsample = seq(1, length(twalk$output[,1]), th)
sample <- tail(twalk$output[subsample,],s_size)

##
# Plot results
##

# energy plot
par(mfrow=c(1,1))
plot(tail(twalk$Us[subsample],s_size),type = 'l')
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
abline(v = boun_p,col=rgb(1,0,0,.8))
# data vs simulated over ages
ages <- Random_bacon(site$depth)
ages <- ages - ages[1]
top_ages <- head(ages , n_data)
plot(top_ages,Data,type='l',col=rgb(1,0,0,.5),ylim=c(0,max(Data*1.1)))
polygon(c(top_ages,rev(top_ages) ), c(Data*(.9), rev(Data*(1.1))), col=rgb(1,0,0,.5),border = NA)
for (s in 1:dim(sample)[1]){
  lines(top_ages,simu(sample[s,]),col=rgb(0,0,0,.15))
  abline(v=approx(site$depth,site$zeroed_age,sample[s,6])$y,col=rgb(0,1,0,.01))
}
abline(v=approx(site$depth,site$zeroed_age,boun_p)$y,col=rgb(1,0,0,.8))
# energy plot
par(mfrow=c(1,1))
plot(tail(twalk$Us[subsample],s_size),type = 'l')


# end of profiling - delete before uploading
Rprof()
summaryRprof(tmp)

state = matrix(c(tail(twalk$output,1),tail(twalk$outputp,1)),ncol=length(tail(twalk$output,1)),byrow = T)

write.csv(state,'state.csv',col.names = F,row.names = F)





