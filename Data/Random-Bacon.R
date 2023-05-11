###################################################
###  Code to get random age-depth models        ###
###           from Bacon                        ###
###  By Marco A. Aquino-Lopez                   ###
###################################################
library(rbacon)
size_sample = 1e+3

setup_RandBacon <- function(core,folder,run){
  Bacon(coredir = folder,core = core,run = run,ask = F, suggest = F)
  agedepth()
  age_depth_models <<- sapply(info$elbows,Bacon.Age.d)
  elbows <<- info$elbows
}



Random_bacon <- function(x){
  ad_model <- age_depth_models[as.integer(runif(1,1,dim(age_depth_models)[1])),]
  approx(elbows,ad_model,x)$y
}


