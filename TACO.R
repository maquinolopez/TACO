###########################################
###          Bayesian Weighted Clymo          ###
###          By Marco A. Aquino-Lopez         ###
### 
### TACO - Transitional Analysis through Clymo Optimization.
###########################################
rm(list = ls())

Taco <- function(core, folder, # Core and name parameters 
                 th = 100, s_size = as.integer(1e+3), burnin = as.integer(1e+3) , # MCMC parameters
                 ACboun_mean_prior= 30, ACboun_sd_prior= 5 ,  # Prior for the boundery
                 delay_mean_prior = 20, delay_sd_prior = 10,   # Delay parameter
                 peat_add_mean = 500 , peat_add_sd = 400.4,   # peat addition
                 alpha_mean = .0004, alpha_sd = .00013 ,            # decay Acre
                 alpha2_mean = .0004, alpha2_sd = .00013 ,    # decay Cat
                 epsilon = 1e-2,                              # tolerance parameter
                 newrun = TRUE){
  #### Read data and prepare variables ####
  
  # data structure and headers
  # depth,	bulk_density,	carbon,	zeroed_age,	biome
  Data <- read.csv(paste0(folder, '/', core,'/', core, '.csv'))
  # Remove rows with NA
  Data <- na.omit(Data)
  
  # Check if the cleaned data frame is empty
  if (nrow(Data) == 0) {
    stop("The data frame is empty after removing NA values. Check input data for NAs.")
  } 
  
  if(any(diff(Data$zeroed_age)<0)){
    stop("Double check your age-depth model")
  }
  
  # Set up the variable which will be used
  n_data <- length(Data$depth)
  epsilon1  <- log( (1-epsilon)/epsilon )# this is 1/ \epsilon
  error_p <- .1
  depths <- Data$depth
  top_ages <- Data$zeroed_age

  bot_ages <- c(Data$zeroed_age[-1], tail(Data$zeroed_age,1) + tail(diff(Data$zeroed_age),1) )
  top_depth <- Data$depth
  bot_depth <- Data$depth + diff(Data$depth)[1] # assuming same thickness of samples
  Data$curr_mass_den <- Data$bulk_density * 1.0e+4 # transform from g/cm^3 to g/m^3
  # Data$curr_mass_den <- Data$carbon
  sd_ll <- Data$curr_mass_den * error_p
  Dat <- Data$curr_mass_den
  
  peat_add_mean <- peat_add_mean #*0.0001 
  peat_add_sd <- peat_add_sd #* 0.0001 
  
  
  last_depth <- tail(bot_depth,1)
  # midpoints
  depths <- tail(Data$depth,n_data) + (Data$depth[2]-Data$depth[1])/2
  
  # variables for mcmc
  m0_upp <- 5e+3
  alpha_lim <- 1e-2
  l_depth <- tail(bot_depth,1)
  b_depth <- bot_depth[1]
  

  #### Functions ####
  
  age_depth <- function(x){
    approx(x = c(0,bot_depth),y = c(0,bot_ages),xout = x)$y
  }

  
  w <- function(param){ # Pondering function (it is a logistic regression with a fix)
    k <- epsilon1 / param[4] # grow parameter
    v <- param[5]
    ws <- 1 / (1 + exp(- k * (depths - v)))
    # old function, should do the same
    # b_1 = log(epsilon-1)/param[4] # param[4] is the delay
    # b_0 = -b_1*param[5] # param[5] is the A/C boundary
    # ws <- 1/ (1 + exp(-b_0-b_1*depths))
    return(ws)
  }

  
  # M_t <- function(m0, alpha, delta_t_x){
  #   # delta_t_x is the time diferential in a sample
  #   # Mi = (m0/alpha) * ( (delta_t_x) + (exp(-alpha * delta_t_x ) - 1 )/alpha  )
  #   Mi = (m0/alpha) * ( 1 - exp(-alpha * delta_t_x ))
  #   return(Mi)
  # }

  M_t <- function(m0, alpha, t_x, t_x_delta){
    # t_x is the bottom age
    # t_x_delta is the top age t(x-delta)
    #(note that t_x > t_x_delta)
    # Mi = (m0/alpha) * ( (t_x - t_x_delta) + (exp(-alpha * t_x ) - exp(-alpha * t_x_delta) )/alpha  )
    Mi = (m0/alpha) * ( exp(-alpha * t_x_delta) - exp(-alpha * t_x ))
    return(Mi)
  }
  ##### unused code ####
  
  # old version of the simu function and m0s function
  # m0s <- function(param){
  #   # param[1]: peat addition
  #   # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
  #   # param[4] is the delay (\Gamma / 2)  nota: revisar que este tomando el intervalo mencionado.
  #   # param[5] is the A/C boundary (\epsilon)
  #   # M_t(m0, alpha, t_x, t_x_delta)
  #   # t_x is the bottom
  #   # t_x_delta is the top age (note that t_x > t_x_delta)
  #   T_a = age_depth(param[5] - param[4] ) # this one is at the top
  #   T_c = age_depth( param[5] + param[4]  )
  # 
  #   # m0_a <- M_t(param[-c(1,2,3,4,5)],param[2], T_a,0)
  #   m0_a <- M_t(param[1],param[2], T_a,0)
  #   
  #   # Note that this middle section will decay half of the time as acretome and the other as catotelm
  #   m0_c <- .5 * (M_t(m0_a,param[2], T_c - T_a, 0) +  
  #                 M_t(m0_a,param[3], T_c - T_a, 0) )
  # 
  #   # m0 <- data.frame(m0 = param[-c(1,2,3,4,5)], m0a=m0_a , m0c=m0_c  )
  #   m0 <- c(param[1], m0_a , m0_c  )
  # 
  #   return(m0)
  # }
  # 

  # simu <- function(param){
  #   # param[1]: peat addition for the overall system
  #   # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
  #   # param[4] is the delay (\Gamma / 2)
  #   # param[5] is the A/C boundary (\epsilon)
  #   # param[-c(1,2,3,4,5)] are the section additions
  #   # Select samples per intervals
  #   # Note that this has to change when the age-depth model is moving
  #   # Calculate the times where the phases changes
  #   T_a = param[5] - param[4]
  #   T_c = param[5] + param[4]
  #   # Calculate the transicion function for each sample
  #   ws <- w(param)
  #   # Calculate which samples are in each phase
  #   indx_phase1 <- which(top_depth <= T_a)
  #   indx_phase2 <- which(top_depth >  T_a & 
  #                        top_depth <= T_c)
  #   indx_phase3 <- which(top_depth >  T_c)
  #   # get the m0s for each phase
  #   M_0s <- m0s(param)
  # 
  #   # Make Ts ages
  #   T_a = age_depth(top_depth[indx_phase2[1]]) #param[5] - param[4] )
  #   T_c = age_depth(top_depth[indx_phase3[1]]) #param[5] + param[4] )
  # 
  #   # get the simulation for each phase
  #   # mass in the first phase
  #   # mi_1 <- M_t(M_0s$m0[indx_phase1], param[2], bot_ages[indx_phase1], top_ages[indx_phase1])
  #   mi_1 <- M_t(M_0s[1], param[2], bot_ages[indx_phase1], top_ages[indx_phase1])
  #   # mass in the second phase
  #   # mi_2 <- ws[indx_phase2] * M_t(M_0s$m0a[indx_phase2], param[2], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a) +
  #   #   (1 - ws[indx_phase2]) * M_t(M_0s$m0a[indx_phase2], param[3], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a)
  #   mi_2 <- ws[indx_phase2] * M_t(M_0s[2], param[2], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a) +
  #     (1 - ws[indx_phase2]) * M_t(M_0s[2], param[3], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a)
  #   # mass in the last phase
  #   # mi_3 <- M_t(M_0s$m0c[indx_phase3],param[3], bot_ages[indx_phase3] - T_c, top_ages[indx_phase3] - T_c )
  #   mi_3 <- M_t(M_0s[3],param[3], bot_ages[indx_phase3] - T_c, top_ages[indx_phase3] - T_c )
  #   #
  #   z <- c(mi_1,mi_2,mi_3)
  #   if ( any( z < 0 ) ){
  #     print(z)
  #   }
  #   # debugging  
  #   # if (any(is.na(z))){
  #   #   print(which(is.na(z)))
  #   #   print(round(param))
  #   # }
  #   return(z)
  # }
  #####
  
  simu <- function(param){
    # param[1]: peat addition for the overall system
    # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
    # param[4] is the delay (\Gamma / 2)
    # param[5] is the A/C boundary (\epsilon)
    # param[-c(1,2,3,4,5)] are the section additions
    # Select samples per intervals
    # Note that this has to change when the age-depth model is moving
    # Calculate the times where the phases changes
    T_a = age_depth(param[5] - param[4])
    T_c = age_depth(param[5] + param[4])
    # Calculate the transicion function for each sample
    ws <- w(param)
    # Calculate which samples are in each phase
    new_bot <- sort(c(bot_ages,c(T_a,T_c)))
    new_top <- sort(c(top_ages,c(T_a,T_c)))
    loca_tmp_a <- which(new_top == T_a)
    loca_tmp_c <- which(new_top == T_c)

    ws <- append(ws, ws[loca_tmp_a] , after = loca_tmp_a)
    ws <- append(ws, ws[loca_tmp_c + 1] , after = loca_tmp_c + 2)
    
    # Pre-allocate the vector with the correct size (n_data + 2 for the loop, +1 for the initial value)
    m_i <- numeric(n_data + 2)

    # Initialize a variable to keep track of the last value
    m0 <- param[1]
    
    # Loop through without growing m_i
    for (i in 1:(n_data + 2)) {
      # each face influence
      # if(new_bot[i] < T_c){
        Ac_inf <-  M_t(m0, param[2],  new_bot[i] , new_top[i])
      # }else{
      #   Ac_inf <- 0 
      # }      
      if(new_bot[i] > T_a){
        Ca_inf <-  M_t(m0, param[3],  new_bot[i] - T_a, new_top[i] - T_a)  
      }else {
        Ca_inf <- 0 
      }
      
      # Sample creation
      mi <- ws[i] * Ac_inf + (1 - ws[i]) * Ca_inf
      # Update the last_value and assign it to the correct position in m_i
      # m0 <- mi/time_win
      m_i[i] <- mi
    }
    # Add the appended values to the original ones
    m_i[loca_tmp_a] <- m_i[loca_tmp_a] + m_i[loca_tmp_a + 1]
    m_i[loca_tmp_c + 2] <- m_i[loca_tmp_c + 2] + m_i[loca_tmp_c + 3]
    
    # Since index_c > index_a, remove index_c first to avoid shifting issues
    m_i <- m_i[-c(loca_tmp_a + 1, loca_tmp_c + 3)]
    return(m_i)
  }
  
  
  # Rcpp::sourceCpp("/Users/ma2060/GitHub/Bayesian_Carbon_Acc/simu_functions.cpp")
  
  log_tdist <- function(X, Mu, sigma, a, b) {
    sigma_squared = sigma^2 # Pre-compute
    
    term1 = ((2 * a + 1) / 2) * log(b + (X - Mu)^2 / (2 * sigma_squared))
    term2 = .5 * log(sigma_squared)
    
    results =  (term1 + term2)

    # Sum the results from all pairs
    total_result = -sum(results)
        
    return(total_result)
  }
  
  ll <- function(param){
    # z = simu(param, bot_depth, bot_ages, top_ages)
    # z = simu(param)
    # print(z)
    l = dnorm(simu(param), mean = Dat,sd = sd_ll, log = TRUE)
    # l = log_tdist(Dat,simu(param),sd_ll,a = 3,b = 4)
    # l = dt(z/sd_ll ,df = 3/4,log = T)
    # debugging
    # mask = !(l == -Inf)
    # l <- sum(l[mask])
    # if (any(is.infinite(l)))(l = -.Machine$double.xmax )
  
    return(sum(l))
  }  
  
  prior_d = function(param){
    # param[1]: peat addition
    d1 = dnorm(param[1],mean = peat_add_mean ,sd =  peat_add_sd , log=TRUE)
    # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
    # d1 = d1 + dnorm(param[2],mean = alpha_mean ,sd =  alpha_sd , log=TRUE) 
    # d1 = d1 + dnorm(param[3],mean = alpha2_mean  ,sd =  alpha2_sd , log=TRUE)
    # param[4] is the delay (\Gamma / 2)  
    d1 = d1 + dnorm(param[4], mean = delay_mean_prior, sd = delay_sd_prior, log =TRUE)
    # param[5] is the A/C boundary (\epsilon)
    d1 = d1 + dnorm(param[5], mean= ACboun_mean_prior, sd = ACboun_sd_prior, log =TRUE)
    # param[-c(1,2,3,4,5)]  all other m0s
    # d1 = d1 + sum(dnorm(param[-c(1,2,3,4,5)] , mean=param[1] , sd = .1 * peat_add_sd, log =TRUE) ) / n_data
    return( d1 )
  }
  
  obj <- function(param){
    ob <- loglike +  prior_d(param) 
    return( -ob )
  }
  
  ini = function(){
      m0 = abs(rnorm(1,mean = peat_add_mean ,sd =  peat_add_sd ))
      alphas = rev(sort(abs(rnorm(2,mean = c(alpha_mean,.01*alpha_mean) ,sd =  alpha_sd )  )))
      delay = abs(rnorm(1, mean= 3 * mean(diff(bot_depth)), sd = 3))
      AC_bou = runif(1, min= bot_depth[2], max = tail(bot_depth,1)  )
      # AC_bou = abs(rnorm(1, mean= mean(bot_depth), sd = 2))
      # allothers = abs(rnorm(n_data,mean = m0 ,sd = .01 * peat_add_sd ))
      # ini0 <- c(m0,alphas,delay,AC_bou,allothers)

    
    ini0 <- c(m0,alphas,delay,AC_bou)
    return(ini0)
  }


  
  supp <- function (param){
    # depths at which the phases changes
    d_a = param[5] - param[4]
    d_c = param[5] + param[4]
    
    upper = c( m0_upp, # peat addition (m0) top lim
               alpha_lim, # alpha acretelm 
               param[2], # alpha2 catotelm
               delay_mean_prior * 4, # delay should not be bigger than 4 times the prior one
               l_depth -  2 * param[4] #, # the AC boundery should be smaller than the last depth
               # rep(m0_upp,n_data)  # this line is to add multiple m0s
               ) 
    

    # if(!any(is.na(param)) ){
      if (all(param > 0) & d_a > b_depth & param[4] > b_depth & all(param < upper) & d_c < l_depth){
        loglike <<- ll(param)
        if(  !(is.na(loglike) || is.infinite(loglike))  ){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else{
        return(FALSE)
      }
    # }else{return(FALSE)}
  }
  
  #### Running Twalk ####
  message('Sale una orden con todo!!!')
  
  # load the twalk
  source_twalk('~/GitHub/Bayesian_Carbon_Acc/')
  
  # Check if state file exists
  state_exists <- file.exists( paste0(folder,'/',core,'/',core,'_state.csv'))
  
  if(state_exists & !newrun){
    # State file exists, load it
    state <- as.matrix(read.csv(paste0(folder,'/',core,'/',core,'_state.csv')))
    x0 = state[1,]
    x1 = state[2,]
    
  } else {
    # No state file or rerun flag is TRUE, initialize normally
    x0 = ini()
    while (!supp(x0)){  x0 = ini()      } 
    x1 = ini()
    while (!supp(x1))(    x1 = ini())
  }
  
  th = length(x1) * th
  cat(paste0('initial points are: \n'))
  cat(x0)
  cat('\n')
  cat(x1)
  cat('\n')
  twalk <- Runtwalk(dim = length(x1),  
                    Tr = s_size,thinning = th, burnin= burnin *th,
                    Obj = obj, Supp = supp,
                    x0 = x0, xp0 = x1)
  
  
  
  # Save state
  state = matrix(c(tail(twalk$output,1),tail(twalk$outputp,1)),
                 ncol=length(tail(twalk$output,1)),byrow = T)
  write.csv(state,paste0(folder,'/',core,'/',core,'_state.csv'),row.names = F)
  
  # Subsample
  sample <- tail(twalk$output,s_size)
  
  # Save results
  write.csv(sample, paste0(folder,'/',core,'/',core,'_output.csv'),row.names = FALSE)
  
  
  
  #### generate samples for plotting  ####
  # Initialize matrix to store assembles
  # transition
  w_mat <- matrix(NA, nrow = nrow(sample), ncol = n_data)
  
  # bulk density simulations
  Bd_mat <- matrix(NA, nrow = nrow(sample), ncol = n_data)
  
  # Loop through each row of c
  for(i in 1:nrow(sample)) {
    # Apply tau to breaks with this row of c as params
    Bd_row <- simu( sample[i,])
    # Store result in matrix
    Bd_mat[i,] <- Bd_row
    # Apply to w function 
    w_row <- w( sample[i,])
    # Store result in matrix
    w_mat[i,] <- w_row
  }  

  
  #### Plot the results ####
  message("Commencing the plotting process...")
  
  plots <- function(saveplot = TRUE){
    
    if(saveplot){pdf(paste0(folder,core,'/Taco_',core,'.pdf'))}
    layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,8,8,8,
                    1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,8,8,8,
                    1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,8,8,8,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6), 9, 18, byrow = TRUE))
    # Plot the energy
    {
      plot(twalk$Us[-1],col='gray40',type='l',yaxt="n",ylab='',xlab = 'iterations',main = " log objective" )
      # lines(twalk$Ups[-1],col='gray70')
    }
    
    # Plot M0
    {
      d <- density(sample[,1])
      plot(d,col='blue',type='l',yaxt="n",ylab='',xlab = 'g/m^2' , main = 'Global m0' )
    }
    # Plot alpha acretelm
    {
      d <- density(sample[,2])
      plot(d,col='blue',type='l',yaxt="n",ylab='',xlab = '1/yr',main = 'Decomp Acrotelm')
    }
    # Plot alpha catotelm
    {
      d <- density(sample[,3])
      plot(d,col='blue',type='l',yaxt="n",ylab='',xlab = '1/yr',main = 'Decomp Catotelm')
    }
    # Plot age-depth model
    {
      plot(bot_ages,bot_depth,type = 'l',col = 'gray10', ylim = rev(range(c(bot_depth, top_depth))),
           ylab='Depth',xlab='Age')
      abline(h=sample[,5],col=rgb(0,.2,1,.01))
      abline(h=sample[,5]-sample[,4],col=rgb(0,1,.2,.01))
      abline(h=sample[,5]+sample[,4],col=rgb(0,1,.2,.015))
    }
    # Transition and agreement
    {
      plot(1, 1, type = "n", xlim = range(Data$curr_mass_den) + c(-error_p, error_p)*range(Data$curr_mass_den)*1.01,
           ylim = rev(range(c(bot_depth, top_depth))), xlab = "g/m^2", ylab = "Depth (cm)",
           main = "")
      # plot(colMeans(tar_mt),tar_ages , type='l', xlab=ylabel, ylab=xlabel, ylim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1)) 
      
      # This plots the gray scale
      for(i in 1:ncol(Bd_mat) ){
        h <- hist(Bd_mat[,i], plot = FALSE,breaks=150)
        cols <- gray(1-h$counts/max(h$counts),alpha = .4)
        # Plot non-zero rects 
        rect(ybottom = bot_depth[i], ytop = top_depth[i],
               xleft = h$breaks[-1], #head(h$breaks, -1),
               xright = head(h$breaks, -1),# h$breaks[-1],
               col = cols, border = NA)
        points(mean(Bd_mat[,i]),bot_depth[i]-.5,col=rgb(0,1,.1,.4),pch=16,cex=.5)
      }
      
      # Adding simulation lines
      # for (i in 1:nrow(Data)) {
        # for (j in 1:nrow(Bd_mat)) {
        #   lines(c(Bd_mat[ j,i], Bd_mat[ j,i]), c(bot_depth[i], top_depth[i]), 
        #         col = rgb(0, 0, 1, 0.1))  # Blue with transparency
        #   
        # }
      # }
      # Adding the data squares
      for (i in 1:nrow(Data)) {
        rect(xleft = Data$curr_mass_den[i] - error_p * Data$curr_mass_den[i],
             xright = Data$curr_mass_den[i] + error_p * Data$curr_mass_den[i],
             ybottom = bot_depth[i],
             ytop = top_depth[i],
             border = NA,
             col = rgb(1,0,0,.3))  
      } 
  

      
      
      # Superimpose the second plot
      par(new=TRUE)
      plot(w_mat[1,], depths, col = rgb(0, 0, 0, alpha=0.01),axes=FALSE, 
           xlab="", ylab="",ylim = rev(range(c(bot_depth, top_depth))))
      # Adding simulation lines for each column in w_mat
      for (j in 2:nrow(w_mat)) {
        lines(w_mat[j,], depths, col = rgb(0, 0, 0, alpha=0.01))  # Blue with transparency
      }
      
      # Add the second x-axis at the top
      axis(3)  # 3 refers to the top position
      title("Transition", line=1, outer=TRUE, adj=0)  # Label for the new x-axis
    }
    
    # Plot delay
    {
      d <- density(sample[,4])
      plot(d,col='blue',type='l',yaxt="n",ylab='',xlab = 'cm' , main = 'Delay' )
    }
    # Plot A/C Boundery

    {
      d <- density(sample[,5])
      plot(d,col='blue',type='l',yaxt="n",ylab='',xlab = 'cm'  , main = 'A/C' )
    }

    
    if(saveplot){dev.off()}
  }
  
  plots()
  plots(saveplot = F)
  
  ###### Final message #######
  # message('The IAT and acceptance ratios is:\n')
  twalk$w_matrix <- w_mat
  twalk$Bd_mat <- Bd_mat
  twalk$Data <- Data
  twalk$M0 <- sample[,1]
  twalk$Deco_Acre <- sample[,2]
  twalk$Deco_Cato <- sample[,3]
  twalk$AC <- sample[,5]
  
  
  
  iat <- IAT(twalk,to=twalk$Tr)

  # 
  cat("\n================== IAT DIAGNOSTIC ==================\n")
  cat("IAT Value:", iat, "\n")
  if (iat < 5) {
    cat("Interpretation: The chain exhibits low correlation among successive samples.\n")
    cat("Recommendation: Current settings appear satisfactory.\n")
  } else {
    cat("Interpretation: The chain exhibits high correlation among successive samples.\n")
    cat("Recommendation: Consider increasing the thinning value and rerunning the chain.\n")
  }
  
  
  
  
  
  message('Servido joven!!!')
  
  return(twalk)

  }
   



#####
# Load twalk
source_twalk <- function(folder) {
  # Construct the path to the twalk.R file in the specified folder
  twalk_path <- file.path(folder, "twalk.R")
  # Check if the twalk.R file exists in the specified folder
  if (file.exists(twalk_path)) {
    source(twalk_path)
    message("Successfully loaded 'twalk.R' from", folder, "directory.\n")
  } else {
    warning("File twalk.R was not found in the specified folder.")
  }
}






# Nota: revisar unidades

# 186
# "Tasiusaq"

# 199
# "Wylde Lake bog"

# 140
# "PAT-HB-2010"
