
age_depth <- function(x){
  approx(x = c(0,bot_depth),y = c(0,bot_ages),xout = x)$y
}


w <- function(param){ # Pondering function (it is a logistic regression with a fix)
  k <- epsilon1 / param[4] # grow parameter
  v <- param[5]
  ws <- 1/(1 + exp(- k * (depths-v)) )
  # old function, should do the same
  # b_1 = log(epsilon-1)/param[4] # param[4] is the delay
  # b_0 = -b_1*param[5] # param[5] is the A/C boundary
  # ws <- 1/ (1 + exp(-b_0-b_1*depths))
  return(ws)
}

M_t <- function(m0, alpha, t_x, t_x_delta){
  # t_x is the bottom age
  # t_x_delta is the top age t(x-delta) (note that t_x > t_x_delta)
  # Mi = (m0/alpha) * ( (t_x - t_x_delta) + (exp(-alpha * t_x ) - exp(-alpha * t_x_delta) )/alpha  )
  Mi = (-m0/alpha) * (exp(-alpha * t_x ) - exp(-alpha * t_x_delta) )
  return(Mi)
}


len_data <- length(bot_depth)
m0s <- function(param){
  # param[1]: peat addition
  # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
  # param[4] is the delay (\Gamma / 2)  nota: revisar que este tomando el intervalo mencionado.
  # param[5] is the A/C boundary (\epsilon)
  # M_t(m0, alpha, t_x, t_x_delta)
  # t_x is the bottom
  # t_x_delta is the top age (note that t_x > t_x_delta)
  T_a = age_depth(param[5] - param[4] ) # this one is at the top
  T_c = age_depth( param[5] + param[4]  )

  # m0_a <- M_t(param[-c(1,2,3,4,5)],param[2], T_a,0)
  m0_a <- M_t(rep(param[1],len_data),param[2], T_a,0)
  m0_c <- M_t(m0_a,param[3], T_c - T_a, 0)
  # m0 <- data.frame(m0 = param[-c(1,2,3,4,5)], m0a=m0_a , m0c=m0_c  )
  m0 <- data.frame(m0 = rep(param[1],len_data), m0a=m0_a , m0c=m0_c  )

  return(m0)
}

simu <- function(param){
  # param[1]: peat addition for the overall system
  # param[2] and param[3]: are decomposition  (acrotelm - catotelm)
  # param[4] is the delay (\Gamma / 2)
  # param[5] is the A/C boundary (\epsilon)
  # param[-c(1,2,3,4,5)] are the section additions
  # Select samples per intervals
  # Note that this has to change when the age-depth model is moving
  # Calculate the times where the phases changes
  T_a = param[5] - param[4]
  T_c = param[5] + param[4]
  # Calculate the transicion function for each sample
  ws <- w(param)
  # Calculate which samples are in each phase
  indx_phase1 <- which(bot_depth <= T_a)
  indx_phase2 <- which(bot_depth > T_a & bot_depth <= T_c)
  indx_phase3 <- which(bot_depth > T_c)
  # get the m0s for each phase
  M_0s <- m0s(param)
  # print(M_0s)
  # Make Ts ages
  T_a = age_depth(param[5] - param[4] )
  T_c = age_depth(param[5] + param[4] )
  # get the simulation for each phase
  # mass in the first phase
  mi_1 <- M_t(M_0s$m0[indx_phase1], param[2], bot_ages[indx_phase1], top_ages[indx_phase1])
  # mass in the second phase
  mi_2 <- ws[indx_phase2] * M_t(M_0s$m0a[indx_phase2], param[2], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a) +
    (1 - ws[indx_phase2]) * M_t(M_0s$m0a[indx_phase2], param[3], bot_ages[indx_phase2] - T_a, top_ages[indx_phase2] - T_a)
  # mass in the last phase
  mi_3 <- M_t(M_0s$m0c[indx_phase3],param[3], bot_ages[indx_phase3] - T_c, top_ages[indx_phase3] - T_c )
  #
  z <- c(mi_1,mi_2,mi_3)

  return(z)
}