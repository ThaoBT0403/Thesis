library(Rprobit)
# Input the R file "redraw_choices_from_data"
source("/Users/thanhthaobui/Documents/Abschlussarbeit Bauer/redraw_choices_from_data.R") 

# Input the R file "build_par_from_mod"
source("/Users/thanhthaobui/Documents/Abschlussarbeit Bauer/build_par_from_mod.R")

# Input the R file "substract_choice_regressor_from_data"
source("/Users/thanhthaobui/Documents/Abschlussarbeit Bauer/substract_choice_regressor_from_data.R")

# Model structure for non-parametric analyses
#Thao: Estimate one beta, mixing variance fixed to 0, 3 choices - variance is diag(3)
mode_np <- mod_nonpara_cl$new(alt  = 3,
                              Hb   = matrix(1,1,1),
                              fb   = matrix(0,1,1),
                              HO   = matrix(0,0,0), 
                              fO   =  matrix(0,0,0),
                              HL   = matrix(0,6,0),
                              fL   = matrix(0,6,1),
                              ordered = FALSE
)

mode_np$fL[1] = 1
mode_np$fL[4] = 1
mode_np$fL[6] = 1

params = matrix(seq(from=-4,to=6,by= 0.1),nrow = 1)
length(params)
mode_np$set_grid_points(params)



# Define the object from which the data is drawn
#Thao: Gaussian mixed MNP with 1 beta estimated, variance of mixing estimated, 3 options with variance fixed
mod0 <- mod_cl$new(alt  = 3,
                   Hb   = matrix(1,1,1),
                   fb   = matrix(0,1,1),
                   HO   = matrix(1,1,1),
                   fO   =  matrix(0,1,1),
                   HL   = matrix(0,6,0),
                   fL   = matrix(0,6,1),
                   ordered = FALSE
)

mod0$fL[1] = 1
mod0$fL[4] = 1
mod0$fL[6] = 1

theta0 <- c(1,.5)

# Set up Rprobit object
#Thao: random a dataset with N= 500 decision makers, each has 5 choice occasions, each occasion has 3 alts
N= 500
control_simul = list(Tp = rep(5,N))

form <- choice ~ V1 | 0
re <- c("V1")

Rprobit_obj <- setup_Rprobit(
  form = form,
  mod = mod0,
  re = re,
  seed = 17,
  theta_0 = theta0,
  control = control_simul 
)
Rprobit_obj$data$data[1] 
Rprobit_obj$data_raw$df[c(1:5), ] #Thao: data is just comparable with data_raw but is subtracted the value of the first choice.
dim(Rprobit_obj$data_raw$df)


Rprobit_obj$control$probit <- FALSE
Rprobit_obj$control$control_weights$cml_pair_type <- 0 # Full pairwise CML function used.



# Redraw from different params
mu0 = 2
sd0 = 1
thetas <- matrix(rnorm(2*N,mean= mu0,sd = sd0),ncol =2)
thetas[,2] <- 0
seed = 1

Rprobit_rd <- redraw_choices_from_data( Rprobit_obj, seed = seed, thetas) #Thao: New y is put into dataset to reflex the effect of random error term
# The 2 forms of the dataset used to estimate models
Rprobit_rd$data$data[[1]]
Rprobit_rd$data_raw$df[1:10,]
dim(Rprobit_rd$data_raw$df)

# Number of replications
#M = 10

# Thao: Grid points to determine the CDF
cdf_knots <- seq(from = -4, to = 6, by = 0.01) 
length(cdf_knots)


#######################################################
### Normal estimation with normal mixed parameter 
#######################################################
#Thao: cr_MMNP <- eval_MMNP( Rprobit_rd, mod0, cdf_knots)
#Thao: put the data to the variables in the function to test the function.

Rprobit_rd$mod <- mod0
Rprobit_rd$theta_0  = c(1,.5)
Rprobit_obj <- Rprobit_rd
mode <- mod0

######################

eval_MMNP <- function( Rprobit_obj, mode, cdf_knots){
  
  N <- Rprobit_obj$data$N
  Nest <- ceiling(N/2)
  Nval <- N-Nest 
  
  # out of sample data
  Rpro_out <- Rprobit_obj$clone(deep = TRUE)
  
  # split data 
  data_est <- Rpro_out$data
  data_val <- data_est$clone(deep=TRUE) #Thao: When without clone and change then the origianl will be affected
  data_est$set_data(data_est$data[c(1:Nest)])
  data_val$set_data(data_val$data[c((Nest+1):N)])
  
  data_est_raw <- Rpro_out$data_raw
  data_val_raw <- data_est_raw$clone(deep=TRUE)
  Nmax <- sum(data_est_raw$Tp[1:Nest])
  data_est_raw$set_df(data_est_raw$df[1:Nmax,])
  data_val_raw$set_df(data_val_raw$df[c((Nmax+1):sum(data_val_raw$Tp)),]) 
  
  Rpro_out$data <- data_est
  Rpro_out$data_raw <- data_est_raw
  Rpro <- fit_Rprobit(Rpro_out, init_method = "random", control_nlm = NULL)
  summary(Rpro)
  beta = Rpro$theta[1]
  seconds <- as.numeric(as.difftime(Rpro$info$estimation_time, units = "secs"))
  
  # Calculate CDF for estimation on points cdf_knots
  par <- build_par_from_mod(Rpro$theta, Rpro$mod) #Thao: the difference between Rpro$theta: estimated one, ^2 equal Rpro$theta[2] and Rpro$theta0: initial value, 
  mu <- par$b 
  sd <- sqrt(par$Omega) 
  
  cdf_est <- cdf_knots*0 
  cdf_est <- pnorm(cdf_knots,mean=mu,sd=sd) 
  
  sum_est <- summary(Rpro)
  
  Rpro$data <- data_val
  Rpro$data_raw <- data_val_raw
  
  sum_val <- summary(Rpro) 
  
  # evaluate criterion 
  ll <- sum_est$ll
  ll_in <- sum(log(sum_est$predictions[,"choice_prob"]))
  ll_out <- sum(log(sum_val$predictions[,"choice_prob"]))
  
  # return results
  structure(
    list(
      ll = ll,
      ll_out = ll_out, 
      cdf = cdf_est,
      beta = beta,
      seconds = seconds
    )
  )
}

Rprobit_rd$mod <- mod0
Rprobit_rd$theta_0  = c(1,.5)


cr_MMNP <- eval_MMNP( Rprobit_rd, mod0, cdf_knots)


#######################################################
### estimation with LC MMNP model  
#######################################################

form <- choice ~ V1 | 0
re <- c("V1")
mod_LC <- mod_latclass_cl$new(  #Thao: estimate one beta, mixing variance is fixed at 0.5^2, variance of error is fixed diag(3)
  Hb   = as.matrix(c(1),ncol=1),
  fb   = as.matrix(c(0),ncol=1),
  HO   = matrix(0,1,0),
  fO   = matrix(.5,1,1),
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE,
)
mod_LC$fL[c(1,4,6),1]=1
mod_LC$num_class <- 3
theta_0 = c(0,1,2,0,0) #Thao: first 3: beta in 3 classes, the last 2: their weight with the first weight fixed at 0

Rprobit_rd$theta_0 <- theta_0
Rprobit_rd$mod <- mod_LC

#Thao:
Rprobit_obj <- Rprobit_rd
Rprobit_rd$data$data
Rprobit_rd$data_raw$df

###################################

eval_LC_MMNP <- function( Rprobit_obj, mod_LC, cdf_knots){
  
  N <- Rprobit_obj$data$N
  Nest <- ceiling(N/2)
  Nval <- N-Nest 
  
  # out of sample data
  Rpro_out <- Rprobit_obj$clone(deep = TRUE)
  
  # split data 
  data_est <- Rpro_out$data
  data_val <- data_est$clone(deep=TRUE)
  data_est$set_data(data_est$data[c(1:Nest)])
  data_val$set_data(data_val$data[c((Nest+1):N)])
  
  data_est_raw <- Rpro_out$data_raw #Thao: The difference between $data_raw and $data
  data_val_raw <- data_est_raw$clone(deep=TRUE)
  Nmax <- sum(data_est_raw$Tp[1:Nest])
  data_est_raw$set_df(data_est_raw$df[1:Nmax,])
  data_val_raw$set_df(data_val_raw$df[c((Nmax+1):sum(data_val_raw$Tp)),]) 
  
  Rpro_out$data <- data_est
  Rpro_out$data_raw <- data_est_raw
  Rpro_out$mod <- mod_LC
  
  Rpro_out$control$probit <- FALSE
  Rpro_out$control$el = 1
  Rpro_out$control$control_weights$cml_pair_type <- 0
  
  Rpro <- fit_LC_Rprobit(Rpro_out,init_method = "kmeans")
  summary(Rpro)
  
  # calculate time 
  seconds <- as.numeric(as.difftime(Rpro$info$estimation_time, units = "secs"))
  
  # calculate CDF for estimation on points cdf_knots
  tot_params <- Rpro$mod$lthb + Rpro$mod$lthO +Rpro$mod$lthL
  
  num_class <- Rpro$mod$num_class
  pi_params <- Rpro$theta[(num_class*tot_params+1):length(Rpro$theta)]
  eta <- c(0,pi_params)
  pi <- exp(eta)
  pi <- pi/sum(pi)
  beta_esti <- sum(pi*Rpro$theta[1:(tot_params*num_class)])
  cdf_est <- cdf_knots*0 
  #Thao: jb <- 1
  for (jb in 1:num_class){
    theta_jb <- Rpro$theta[(jb-1)*tot_params + c(1:tot_params)]
    par <- build_par_from_mod(theta_jb, Rpro$mod)
    mu <- par$b 
    sd <- sqrt(par$Omega)
    
    cdf_est <- cdf_est+ pnorm(cdf_knots,mean=mu,sd=sd)*pi[jb]
    
  } 
  sum_est <- summary(Rpro)
  
  Rpro$data <- data_val
  Rpro$data_raw <- data_val_raw
  
  sum_val <- summary(Rpro)
  
  # evaluate criterion 
  ll <- sum_est$ll
  ll_in <- sum(log(sum_est$predictions[,"choice_prob"]))
  ll_out <- sum(log(sum_val$predictions[,"choice_prob"]))
  
  # return results
  structure(
    list(
      ll = ll,
      ll_out = ll_out, 
      cdf = cdf_est,
      beta = beta_esti,
      seconds = seconds
    )
  )
}


cr_LC_MMNP <- eval_LC_MMNP( Rprobit_rd, mod_LC, cdf_knots)



#########################################
# estimation with non-parametric method  
#########################################
#Thao: 
Rprobit_obj <- Rprobit_rd
###############################

eval_nonpara <- function( Rprobit_obj, mode_np, cdf_knots, cml_pair_type = 0){
  
  N <- Rprobit_obj$data$N
  Nest <- ceiling(N/2)
  Nval <- N-Nest 
  
  Rprobit_obj$mod <- mode_np
  
  # out of sample data
  Rpro_out <- Rprobit_obj$clone(deep = TRUE)
  
  # split data 
  data_est <- Rpro_out$data
  data_val <- data_est$clone(deep=TRUE)
  data_est$set_data(data_est$data[c(1:Nest)])
  data_val$set_data(data_val$data[c((Nest+1):N)])
  
  data_est_raw <- Rpro_out$data_raw
  data_val_raw <- data_est_raw$clone(deep=TRUE)
  Nmax <- sum(data_est_raw$Tp[1:Nest])
  data_est_raw$set_df(data_est_raw$df[1:Nmax,])
  data_val_raw$set_df(data_val_raw$df[c((Nmax+1):sum(data_val_raw$Tp)),]) 
  
  Rpro_out$data <- data_est
  Rpro_out$data_raw <- data_est_raw
  
  Rpro_out$mod <- mode_np
  Rpro_out$theta_0 <- rep(1/mode_np$num_grid_points,mode_np$num_grid_points)
  Rpro_np <- fit_nonpara_Rprobit(Rpro_out, init_method = "theta", control_nlm = NULL, cml_pair_type = cml_pair_type)
  probs_np =  Rpro_np$theta # calculate_grid_probs(mode_np,Rpro_np$theta) #Thao: Prob of different b fixed given = params
  
  seconds <- as.numeric(as.difftime(Rpro_np$info$estimation_time, units = "secs"))
  
  betas <- Rpro_np$mod$params
  beta_esti <- sum(betas*probs_np)
  cdf_est <- cdf_knots*0 #Thao: jj <- 1 length(cdf_knots)
  for (jj in 1:length(cdf_knots)){
    cdf_est[jj] <- sum(probs_np[(betas< cdf_knots[jj])]) #Thao: sum the prob estimated from the mod of all points in cdf_knots. betas: theta0,prob_np: prob estimated of theta0, cdf_knots: points to calculate cdf
  }
  
  data_tr <- data_est$clone()
  data_tr$data <- substract_choice_regressor_from_data(data_tr$data)
  probs_new <- choice_probs_nonpara(data_tr, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1) #Thao?
  dim(probs_new)
  probs_new <- probs_new[, -dim(probs_new)[2]]
  probs_choices <- probs_new %*% probs_np
  ll <- sum(log(probs_choices))
  ll_in <- ll
  
  data_tr <- data_val$clone()
  data_tr$data <- substract_choice_regressor_from_data(data_tr$data)
  probs_new <- choice_probs_nonpara(data_tr, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1)
  probs_new <- probs_new[, -dim(probs_new)[2]]
  dim(probs_new)
  probs_choices <- probs_new %*% probs_np
  ll_out <- sum(log(probs_choices))
  
  # return results
  structure(
    list(
      ll = ll,
      ll_out = ll_out, 
      cdf = cdf_est,
      beta = beta_esti,
      seconds = seconds
    )
  )
}

# test the proc
cr_nonpara <- eval_nonpara( Rprobit_rd, mode_np, cdf_knots)


##############################
### now start the engine ####
##############################

cdf_true <- cdf_knots*0 
cdf_true <- pnorm(cdf_knots,mean = mu0, sd = sd0)

library(ggplot2)

ggplot_data <- rbind(
  data.frame(x = cdf_knots, y = cdf_true, method = "true"),
  data.frame(x = cdf_knots, y = cr_MMNP$cdf, method = "MMNP"),
  data.frame(x = cdf_knots, y = cr_LC_MMNP$cdf, method = "LC_freq"),
  data.frame(x = cdf_knots, y = cr_nonpara$cdf, method = "nonpara")
)

ggplot(ggplot_data, aes(x = x, y = y, color = method)) +
  geom_line() 


simul_nonpara <- function(form, mod0, mode, mod_LC, mode_np, mod_LC_B, re, cdf_knots, cdf_true, latent_class, latent_class2, cml_pair_type = 0, control_simul, thetas, M= 100){
  # define outputs
  timing = matrix(0,M,3)
  ll_ins = timing
  ll_outs = timing 
  cdf_diff = timing
  beta_M <- timing
  
  cdf_MMNP = matrix(0,M,length(cdf_knots))
  cdf_LC_MMNP = matrix(0,M,length(cdf_knots))
  cdf_np = matrix(0,M,length(cdf_knots))
  
  
  labs = c('MMNP','LC_MMNP','nonpara')
  colnames(ll_ins) <- labs
  colnames(ll_outs) <- labs
  colnames(cdf_diff) <- labs
  colnames(timing) <- labs
  colnames(beta_M) <- labs
  
  # run replications 
  pb <- txtProgressBar(max = M,style = 3)
  progress <- function(n) setTxtProgressBar(pb,n) 
  opts <- list(progress = progress) 
  
  for (m in 1:M){ #Thao: m <- 1
    setTxtProgressBar(pb, m)
    # generate data 
    Rprobit_obj <- setup_Rprobit(form = form,mod = mod0,re = re,seed = m,theta_0 = thetas[1,],control = control_simul) 
    Rprobit_obj$control$probit <- FALSE
    Rprobit_obj$control$control_weights$cml_pair_type <- 0 # full pairwise CML function to be used.
    Rprobit_rd <- redraw_choices_from_data( Rprobit_obj, seed = 2*m, thetas)
    
    # now calculate the various methods 
    # MMNP 
    Rprobit_rd$mod <- mod0
    Rprobit_rd$theta_0  = c(1,.5)
    
    #time_1 <- Sys.time()
    cr_MMNP <- eval_MMNP( Rprobit_rd, mod0, cdf_knots)
    #time_2 <- Sys.time()
    timing[m,1] <- cr_MMNP$seconds
    #beta
    beta_M[m,1] <- cr_MMNP$beta
    
    # LC 
    theta_0 = c(0,1,2,0,0)
    
    Rprobit_rd$theta_0 <- theta_0
    Rprobit_rd$mod <- mod_LC
    #time_1 <- Sys.time()
    cr_LC_MMNP <- eval_LC_MMNP( Rprobit_rd, mod_LC, cdf_knots)
    #time_2 <- Sys.time()
    timing[m,2] <- cr_LC_MMNP$seconds
    #beta
    beta_M[m,2] <- cr_LC_MMNP$beta
    
    # non parametric 
    #time_1 <- Sys.time()
    cr_nonpara <- eval_nonpara( Rprobit_rd, mode_np, cdf_knots, cml_pair_type)
    #time_2 <- Sys.time()
    timing[m,3] <- cr_nonpara$seconds
    #beta
    beta_M[m,3] <- cr_nonpara$beta
    
    # collect all results 
    ll_ins[m,] <- c(cr_MMNP$ll,cr_LC_MMNP$ll,cr_nonpara$ll)
    ll_outs[m,] <- c(cr_MMNP$ll_out,cr_LC_MMNP$ll_out,cr_nonpara$ll_out)
    cdf_diff[m,] <- c(sum(abs(cdf_true-cr_MMNP$cdf)),sum(abs(cdf_true-cr_LC_MMNP$cdf)),sum(abs(cdf_true-cr_nonpara$cdf)))
    
    # cdfs 
    cdf_MMNP[m,] <- cr_MMNP$cdf
    cdf_LC_MMNP[m,] <- cr_LC_MMNP$cdf
    cdf_np[m,] <- cr_nonpara$cdf
    
  } 
  
  close(pb) 
  
  # set up output 
  cdfs <- list(cdf_MMNP, cdf_LC_MMNP, cdf_np)
  
  out <- list(timing = timing, ll_ins = ll_ins, ll_outs = ll_outs, cdf_diff = cdf_diff,cdfs =cdfs, beta_M = beta_M)
  return(out)
} 

N=50
control_simul = list(Tp = rep(5,N))

M =200
cml_pair_type = 0

out_norm <- simul_nonpara(form, mod0, mode, mod_LC, mode_np, mod_LC_B, re, cdf_knots, cdf_true, latent_class, latent_class2, cml_pair_type, control_simul, thetas, M)


apply(out_norm$timing,2,mean)

apply(out_norm$ll_outs,2,mean)
apply(0.01*out_norm$cdf_diff,2,mean)

# Overview about the distribution of estimated beta
summary(out_norm$beta_M)

# deal with estimated cdfs 
cdf_MMNP <- out_norm$cdfs[[1]]
cdf_LC_MMNP <- out_norm$cdfs[[2]]
cdf_nonpara <- out_norm$cdfs[[3]]

#Thao: Draw the graphs after 200 simulation times.

# calculate mean cdfs for the three methods 
cdf_MMNP_mean <- apply(cdf_MMNP,2,mean)
cdf_LC_MMNP_mean <- apply(cdf_LC_MMNP,2,mean)
cdf_nonpara_mean <- apply(cdf_nonpara,2,mean)

library(ggplot2)

ggplot_data <- rbind(
  data.frame(x = cdf_knots, y = cdf_true, label = "true"),
  data.frame(x = cdf_knots, y = cdf_MMNP_mean, label = "MMNP"),
  data.frame(x = cdf_knots, y = cdf_LC_MMNP_mean, label = "LC_freq"),
  data.frame(x = cdf_knots, y = cdf_nonpara_mean, label = "nonpara")
)

ggplot(ggplot_data, aes(x = x, y = y, color = label)) +
  geom_line() +
  scale_color_manual(values = c("MMNP" = "green", "true" = "red", "LC_freq" = "blue", "nonpara" = "deepskyblue")) +
  labs(title = "The mixing distribution from three mixed models in average and its true distribution",
       x = "Grid points",
       y = "CDF") +
  theme_minimal()
########################################
######### Thao: Draw the distribution of mixing after 200 simulation times 
######### From MMNP models
install.packages("tidyverse")
library(tidyverse)

data_wide <- data.frame(x = cdf_knots)
for (i in 1:M) {
  data_wide[[paste0("cdf", i)]] <- cdf_MMNP[i,]  # Random cumulative sums
}

# Step 2: Reshape to long format
data_long <- data_wide %>%
  pivot_longer(cols = starts_with("cdf"), names_to = "cdf", values_to = "value") %>%
  mutate(label = "MMNP", color = "green", size = 0.3)  # Assign same label and color to 200 variables

# Step 3: Add the extra variable with a different label and color
extra_data <- data.frame(x = cdf_knots, value = cdf_true, 
                         cdf = "cdf_true", label = "true", color = "red", size = 0.7)

# Combine the two datasets
combined_data <- bind_rows(data_long, extra_data)

# Step 4: Plot the data
ggplot(combined_data, aes(x = x, y = value, group = cdf)) +
  geom_line(aes(color = label, size = size)) +  # Line plot with different colors and labels
  geom_line(data = extra_data, aes(x = x, y = value, color = label, size = size)) +
  scale_color_manual(values = c("MMNP" = "green", "true" = "red")) +
  scale_size_continuous(range = c(0.3, 0.7)) + 
  labs(title = "The mixing distribution from the MMNP model and its true distribution",
       x = "Grid points",
       y = "CDF") +
  theme_minimal()

######### From LC-MMNP models

data_wide <- data.frame(x = cdf_knots)
for (i in 1:M) {
  data_wide[[paste0("cdf", i)]] <- cdf_LC_MMNP[i,]  # Random cumulative sums
}

# Step 2: Reshape to long format
data_long <- data_wide %>%
  pivot_longer(cols = starts_with("cdf"), names_to = "cdf", values_to = "value") %>%
  mutate(label = "LC-MMNP", color = "blue", size = 0.3)  # Assign same label and color to 200 variables

# Step 3: Add the extra variable with a different label and color
extra_data <- data.frame(x = cdf_knots, value = cdf_true, 
                         cdf = "cdf_true", label = "true", color = "red", size = 0.7)

# Combine the two datasets
combined_data <- bind_rows(data_long, extra_data)

# Step 4: Plot the data
ggplot(combined_data, aes(x = x, y = value, group = cdf)) +
  geom_line(aes(color = label, size = size)) +  # Line plot with different colors and labels
  geom_line(data = extra_data, aes(x = x, y = value, color = label, size = size)) +
  scale_color_manual(values = c("LC-MMNP" = "blue", "true" = "red")) +
  scale_size_continuous(range = c(0.3, 0.7)) + 
  labs(title = "The mixing distribution from the LC-MMNP model and its true distribution",
       x = "Grid points",
       y = "CDF") +
  theme_minimal()

######### From non-parametric MMNP models

data_wide <- data.frame(x = cdf_knots)
for (i in 1:M) {
  data_wide[[paste0("cdf", i)]] <- cdf_nonpara[i,]  # Random cumulative sums
}

# Step 2: Reshape to long format
data_long <- data_wide %>%
  pivot_longer(cols = starts_with("cdf"), names_to = "cdf", values_to = "value") %>%
  mutate(label = "non-parametric MMNP", color = "deepskyblue", size = 0.3)  # Assign same label and color to 200 variables

# Step 3: Add the extra variable with a different label and color
extra_data <- data.frame(x = cdf_knots, value = cdf_true, 
                         cdf = "cdf_true", label = "true", color = "red", size = 0.7)

# Combine the two datasets
combined_data <- bind_rows(data_long, extra_data)

# Step 4: Plot the data
ggplot(combined_data, aes(x = x, y = value, group = cdf)) +
  geom_line(aes(color = label, size = size)) +  # Line plot with different colors and labels
  geom_line(data = extra_data, aes(x = x, y = value, color = label, size = size)) +
  scale_color_manual(values = c("non-parametric MMNP" = "deepskyblue", "true" = "red")) +
  scale_size_continuous(range = c(0.3, 0.7)) + 
  labs(title = "The mixing distribution from the non-parametric MMNP model and its true distribution",
       x = "Grid points",
       y = "CDF") +
  theme_minimal()
  

