#####################################################
#### Helveston et al. example                  ######
#####################################################

#########################
## data preparations 
#########################

library("logitr")
library(oeli) 
source("/Users/thanhthaobui/Documents/Abschlussarbeit Bauer/substract_choice_regressor_from_data.R")

## load data directly 
data <- load(file = "/Users/thanhthaobui/Documents/Abschlussarbeit Bauer/Data set USA China/helveston_data_wide_china.RData")

# Separate the dataset into Calibration and evaluation subset with 360 and 88 deciders in 2 sets

data_est <- data_wide[1:(360*15),]
data_val <- data_wide[(360*15 + 1):nrow(data_wide),]

#######################################
# estimate MNP model without mixing ----------------------------------------------------------
#######################################

make_Hb <- function(P) {
  diag(P)
}

make_fb <- function(P) {
  matrix(0, P, 1)
}

make_HO <- function(P_re) {
  if (P_re > 0) {
    ltho <- P_re * (P_re + 1) / 2
    if (P_re>1){
      ind <- fdiag(P_re)+1
    } else {
      ind <- 1
    }
    return(diag(ltho)[, ind, drop = FALSE])
  } else {
    return(matrix(0, 0, 0))
  }
}

make_fO <- function(P_re) {
  if (P_re > 0) {
    ltho <- P_re * (P_re + 1) / 2
    matrix(0, ltho, 1)
  } else {
    return(matrix(0, 0, 0))
  }
}

form <- as.formula(paste("choice ~",  paste(vars, collapse = " + ") ,  "| 0"))
P <- length(vars)

# re = NULL 
re = NULL
P_re <- length(re)

# load package 
library(Rprobit)

mod_hel <- mod_cl$new(
  alt = 3,
  Hb = make_Hb(P),
  fb = make_fb(P),
  HO = make_HO(P_re),
  fO = make_fO(P_re),
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mod_hel$fL[1] = 1
mod_hel$fL[4] = 1
mod_hel$fL[6] = 1

hel_obj_est <- ?setup_Rprobit(
  form = form,
  data_raw = data_est,
  re = NULL,
  id = "id",
  mod = mod_hel
)

hel_obj_val <- setup_Rprobit(
  form = form,
  data_raw = data_val,
  re = NULL,
  id = "id",
  mod = mod_hel
)

hel_obj_est$control$probit <- FALSE

theta_0 <- c(true_b)
hel_obj_est$set_theta(theta_0)

### Fit the MNP model to the data
hel_fit <- fit_Rprobit(hel_obj_est)
sum_est <- summary(hel_fit)

## Thao: Doing the t-test t check whether the parameters are statistically siginificant.
t <- sum_est$par$b/sum_est$par$b_sd

df <- 360*15 - length(theta_0)
p_value <- 2 * pt(-abs(t), df)
s <- rep(NULL,length(sum_est$par$b))

for (i in 1:length(sum_est$par$b)){
  if (p_value[i]<0.001) {
    s[i] <- "***"
  } else if (p_value[i]<0.01) {
    s[i] <- "**"
  } else if (p_value[i]<0.05) {
    s[i] <- "*"
  } else {
    s[i] <- ""
  }
}

result <- cbind(as.data.frame(hel_fit$vars), as.data.frame(paste0(round(sum_est$par$b,5), "(", round(sum_est$par$b_sd,5), ")")),s)
colnames(result) <- c("Variable", "MNP", "Significant level")

install.packages("knitr")
library(knitr)

kable(result, "latex")

#Thao: Fit the model to the evaluation dataset to check the prediction ability

hel_fit$data_raw <- hel_obj_val$data_raw
hel_fit$data <- hel_obj_val$data

sum_val <- summary(hel_fit)
ll_est <- sum(log(sum_est$predictions[,"choice_prob"]))
ll_val <- sum(log(sum_val$predictions[,"choice_prob"]))

# Calculate the confusion matrix for the validation dataset with rows as the true choice and columns as the predicted choice
tb_no <- table(sum_val$predictions[,4], c("A","B","C")[apply(sum_val$predictions[,1:3], 1, which.max)])
#The proportion of accuracy
sum(diag(tb_no))/sum(tb_no)
kable(tb_no, "latex")
######################################################################### now make opCost random 
re = "opCost"
P_re <- length(re)

mod_hel <- mod_cl$new(
  alt = 3,
  Hb = make_Hb(P),
  fb = make_fb(P),
  HO = make_HO(P_re),
  fO = make_fO(P_re),
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mod_hel$fL[1] = 1
mod_hel$fL[4] = 1
mod_hel$fL[6] = 1

hel_obj <- setup_Rprobit(
  form = form,
  data_raw = data_est,
  re = re,
  id = "id",
  mod = mod_hel
)

hel_obj$control$probit <- FALSE
theta_0 <- c(hel_fit$theta, rep(0.05, P_re))
hel_obj$set_theta(theta_0)

### Fit the Gaussian MMNP model to the calibration subset
hel_fit0 <- fit_Rprobit(hel_obj, cml_pair_type = 0)
sum_est0 <- summary(hel_fit0)

## Thao: Doing the t-test t check whether the parameters are statistically siginificant.
t <- sum_est0$par$b/sum_est0$par$b_sd

df <- 360*15 - length(theta_0)
p_value <- 2 * pt(-abs(t), df)
s <- rep(NULL,length(sum_est0$par$b))

for (i in 1:length(sum_est0$par$b)){
  if (p_value[i]<0.001) {
    s[i] <- "***"
  } else if (p_value[i]<0.01) {
    s[i] <- "**"
  } else if (p_value[i]<0.05) {
    s[i] <- "*"
  } else {
    s[i] <- ""
  }
}

result0 <- cbind(as.data.frame(hel_fit0$vars), as.data.frame(paste0(round(sum_est0$par$b,5), "(", round(sum_est0$par$b_sd,5), ")")),s)
colnames(result0) <- c("Variable", "MMNP", "Significant level")

kable(result0, "latex")

#Thao: Fit the model to the evaluation dataset to check the prediction ability

hel_fit0$data_raw <- hel_obj_val$data_raw
hel_fit0$data <- hel_obj_val$data

sum_val0 <- summary(hel_fit0)
sum(log(sum_est0$predictions[,"choice_prob"]))
sum(log(sum_val0$predictions[,"choice_prob"]))

# Calculate the confusion matrix for the validation dataset with rows as the true choice and columns as the predicted choice
tb_Gau <- table(sum_val0$predictions[,4], c("A","B","C")[apply(sum_val0$predictions[,1:3], 1, which.max)])
sum(diag(tb_Gau))/sum(tb_Gau)
kable(tb_Gau, "latex")
##########################
## Latent class model 
##########################

mode_LC <- mod_latclass_cl$new(
  Hb   = matrix(0,16,1),
  fb   = matrix(0,16,1),
  HO   = matrix(1,1,0),
  fO   =  matrix(1,1,1)*0.01,
  HL   = matrix(0,6,0),
  fL   = matrix(0,6,1),
  ordered = FALSE
)

mode_LC$fL[1] = 1
mode_LC$fL[4] = 1
mode_LC$fL[6] = 1

mode_LC$Hb[1,1] <- 1
mode_LC$fb <- as.matrix(hel_fit0$theta[1:16],ncol=1)
mode_LC$fb[1] <- 0 

mode_LC$num_class <- 3

theta_0LC <- c(-0.5,-0.3,-0.1,0,0)

hel_obj$mod <- mode_LC

hel_obj$control$probit <- FALSE
hel_obj$control$el = 1
hel_obj$control$control_weights$cml_pair_type <- 0
hel_obj$theta_0 <- theta_0LC
hel_obj$theta <- theta_0LC

### Fit the latent class model to the calibration subset:
hel_LC <- fit_LC_Rprobit(hel_obj, init_method = "theta")
sum_est_LC <- summary(hel_LC)

## Thao: Doing the t-test t check whether the parameters are statistically siginificant.
sum(exp(c(0, hel_LC$theta[4:5]))/sum(exp(c(0, hel_LC$theta[4:5])))*para1[1:3])
para_sd1 <- sqrt(diag(hel_LC$vv)) #Thao: Check with the prob about the last var is var of mixing
para1 <- hel_LC$theta
t1 <- para1/para_sd1
df1 <- 360*15 - length(hel_LC$theta)
p_value1 <- 2 * pt(-abs(t1), df1)

s1 <- rep(NULL,length(para1 ))
for (i in 1:length(para1)){
  if (p_value1[i]<0.001) {
    s1[i] <- "***"
  } else if (p_value1[i]<0.01) {
    s1[i] <- "**"
  } else if (p_value1[i]<0.05) {
    s1[i] <- "*"
  } else {
    s1[i] <- ""
  }
}

result1 <- cbind(as.data.frame(paste0(round(para1,5), "(", round(para_sd1,5), ")")),s1)
colnames(result1) <-  c("LC_MMNP","SL")

kable(result1, "latex")


#Thao: Fit the model to the evaluation dataset to check the prediction ability

hel_LC$data_raw <- hel_obj_val$data_raw
hel_LC$data <- hel_obj_val$data

sum_val_LC <- summary(hel_LC)
sum(log(sum_est_LC$predictions[,"choice_prob"]))
sum(log(sum_val_LC$predictions[,"choice_prob"]))

# Calculate the confusion matrix for the validation dataset with rows as the true choice and columns as the predicted choice
tb_LC <- table(sum_val_LC$predictions[,4], c("A","B","C")[apply(sum_val_LC$predictions[,1:3], 1, which.max)])
sum(diag(tb_LC))/sum(tb_LC)
kable(tb_LC, "latex")
tb_nonpara
#########################################
## non-parametric estimation 
#########################################
params = matrix(seq(from=-0.9,to=.5,by= 0.02),nrow=1)

mode_np <- mod_nonpara_cl$new(alt  = 3,
                              Hb   = matrix(0,16,1),
                              fb   = matrix(0,16,1),
                              HO   = matrix(0,1,0),
                              fO   =  matrix(1,1,1)*0.01,
                              HL   = matrix(0,6,0),
                              fL   = matrix(0,6,1),
                              ordered = FALSE
)

mode_np$fL[1] = 1
mode_np$fL[4] = 1
mode_np$fL[6] = 1

mode_np$Hb[1,1] <- 1
mode_np$fb <- as.matrix(hel_fit0$theta[1:16],ncol=1)
mode_np$fb[1] <- 0 

mode_np$set_grid_points(params)

hel_obj_np <- setup_Rprobit(
  form = form,
  data_raw = data_est,
  re = re,
  id = "id",
  mod = mode_np
)

hel_obj_np$control$probit <- FALSE
hel_obj_np$control$control_weights$cml_pair_type <- 0 # full pairwise CML function to be used.
hel_obj_np$mod <- mode_np

### Fit the non-parametric MMNP model to the data
Rpro_np <- fit_nonpara_Rprobit(hel_obj_np, init_method = "random", control_nlm = NULL, cml_pair_type = 0)

### Calculate some summary about the parameter value from the non-parametric model
probs_np <- Rpro_np$theta
plot(params,probs_np,col="red")
which.max(probs_np)
params[42]
probs_np[42]
sum(probs_np*params)

### Calculate the sum log-likelood for the calibration subset
data_est <- hel_obj_est$data
data_tr <- data_est$clone()
data_tr$data <- substract_choice_regressor_from_data(data_tr$data)
probs_new <- choice_probs_nonpara(data_tr, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1) #Thao?
dim(probs_new)
probs_new <- probs_new[, -dim(probs_new)[2]]
probs_choices <- probs_new %*% probs_np
dim(probs_new)
length(probs_np)
ll <- sum(log(probs_choices))
ll_in <- ll

### Calculate the confusion matrix for the non-parametric model.
### The steps are more complex than the latent class and Gaussian MMNP models as it is not available to use predict function
### to calculated the choice probability in the validation subset but need to create some other functions to compute the choice probabilities
data_val <- hel_obj_val$data
data_tr <- data_val$clone()
data_trA <- data_val$clone()
data_trB <- data_val$clone()
data_trC <- data_val$clone()

#Thao: Calculate the choice probablity of different choice by utilizing the function choice_probs_nonpara
for (i in 1:length(data_trA$data)){
  data_trA$data[[i]]$y <- as.factor(rep("A", 15))
}

for (i in 1:length(data_trB$data)){
  data_trB$data[[i]]$y <- as.factor(rep("B", 15))
}

for (i in 1:length(data_trC$data)){
  data_trC$data[[i]]$y <- as.factor(rep("C", 15))
}

data_tr$data <- substract_choice_regressor_from_data(data_tr$data)
data_trA$data <- substract_choice_regressor_from_data(data_trA$data)
#Thao: Cannot use the function for choice B and C because it just give the same result as subtract the choice A
# then need to do it manually below
#data_trB$data <- substract_choice_regressor_from_data(data_trB$data)
#data_trC$data <- substract_choice_regressor_from_data(data_trC$data)

for (i in 1:length(data_trB$data)){
  for (j in 1:15){
    data_trB$data[[i]]$X[[j]][3,] <- data_trB$data[[i]]$X[[j]][3,] - data_trB$data[[i]]$X[[j]][2,]
    data_trB$data[[i]]$X[[j]][1,] <- data_trB$data[[i]]$X[[j]][1,] - data_trB$data[[i]]$X[[j]][2,]
    data_trB$data[[i]]$X[[j]] <- data_trB$data[[i]]$X[[j]][-2,]
  }
}

for (i in 1:length(data_trC$data)){
  for (j in 1:15){
    data_trC$data[[i]]$X[[j]][2,] <- data_trC$data[[i]]$X[[j]][2,] - data_trC$data[[i]]$X[[j]][3,]
    data_trC$data[[i]]$X[[j]][1,] <- data_trC$data[[i]]$X[[j]][1,] - data_trC$data[[i]]$X[[j]][3,]
    data_trC$data[[i]]$X[[j]] <- data_trC$data[[i]]$X[[j]][-3,]
  }
}

probs_new <- choice_probs_nonpara(data_tr, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1)
probs_newA <- choice_probs_nonpara(data_trA, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1)
probs_newB <- choice_probs_nonpara(data_trB, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1)
probs_newC <- choice_probs_nonpara(data_trC, Rpro_np$mod, Rpro_np$control, cml_pair_type=-1)

probs_new <- probs_new[, -dim(probs_new)[2]]
probs_newA <- probs_newA[, -dim(probs_newA)[2]]
probs_newB <- probs_newB[, -dim(probs_newB)[2]]
probs_newC <- probs_newC[, -dim(probs_newC)[2]]


probs_choices <- probs_new %*% probs_np
probs_choicesA <- probs_newA %*% probs_np
probs_choicesB <- probs_newB %*% probs_np
probs_choicesC <- probs_newC %*% probs_np

probs_choicesA + probs_choicesB + probs_choicesC

### Thao: The value of the sum log-likelihood for the validation subset.
ll_out <- sum(log(probs_choices))

prob_nonprara <- as.data.frame(cbind(probs_choicesA, probs_choicesB, probs_choicesC))

colnames(prob_nonprara) <- c("A", "B", "C")

prob_nonprara$choice <- as.factor(c("A", "B", "C"))[apply(prob_nonprara, 1, which.max)]

### Thao: The confusion matrix
tb_nonpara <- table(sum_val_LC$predictions[,4], prob_nonprara[,4])

sum(diag(tb_nonpara))/sum(tb_nonpara)

kable(tb_nonpara, "latex")
