# ==============================================================================
# ================================= Packages ===================================
# ==============================================================================
#aCGH - Diffuse prior

library("wavethresh")
library("EbayesThresh")
library('waveslim')
library('changepoint')
library('ggplot2')
library('latex2exp')
library('coda')
library('ggtext')
library('bfw')

set.seed(123)
# ==============================================================================
# ============================= Preparing the data =============================
# ==============================================================================

data(Lai2005fig4)#data(Lai2005fig3)
y <- Lai2005fig4$GBM29#y <- Lai2005fig3$GBM31
data_size <- length(y)
aux_data_size <- 2^ceiling(log2(data_size ))
wavelet_adapted_size <- aux_data_size - data_size
serie <- seq_len(data_size)/data_size

# ============================ Global variables ================================

#Gibbs sampler: Setting down the size of the chain, burn-in and lag
burn <- 1000
lags <- 5
nchain <- 1000
BB <- burn + lags*nchain
n.verbose <- 1000
verbose <- TRUE

#Family choice
family_choice <- "Coiflets" #"DaubLeAsymm", "DaubExPhase", "Coiflets"

#Number of vanishing moments
number_vm <- 5

# ==============================================================================
# =========================== Important Functions ==============================
# ==============================================================================

vectorize_coeff <- function(data, vm, fc){
  wavelet1 <- wd(data, filter.number = vm, family = fc)
  vector1 <- c(1:length(data))
  vector1[1] <- wavelet1$C[2^(log2(length(data))+1)-1]
  vector1[2:length(data)] <- wavelet1$D
  return(vector1)
}

inverse_coeff <- function(coeff, vm, fc){
  zero <- rep(0, length(coeff))
  zero_wd <- wd(zero, filter.number = vm, family = fc)
  zero_wd$C[2^(log2(length(coeff))+1)-1] <- coeff[1]
  zero_wd$D <- coeff[2:length(coeff)]
  return(wr(zero_wd))    
}

bayesthresh_coeff <- function(coeff, vm, fc){
  zero <- rep(0, length(coeff))
  zero_wd <- wd(zero, filter.number = vm, family = fc)
  zero_wd$C[2^(log2(length(coeff))+1)-1] <- coeff[1]
  zero_wd$D <- coeff[2:length(coeff)]
  wavelet1 <- threshold(zero_wd,levels= 3:(nlevelsWT(zero_wd)-1),
                        type = "soft", policy = "BayesThresh",
                        by.level = FALSE, dev = madmad,
                        boundary = FALSE, alpha = 0.5,
                        beta = 1, C1 = NA, C2 = NA, 
                        C1.start = 100)
  zero[1] <- wavelet1$C[2^(log2(length(coeff))+1)-1]
  zero[2:length(coeff)] <- wavelet1$D
  return(zero)
}

# ==============================================================================
# ============================== Gibbs Sampler =================================
# ==============================================================================

mu <- matrix(nrow=BB,ncol = 2)
phi <- matrix(nrow=BB,ncol = 2)
coef_wavelet <- matrix(nrow=BB,ncol=aux_data_size)
alpha <- matrix(nrow=BB,ncol=data_size) 
inv_coeff <- matrix(nrow=BB,ncol=aux_data_size)

prior_mean1 <- quantile(y, probs = 0.25)
prior_mean2 <- quantile(y, probs = 0.75)
prior_variance1 <- var(y)
prior_variance2 <- var(y)
prior_alpha_gamma <- 0.01
prior_beta_gamma <- 0.01

#Defining the start values
mu[1,]<-c(prior_mean1,prior_mean2)
phi[1,]<- c(1/prior_variance1, 1/prior_variance2)
coef_wavelet[1,] <- vector("double", aux_data_size) 
alpha[1,] <- rep_len(1, data_size)  
prob0start <- (1-alpha[1,])*dnorm(y, mean = mu[1,1], 
                                  sd = 1/sqrt(phi[1,1]))
prob1start <- alpha[1,]*dnorm(y, mean = mu[1,2], 
                              sd = 1/sqrt(phi[1,2]))
probvalstart <- prob1start/(prob0start + prob1start)

pred_label <- rbinom(data_size, 1, prob = probvalstart)
inv_coeff[1,] <- inverse_coeff(coef_wavelet[1,], vm = number_vm, fc = family_choice) 
lat_var <- rnorm(aux_data_size, mean = inv_coeff[1,], 1)

# ================================== Gibbs =====================================

for(ii in 2:BB){
  #Creating the auxiliary variables for calculating the conditional posterior 
  n1 <- sum(pred_label)
  n0 <- data_size - n1
  y1 <- y[pred_label==1]
  y0 <- y[pred_label==0]
  sy1 <- sum(y1)
  sy0 <- sum(y0)
  
# ====================== Generating mean.lower =================================
  var_mu_lower <- 1/(n0*phi[ii-1,1]+ 1/prior_variance1)
  mean_mu_lower <- var_mu_lower *(sy0*phi[ii-1,1] +
                                    prior_mean1/prior_variance1)
  mu[ii,1] <- rnorm(n = 1, mean = mean_mu_lower, sd = sqrt(var_mu_lower))
  
# ====================== Generating prec.lower =================================  
  alpha_gamma_lower <- prior_alpha_gamma + n0/2
  beta_gamma_lower <- prior_beta_gamma + sum((y0 - mu[ii,1])^2)/2
  phi[ii,1] <- rgamma(1, shape = alpha_gamma_lower, rate = beta_gamma_lower)
  
# ====================== Generating mean.higher ================================  
  var_mu_higher <- 1/(n1*phi[ii-1,2]+ 1/prior_variance2)
  mean_mu_higher <- var_mu_higher *(sy1*phi[ii-1,2] +
                                      prior_mean2/prior_variance2)
  mu[ii,2] <- rnorm(n = 1, mean = mean_mu_higher, sd = sqrt(var_mu_higher))
  
# ====================== Generating prec.higher ================================ 
  alpha_gamma_higher <- prior_alpha_gamma + n1/2
  beta_gamma_higher <- prior_beta_gamma + sum((y1 - mu[ii,2])^2)/2
  phi[ii,2] <- rgamma(1, shape = alpha_gamma_higher, rate = beta_gamma_higher)
  
  
# ============================== Constraint ====================================
  #To dealing with label switching, we stablish that the pairs (\mu_k, \phi_k) 
  #are ordered under the constraint \mu_1 < \mu_2 and \phi_1 < \phi_2
  
  if(mu[ii,2] < mu[ii,1]){
    mu[ii,]<- mu[ii,2:1]
    phi[ii,]<- phi[ii,2:1]
  }
  
# ====================== Generating pred_label =================================
  prob0 <- (1-alpha[ii-1,])*dnorm(x = y, mean = mu[ii,1],
                                  sd = 1/sqrt(phi[ii,1]))
  prob1 <- alpha[ii-1,]*dnorm(x = y, mean = mu[ii,2],
                              sd = 1/sqrt(phi[ii,2]))
  probval <- prob1/(prob0 + prob1)
  
  pred_label <- rbinom(n = data_size, size = 1, prob = probval)
  
# ========================= Generating lat_var =================================
  pred_label_aux <- c(pred_label,rev(pred_label)[1:wavelet_adapted_size])
  U <- runif(n = aux_data_size)
  upper_bound <- (pnorm(q = 0, mean = inv_coeff[ii-1,], sd = 1))^(1-pred_label_aux)
  lower_bound <- (pnorm(q = 0, mean = inv_coeff[ii-1,], sd = 1))*pred_label_aux
  prob_l <- (U*(upper_bound - lower_bound))+lower_bound
  lat_var <- qnorm(prob_l, mean = inv_coeff[ii-1,], sd = 1)
  
# ============================ Generating theta ================================  
  coef_wavelet[ii,] <- rnorm(n = aux_data_size, 
                             mean = vectorize_coeff(data = lat_var, 
                                                    vm = number_vm, 
                                                    fc = family_choice), 
                             sd = 1)
  
  #Regularizing the coefficients
  coef_wavelet[ii,] <- bayesthresh_coeff(coeff = coef_wavelet[ii,], 
                                         vm = number_vm, fc = family_choice)
  
# =========================== Generating inv_coeff =============================
  # inv_coeff[ii,] <- qnorm(alpha_behavior, mean = 0, sd = 1)
  inv_coeff[ii,] <- inverse_coeff(coef_wavelet[ii,], vm = number_vm, 
                                  fc = family_choice)
  
# ========================== Generating alpha ==================================
  #alpha[ii,] <- alpha_behavior
  alpha[ii,] <- pnorm(inv_coeff[ii,], 0, 1)[1:data_size]
  
  if(!ii%%n.verbose & verbose)
    print(paste0(ii, " iterations of ", BB, "."))
}


# =========================== Burn-in and lag ==================================
#Creating the matrices with burn-in and lag

N <- seq(burn+1, BB, lags)

matrix_mu <- mu[N,]
matrix_phi <- phi[N,]
matrix_alpha <- alpha[N,]

#estimated values by median
estimated_mu_median <- apply(matrix_mu,2,median)
estimated_phi_median <- apply(matrix_phi,2,median)
estimated_alpha_median <- apply(matrix_alpha,2,median)

#estimated values by mean
estimated_mu_mean <- apply(matrix_mu,2,mean)
estimated_phi_mean <- apply(matrix_phi,2,mean)
estimated_alpha_mean <- apply(matrix_alpha,2,mean)


plot(estimated_alpha_median,type="l",col="black", lwd = 1, 
     ylab=expression(alpha[t]))
lines(seq(1,data_size,1), estimated_alpha_mean, lty=3, lwd = 3,col = "blue")


plot(y, ylab="", pch = 1, col = "blue")

