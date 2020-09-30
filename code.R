library(tidyverse)
library(rstanarm)
country.code <- 'PA_subset'
## NHPP ##
# ref_dat: data for points in reference set, with relevant covariates
# nonref_dat: data for points in non-reference set, with relevant covariates
ref_pts <- ref_dat %>%
 dplyr::select(elev, lc, protected,urban,
               road_dist, city_dist)%>%
 mutate(lc = factor(lc)) 

nonref_pts <- nonref_dat %>%
 dplyr::select(elev, lc, 
               protected,urban,road_dist, city_dist)

Xobs <- ref_pts %>%
 mutate(cat = "obs")
Xref <- nonref_pts %>%
 mutate(cat = "ref")
temp <- rbind(Xobs, Xref) %>%
 dplyr::mutate_if(is.numeric, scale)
Xref <-  model.matrix(~., temp %>%
                       filter(cat == "obs") %>% 
                       dplyr::select(-cat))
Xnonref <-  model.matrix(~., temp %>%
                          filter(cat == "ref") %>% 
                          dplyr::select(-cat))

nhpp <- function(I, burnin, b_prop_sd, X, Xq, Xstack,wj, thin_amt){
 sd_beta <- 100
 #INITS
 p <- ncol(X)
 N <- nrow(X)
 Q <- nrow(Xq)
 beta <- rep(0,p)
 BETAS <-matrix(NA, nrow = I, ncol = p)
 beta.acc <- matrix(0, nrow = I+burnin, ncol = p)
 
 for (i in 1:(I+burnin)){
  # update beta
  for (s in sample(1:p)){
   beta_star <- beta
   beta_star[s] <- rnorm(1, beta[s],b_prop_sd[s])
   log_lhr <- (sum(X[,s] * beta_star[s]) - sum(wj* exp(Xstack %*% beta_star))) -
    (sum(X[,s] *beta[s]) - sum(wj* exp(Xstack %*% beta))) + 
    dnorm(beta_star[s], 0, sd_beta, log = T) - 
    dnorm(beta[s], 0, sd_beta, log = T)
   if (log_lhr > log(runif(1))){
    beta[s] <- beta_star[s]
    beta.acc[i, s] <- 1
   }
  }
  
  
  if (i > burnin){
   store <- i-burnin
   BETAS[store, ] <- beta
  }
 }
 thin <- seq(1,I, thin_amt)
 return(list(BETAS[thin,], colMeans(beta.acc)))
}

b_prop_sd <- c(0.125, 0.125, 0.175, 0.15, 0.35, 0.25, 0.15, 0.15, 0.15, 0.15)

X <- as.matrix(Xref)
Xq <-as.matrix(Xnonref)
Xstack<-rbind(X,Xq)
p <- ncol(X)
N_nonref <- nrow(Xq)
N_ref <- nrow(X)
N <- N_ref + N_nonref
D <- readRDS(paste0(country.code, "_area.Rda"))
wj <- rep(D/N, N)
I <- 10000
burnin <- 10000
thin <- 10

## NNGP
# ref_nbs: list of neighbors for each point in reference set
# ref_dists: list of distances for each point in reference set to its neighbors 
# nonref_nbs: list of neighbors for each point not in reference set, where neighbors come from reference set
# nonref_dists: list of distances for each point in non-reference set to its neighbors
nbs <- c(ref_nbs, nonref_nbs)
dists <- c(ref_dists, nonref_dists)
N_ref <- length(ref_nbs)
N_nonref <- length(nonref_nbs)
N <- N_ref + N_nonref

# for each point in reference set, create list of points where it is in their neighbor sit
u_ls <- list()
for (i in 1:N_ref){
 temp <- c()
 for(j in 1:N){
  if (i %in% nbs[[j]]){
   temp <- c(temp, j)
  }
 }
 if(length(temp) == 0){
  u_ls[[i]]<- NA
 } else{
  u_ls[[i]] <- temp
 }
}

# conditional means
get_omega_ai <- function(dists, phi, N){
 omega <- list()
 ai <- list()
 for(i in 1:N){
  m <- nrow(dists[[i]])
  if(length(m) > 0){
   C22 <- exp(-(1/phi)*dists[[i]][2:m, 2:m])
   C12 <- exp(-(1/phi)*dists[[i]][1, 2:m])
   omega[[i]] <- 1 - C12 %*% solve(C22) %*% C12
   ai[[i]] <- C12 %*% solve(C22)
  } else{
   omega[[i]] <- 1
   ai[[i]] <- 0
  }
 }
 return(list(omega, ai))
}

run_lgcp <- function(I, burnin, X, Xq, N_ref, N_nonref, N, b_prop_sd, wj, phi, thin_amt) {
 # set up 
 omega_ai <- get_omega_ai(dists, phi, N)
 omega <- omega_ai[[1]]
 ai <- omega_ai[[2]]
 Xstack <- rbind(X, Xq)
 
 ## hyper params
 sd_beta <- 100
 z_ref_prop_sd <- 0.2
 z_nonref_prop_sd <- 0.2
 s2_prop_sd <- 0.075
 a0 <- 2; b0 <- 1
 
 ## inits ##
 beta <-  rep(0,p)
 s2 <- 0.5
 Z <- rnorm(N, 0, sqrt(s2))
 XB <- X%*%beta
 XXB <- Xstack %*% beta
 
 ## store ##
 BETAS <-matrix(NA, nrow = I, ncol = p)
 S2s <- rep(NA, nrow = I)
 ZS <- matrix(NA, nrow = I,ncol = N)
 Z.acc <- matrix(0,nrow = I+burnin, ncol = N)
 beta.acc <- matrix(0, nrow = I+burnin, ncol = p)
 s2.acc <- 0
 
 
 ## GIBBS SAMPLER ##
 for (g in 1:(I+burnin)){
  # update z in ref set
  z_temp <- rnorm(N, 0, z_ref_prop_sd) + Z
  for (i in 1:N_ref){
   Z_star = Z
   Z_star[i] = z_temp[i]
   # first part of prior
   if (i == 1){ 
    mtemp = 0;
    vartemp = s2;
   } else{
    mtemp = ai[[i]]%*%Z[nbs[[i]]]
    vartemp = s2 * omega[[i]]
   }
   
   # second part of prior
   u_nbs <- u_ls[[i]]
   if (length(u_nbs) == 1){
    mtemp2_star = ai[[i]] * Z_star[nbs[[u_nbs]]]
    mtemp2 = ai[[i]] * Z[nbs[[u_nbs]]]
   } else{
    mtemp2_star <- apply(matrix(1:length(u_nbs), nrow = 1), 2, function(x){ai[[u_nbs[x]]]%*%Z_star[nbs[[u_nbs[x]]]]})
    mtemp2 <- apply(matrix(1:length(u_nbs), nrow = 1), 2, function(x){ai[[u_nbs[x]]]%*%Z[nbs[[u_nbs[x]]]]})
   }
   vartemp2 <- unlist(omega[u_nbs]) * s2
   
   
   if (length(u_nbs) == 0){
    lhr <- -wj[i]*exp(XXB[i] - s2/2 +Z_star[i]) - (-wj*exp(XXB[i] - s2/2 + Z[i]))+
     Z_star[i] - Z[i] +
     -0.5*(((Z_star[i] - mtemp)^2/vartemp ) -((Z[i] - mtemp)^2/vartemp) )
   } else{
    lhr <- -wj[i]*(exp(XXB[i] - s2/2 + Z_star[i]) - exp(XXB[i] - s2/2 + Z[i]))+
     Z_star[i] - Z[i] +
     -0.5*(((Z_star[i] - mtemp)^2/vartemp ) -((Z[i] - mtemp)^2/vartemp) ) +
     -0.5*( sum(((Z[u_nbs] - mtemp2_star)^2/vartemp2) - ((Z[u_nbs] - mtemp2)^2/vartemp2))) 
   }
   if(lhr > log(runif(1))){
    Z[i] = Z_star[i]
    Z.acc[g,i] <- 1
   }
  }
  
  # update z outside ref set
  Z_temp = rnorm(N_nonref, 0, z_nonref_prop_sd) + Z[(1+N_ref):N]
  for (i in (1+N_ref):N){
   z_star = Z_temp[i - N_ref]
   mtemp = ai[[i]] %*% Z[ nbs[[i]]]
   vartemp = s2*omega[[i]]
   
   lhr <- -wj[i]*(exp(XXB[i] - s2/2 + z_star) - exp(XXB[i] - s2/2 + Z[i]) ) +
    -0.5/vartemp*( (z_star - mtemp)^2 - (Z[i] - mtemp)^2)
   if(lhr>log(runif(1))){
    Z[i] <- z_star
    Z.acc[g,i] <- 1
   }
  }
  
  # update s2
  ls2 = log(s2)
  ls2_prop = rnorm(1,0, s2_prop_sd) + ls2
  s2_prop = exp(ls2_prop)
  if(s2_prop >0){
   mmtemp <- rep(NA, N)
   mmtemp[1] <- 0#Z[1]
   for (i in 2:N){
    mmtemp[i] <- as.numeric(ai[[i]] %*% Z[nbs[[i]]])
   }
   lh_prev <- sum(XB - s2/2 + Z[1:N_ref]) - sum(wj*exp(XXB - s2/2 + Z))  
   lh_prop <- sum(XB - s2_prop/2 + Z[1:N_ref])- sum(wj*exp(XXB - s2_prop/2 + Z))
   lhr <-  lh_prop - lh_prev +
    -0.5*N*(ls2_prop - ls2) +
    -0.5*(sum((Z - mmtemp)^2 / unlist(omega)))*(1/s2_prop - 1/s2)+
    dgamma(1/s2_prop, shape = a0, scale = b0, log = T) - dgamma(1/s2, shape = a0, scale = b0, log = T)+
    ls2_prop - ls2
   
   if (lhr > log(runif(1))){
    s2 <- s2_prop
    s2.acc <- s2.acc + 1
   }
   
  }
  
  
  ## update beta ##
  for (s in sample(1:p)){
   beta_star <- beta
   beta_star[s] <-  rnorm(1,beta[s] ,b_prop_sd[s])
   
   log_lhr <- (sum(X%*%beta_star + Z[1:N_ref] - s2/2) - sum(wj*exp(Xstack%*%beta_star + Z - s2/2)))-
    (sum(X %*% beta + Z[1:N_ref] - s2/2) - sum(wj*exp(Xstack%*%beta + Z - s2/2))) +
    dnorm(beta_star[s], 0, sd_beta, log = T) -
    dnorm(beta[s], 0, sd_beta, log = T)
   if (log_lhr > log(runif(1))){
    beta <- beta_star
    beta.acc[g, s] <- 1
   }
   
  }
  XB <- X%*%beta
  XXB <- Xstack%*%beta
  
  ## store ##
  if (g > burnin){
   store <- g - burnin
   BETAS[store,] <- beta
   S2s[store] <- s2
   ZS[store,] <- Z
  }
 }
 beta.acc.ratios <- colMeans(beta.acc)
 Z.acc.ratios <- colMeans(Z.acc)
 s2.acc.ratio <- s2.acc/(burnin+I)
 acc.ratios <- c(beta.acc.ratios, s2.acc.ratio, Z.acc.ratios)
 thin <- seq(1, I, thin_amt)
 return(list(BETAS[thin,], S2s[thin], ZS[thin,], acc.ratios))
}

## Effort ##
# effort_df: data containing covariates and metrics for effort
# Dist: distance matrix for points in effort_df
# w: choice of effort

Xw <- model.matrix(~prop_weekend + prop_morn +  lc, effort_df %>%
                    mutate(road_dist = scale(road_dist),
                           city_dist = scale(city_dist),
                           lc = relevel(factor(lc), "forest")))

N <- nrow(Xw)
# psi_s2 <- spatial variance 
# psi_phi <- spatial decay param

calc_Sigma <- function(N, D, s2, phi){
 Sigma <- s2*exp(-(1/phi)*D)
 return(Sigma)
}


# GIBBS SAMPLER #
run_effort <- function(I, burnin, thin_amt, w, Xw, N, D, phi){
 # hyperparms 
 gamma_s2 <- 100
 a <- b <- 2
 phi_prop_s2 <- 0.1
 s2_prop_s2 <- 0.1
 
 # set-up/inits
 D_inv <- solve(D)
 q <- ncol(Xw)
 I_q <- diag(q)
 
 psi <- rnorm(N)
 psi_s2 <- 0.5
 psi_phi <- phi
 psi_Corr <- exp(-(1/phi)*D) #calc_Sigma(N,D, psi_s2,psi_phi)
 psi_Corr_inv <- chol2inv(chol(psi_Corr))
 psi_Sigma <- psi_s2 * psi_Corr
 gamma <- rep(0, q)
 log_w <- log(w + 0.01)
 # STORE #
 S2 <-  rep(NA, I)
 PSI <- matrix(NA, nrow = I, ncol = N)
 GAMMA <- matrix(NA, nrow = I, ncol = q)
 s2.acc <- rep(0, I+burnin)
 
 # sampler #
 for (i in 1:(I+burnin)){
  prec_psi <- (1/psi_s2)*psi_Corr_inv
  gamma_v <- solve((1/gamma_s2)*I_q + t(Xw)%*%prec_psi%*%(Xw))
  gamma_m <- gamma_v %*% t(Xw)%*%prec_psi%*%(log_w - psi)
  gamma <- rmvn(1, gamma_m, gamma_v)
  Xw_gamma <- Xw %*% t(gamma)

  nu <- t(rmvn(1, rep(0,N), psi_Sigma))
  u <- runif(1)
  
  # previous log-like
  logy <- dmvn(log_w, Xw_gamma + psi, psi_Sigma, log = T) + log(u)
  
  # init. proposal + proposal bracket
  theta <- runif(1, min = 0, max = 2*pi)
  theta_min <- theta-2*pi
  theta_max <- theta
  
  flag = 0
  while(flag == 0){
   # propose nu and f
   psi_star <- psi*cos(theta) + nu*sin(theta)
   if (dmvn(log_w, Xw_gamma + psi_star, psi_Sigma, log = T) > logy){
    psi <- psi_star
    flag <- 1 # accept nu_star
    break
   } else{
    # shrink bracket
    if (theta < 0){
     theta_min <- theta
    } else{
     theta_max <- theta
    }
    theta <- runif(1, theta_min, theta_max)
   }
  }
  
  # update s2
  s2_prop <- rnorm(1)*s2_prop_s2 + psi_s2
  if(s2_prop > 0){
   psi_Sigma_prop <- s2_prop*psi_Corr
   lhr <- dmvn(log_w, Xw_gamma + psi, psi_Sigma_prop, log = T) - 
    dmvn(log_w, Xw_gamma + psi, psi_Sigma, log = T) +
    dgamma(1/s2_prop, a, b, log = T) - dgamma(1/psi_s2, a, b, log = T)
   if(lhr > log(runif(1))){
    psi_s2 <- s2_prop
    psi_Sigma <- psi_Sigma_prop
    s2.acc[i] <- 1
   }
  }
  
 
  ## store ## 
  if (i > burnin){
   store <- i - burnin
   GAMMA[store,] <- gamma
   S2[store] <- psi_s2
   PSI[store,] <- psi
  }
 }
 
 thin <- seq(1, I, thin_amt)
 return(list(GAMMA[thin,], S2[thin], PSI[thin,]))
}


## Poisson regression for species abundances ##
# Y: dataframe of species counts and coordinates
# w: vector of effort metric for the observed locations
# ref_dat: data for points in reference set, with relevant covariates
X <- ref_dat %>%
  dplyr::select(elev, lc, protected,urban,
                road_dist, city_dist)%>%
  mutate(lc = factor(lc),
         elev = scale(elev))
X_effort <- cbind(X,w = w)
spec_ls <- c("Junco hyemalis", "Zenaida macroura","Agelaius phoeniceus")

glm1_ls <- glm2_ls <- pred1_ls <- pred2_ls <-  list()
for (i in 1:length(spec_ls)){
  sp <- spec_ls[[i]]
  G <- 5000
    glm1_ls[[i]] <- stan_glm(Y[,sp] ~., data =  X, family = poisson, prior = normal(0, 5, autoscale = F), 
                             iter = G, chains = 2)
    glm2_ls[[i]] <- stan_glm(Y[,sp] ~., data =  X_effort, family = poisson, prior = normal(0, 5, autoscale = F), 
                             iter = G, chains = 2)
}

## if splitting data into train/test set:
# post_beta: posterior samples for coefficients for NNGP
# post_zs2: posterior samples for spatial variance for NNGP
# post_z: posterior samples for spatial random effects for NNGP
# post_gamma: posterior samples for coefficients for effort
# post_ws2: posterior samples for spatial variance for effort 
# post_delta: posterior samples for shared process coefficient for effort
# test_df: dataframe for covariates for held out data
# test_coords: coordinates of locations for held out data
# test_nbs: list of neighbors of held out locations, where neighbors come from the reference set
# test_dist: list of distance matrices for the held out locations and their neighbors
# Y: dataframe of species counts and coordinates for train set
# w: vector of effort metric for the observed locations in train set
# ref_dat: data for points in reference set, with relevant covariates, for locations in train set
# Y_test: species counts for held out locations

X_test <- test_df %>%
  dplyr::mutate_if(is.numeric, scale)  %>%
  dplyr::select(elev, lc, protected, urban, road_dist, city_dist)
Xtest <- model.matrix(~., X_test)
Xw_test <-  model.matrix(~prop_weekend +prop_morn  + lc, test_df %>%
                           mutate(road_dist = scale(road_dist),
                                  city_dist = scale(city_dist)))
N_test <- nrow(X_test)

nngp_lambdas <- matrix(nrow = I, ncol = nrow(X))
for (i in 1:I){
  nngp_lambdas[i,] <- exp(Xtest %*% post_beta[i,] + post_z[i,1:N_ref] - post_z_s2[i]/2)
}
post_lambdas <- colMeans(nngp_lambdas)


omega_ai <- get_omega_ai(test_dists, phi, N_test)
omega <- omega_ai[[1]]
ai <- omega_ai[[2]]
z_krig_test <-  matrix(NA, nrow = I, ncol = N_test)
for (i in 1:I){
  cond_m <- rep(NA, N_test)
  cond_v <- post_z_s2[i]*unlist(omega)
  for (k in 1:N_test){
    cond_m[k] <- ai[[k]]%*%post_z[i, test_nbs[[k]]]
  }
  z_krig_test[i,] <- rnorm(N_test, cond_m, cond_v)
}

# prediction for test point effort
w_pred_test <- matrix(NA, nrow = I, ncol = N_test)
for (i in 1:I){
  m <- Xw_test%*%post_gamma[i,] + post_delta[i]*z_krig_test[i,]
  w_pred_test[i,] <- rnorm(N_test, m, sqrt(post_w_s2[i]))
}
w_pred <- exp(colMeans(w_pred_test))

X <- ref_dat %>%
  dplyr::select(elev, lc, protected,urban,
                road_dist, city_dist)%>%
  mutate(lc = factor(lc),
         elev = scale(elev))
X_effort <- cbind(X,w = w)

X_pred <-test_df %>%
  dplyr::select(elev, lc, protected,urban,
                road_dist, city_dist)%>%
  mutate(lc = factor(lc),
         elev = scale(elev))
X_pred_effort <- cbind(X_pred, w = w_pred)

### fit ###
spec_ls <- c("Junco hyemalis", "Zenaida macroura","Agelaius phoeniceus")
glm1_ls <- glm2_ls <- pred1_ls <- pred2_ls <-  list()

for (i in 1:length(spec_ls)){
  sp <- spec_ls[[i]]
  I <- 5000
  glm1_ls[[i]] <- stan_glm(Y[,sp] ~., data =  X, family = poisson, prior = normal(0, 5, autoscale = F), 
                           iter = I, chains = 2)
  glm2_ls[[i]] <- stan_glm(Y[,sp] ~., data =  X_effort, family = poisson, prior = normal(0, 5, autoscale = F), 
                           iter = I, chains = 2)
  pred1_ls[[i]] <- pred1 <- posterior_predict(glm1_ls[[i]], newdata = X_pred)
  pred2_ls[[i]] <- pred2 <- posterior_predict(glm2_ls[[i]], newdata = X_pred_effort)
  
}




