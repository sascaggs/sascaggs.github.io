

fw_detection_sim = function(
    
  seed = NULL,
  K    = 5,   
  S    = 30,
  
  mu_alpha    = -2.5, # base rate
  sigma_alpha = 1,    # variation 
  sigma_u     = 1.5,  # consumer random effect
  sigma_v     = 1.5,  # resource random effect
  
  mu_p        = -2,   # detection intercept
  mu_beta     = 1,    # mean effect of consumer mass on preferences
  sigma_beta  = 0.5,  # variation in beta
  mu_gamma    = 1,    # mean effect of size ratio on deection 
  sigma_gamma = 0.3,  # variation in gamma
  
  mu_M        = 0,    # mean Mass
  sigma_M     = 1.5)  # variation in mass 

{
  
  # ---- Setup ----
  set.seed(seed)
  inv_logit = function(x) exp(x) / (1 + exp(x))
  
  # consumer mass 
  M = matrix(NA_real_, nrow = S, ncol = K)
  
  # generate varying base rates and detection probabilities, and body mass effects 
  alpha_k = rnorm(K, mu_alpha, sigma_alpha)
  beta_k = rnorm(K, mu_beta, sigma_beta)
  gamma_k = rnorm(K, mu_gamma, sigma_gamma)
  
  # ---- Outputs ----
  A_Mtrue = array(0, dim = c(S,S,K))
  A_Mobs  = array(0, dim = c(S,S,K))
  
  
  for(k in 1:K) {
    # Random effects
    u_i = rnorm(S, mean = 0, sd = sigma_u)
    v_j = rnorm(S, mean = 0, sd = sigma_v)
    
    # Mass 
    mass_k = exp(rnorm(S, mu_M, sigma_M))
    M[,k] = log10(mass_k)
    
    # Loop 
    for(i in 1:S) for(j in 1:S) {
      
      # true links 
      pi_ijk = inv_logit(alpha_k[k] + u_i[i] + v_j[j] + beta_k[k] * M[i,k])
      
      # detection probability 
      ratio_ijk = (M[i,k] - M[j,k])
      p_ijk = inv_logit(mu_p + gamma_k[k] * ratio_ijk)
      
      # true and observed links 
      A_Mtrue[i,j,k] = rbinom(1, 1, pi_ijk)
      A_Mobs[i,j,k] = rbinom(1, 1, p_ijk * pi_ijk)
      
    }
  }
  return(list(
    true = A_Mtrue, 
    obs  = A_Mobs, 
    params = data.frame(alpha_k, beta_k, gamma_k), 
    data = M
  ))
}

# run sims 
mu_beta_seq = c(0,1,3)
mu_gamma_seq = c(0,1,3)
params = expand.grid(mu_beta_seq, mu_gamma_seq)

output = list()
for( s in 1:nrow(params)) {
  beta = params[s,1]
  gamma = params[s,2]
  output[[s]] = fw_detection_sim(mu_beta = beta, mu_gamma = gamma)
}

# Create stan data 

make_stan_data <- function(sim) {
  
  A <- sim$obs   # array [S × S × K]
  M <- sim$data  # matrix [S × K]
  
  list(
    A      = A,
    N_cons = dim(A)[1],
    N_res  = dim(A)[2],
    K      = dim(A)[3],
    M      = M
  )
}

stan_datasets <- lapply(output, make_stan_data)

# Fit first model 
library(cmdstanr)
det_model = cmdstan_model("C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-20-food-web-measurement-2/detection_model_biomass_det.stan")


fit0_0 = det_model$sample(
  data = stan_datasets[[1]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.95
)

fit1_0 = det_model$sample(
  data = stan_datasets[[2]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.95
)

fit3_0 = det_model$sample(
  data = stan_datasets[[3]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.95
)

fit0_1 = det_model$sample(
  data = stan_datasets[[4]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.95
)

fit1_1 = det_model$sample(
  data = stan_datasets[[5]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.98
)

fit3_1 = det_model$sample(
  data = stan_datasets[[6]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.98
)

fit0_3 = det_model$sample(
  data = stan_datasets[[7]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.98
)

fit1_3 = det_model$sample(
  data = stan_datasets[[8]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.98
)

fit3_3 = det_model$sample(
  data = stan_datasets[[9]], 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4, 
  adapt_delta = 0.98
)


# get params for comparison 
param_df = bind_rows(lapply(output, `[[`, 'params'))
param_df$scenario = rep(1:9, each = 5)
param_df$scenario_params = rep(paste(params[,1], params[,2], sep = '_'), each=5)
param_df$k = 1:5
colnames(param_df)[c(1,2,3)] = c('.alpha_k','.beta_k','.gamma_k')

library(tidybayes)


fits = list(fit0_0, fit0_1, fit0_3, fit1_0, fit1_1, fit1_3, fit3_0, fit3_1, fit3_3)

summarise_fits = function(fit, scenario_id) {
  fit |> 
    spread_draws(beta_k[k], gamma_k[k]) |> 
    mean_qi() |> 
    mutate(scenario = scenario_id)
} 

post = map_df(
  .x = seq_along(fits), 
  .f = ~ summarise_fits(fits[[.x]], scenario_id = .x)
)

post_df = merge(param_df, post, by = c('scenario','k'))

post_df |> 
  ggplot(aes(x=beta_k, y=k)) + 
  geom_errorbar(aes(xmin = beta_k.lower, xmax = beta_k.upper),
                width = 0) + 
  geom_point() + 
  geom_point(aes(x=.beta_k, y=k), color='red') + 
  facet_wrap(~scenario_params)


params |> 
  mutate(label = TeX(paste0("$\\beta_k = ", Var1,"$")))

TeX(paste0("$\\beta_k = ", params$Var1,"$"))







#