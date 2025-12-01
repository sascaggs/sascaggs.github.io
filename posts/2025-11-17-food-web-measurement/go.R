library(tidyverse)
library(latex2exp)

# ---- Setup ----
set.seed(777)
inv_logit = function(x) exp(x) / (1 + exp(x))

# ---- Parameters ----
K = 3   
S = 30
mu_alpha = -2
sigma_alpha = 1
sigma_u = 2
sigma_v = 1
mu_p = -1
sigma_p = 0.5

# generate varying base rates and detection probabilities 
alpha_k = rnorm(K, mu_alpha, sigma_alpha)
p_k = inv_logit(rnorm(K, mu_p, sigma_p))

# ---- Outputs ----
A_true = array(0, dim = c(S,S,K))
A_obs  = array(0, dim = c(S,S,K))

for(k in 1:K) {
  
  # Random effects
  u_i = rnorm(S, mean = 0, sd = sigma_u)
  v_j = rnorm(S, mean = 0, sd = sigma_v)
  
  # Loop 
  for(i in 1:S) for(j in 1:S) {
    
    pi_ijk = inv_logit(alpha_k[k] + u_i[i] + v_j[j])
    
    A_true[i,j,k] = rbinom(1, 1, pi_ijk)
    A_obs[i,j,k] = rbinom(1, 1, p_k[k] * pi_ijk)
  }
}

data_list = list(
  N_cons = S, 
  N_res  = S, 
  K      = K, 
  A      = A_obs
)

library(cmdstanr)

run = F
if(run) {
mod <- cmdstan_model("C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/detection_model.stan")

fit1 <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.95
)

fit1$save_output_files(dir = 'C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/fits')

saveRDS(fit1, 'C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/fits/detection_model1.RDS')
}

fit1 = readRDS('C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/fits/detection_model1.RDS')


library(tidybayes)

focal_pars = data.frame(k = factor(1:k), alpha_k, p_k)
focal_pars

alpha_k_draws = fit1 |> 
  spread_draws(mu_alpha, sigma_alpha, alpha_z[k]) |>
  group_by(k) |> 
  mutate(alpha_k_hat = mu_alpha + sigma_alpha * alpha_z) 

alpha_k_draws |> 
  ggplot(aes(x=alpha_k_hat, y=factor(k))) + 
  stat_pointinterval(color='black', .width = 0.9) + 
  labs(x = TeX("$\\hat{\\alpha}_k$"), y = 'k') + 
  geom_point(data=focal_pars, 
             aes(x = alpha_k, y=k), 
             color = 'red', pch=8, size=2, stroke=2)
  

fit |> 
  spread_draws(p_z[web_k], sigma_p, mu_p) |> 
  group_by(web_k) |> 
  mutate(prob_k = inv_logit(mu_p + sigma_p*p_z)) |>  
  ggplot(aes(x=prob_k, y=factor(web_k))) + 
  stat_slab()

