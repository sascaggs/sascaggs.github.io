library(tidyverse)
library(latex2exp)

# graphs
graph_theme = theme(
  # panels
  panel.background = element_rect(
    color  = 'black', 
    fill   = '#ffffffff', 
    size   = 1 ), 
  panel.grid = element_blank( ), 
  panel.spacing = unit(15, 'pt'), 
  # axes
  axis.ticks  = element_line(
    color = 'black', 
    size  = 0.5 ), 
  axis.ticks.length = unit(2, 'mm'), 
  axis.text = element_text(color='black'), 
  # strips
  strip.background = element_rect(
    color = '#ffffffff', 
    fill  = '#ffffffff',), 
  strip.text = element_text(
    color  = 'black', 
    vjust  = 1.2,
    hjust  = 0,
    size   = 10, 
    margin = unit( c(4,0,4,0), 'mm') ), 
  # title
  plot.title = element_text(
    size  = 14, 
    hjust = 0, 
    vjust = 3),
  # legend 
  legend.key = element_blank()
)

theme_set(graph_theme)

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
  stat_pointinterval(color='black', .width=c(0.67,0.9)) + 
  labs(x = TeX("$\\hat{\\alpha}_k$"), y = 'k') + 
  geom_point(data=focal_pars, 
             aes(x = alpha_k, y=k), 
             color = 'magenta', pch=8, size=2, stroke=2)
  

p_k_draws = fit1 |> 
  spread_draws(p_z[k], sigma_p, mu_p) |> 
  group_by(k) |> 
  mutate(p_k_hat = inv_logit(mu_p + sigma_p*p_z))

p_k_draws |> 
  ggplot(aes(x=p_k_hat, y=factor(k))) + 
  stat_pointinterval(color='black', .width=c(0.67,0.9)) + 
  labs(x = TeX("$\\hat{p}_k$"), y = 'k') + 
  geom_point(data=focal_pars, 
             aes(x = p_k, y=k), 
             color = 'cyan3', pch=8, size=2, stroke=2)







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
mu_b = 1
sigma_b = 0.3

# consumer mass 
M = matrix(NA_real_, nrow = S, ncol = K)

# generate varying base rates and detection probabilities, and body mass effects 
alpha_k = rnorm(K, mu_alpha, sigma_alpha)
p_k = inv_logit(rnorm(K, mu_p, sigma_p))
b_k = rnorm(K, mu_b, sigma_b)

# ---- Outputs ----
A_true = array(0, dim = c(S,S,K))
A_obs  = array(0, dim = c(S,S,K))
A_Mtrue = array(0, dim = c(S,S,K))
A_Mobs  = array(0, dim = c(S,S,K))


for(k in 1:K) {
  # Random effects
  u_i = rnorm(S, mean = 0, sd = sigma_u)
  v_j = rnorm(S, mean = 0, sd = sigma_v)
  
  
  # Mass 
  mass_k = exp(rnorm(S, mu_b, sigma_b))
  M[,k] = log(mass_k)
  
  # Loop 
  for(i in 1:S) for(j in 1:S) {
    
    pi_ijk = inv_logit(alpha_k[k] + u_i[i] + v_j[j])
    A_true[i,j,k] = rbinom(1, 1, pi_ijk)
    A_obs[i,j,k] = rbinom(1, 1, p_k[k] * pi_ijk)
    
    # body mass model 
    .pi_ijk = inv_logit(alpha_k[k] + (u_i[i] + b_k[k] * M[i,k]) + v_j[j])
    A_Mtrue[i,j,k] = rbinom(1, 1, .pi_ijk)
    A_Mobs[i,j,k] = rbinom(1, 1, p_k[k] * .pi_ijk)
    
  }
}

data_list = list(
  N_cons = S, 
  N_res  = S, 
  K      = K, 
  M      = M, 
  A      = A_Mobs
)

run = T
if(run) {
  mod2 <- cmdstan_model("C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/detection_model_bodymass.stan")
  
  fit2 <- mod2$sample(
    data = data_list,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    adapt_delta = 0.98
  )
  
  fit2$save_output_files(dir = 'C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/fits')
  
  saveRDS(fit2, 'C:/Users/scaggs.32/OneDrive - The Ohio State University/Professional/sascaggs.github.io/posts/2025-11-17-food-web-measurement/fits/detection_model_bodymass.RDS')
}

focal_pars2 = data.frame(k = factor(1:k), 
                        alpha_k, p_k, b_k)
focal_pars2

alpha_k_draws2 = fit2 |> 
  spread_draws(mu_alpha, sigma_alpha, alpha_z[k]) |>
  group_by(k) |> 
  mutate(alpha_k_hat = mu_alpha + sigma_alpha * alpha_z) 

alpha_k_draws2 |> 
  ggplot(aes(x=alpha_k_hat, y=factor(k))) + 
  stat_pointinterval(color='black', .width=c(0.67,0.9)) + 
  labs(x = TeX("$\\hat{\\alpha}_k$"), y = 'k') + 
  geom_point(data=focal_pars2, 
             aes(x = alpha_k, y=k), 
             color = 'white', pch=8, size=4, stroke=2) +
  geom_point(data=focal_pars2, 
             aes(x = alpha_k, y=k), 
             color = 'magenta', pch=8, size=2, stroke=2)


p_k_draws2 = fit2 |> 
  spread_draws(p_z[k], sigma_p, mu_p) |> 
  group_by(k) |> 
  mutate(p_k_hat = inv_logit(mu_p + sigma_p*p_z))

p_k_draws2 |> 
  ggplot(aes(x=p_k_hat, y=factor(k))) + 
  stat_pointinterval(color='black', .width=c(0.67,0.9)) + 
  labs(x = TeX("$\\hat{p}_k$"), y = 'k') + 
  geom_point(data=focal_pars2, 
             aes(x = p_k, y=k), 
             color = 'white', pch=8, size=4, stroke=2) +
  geom_point(data=focal_pars, 
             aes(x = p_k, y=k), 
             color = 'cyan3', pch=8, size=2, stroke=2)


b_k_draws2 = fit2 |> 
  spread_draws(b_k[k]) |> 
  group_by(k) 

b_k_draws2 |> 
  ggplot(aes(x=b_k, y=factor(k))) + 
  stat_pointinterval(color='black', .width=c(0.67,0.9)) + 
  labs(x = TeX("$\\hat{b}_k$"), y = 'k') + 
  geom_point(data=focal_pars2, 
             aes(x = b_k, y=k), 
             color = 'white', pch=8, size=4, stroke=2) +
  geom_point(data=focal_pars, 
             aes(x = b_k, y=k), 
             color = '#3300ff', pch=8, size=2, stroke=2)














