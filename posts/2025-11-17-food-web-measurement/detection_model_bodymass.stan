
data {
  int<lower=1> N_cons;
  int<lower=1> N_res;
  int<lower=1> K;
  matrix[N_cons,K] M;
  // observed adjacency
  array[N_cons, N_res, K] int<lower=0, upper=1> A; 
}


parameters {
  // body mass slopes 
  vector[K] b_z;
  real mu_b;
  real<lower=0> sigma_b;
  
  // non-centered latent effects
  vector[K] p_z;
  vector[K] alpha_z;
  matrix[N_cons, K] u_z;
  matrix[N_res, K] v_z;
  
  // hyperparameters
  real mu_p; 
  real<lower=0> sigma_p; 
  real mu_alpha; 
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_u; 
  real<lower=0> sigma_v;
}

model {
  // b priors
  b_z ~ normal(0,1);
  mu_b ~ normal(0,1); 
  sigma_b ~ exponential(1);
  
  // priors for latent z-parameters 
  p_z ~ normal(0,1); 
  alpha_z ~ normal(0,1);
  to_vector(u_z) ~ normal(0,1); 
  to_vector(v_z) ~ normal(0,1); 
  
  //priors for hyperparameters 
  mu_p ~ normal(0,1);
  sigma_p ~ exponential(1);
  mu_alpha ~ normal(0,1);
  sigma_alpha ~ exponential(1); 
  sigma_u ~ exponential(1);
  sigma_v ~ exponential(1);
  
  // likelihood 
  for(k in 1:K) {
    // detection random effect
    real p_k = inv_logit(mu_p + sigma_p * p_z[k]);
    
    // link random effect
    real alpha_k = mu_alpha + sigma_alpha * alpha_z[k];
    
    // body mass slopes 
    real b_k = mu_b + sigma_b * b_z[k];
    
    for(i in 1:N_cons) {
      // consumer random effects 
      real u_ik = sigma_u * u_z[i,k];
      
      for(j in 1:N_res) {
        // resource random effect 
        real v_jk = sigma_v * v_z[j,k];
        
        // link probability 
        real pi_ijk = inv_logit(alpha_k + (u_ik + b_k*M[i,k]) + v_jk);
        
        // detection adjusted probability 
        real q_ijk = p_k * pi_ijk;
        
        // likelihood
        A[i,j,k] ~ bernoulli(q_ijk);
      }
    }
  }
}

