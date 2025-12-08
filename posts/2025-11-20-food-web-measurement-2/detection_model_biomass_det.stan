
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
  vector[K] beta_k;
  
  // body mass detection
  vector[K] gamma_k;
  
  // non-centered latent effects
  vector[K] alpha_z;
  matrix[N_cons, K] u_z;
  matrix[N_res, K] v_z;
  
  // hyperparameters
  real mu_p; 
  real mu_alpha; 
  real<lower=0> sigma_alpha; 
  real<lower=0> sigma_u; 
  real<lower=0> sigma_v;
  
}

model {
  
  // priors 
  beta_k ~ normal(0,1);
  gamma_k ~ normal(0,1);
  
  // priors for latent z-parameters 
  alpha_z ~ normal(0,1);
  
  to_vector(u_z) ~ normal(0,1); 
  to_vector(v_z) ~ normal(0,1); 
  
  //priors for hyperparameters 
  mu_p ~ normal(0,1);
  mu_alpha ~ normal(0,1);
  sigma_alpha ~ exponential(1); 

  sigma_u ~ exponential(1);
  sigma_v ~ exponential(1);
  
  // likelihood 
  for(k in 1:K) {
    
    // link random effect
    real alpha_k = mu_alpha + sigma_alpha * alpha_z[k];
    
    for(i in 1:N_cons) {
      
      // consumer random effects 
      real u_ik = sigma_u * u_z[i,k];
      
      for(j in 1:N_res) {
        
        // resource random effect 
        real v_jk = sigma_v * v_z[j,k];
        
        // detection random effect
        real p_ijk = inv_logit(mu_p + gamma_k[k] * (M[i,k] - M[j,k]) );
        
        // link probability 
        real pi_ijk = inv_logit(alpha_k + u_ik + v_jk + beta_k[k] * M[i,k]);
        
        // detection adjusted probability 
        real q_ijk = p_ijk * pi_ijk;
        
        // likelihood
        A[i,j,k] ~ bernoulli(q_ijk);
      }
    }
  }
}

