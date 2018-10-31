//mixture model for noise+data+0spikes
//Updats on 24 Jul 2018
//hierarchical model for the noise part

// no noise part for cell cycle

//in this version we add dirichlet process for predicting number of clusters
functions{
  
  // Returns a vector with zinb evaluated for each value of x
  vector log_zinb_vec(int[] x, real mu, real phi, real theta0) {
    vector[size(x)] out_val;
    real phi_min_1;
    real mu_p_phi;
    real l_aux_ratio1;
    real l_aux_ratio2;
    real ldens_xeq0;
    real l_1_m_theta;
    
    // Auxiliary vars to speed up
    mu_p_phi     = mu + phi;
    phi_min_1    =  phi - 1;
    l_aux_ratio1 = log(mu/mu_p_phi);
    l_aux_ratio2 = phi*log(phi/mu_p_phi);
    l_1_m_theta  = bernoulli_lpmf(0 | theta0); // log(1-theta0);
    
    // ZINB log Desity for x = 0
    ldens_xeq0 = log_sum_exp(bernoulli_lpmf(1 | theta0), l_1_m_theta + l_aux_ratio2);

    for(cell_i in 1:size(x)){
      if(x[cell_i] == 0){ out_val[cell_i] = ldens_xeq0;}
      else{
        out_val[cell_i] = l_1_m_theta + lchoose(x[cell_i] + phi_min_1, x[cell_i]) + x[cell_i] * l_aux_ratio1 + l_aux_ratio2;
      }
      // if(is_inf(out_val[cell_i]) || is_nan(out_val[cell_i])){
      //   print("=======================");
      //   print(l_1_m_theta);
      //   print(lchoose(x[cell_i] + phi_min_1, x[cell_i]));
      //   print(x[cell_i] * l_aux_ratio1);
      //   print(l_aux_ratio2);
      //   print(out_val[cell_i]);
      //   print(mu);
      //   print(phi);
      //   print(theta0);
      //   print(x[cell_i] );
      //   print("=======================");
      // }
    }
    return(out_val);
  }
  
  vector log_zinb_vec2(int[] x_indx, int[] x_unique, real mu, real phi, real theta0) {
    vector[size(x_indx)] out_val;
    vector[size(x_unique)] ll_unique;
    // Calculate for unique values
    ll_unique = log_zinb_vec(x_unique, mu, phi, theta0);
    out_val   = ll_unique[x_indx];
    return(out_val);
  }
  
  
  real log_zinb(int x, real mu, real phi, real theta0) {
    real lpmf_zinb;
    if(x == 0) lpmf_zinb = log_sum_exp(bernoulli_lpmf(1 | theta0), bernoulli_lpmf(0 | theta0) + neg_binomial_2_lpmf(0 | mu, phi));
    else lpmf_zinb = bernoulli_lpmf(0 | theta0) + neg_binomial_2_lpmf(x | mu, phi);
    return(lpmf_zinb);
  }
  
}

data {
  int<lower=0> c; //  no. of clusters
  int<lower=0> N; //  no. of cells
  int<lower=1> k; // no. of proteins
  int x[N, k];    // Value for each sample on each dimension
  int<lower=0> d; // Day Index
  int day_index[N];
  int indx_tp[d,2]; // Index from and length to for each time point
  // Data to only calculate ZINB for each count value once 
  int num_unique_counts; // Number of unique count values of any protein
  int unique_prot_x[num_unique_counts];
  int indx_prot[k,2]; // Index from and length to for each protein. For
  int x_indx[N, k];  // Index of proteins
}

parameters {
  
  //real<lower = 0.1, upper = 1000>  mu_bg[k]; //, upper = 40
  real<lower = 0.1>  mu_s_cp[c, k];
  real<lower = 0.1>  phi_s_cp[c, k]; // Numerical stability of the negative binomial
  real<lower=0, upper =1>   theta0_s_cp[c, k];
  
  //adding a scaling factor
  //real<lower=0, upper = 20> scaling_t[d-1];
  
  //adding the dirichlet parameter
  // real<lower=0, upper = 1> v[d, c];
  real<lower=0, upper = 1> v[c];
  
  real<lower=0> alpha; // hyper prior DP(alpha,base)
}

transformed parameters{
  
  simplex[c] theta[d];
  
  // for(i in 1:d) theta[i, 1] = v[i, 1];
  for(i in 1:d) theta[i, 1] = v[1];
  
  // stick-break process based on The BUGS book Chapter 11 (p.294)
  for(i in 1:d) {
    
    for(j in 2:(c-1)){
      
      // theta[i, j]= v[i, j]*(1-v[i, j-1])*theta[i, j-1]/v[i, j-1]; 
      theta[i, j]= v[j]*(1-v[j-1])*theta[i, j-1]/v[j-1]; 
  }
  theta[i, c]= 1-sum(theta[i, 1:(c-1)]); // to make a simplex.
  }
  
}

model {
  real unsc_log_lik2;
  matrix[N, c] log_lik_clust_mat;
//  vector[N] log_lik_clust_mat[c];
  vector[c] log_like_i;
  real scale_use;
  vector[k] log_log_i_p; 
  vector[N] log_log_i_c; 
  vector[c] log_theta[d];
  int pos1;
  int pos2;
  
 // real 
  alpha ~ gamma(6,1);
  
 // mu_bg    ~ normal(0, 20);

  for(cluster_i in 1:c) {
    mu_s_cp[cluster_i,]     ~ normal(0, 2000);
    phi_s_cp[cluster_i,]    ~ normal(0, 2000);
    theta0_s_cp[cluster_i,] ~ beta(2, 2);
  }
  
  for(cluster_i in 1:c) v[cluster_i] ~ beta(1, alpha);
  
  for(day_i in 1:d){
    //theta[day_i] ~ dirichlet(rep_vector(2, c));
    //Calculate log theta
    log_theta[day_i] = log(theta[day_i]);
  }
  
  //scaling_t    ~ normal(0, 5);
  
 // Fill mat, need to do this for ragged array
 log_lik_clust_mat = rep_matrix(0, N, c);
 
 for(tp in 1:d) {
    // if(tp==1) scale_use = 1;
    // else scale_use = scaling_t[tp-1];
    pos1 = indx_tp[tp, 1];
    pos2 = indx_tp[tp, 2];
    for(cluster_i in 1:c) {
      for(prot_i in 1:k) {
        if(is_inf(mu_s_cp[cluster_i, prot_i]) || is_nan(mu_s_cp[cluster_i, prot_i])){
        print("=======================");
        // print(mu_bg[prot_i]*scale_use + mu_s_cp[cluster_i, prot_i]);
        // print(mu_bg[prot_i]);
        print(prot_i);
        // print(scale_use);
        print(mu_s_cp[cluster_i, prot_i]);
        print("=======================");
      }

        log_lik_clust_mat[pos1:pos2, cluster_i] = log_lik_clust_mat[pos1:pos2, cluster_i] + 
          log_zinb_vec2(x_indx[pos1:pos2, prot_i], unique_prot_x[indx_prot[prot_i,1]:indx_prot[prot_i,2]], 
                       mu_s_cp[cluster_i, prot_i], phi_s_cp[cluster_i, prot_i], theta0_s_cp[cluster_i, prot_i]);
      }
    }
  }

  for(cell_i in 1:N) target += log_sum_exp(log_theta[day_index[cell_i]] + log_lik_clust_mat[cell_i]');
    
}
