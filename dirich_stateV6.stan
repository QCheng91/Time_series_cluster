//State transition Dirichlet distribution
//cell death state
functions {
  
  matrix lambda_to_alpha(int[,] knn_weight, int n_state, vector lambda) {
    int len_knn;
    int indx_start;
    int indx_end;
    matrix[n_state + 1, n_state + 1] alpha;

    len_knn = dims(knn_weight)[1];

    alpha = rep_matrix(0., n_state + 1, n_state + 1);
  
    for(j in 1:len_knn){
     //for(i in 1:n_rate){
         int pos1 = knn_weight[j,1];
         int pos2 = knn_weight[j,2];
         
         alpha[pos1, pos2] = lambda[j];
     }
     indx_start = len_knn+1;
     indx_end   = num_elements(lambda);
     alpha[n_state+1, 1:n_state] = lambda[len_knn+1:indx_end]';
     for(i in 1:n_state + 1) alpha[i,i] = -sum(alpha[1:n_state+1, i]);
   
    return(alpha); 
  }
  
  
  matrix lambda_to_alpha_2(int[,] knn_weight, int n_state, vector lambda) {
    
    matrix[n_state, n_state] alpha;
  
    for(i in 1:n_state){
      for(j in 1:n_state){
        
        alpha[i,j] = 0;
      }
    }
    
    for(k in 1:n_state) {
      
      if(k==1){
        
        alpha[k, k] = -lambda[k];
        
      }
      
      else if(k==n_state) {
        
        alpha[k, k-1] = lambda[k-1];
      }
      
      else {
        
        alpha[k, k-1] = lambda[k-1];
        alpha[k, k] = -lambda[k];
      }
      
    }
    return(alpha); 
  }
}

data {

  int<lower=0> n_state; //number of states
  int<lower=0> n_time; //number of time points
  vector[n_state] x_t_in[n_time]; //data at t
  //vector[n_state] y[n_state]; //Hamming distance
  // int<lower=0> k; //no. of nearest neighbours
  vector[n_time] time_points; //time measured
  int<lower=0> n_rate_in;
  int knn_weight[n_rate_in, 2]; //knn weight for clusters
}
// make knn symmetric and make list yes
// make alpha transformed param yes
// lambda of length list knn yes 
// make jump point in one yes
// Compare theta to predicted number of cells in each cluster 

transformed data {
  vector[n_state] x_t[n_time]; //data at t
  int n_rate;
  
  for(t in 1:n_time) {
   x_t[t] = x_t_in[t] + 1e-6;
   x_t[t] = x_t[t]/sum(x_t[t]);
  }
  
  n_rate = n_rate_in + n_state;
}


parameters {
  vector<lower=0>[n_rate] lambda1;
  vector<lower=0>[n_rate] lambda2;
  // real<lower = 0> Sigma;
  vector<lower = 1000>[n_state] Sigma_vec1;
  vector<lower = 1000>[n_state] Sigma_vec2;
  // real<lower = 0> kk;
  
  vector<lower = 0>[n_rate] lambda_h;
  real<lower = 0> tau;
}

model {

  vector[n_state ] x_t_1;
  vector[n_state + 1] x_t_1_long;
  vector[n_state + 1] x_t_0_aux;
  
  // vector[n_state] aux_vec;
  matrix[n_state + 1, n_state + 1] alpha1;
  matrix[n_state + 1, n_state + 1] alpha2;

  vector[n_state + 1] x_0;
  vector[n_state + 1] x_11;
  // kk ~ normal(0, 10);
  
  //lambda ~ double_exponential(0, kk);
  lambda_h ~ cauchy(0, 1);
  tau ~ cauchy(0, 1);

  for(i in 1:n_rate) lambda1[i] ~ normal(0, lambda_h[i] * tau);
  for(i in 1:n_rate) lambda2[i] ~ normal(0, lambda_h[i] * tau);

  //Sigma  ~ normal(0, 1000);
  Sigma_vec1  ~ normal(0, 3000);
  Sigma_vec2  ~ normal(0, 3000);
  
  alpha1 = lambda_to_alpha(knn_weight, n_state, lambda1);
  alpha2 = lambda_to_alpha(knn_weight, n_state, lambda2);
  
  x_0            = rep_vector(0, n_state + 1);
  x_0[1:n_state] = x_t[1];
  x_11           = matrix_exp((11-time_points[1]) * alpha1) * x_0;
  
  for(t in 2:n_time) {
    
    if(time_points[t] <= 11) {
      
      x_t_1_long = matrix_exp((time_points[t]-time_points[1]) * alpha1) * x_0;
      x_t_1 = x_t_1_long[1:n_state] + 1e-6;
      x_t_1 = x_t_1/sum(x_t_1);
      x_t[t] ~ beta(Sigma_vec1 .* x_t_1, Sigma_vec1 .* (1-x_t_1));
      // x_t[t] ~ dirichlet(x_t_1*Sigma + 1e-4);
      // x_t[t] ~ normal(x_t_1, Sigma_vec);
    
    }
 
    else if(time_points[t] > 11) {
      
      x_t_1_long = matrix_exp((time_points[t]-11) * alpha2) * x_11;
      x_t_1 = x_t_1_long[1:n_state] + 1e-6;
      x_t_1 = x_t_1/sum(x_t_1);
      x_t[t] ~ beta(Sigma_vec2 .* x_t_1, Sigma_vec2 .* (1-x_t_1));
      // aux_vec = matrix_exp((11-time_points[1]) * alpha) * x_t[1];
      // x_t_1   = matrix_exp((time_points[t]-11) * alpha) * aux_vec;
      // x_t[t] ~ dirichlet(x_t_1[t]*Sigma + 1e-4);

    }

    // else{
    //   x_t_1[t] = matrix_exp((time_points[t]-14) * alpha) * x_t_1[9];
    //   x_t[t] ~ dirichlet(x_t_1[t]*Sigma + 1e-4);
    //   
    // }
   
  }
}

generated quantities{
  vector[n_state + 1] x_hat;
  vector[n_state + 1] x_0;
  vector[n_state + 1] x_0_11;
  int t = 2;
  matrix[n_state + 1, n_state + 1] alpha1;
  matrix[n_state + 1, n_state + 1] alpha2;
  // simplex[n_state] x_t_1[n_time];
  vector[n_state] x_t_1[n_time];
  x_t_1[1] = x_t[1];

  alpha1 = lambda_to_alpha(knn_weight, n_state, lambda1);
  alpha2 = lambda_to_alpha(knn_weight, n_state, lambda2);
  
  x_0            = rep_vector(0, n_state + 1);
  x_0[1:n_state] = x_t[1];
  x_0_11           = matrix_exp((11-time_points[1]) * alpha1) * x_0;
  
  while(t <= n_time){
    
    if(time_points[t] <= 11) {
      
      x_hat    = matrix_exp((time_points[t]-time_points[1]) * alpha1) * x_0;
      x_hat[1:n_state] = x_hat[1:n_state]/sum(x_hat[1:n_state]);
 
      // x_t_1[t] = dirichlet_rng(x_hat*Sigma);
      
      for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec1[i] * x_hat[i], Sigma_vec1[i] * (1-x_hat[i]));
      
    }
 
    else if(time_points[t] > 11) {

      x_hat    = matrix_exp((time_points[t]-11) * alpha2) * x_0_11;
      x_hat[1:n_state] = x_hat[1:n_state]/sum(x_hat[1:n_state]);

      for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec2[i] * x_hat[i], Sigma_vec2[i] * (1-x_hat[i]));

    }
    // 
    // else{
    //   
    //   x_hat    = matrix_exp((time_points[t]-time_points[9]) * alpha) * x_t[9];
    //   x_t_1[t] = dirichlet_rng(x_hat*Sigma);
    //   
    // }
    
    t = t + 1;
  }

}
