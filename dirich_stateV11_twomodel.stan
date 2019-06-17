//a new version including cell cycle
//25/03/19
//State transition Dirichlet distribution
//fit both cell sum and proportion data

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
  
  matrix lambda_to_alpha_3(int[,] knn_weight, int n_state, vector lambda, vector lambda_net) {
    int len_knn;
    int indx_start;
    int indx_end;
    matrix[n_state, n_state] alpha; //cell cycle - cell death

    len_knn = dims(knn_weight)[1];

    alpha = rep_matrix(0., n_state, n_state);
  
    for(j in 1:len_knn){
         int pos1 = knn_weight[j,1];
         int pos2 = knn_weight[j,2];
         
         alpha[pos1, pos2] = lambda[j];
     }
     
     indx_start = len_knn+1;
     indx_end   = num_elements(lambda);
     
     //cell cycle - cell death - transitions
     for(i in 1:n_state) alpha[i,i] = lambda_net[i] - sum(alpha[1:n_state, i]);
     
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
  int y_t_in[n_time]; //total number of viable cells per time point
  real y_t_sd_in[n_time]; //total number of viable cells per time point
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
  real y_t[n_time];
  
  int n_rate;
  
  for(t in 1:n_time) {
    
   x_t[t] = x_t_in[t] + 1e-6;
   x_t[t] = x_t[t]/sum(x_t[t]);
   
   y_t[t] = y_t_in[t]+1e-6;
  }
  
  n_rate = n_rate_in;
}

parameters {
  vector<lower = 0>[n_rate] lambda1;
  vector<upper = 2.5>[n_state] lambda_net;
  vector<lower = 0>[n_state] Sigma_vec1;
  real<lower = 0> Sigma_vec2;
  vector<lower = 0>[n_time] true_y;
}

model {

  vector[n_state] x_t_1;
  vector[n_state] x_t_1_long;
  matrix[n_state, n_state] alpha1;
  vector[n_state] x_0;
 
  lambda1 ~ normal(0, 2);
  lambda_net ~ normal(2, 1);
  Sigma_vec1  ~ normal(500, 1000);
  Sigma_vec2  ~ normal(0,10);
  
  true_y ~ normal(y_t, y_t_sd_in);
  
  alpha1 = lambda_to_alpha_3(knn_weight, n_state, lambda1, lambda_net);
  
  x_0 = x_t[1]*y_t_in[1];
  
  for(t in 2:n_time) {
    
      x_t_1_long = matrix_exp((time_points[t]-time_points[1]) * alpha1) * x_0;
      x_t_1 = x_t_1_long[1:n_state]+1e-6;
      x_t_1 = x_t_1/sum(x_t_1);
      
      x_t[t] ~ beta(Sigma_vec1 .* x_t_1, Sigma_vec1 .* (1-x_t_1));
      // for(i in 1:n_state) x_t[t, i] ~ neg_binomial_2(x_t_1[i], Sigma_vec1[i]);
      y_t_in[t] ~ neg_binomial_2(sum(x_t_1_long), Sigma_vec2);
      //true_y[t] ~ normal(sum(x_t_1_long), Sigma_vec2);
      //y_t_in[t] ~ normal(sum(x_t_1_long), Sigma_vec2);
   
  }
}

// generated quantities{
//   vector[n_state] x_hat;
//   vector[n_state] x_0;
//   vector[n_time] y_t;
//   
//   int t = 2;
//   matrix[n_state, n_state] alpha1;
//  
//   vector[n_state] x_t_1[n_time];
//   real x_sum;
//   
//   x_t_1[1] = to_vector(x_t[1]);
//   y_t[1] = y_t_in[1];
// 
//   alpha1 = lambda_to_alpha_3(knn_weight, n_state, lambda1, lambda_net);
// 
//   x_0 = x_t_1[1]*y_t_in[1];
//  
//   while(t <= n_time){
//     
//     x_hat    = matrix_exp((time_points[t]-time_points[1]) * alpha1) * x_0;
//     x_sum = sum(x_hat);
//     x_hat = (x_hat + 1e-6)/sum(x_hat);
//     
//     // for(i in 1:n_state) x_t_1[t, i] = neg_binomial_2_rng(x_hat[i], Sigma_vec1[i]);
//     for(i in 1:n_state) x_t_1[t, i] = beta_rng(Sigma_vec1[i]* x_hat[i], Sigma_vec1[i]* (1-x_hat[i]));
//     
//     y_t[t] = neg_binomial_2_rng(x_sum, Sigma_vec2); // (x_sum, Sigma_vec2);neg_binomial_2_rng normal_rng
//     
//     t = t + 1;
//   }
// 
// }
