#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

//log-sum-exp trick

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}


// [[Rcpp::export]]
double log_zinb(int x, double mu, double phi, double theta0) {
  double lpmf_zinb;
  vec log_sum_exp_aux(2);
  if(x == 0){
    log_sum_exp_aux[0] = log(theta0);
    log_sum_exp_aux[1] = log(1-theta0) + R::dnbinom_mu(0, phi, mu, 1);
    lpmf_zinb = logSumExp(log_sum_exp_aux);
  }     
  else lpmf_zinb = log(1-theta0) + R::dnbinom_mu(x, phi, mu, 1);
  return(lpmf_zinb);
}

// [[Rcpp::export]]
mat predict_cluster_prob(umat data_x1, mat theta_fit, mat theta0_c_fit, mat mu_c_fit, mat phi_c_fit, vec mu_bg_fit, vec st_fit, vec day_indx){
  int num_cells = data_x1.n_rows;
  int num_clust = theta_fit.n_rows;
  int num_prot  = mu_bg_fit.n_elem;

  double st_i, theta_c_t, log_lik_i;
  double mu_c_p, theta_0, phi_c_p, mu_bg, mu_i;
  int cell_indx, clust_i, prot_i;
  
  mat log_theta_mat = log(theta_fit + 1e-10);
  urowvec cell_i;
  arma::mat all_cells(num_cells, num_clust);
  vec clust_lik(num_clust);
  

  for(cell_indx = 0; cell_indx < num_cells; cell_indx++){
   
    cell_i = data_x1.row(cell_indx);
    
    st_i   = st_fit(day_indx(cell_indx)-1);

    clust_lik.fill(0);

    for(clust_i  = 0; clust_i < num_clust; clust_i++) {
      theta_c_t = log_theta_mat(clust_i, day_indx(cell_indx)-1);
      log_lik_i = 0;
      for(prot_i = 0; prot_i < num_prot; prot_i++){
        theta_0 = theta0_c_fit(prot_i, clust_i);
        mu_c_p  = mu_c_fit(prot_i, clust_i);
        phi_c_p = phi_c_fit(prot_i, clust_i);
        mu_bg   = mu_bg_fit(prot_i);
        mu_i    = mu_bg*st_i + mu_c_p;
        log_lik_i = log_lik_i + log_zinb(cell_i(prot_i), mu_i, phi_c_p, theta_0);
      }
      clust_lik(clust_i) = log_lik_i + theta_c_t;

      all_cells(cell_indx, clust_i) = clust_lik(clust_i);
    }

   all_cells.row(cell_indx) = exp(all_cells.row(cell_indx) - all_cells.row(cell_indx).max());
   all_cells.row(cell_indx) = all_cells.row(cell_indx) / sum(all_cells.row(cell_indx));
      
     //.index_max() + 1;
    //all_cells(cell_indx) = clust_lik.index_max() + 1;
    // cout << cell_indx;
 }
  return(all_cells);
}


// [[Rcpp::export]]
vec predict_cluster(umat data_x1, mat theta_fit, mat theta0_c_fit, mat mu_c_fit, mat phi_c_fit, vec mu_bg_fit, vec st_fit, vec day_indx){
  int num_cells = data_x1.n_rows;
  int num_clust = theta_fit.n_rows;
  int num_prot  = mu_bg_fit.n_elem;
  double st_i, theta_c_t, log_lik_i;
  double mu_c_p, theta_0, phi_c_p, mu_bg, mu_i;
  int cell_indx, clust_i, prot_i;
  mat log_theta_mat = log(theta_fit + 1e-10);
  urowvec cell_i;
  arma::mat all_cells(num_cells, num_clust);
  vec clust_lik(num_clust);
  
  vec cell_cluster_max(num_cells);
  

  mat prob_matx = predict_cluster_prob(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx);
  
  for(cell_indx = 0; cell_indx < num_cells; cell_indx++){
    
    
    cell_cluster_max(cell_indx) = prob_matx.row(cell_indx).index_max() + 1;

  }
  return(cell_cluster_max);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


// /*** R
// timesTwo(42)
// */
