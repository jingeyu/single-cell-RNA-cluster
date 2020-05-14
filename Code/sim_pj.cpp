#include <RcppArmadillo.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
// find indices of R_lw = s
uvec matches( mat R_t, int s){
  uvec index = find(R_t == s);
  return index;
}

// the following function returns the coordinates of i in matrix R
// (R was expanded to a vector by column)
// [[Rcpp::export]]
vec r_tran(int i, int L = 10){
  vec ind;
  if((i+1) % L == 0){
    ind <<  L-1  <<  (i+1)/L - 1;
  }else{
    ind << (i+1) % L - 1 << (i+1)/L;
  }
  return ind;
}

// [[Rcpp::export]]
mat mu_update(mat X, vec sgm_sq_star_t, mat R_t,
              double tau_mu, double eta_mu, int S, int G){
  mat mu(G, S); mu.fill(1.0);
  vec X_hat(G);
  for(int s = 0; s < S; s++){
    uvec matchx = matches(R_t, s+1);
    int N_s = matchx.n_elem;
    if(N_s == 0)
      X_hat = vec(G, fill::zeros);
    else
      X_hat = sum(X.cols(matchx), 1);
    // vec tmp = sum(h_t.each_row() % p_t.col(s).t(), 1);
    vec sgm_sq_hat, mu_hat;
    sgm_sq_hat = 1/(N_s / sgm_sq_star_t  + 1 / pow(tau_mu,2));
    mu_hat = (X_hat / sgm_sq_star_t + eta_mu / pow(tau_mu,2)) % sgm_sq_hat;
    mu.col(s) = mu_hat + randn<vec>(G) % sqrt(sgm_sq_hat);
  }
  return mu;
}


// [[Rcpp::export]]
vec sgm_sq_star_update(mat X, mat R_t, mat mu_t,
                       int S, int G, int N,
                       double alpha, double beta){
  
  vec sgm_sq_t(G,fill::zeros);
  vec s(G); s.fill(beta);
  double shape = N/2.0 + alpha;
  for(int i = 0; i < N; i++){
    s = s + 0.5 * pow(X.col(i) - mu_t.col(R_t(i)-1), 2);
  }
  vec scale = 1/s;
  for(int g = 0; g < G; g++){
    sgm_sq_t.at(g) = 1/randg(distr_param(shape, scale.at(g)));
  }
  
  return sgm_sq_t;
}



// [[Rcpp::export]]
Rcpp::List neigh_index(int i, int L = 10, int W = 10){
  int l = r_tran(i, L).at(0);
  int w = r_tran(i, L).at(1);
  vec tmp;
  tmp << l << l-1 << l+1 << l <<  w-1<< w << w << w+1;
  mat Tmp = mat(tmp);
  Tmp.set_size(4,2);
  umat a = Tmp >= L || Tmp>= W || Tmp < 0 ;
  arma::mat N = arma::conv_to<arma::mat>::from(a); //convert type umat to mat .
  vec rows = sum(N,1);
  uvec b = find(rows == 0);//length of b is larger than 1;
  mat ind_r =   Tmp.rows(b);
  vec ind_x(ind_r.n_rows);
  for(int i = 0; i < ind_r.n_rows; i++){
    ind_x.at(i) = ind_r(i, 0)  + ind_r(i, 1) * L;
  }
  
  return Rcpp::List::create(Rcpp::Named("ind_r") = ind_r,
                            Rcpp::Named("ind_x") = ind_x );
}
//this function is to judge R_lw's neighbor elements.


//update the region vector
// [[Rcpp::export]]
mat  R_update(mat X, mat R_t, mat mu_t, mat theta_t, vec sgm_sq_t,
              int S, int G, int L, int W){
  mat R_hat(L,W,fill::zeros);
  for(int i = 0; i < L * W; i++){
    vec r = r_tran(i, L);
    Rcpp::List neigh = neigh_index(i,L,W);
    mat neigh_r = neigh.at(0);
    vec tmp(S,fill::zeros);
    for(int s = 0; s < S; s++){
      vec nor = log_normpdf(X.col(i), mu_t.col(s), sqrt(sgm_sq_t));
      double a = sum(nor);
      for(int j = 0; j < neigh_r.n_rows; j++){
        if(R_t(neigh_r(j,0),neigh_r(j,1)) != s+1){
          a += theta_t(s,R_t(neigh_r(j,0),neigh_r(j,1)) - 1);
        }
      }
      tmp.at(s) = a;
    }
    tmp = exp(tmp - tmp.max());
    vec prob = tmp/sum(tmp);
    uvec fullvec;
    fullvec = regspace<uvec>(1, S);
    R_hat(r.at(0),r.at(1)) =  Rcpp::RcppArmadillo::sample(fullvec, 1, false, prob).at(0);
  }
  return R_hat;
}


//calculate H(R|theta) function.
// [[Rcpp::export]]
double H_fun(mat R_t, mat theta_t, int N){
  double sum = 0.0;
  for(int i = 0; i < N; i++){
    uvec lw = neigh_index(i).at(1);
    vec r_nei = R_t(lw);
    uvec index = find(r_nei != R_t.at(i)); 
    vec neibor = r_nei(index);  
    // if(neibor.n_elem == 0)
    //   continue;
    for(int j = 0; j < neibor.n_elem; j++)
      sum += theta_t(R_t.at(i)-1, neibor.at(j)-1);
  }
  return - sum;
}


// [[Rcpp::export]]
mat theta_update(mat X, mat R_t, mat mu_t, mat theta_t, vec sgm_sq_t,
                 int S, int G, int L, int W, 
                 double tau_0, double eta_theta, double tau_theta){
  
  //M, interation number of DMC.
  mat theta_star = theta_t;
  //sample theta matrix from proposal distribution (normal)
  //upper triangular matrix.
  for(int i = 1; i < S; i++){
    for(int j = 0; j < i; j++){
      double theta_star_ji = randn<double>() * tau_0 + theta_t(j,i);
      // theta_star(j, i) = randn<double>() * tau_0 + Theta.slice(m-1)(j,i);
      // theta_star(i,j) = theta_star(j,i);
      mat y = R_update(X, R_t, mu_t, theta_star, sgm_sq_t,S,G,L,W);
      double log_ratio = log_normpdf(theta_star(i,j), eta_theta, tau_theta) +
        log_normpdf(theta_t(i, j), theta_star(i,j), tau_0) -
        log_normpdf(theta_t(i,j), eta_theta, tau_theta) - 
        log_normpdf(theta_star(i,j), theta_t(i,j), tau_0) + 
        H_fun(R_t, theta_t, L*W) + H_fun(y, theta_star, L*W) - 
        H_fun(R_t, theta_star, L*W) - H_fun(y, theta_t, L*W);
      
      double logu = log(randu<double>());
      if(log_ratio > logu){
        theta_star(j, i) = theta_star_ji;
        theta_star(i, j) = theta_star_ji;
      }
    }
  }
  
  return theta_star;
}

// [[Rcpp::export]]
mat h_update(mat Z, vec C_t, vec sgm_sq_t, 
             double eta_h, double tau_h, 
             int G, int K){
  mat h_t(G, K, fill::zeros) ;
  vec Z_hat(G,fill::zeros);
  for(int k = 0; k < K; k++){
    uvec matchx = matches(C_t, k+1);
    int N_k = matchx.n_elem;
    if(N_k == 0)
      Z_hat = vec(G, fill::zeros);
    else
      Z_hat = sum(Z.cols(matchx), 1);
    vec sgm_sq_hat, h_hat;
    sgm_sq_hat = 1/(N_k / sgm_sq_t  + 1 / pow(tau_h,2));
    h_hat = (Z_hat / sgm_sq_t + eta_h / pow(tau_h,2)) % sgm_sq_hat;
    h_t.col(k) = h_hat + randn<vec>(G) % sqrt(sgm_sq_hat);
  }
  return h_t;
}

// [[Rcpp::export]]
vec sgm_sq_update(mat Z, mat h_t, vec C_t, 
                  double alpha_1, double beta_1, 
                  int G, int M){
  vec sgm_sq_t(G, fill::zeros);
  vec s(G); s.fill(beta_1);
  double shape = M/2.0 + alpha_1;
  for(int i = 0; i < M; i++)
    s = s + 0.5 * pow(Z.col(i) - h_t.col(C_t.at(i)-1), 2);
  vec scale = 1/s;
  for(int g = 0; g < G; g++)
    sgm_sq_t.at(g) = 1/randg(distr_param(shape, scale.at(g)));
  return sgm_sq_t;
}

// [[Rcpp::export]]
vec rDir(vec alpha){
  int l = alpha.n_elem;
  vec x(l,fill::zeros);
  for(int i = 0; i < l; i++){
    arma::vec b = arma::randg<arma::vec>(1, distr_param(alpha.at(i), 1.0));
    x.at(i) = b.at(0);
  }
  vec p = x/sum(x);
  return p;
}


// [[Rcpp::export]]
vec C_update(mat Z, mat h_t, vec  C_t, vec sgm_sq_t,
             vec gam_1, int M, int K){
  for(int i = 0 ; i < M; i++){
    vec N_k = zeros<vec>(K);
    for(int k = 0; k < K; k++){
      uvec index = find(C_t == k+1);
      N_k.at(k) = index.n_elem;
    }
    vec log_pai_i = log(rDir(gam_1 + N_k));
    vec tmp(K,fill::zeros);
    for(int k =0; k < K; k++){
      vec nor = log_normpdf(Z.col(i), h_t.col(k), sqrt(sgm_sq_t));
      double a = sum(nor);
      tmp.at(k) = a;
    }
    tmp = tmp + log_pai_i;
    tmp = exp(tmp - max(tmp));
    vec prob = tmp/sum(tmp);
    uvec fullvec;
    fullvec = regspace<uvec>(1, K);
    C_t.at(i) =  Rcpp::RcppArmadillo::sample(fullvec, 1, false, prob).at(0);
  }
  return C_t;
}

// [[Rcpp::export]]
double BIC_k(mat X, mat mu, vec C, vec sgm_sq){
  int M = X.n_cols;
  int G = X.n_rows;
  int k = mu.n_cols;
  double tmp = 0.0;
  for(int i = 0; i < M; i++){
    vec temp = log_normpdf(X.col(i), mu.col(C.at(i)-1), sqrt(sgm_sq));
    tmp = tmp + sum(temp);
  }
  double BIC_re = -2 * tmp + k * log(G*M);
  return BIC_re;
}

  
  
  
  
  