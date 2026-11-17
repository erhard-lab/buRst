#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool isInteger(double x) {
  if (ISNAN(x)) return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) return false;
  return true;
}

// [[Rcpp::export]]
NumericVector cpp_gmRNA_new(double n, double km, double B, double d, double tend) {
  Rcpp::RNGScope scope;
  
  if (!isInteger(n) || n < 0) return NumericVector(0);
  int nn = (int)n;
  
  NumericVector res(nn);
  
  for (int i = 0; i < nn; i++) {
    double x = 0.0;
    double t = 0.0;
    
    while (t < tend) {
      double lambda1 = km;
      double lambda2 = d * x;
      double rate = lambda1 + lambda2;
      
      if (rate <= 0.0) break;
      
      // R::rexp uses scale = 1/rate
      double tau = R::rexp(1.0 / rate);
      if (t + tau > tend) break;
      t += tau;
      
      double u = R::runif(0.0, 1.0);
      if (u * rate < lambda1) {
        // burst
        int burst = R::rgeom(1.0 / (1.0 + B));
        x += burst;
      } else {
        // degradation (only if x>0)
        if (x > 0) x -= 1.0;
      }
    }
    
    res[i] = x;
  }
  
  return res;
}

// [[Rcpp::export]]
NumericVector cpp_gmRNA_old(double n, double d, double tend, std::vector<int> init) {
  Rcpp::RNGScope scope;
  
  if (!isInteger(n) || n < 0) return NumericVector(0);
  int nn = (int)n;
  
  if ((int)init.size() < nn) stop("init must have length >= n.");
  
  NumericVector res(nn);
  
  for (int i = 0; i < nn; i++) {
    double x = init[i];
    double t = 0.0;
    
    while (t < tend) {
      double rate = d * x;
      if (rate <= 0.0) break;
      
      double tau = R::rexp(1.0 / rate);
      if (t + tau > tend) break;
      t += tau;
      
      if (x > 0) x -= 1.0;
    }
    
    res[i] = x;
  }
  
  return res;
}

// [[Rcpp::export]]
NumericMatrix cpp_gmRNA_joint(int n, double km, double B, double d, double t_end, NumericVector init) {
  Rcpp::RNGScope scope;
  
  if (init.size() != n) stop("Length of init vector must equal n.");
  
  NumericMatrix out(n, 3);
  
  for (int i = 0; i < n; i++) {
    int n_old = init[i];
    int n_new = 0;
    double t = 0.0;
    
    while (t < t_end) {
      int n_total = n_old + n_new;
      
      double rate_synth = km;
      double rate_deg   = d * n_total;
      double rate_total = rate_synth + rate_deg;
      
      if (rate_total <= 0.0) break;
      
      double tau = R::rexp(1.0 / rate_total);
      if (t + tau > t_end) break;
      t += tau;
      
      double u = R::runif(0.0, 1.0);
      if (u < rate_synth / rate_total) {
        int burst = R::rgeom(1.0 / (1.0 + B));
        n_new += burst;
      } else {
        if (n_total > 0) {
          double v = R::runif(0.0, 1.0);
          if (n_new > 0 && v < (double)n_new / (double)n_total) {
            n_new -= 1;
          } else if (n_old > 0) {
            n_old -= 1;
          } else if (n_new > 0) {
            n_new -= 1;
          }
        }
      }
    }
    
    out(i, 0) = n_old + n_new;
    out(i, 1) = n_new;
    out(i, 2) = n_old;
  }
  
  return out;
}