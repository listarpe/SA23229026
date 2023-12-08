#include <Rcpp.h>
using namespace Rcpp;

//' @title pagerank using Rcpp
//' @description pagerank using Rcpp
//' @param adj_matrix adjacency matrix
//' @param alpha Correction probability of transition probability matrix
//' @param max_iter Maximum Number Of Iterations
//' @return Pagerank values for each node
//' @examples
//' \dontrun{
//' r <- pagerankc(adj_matrix, alpha=0.9, max_iter=5000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector pagerankc(NumericMatrix adj_matrix, double alpha=0.85, int max_iter=10000) {
  int n = adj_matrix.nrow();
  NumericMatrix p(n, n);
  for (int i = 0; i < n; i++) {
    double rowsum = 0.0;
    for (int j = 0; j < n; j++) {
      rowsum += adj_matrix(i, j);
    }
    for (int j = 0; j < n; j++) {
      p(i, j) = alpha * adj_matrix(i, j) / rowsum + (1-alpha) / n;
    }
  }
  
  NumericVector r(n, 1.0/double(n));
  NumericVector r0(n);
  
  double err = 1.0;
  int iter = 0;
  while((err > 1e-6) & (iter < max_iter)){
    iter++;
    for (int i = 0; i < n; i++) {
      r0[i] = r[i];
    }
    
    double rsum = 0.0;
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) {
        sum += r0[j]*p(i, j);
      }
      r[j] = sum;
      rsum += sum;
    }
    
    for (int i = 0; i < n; i++) {
      r[i] = r[i] / rsum;
    }
    
    double normrr0 = 0.0;
    double normr = 0.0;
    for (int i = 0; i < n; i++) {
      normrr0 += pow(r[i] - r0[i], 2);
      normr += pow(r[i], 2);
    }
    err= pow(normrr0, 0.5) / pow(normr, 0.5);
  }
  return r;
}

