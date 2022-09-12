#include <Rcpp.h>
using namespace Rcpp;
// #include <RcppEigen.h>
// using namespace Eigen;
// // [[Rcpp::depends(RcppParallel)]]
// #include <RcppParallel.h>
// using namespace RcppParallel;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector ES_score_matrix(NumericVector rw, 
                              NumericVector sum_rset_w,
                              LogicalMatrix sets_j,
                              NumericVector w,
                              NumericVector dec
                              ) {
  int s = sets_j.cols();
  int p = rw.size();
  NumericVector cum_sum(s,0.0);
  NumericVector mx_pos(s,0.0);
  NumericVector mx_neg(s,0.0);
  
  for(int k = 0; k < s; k++) {
    for(int i = 0; i < p; i++) {
      if(sets_j(i,k) == true) {
        cum_sum[k] += abs(rw[i])/sum_rset_w[k];
      } else{
        cum_sum[k] -= w[i]*dec[k];
      }
      if(cum_sum[k] > mx_pos[k]){ mx_pos[k] = cum_sum[k]; }
      if(cum_sum[k] < mx_neg[k]){ mx_neg[k] = cum_sum[k]; }
    }
  }
  // return(pmax(mx_pos,-1*mx_neg)); // two-sided KS stat
  return( -1.0 * (mx_pos+mx_neg));  // GSEA/GSVA stat 
}

// [[Rcpp::export]]
List get_set_stats(NumericVector w, NumericVector r, LogicalMatrix sets_j) {
  int s = sets_j.cols();
  double sum_w = sum(w);
  NumericVector dec(s,0.0);
  NumericVector sum_rw_sets(s,0.0);
  NumericVector norm_factor(s,0.0);
  IntegerVector n_per_set = colSums(sets_j);
  
  for(int k = 0; k < s; k++) {
    NumericVector w_set = w[sets_j(_,k)];
    dec(k) = sum(w_set);
    NumericVector r_set = r[sets_j(_,k)];
    double sum_rw_set = sum(w_set*r_set);
    double rw2_set = sum(pow(w_set*r_set,2));
    double rw_set_norm = pow(sum_rw_set,2.0)/rw2_set;
    NumericVector w_nset = w[!sets_j(_,k)];
    double sum_w_nset = sum(w_nset);
    double w2_nset = sum(pow(w_nset,2));
    double w_nset_norm = pow(sum_w_nset,2.0)/w2_nset;
    sum_rw_sets(k) = sum_rw_set;
    norm_factor(k) = 1.0/sqrt(rw_set_norm*w_nset_norm/(rw_set_norm+w_nset_norm));
  }
  return(List::create(Named("dec") = 1.0 / (sum_w - dec),
                      Named("norm_factor") = norm_factor,
                      Named("sum_rw_sets") = sum_rw_sets));
}

// void get_ES_matrix(RVector<double> rw,
//                               RVector<double> sum_rset_w,
//                               RMatrix<int> sets_j,
//                               RVector<double> w,
//                               RVector<double> dec,
//                               RVector<double> res) {
//   int s = sets_j.ncol();
//   int p = rw.size();
//   // RVector<double> cum_sum(s);
//   // RVector<double> mx_pos(s);
//   // RVector<double> mx_neg(s);
//   // RVector<double> res(s);
//   // double cum_sum = 0;
//   // double mx_pos = 0;
//   // double mx_neg = 0;
//   // double res
//   
//   for(int k = 0; k < s; k++) {
//     double cum_sum = 0;
//     double mx_pos = 0;
//     double mx_neg = 0;
//     for(int i = 0; i < p; i++) {
//       if(sets_j(i,k) == 1) {
//         cum_sum += abs(rw[i])/sum_rset_w[k];
//       } else{
//         cum_sum -= w[i]*dec[k];
//       }
//       if(cum_sum > mx_pos){ mx_pos = cum_sum; }
//       if(cum_sum < mx_neg){ mx_neg = cum_sum; }
//     }
//     res[k] = -1.0 * (mx_pos + mx_neg);
//   }
//   // for(int k=0;k<s;k++) res[k] = -1.0 * (mx_pos[k]+mx_neg[k]);
//   // return(res);
//   // return( -1.0 * (mx_pos+mx_neg));
// }
// 
// // [[Rcpp::export]]
// NumericVector ES_score_matrix2(NumericVector rw_, 
//                               NumericVector sum_rset_w_,
//                               IntegerMatrix sets_j_,
//                               NumericVector w_,
//                               NumericVector dec_
// ) {
//   
//   NumericVector res_(sets_j_.ncol());
//   RVector<double> rw(rw_);
//   RVector<double> sum_rset_w(sum_rset_w_);
//   RMatrix<int> sets_j(sets_j_);
//   RVector<double> w(w_);
//   RVector<double> dec(dec_);
//   RVector<double> res(res_);
//   
//   get_ES_matrix(rw,sum_rset_w,sets_j,w,dec,res);
//   // NumericVector res = );
//   return(wrap(res));
// }
// 
// // NumericMatrix do_ES2(NumericMatrix expr_kernel,NumericMatrix weights,IntegerMatrix sets) {
// //   
// // }



// // [[Rcpp::export]]
// List get_set_stats2(NumericVector w, NumericVector r, LogicalMatrix sets_j) {
//   int s = sets_j.cols();
//   int p = w.size();
//   double sum_w = sum(w);
//   NumericVector dec(s,0.0);
//   NumericVector sum_rset_w(s,0.0);
//   NumericVector sum_set_w(s,0.0);
//   NumericVector norm_factor(s,0.0);
//   IntegerVector n_per_set = colSums(sets_j);
//   
//   for(int k = 0; k < s; k++) {
//     for(int i = 0; i < p; i++) {
//       if(sets_j(i,k)) {
//         dec(k) += w[i];
//         sum_set_w(k) += w[i];
//         sum_rset_w(k) += r[i]*w[i];
//       }
//     }
//     norm_factor[k] = sqrt(sum_set_w[k])/n_per_set[k];
//   }
//   // dec = 1.0 / (sum_w - dec);
//   return(List::create(Named("dec") = 1.0 / (sum_w - dec),
//                       Named("norm_factor") = norm_factor,
//                       Named("sum_rset_w") = sum_rset_w));
// }

// // [[Rcpp::export]]
// NumericVector get_dec(NumericVector w, LogicalMatrix sets_j) {
//   int s = sets_j.cols();
//   int p = w.size();
//   double sum_w = sum(w);
//   NumericVector dec(s,0.0);
//   for(int k = 0; k < s; k++) {
//     for(int i = 0; i < p; i++) {
//       if(sets_j(i,k)) dec(k) += w[i];
//     }
//   }
//   return(1.0/(sum_w - dec));
// }
// 
// 
// // [[Rcpp::export]]
// NumericVector get_sum_rset_w(NumericVector rw, LogicalMatrix sets_j) {
//   int s = sets_j.cols();
//   int p = rw.size();
//   NumericVector sum_rset_w(s,0.0);
//   for(int k = 0; k < s; k++) {
//     for(int i = 0; i < p; i++) {
//       if(sets_j(i,k)) sum_rset_w(k) += abs(rw[i]);
//     }
//   }
//   return(sum_rset_w);
// }
