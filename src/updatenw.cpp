#include <Rcpp.h>
using namespace Rcpp;

// assumes x and y are sorted in increasing order

// [[Rcpp::export]]
IntegerVector shiftVec(IntegerVector x, IntegerVector y) {

  int x_len = x.size();
  int y_len = y.size();
  
  IntegerVector z(x_len);

  int i = 0;
  
  for(int j = 0; j < y_len; j++) {
    while(i < x_len && x[i] < y[j]) {
      z[i] = x[i] - j;
      i++;
    }
  }

  // handle any trailing part of x after getting to the end of y
  while(i < x_len) {
    z[i] = x[i] - y_len;
    i++;
  }

  return(z);
}
