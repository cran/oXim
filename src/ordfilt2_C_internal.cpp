#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib oXim
// [[Rcpp::export]]
NumericMatrix ordfiltInC(NumericMatrix data, double x, NumericMatrix weightedMatrix){
  double nrows = data.nrow();
  double ncols = data.ncol();

  double wmrows = weightedMatrix.nrow();
  double wmcols = weightedMatrix.ncol();

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(wmrows*wmcols);

  for(double j = 1; j < ncols - std::abs(floor(wmcols/2)); j++){
    for(double i = 1; i < nrows - std::abs(floor(wmrows/2)); i++){
      for(double n = 0; n < wmcols; n++){
        for(double m = 0; m < wmrows; m++){
          double index = m*3 + n;
          double a = i + m - 1;
          double b = j + n - 1;
          miniMatrix[index] = data(a, b)*weightedMatrix(m, n);
        }
      }

      std::sort(miniMatrix.begin(), miniMatrix.end());

      emptyData(i, j) = miniMatrix[std::abs(x) - 1];
    }
  }

  return emptyData;
}
