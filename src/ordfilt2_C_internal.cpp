#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ordfiltInC(NumericMatrix data, int x, NumericVector weightedMatrix){
  int nrows = data.nrow();
  int ncols = data.ncol();
  int radius = sqrt(weightedMatrix.size());

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(radius*radius);

  for(int j = 1; j < ncols - floor(radius/2); j++){
    for(int i = 1; i < nrows - floor(radius/2); i++){
      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){
          int index = radius*m + n;
          int a = i + m - 1;
          int b = j + n - 1;
          miniMatrix[index] = data(a, b)*weightedMatrix[index];
        }
      }

      std::sort(miniMatrix.begin(), miniMatrix.end());

      emptyData(i, j) = miniMatrix[x - 1];
    }
  }

  return emptyData;
}
