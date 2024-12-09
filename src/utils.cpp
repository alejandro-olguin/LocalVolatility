#include <Rcpp.h>
using namespace Rcpp;

double linear_interpolation(double s_0, double ds, Rcpp::NumericVector prices){

  int i_star;
  i_star = s_0/ds;
  double sum = 0.;
  // run 2 point Lagrange polynomial interpolation
  sum = sum + (s_0 - (i_star + 1) * ds)/(i_star * ds - (i_star + 1) * ds) * prices(i_star);
  sum = sum + (s_0 - i_star * ds) / ((i_star + 1) * ds - i_star * ds) * prices(i_star + 1);
  return sum;
}

double bilinear_interpolation(NumericVector& x,
                              NumericVector& y,
                              NumericMatrix& z,
                              double& x0,
                              double& y0) {
  int n = x.size();
  int m = y.size();

  // Check if x0, y0 are within bounds
  if (x0 < x[0] || x0 > x[n - 1] || y0 < y[0] || y0 > y[m - 1]) {
    stop("The point (x0, y0) is outside the grid bounds.");
  }

  // Find the indices for the surrounding points
  int i = 0, j = 0;
  while (i < n - 1 && x[i + 1] <= x0) i++;
  while (j < m - 1 && y[j + 1] <= y0) j++;

  // Get the coordinates and values of the surrounding points
  double x1 = x[i], x2 = x[i + 1];
  double y1 = y[j], y2 = y[j + 1];
  double q11 = z(i, j), q21 = z(i + 1, j);
  double q12 = z(i, j + 1), q22 = z(i + 1, j + 1);

  // Compute the interpolation weights
  double denom = (x2 - x1) * (y2 - y1);
  double w1 = ((x2 - x0) * (y2 - y0)) / denom;
  double w2 = ((x0 - x1) * (y2 - y0)) / denom;
  double w3 = ((x2 - x0) * (y0 - y1)) / denom;
  double w4 = ((x0 - x1) * (y0 - y1)) / denom;

  // Return the interpolated value
  return w1 * q11 + w2 * q21 + w3 * q12 + w4 * q22;
}

NumericVector thomas_algorithm(NumericVector& a,
                               NumericVector& b,
                               NumericVector& c,
                               NumericVector& d) {
  size_t n = d.size() - 1;

  NumericVector c_prime = clone(c);
  NumericVector d_prime = clone(d);

  // This updates the coefficients in the first row
  // Note that we should be checking for division by zero here
  c_prime[0] /= b[0];
  d_prime[0] /= b[0];

  // Compute coefficients in the forward sweep
  for (int i = 1; i < n; i++) {
    c_prime[i] /= b[i] - a[i]*c_prime[i-1];
    d_prime[i] = (d_prime[i] - a[i]*d_prime[i-1]) / (b[i] - a[i]*c_prime[i-1]);
  }

  d_prime[n] = (d_prime[n] - a[n]*d_prime[n-1]) / (b[n] - a[n]*c_prime[n-1]);

  // This is the reverse sweep, used to update the solution vector
  for (int i = n; i-- > 0;) {
    d_prime[i] -= c_prime[i]*d_prime[i+1];
  }

  return(d_prime);
}

double isinh(double x) {
  return log(x + std::sqrt(x * x + 1.0));
}

NumericVector tavella_randall(double x0, double alpha, double x_min, double x_max, int n) {
  double c1 = isinh((x_min - x0) / alpha);
  double c2 = isinh((x_max - x0) / alpha);

  Rcpp::NumericVector result(n);

  // Create the sequence from 0 to (n-1)/n
  for (int j = 0; j < n; j++) {
    double i = static_cast<double>(j) / n;
    result[j] = x0 + alpha * std::sinh(c2 * i + c1 * (1.0 - i));
  }

  return result;
}

