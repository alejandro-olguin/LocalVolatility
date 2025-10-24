#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>

// -------- Linear interpolation (1D) --------
// Safer overload when grid does not start at zero (preferred)
double linear_interpolation(double s0,
                            double s_min,
                            double ds,
                            const Rcpp::NumericVector& prices);

// Backward-compatible overload assuming s_min = 0
double linear_interpolation(double s0,
                            double ds,
                            const Rcpp::NumericVector& prices);

// -------- Bilinear interpolation (2D) --------
double bilinear_interpolation(const Rcpp::NumericVector& x,
                              const Rcpp::NumericVector& y,
                              const Rcpp::NumericMatrix& z,
                              const double& x0,
                              const double& y0);

// -------- Thomas tridiagonal solver --------
// a: subdiag (a[0]=0), b: diag, c: superdiag (c[N-1]=0), d: RHS
Rcpp::NumericVector thomas_algorithm(const Rcpp::NumericVector& a,
                                     const Rcpp::NumericVector& b,
                                     const Rcpp::NumericVector& c,
                                     const Rcpp::NumericVector& d);

// -------- Misc --------
double isinh(double x);

// Returns n+1 nodes (0..n), including x_min and x_max
Rcpp::NumericVector tavella_randall(double x0,
                                    double alpha,
                                    double x_min,
                                    double x_max,
                                    int n);

#endif // UTILS_H
