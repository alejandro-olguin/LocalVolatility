#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>

// Declare the function you want to use in another file
double linear_interpolation(double s_0, double ds, Rcpp::NumericVector prices);

double bilinear_interpolation(Rcpp::NumericVector& x,
                              Rcpp::NumericVector& y,
                              Rcpp::NumericMatrix& z,
                              double& x0,
                              double& y0);

Rcpp::NumericVector thomas_algorithm(Rcpp::NumericVector& a,
                                     Rcpp::NumericVector& b,
                                     Rcpp::NumericVector& c,
                                     Rcpp::NumericVector& d);


double isinh(double x);

Rcpp::NumericVector tavella_randall(double x0, double alpha, double x_min, double x_max, int n);

#endif
