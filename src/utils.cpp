// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
using namespace Rcpp;

// ---------- 1D linear interpolation ----------
// Safer overload that knows s_min (use this when s_min != 0)
double linear_interpolation(double s0,
                            double s_min,
                            double ds,
                            const Rcpp::NumericVector& prices) {
  const int N = prices.size();                  // nodes 0..N-1
  if (N < 2) stop("prices must have at least 2 points");

  // Clamp to domain just in case (callers usually handle this already)
  if (s0 <= s_min) return prices[0];
  const double s_max = s_min + ds * (N - 1);
  if (s0 >= s_max)  return prices[N - 1];

  // Cell index
  const double pos = (s0 - s_min) / ds;
  int i = static_cast<int>(std::floor(pos));    // 0..N-2
  double w = pos - i;                           // in [0,1)

  return (1.0 - w) * prices[i] + w * prices[i + 1];
}

// Backward-compatible version (assumes s_min = 0)
double linear_interpolation(double s0,
                            double ds,
                            const Rcpp::NumericVector& prices) {
  return linear_interpolation(s0, 0.0, ds, prices);
}


// ---------- 2D bilinear interpolation ----------
double bilinear_interpolation(const NumericVector& x,
                              const NumericVector& y,
                              const NumericMatrix& z,
                              const double& x0,
                              const double& y0) {
  const int nx = x.size();
  const int ny = y.size();
  if (nx < 2 || ny < 2) stop("x and y must have at least 2 points");
  if (z.nrow() != nx || z.ncol() != ny) stop("z size must match x (rows) and y (cols)");

  // Monotone increasing assumption for grids
  if (x0 < x[0] || x0 > x[nx - 1] || y0 < y[0] || y0 > y[ny - 1]) {
    stop("The point (x0, y0) is outside the grid bounds.");
  }

  // Find cell with upper_bound then step back one (clamped)
  int i_hi = std::upper_bound(x.begin(), x.end(), x0) - x.begin();
  int j_hi = std::upper_bound(y.begin(), y.end(), y0) - y.begin();
  int i = std::max(1, std::min(i_hi, nx - 1)) - 1; // 0..nx-2
  int j = std::max(1, std::min(j_hi, ny - 1)) - 1; // 0..ny-2

  const double x1 = x[i],     x2 = x[i + 1];
  const double y1 = y[j],     y2 = y[j + 1];
  const double q11 = z(i, j), q21 = z(i + 1, j);
  const double q12 = z(i, j + 1), q22 = z(i + 1, j + 1);

  const double dx = x2 - x1,  dy = y2 - y1;
  if (dx == 0.0 || dy == 0.0) stop("Degenerate cell in bilinear_interpolation");

  const double tx = (x0 - x1) / dx;
  const double ty = (y0 - y1) / dy;

  // Standard bilinear blend
  const double r1 = (1.0 - tx) * q11 + tx * q21;
  const double r2 = (1.0 - tx) * q12 + tx * q22;
  return (1.0 - ty) * r1 + ty * r2;
}

// ---------- Thomas algorithm for tridiagonal systems ----------
NumericVector thomas_algorithm(const NumericVector& a,  // subdiag (a[0] unused/0)
                               const NumericVector& b,  // diag
                               const NumericVector& c,  // superdiag (c[N-1] unused/0)
                               const NumericVector& d) { // RHS
  const int N = d.size();
  if (b.size() != N || a.size() != N || c.size() != N) {
    stop("Thomas: a, b, c, d must have the same length (with a[0]=0, c[N-1]=0).");
  }
  if (N < 2) stop("Thomas: system too small");

  NumericVector cp = clone(c);
  NumericVector dp = clone(d);

  // Guard against zero/near-zero pivots
  const double eps = 1e-16;

  // First row
  double denom = b[0];
  if (std::fabs(denom) < eps) stop("Thomas: zero pivot at row 0");
  cp[0] = cp[0] / denom;
  dp[0] = dp[0] / denom;

  // Forward sweep
  for (int i = 1; i < N - 1; ++i) {
    denom = b[i] - a[i] * cp[i - 1];
    if (std::fabs(denom) < eps) stop("Thomas: zero pivot at row " + std::to_string(i));
    cp[i] = cp[i] / denom;
    dp[i] = (dp[i] - a[i] * dp[i - 1]) / denom;
  }

  // Last row
  denom = b[N - 1] - a[N - 1] * cp[N - 2];
  if (std::fabs(denom) < eps) stop("Thomas: zero pivot at last row");
  dp[N - 1] = (dp[N - 1] - a[N - 1] * dp[N - 2]) / denom;

  // Back substitution
  for (int i = N - 2; i >= 0; --i) {
    dp[i] = dp[i] - cp[i] * dp[i + 1];
  }
  return dp; // solution
}



// ---------- Smooth nonuniform grid (Tavella–Randall) ----------
// Returns n+1 nodes, including endpoints x_min and x_max.
// Matches the rest of your PDE code that indexes 0..n.
double isinh(double x) {
  return std::log(x + std::sqrt(x * x + 1.0));
}

//' Generate a sequence based on the Tavella and Randall method
//'
//' This function generates a sequence of numbers following the Tavella and Randall method,
//' which applies a combination of hyperbolic sine transformations to create a smooth sequence
//' between a specified minimum and maximum value.
//'
//' @param x0 Numeric. The central value around which the sequence is generated.
//' @param alpha Numeric. A scaling parameter that controls the shape of the sequence.
//' @param x_min Numeric. The minimum value of the sequence.
//' @param x_max Numeric. The maximum value of the sequence.
//' @param n Integer. The number of values in the generated sequence.
//'
//' @return A numeric vector of length `n` containing the generated sequence.
//'
//' @details The function calculates the sequence by dividing the interval `[0, 1]` into `n`
//' equally spaced points. At each point, it computes the corresponding value based on the
//' Tavella and Randall method, which uses the hyperbolic sine function to ensure a smooth
//' transition between `x_min` and `x_max`.
//'
//' @export
// [[Rcpp::export]]
NumericVector tavella_randall(double x0,
                              double alpha,
                              double x_min,
                              double x_max,
                              int n) {
  if (n < 1) stop("tavella_randall: n must be >= 1");
  if (x_max <= x_min) stop("tavella_randall: x_max must be > x_min");
  if (alpha <= 0.0) stop("tavella_randall: alpha must be > 0");

  const double c1 = isinh((x_min - x0) / alpha);
  const double c2 = isinh((x_max - x0) / alpha);

  // IMPORTANT: n+1 points (0..n)
  Rcpp::NumericVector result(n + 1);
  for (int j = 0; j <= n; ++j) {
    const double t = static_cast<double>(j) / n;   // 0..1 inclusive
    result[j] = x0 + alpha * std::sinh(c2 * t + c1 * (1.0 - t));
  }

  // Make endpoints exact (avoid tiny round-off drift)
  result[0]  = x_min;
  result[n]  = x_max;
  return result;
}
