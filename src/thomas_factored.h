#ifndef THOMAS_FACTORED_H
#define THOMAS_FACTORED_H

// Zero-allocation Thomas algorithm routines for tridiagonal systems.
// All operate on raw double* arrays — no Rcpp, no R API, no heap allocations.
//
// Convention: a[0] = 0 (sub-diagonal), c[N-1] = 0 (super-diagonal).

#include <cmath>
#include <stdexcept>
#include <string>

namespace thomas {

// Full in-place Thomas solve.
// Inputs:  a[N] sub-diag, b[N] diag, c[N] super-diag, d[N] RHS
// Outputs: x[N] solution
// Workspace: cp[N] (caller-owned scratch buffer)
// a, b, c are NOT modified. d is NOT modified. cp is overwritten.
inline void solve(const double* a, const double* b, const double* c,
                  const double* d, double* x, double* cp, int N) {
  static constexpr double eps = 1e-16;
  // Forward sweep
  double denom = b[0];
  if (std::fabs(denom) < eps)
    throw std::runtime_error("thomas::solve: zero pivot at row 0");
  cp[0] = c[0] / denom;
  x[0]  = d[0] / denom;

  for (int i = 1; i < N; ++i) {
    denom = b[i] - a[i] * cp[i - 1];
    if (std::fabs(denom) < eps)
      throw std::runtime_error("thomas::solve: zero pivot at row " + std::to_string(i));
    cp[i] = c[i] / denom;
    x[i]  = (d[i] - a[i] * x[i - 1]) / denom;
  }

  // Back substitution
  for (int i = N - 2; i >= 0; --i) {
    x[i] -= cp[i] * x[i + 1];
  }
}

// Factorize the tridiagonal LHS: compute cp[i] and inv_denom[i].
// This depends only on (a, b, c) — NOT on the RHS.
// After factorization, multiple RHS vectors can be solved cheaply.
//
// cp[i]        = c[i] / denom[i]
// inv_denom[i] = 1.0 / denom[i]
// where denom[0] = b[0], denom[i] = b[i] - a[i]*cp[i-1]
inline void factor(const double* a, const double* b, const double* c,
                   double* cp, double* inv_denom, int N) {
  static constexpr double eps = 1e-16;
  double denom = b[0];
  if (std::fabs(denom) < eps)
    throw std::runtime_error("thomas::factor: zero pivot at row 0");
  inv_denom[0] = 1.0 / denom;
  cp[0] = c[0] * inv_denom[0];

  for (int i = 1; i < N; ++i) {
    denom = b[i] - a[i] * cp[i - 1];
    if (std::fabs(denom) < eps)
      throw std::runtime_error("thomas::factor: zero pivot at row " + std::to_string(i));
    inv_denom[i] = 1.0 / denom;
    cp[i] = c[i] * inv_denom[i];
  }
}

// Solve given a pre-factored LHS. Only touches the RHS vector.
// Inputs:  a[N], cp[N], inv_denom[N] from factor(), d[N] RHS
// Outputs: x[N] solution
// d is NOT modified.
inline void solve_factored(const double* a, const double* cp,
                           const double* inv_denom,
                           const double* d, double* x, int N) {
  // Forward substitution
  x[0] = d[0] * inv_denom[0];

  for (int i = 1; i < N; ++i) {
    x[i] = (d[i] - a[i] * x[i - 1]) * inv_denom[i];
  }

  // Back substitution
  for (int i = N - 2; i >= 0; --i) {
    x[i] -= cp[i] * x[i + 1];
  }
}

} // namespace thomas

#endif // THOMAS_FACTORED_H
