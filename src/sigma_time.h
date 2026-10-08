#ifndef SIGMA_TIME_H
#define SIGMA_TIME_H

#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Calendar-time lookup into a local-volatility matrix.
//
// Convention (package-wide): the columns of a sigma matrix are uniformly
// spaced in *calendar time* t on [0, t_max]; column 0 is the valuation date
// (t = 0) and the last column is t = t_max. Every solver maps its own time
// level to a fractional column position and linearly interpolates between
// the two neighbouring columns, so an option of any maturity tau <= t_max
// reads sigma(S, t) at the correct calendar time, independently of its own
// step size. Positions outside [0, t_max] are clamped (flat extension).
// ---------------------------------------------------------------------------
namespace lvtime {

struct ColLerp {
  int    lo;   // lower column index
  int    hi;   // upper column index (== lo when on a node or clamped)
  double w;    // weight on the upper column, in [0, 1)
};

inline ColLerp col_lerp(double t, double t_max, int ncol) {
  ColLerp c{0, 0, 0.0};
  if (ncol <= 1 || t_max <= 0.0) return c;
  const double p = (t / t_max) * static_cast<double>(ncol - 1);
  if (p <= 0.0) return c;
  if (p >= static_cast<double>(ncol - 1)) { c.lo = c.hi = ncol - 1; return c; }
  c.lo = static_cast<int>(std::floor(p));
  c.hi = c.lo + 1;
  c.w  = p - static_cast<double>(c.lo);
  return c;
}

// Interpolated value from a column-major raw buffer (rows x ncol).
inline double at_raw(const double* sig, int rows, int i, const ColLerp& c) {
  const double lo = sig[i + c.lo * rows];
  if (c.hi == c.lo) return lo;
  const double hi = sig[i + c.hi * rows];
  return (1.0 - c.w) * lo + c.w * hi;
}

// Interpolated value from any matrix-like object with operator()(i, j).
template <typename M>
inline double at(const M& sig, int i, const ColLerp& c) {
  const double lo = sig(i, c.lo);
  if (c.hi == c.lo) return lo;
  return (1.0 - c.w) * lo + c.w * sig(i, c.hi);
}

} // namespace lvtime

#endif // SIGMA_TIME_H
