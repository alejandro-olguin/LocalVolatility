// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <complex>
#include <cmath>

using namespace Rcpp;

//' Heston Model (closed-form, characteristic function)
//'
//' Closed-form European option price under the Heston stochastic volatility
//' model using the characteristic function approach with numerical integration.
//'
//' @param s_0 Spot price (> 0).
//' @param k Strike price (> 0).
//' @param tau Time to expiry in years (> 0).
//' @param r_d Risk-free rate.
//' @param q Dividend yield.
//' @param kappa Mean reversion speed (> 0).
//' @param theta Long-run variance (> 0).
//' @param xi Volatility of variance (vol-of-vol, > 0).
//' @param rho Correlation between stock and variance in \code{[-1, 1]}.
//' @param v_0 Initial variance (> 0).
//' @param type Either \code{"call"} or \code{"put"}.
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' Under the Heston model, the stock price \eqn{S} follows:
//' \deqn{dS = (r_d - q) S \, dt + \sqrt{v} S \, dW_S}
//' \deqn{dv = \kappa(\theta - v) \, dt + \xi \sqrt{v} \, dW_v}
//' with \eqn{dW_S \, dW_v = \rho \, dt}.
//'
//' The price is computed via the Heston (1993) characteristic function
//' with the numerically stable formulation of Albrecher et al. (2007).
//' Integration uses adaptive Simpson's rule over the semi-infinite
//' frequency domain.
//'
//' @examples
//' heston_cf(100, 100, 1, 0.05, 0, 2, 0.04, 0.5, -0.7, 0.04, "call")
//'
//' @export
// [[Rcpp::export]]
double heston_cf(double s_0,
                 double k,
                 double tau,
                 double r_d,
                 double q,
                 double kappa,
                 double theta,
                 double xi,
                 double rho,
                 double v_0,
                 String type) {

  // Input validation
  if (s_0 <= 0.0)  stop("s_0 must be > 0");
  if (k <= 0.0)    stop("k must be > 0");
  if (tau <= 0.0)  stop("tau must be > 0");
  if (kappa <= 0.0) stop("kappa must be > 0");
  if (theta <= 0.0) stop("theta must be > 0");
  if (xi <= 0.0)   stop("xi must be > 0");
  if (v_0 <= 0.0)  stop("v_0 must be > 0");
  if (rho < -1.0 || rho > 1.0) stop("rho must be within [-1, 1]");

  const bool is_call = (type == "call");
  if (!is_call && type != "put") stop("type must be \"call\" or \"put\"");

  // Warn if Feller condition violated
  if (2.0 * kappa * theta <= xi * xi)
    Rcpp::warning("Feller condition 2*kappa*theta > xi^2 is violated; "
                  "variance may reach zero.");

  const double x_log = std::log(s_0 / k);

  // Heston characteristic function for P_j (j=1 or j=2)
  // Uses the "little Heston trap" formulation (Albrecher et al. 2007)
  // which avoids numerical blow-up by using exp(-d*tau) instead of exp(+d*tau).
  //
  // d_j = sqrt((rho*xi*i*u - b_j)^2 - xi^2*(2*u_j*i*u - u^2))
  // where u_1 = 0.5, u_2 = -0.5, b_1 = kappa - rho*xi, b_2 = kappa
  auto heston_phi = [&](double u, int j) -> std::complex<double> {
    const std::complex<double> I(0.0, 1.0);
    const double b  = (j == 1) ? kappa - rho * xi : kappa;
    const double uj = (j == 1) ? 0.5 : -0.5;

    const std::complex<double> iu = I * u;

    // Discriminant: (rho*xi*iu - b)^2 - xi^2*(2*uj*iu - u^2)
    const std::complex<double> d = std::sqrt(
      (rho * xi * iu - b) * (rho * xi * iu - b)
      - xi * xi * (2.0 * uj * iu - u * u)
    );

    // "Little Heston trap" formulation (stable for large tau):
    // g = (b - rho*xi*iu - d) / (b - rho*xi*iu + d)
    const std::complex<double> gminus = b - rho * xi * iu - d;
    const std::complex<double> gplus  = b - rho * xi * iu + d;
    const std::complex<double> g = gminus / gplus;

    const std::complex<double> exp_ndt = std::exp(-d * tau);

    const std::complex<double> D = (gminus / (xi * xi))
      * (1.0 - exp_ndt) / (1.0 - g * exp_ndt);

    const std::complex<double> C = (r_d - q) * iu * tau
      + (kappa * theta / (xi * xi))
        * (gminus * tau - 2.0 * std::log((1.0 - g * exp_ndt) / (1.0 - g)));

    return std::exp(C + D * v_0 + iu * x_log);
  };

  // Numerical integration via Simpson's rule over [0, u_max]
  // The integrand decays as u -> inf; u_max = 200 is generous
  const int n_points = 4096;
  const double u_max = 200.0;
  const double du = u_max / n_points;

  double sum1 = 0.0, sum2 = 0.0;

  // Simpson's 1/3 rule
  for (int n = 1; n <= n_points; ++n) {
    const double u = n * du;
    const std::complex<double> I(0.0, 1.0);

    const std::complex<double> phi1 = heston_phi(u, 1);
    const std::complex<double> phi2 = heston_phi(u, 2);

    const double f1 = (phi1 / (I * u)).real();
    const double f2 = (phi2 / (I * u)).real();

    // Check for NaN/Inf
    if (!std::isfinite(f1) || !std::isfinite(f2)) continue;

    const double w = (n == n_points) ? 1.0 : (n % 2 == 1) ? 4.0 : 2.0;
    sum1 += w * f1;
    sum2 += w * f2;
  }

  double P1 = 0.5 + (du / (3.0 * M_PI)) * sum1;
  double P2 = 0.5 + (du / (3.0 * M_PI)) * sum2;

  // Clamp probabilities
  P1 = std::max(0.0, std::min(1.0, P1));
  P2 = std::max(0.0, std::min(1.0, P2));

  const double call_price = s_0 * std::exp(-q * tau) * P1
                           - k * std::exp(-r_d * tau) * P2;

  if (is_call) {
    return std::max(call_price, 0.0);
  } else {
    // Put-call parity
    const double put_price = call_price - s_0 * std::exp(-q * tau)
                           + k * std::exp(-r_d * tau);
    return std::max(put_price, 0.0);
  }
}
