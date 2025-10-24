// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' European Option 2D (closed-form, ADR)
//'
//' Closed-form Black–Scholes price for a European option on an ADR whose
//' underlying is the product S·X scaled by the ADR ratio n: ADR = (S·X)/n.
//' Assumes constant volatilities for S and X and correlation rho. The
//' effective ADR log-volatility is
//'   sigma_adr = sqrt(sigma_s^2 + sigma_x^2 + 2 * rho * sigma_s * sigma_x).
//'
//' The domestic-measure discounting is used: SX leg discounted by exp(-q * tau),
//' strike leg by exp(-r_d * tau). Note r_f does not enter the closed form directly
//' but is kept for API symmetry with PDE solvers.
//'
//' @param s_0 Stock spot price S(0).
//' @param x_0 FX spot price X(0) (domestic per foreign).
//' @param k Strike price (domestic currency).
//' @param tau Time to expiry (in years, > 0).
//' @param r_d Risk-free rate (domestic).
//' @param r_f Risk-free rate (foreign) — kept for signature symmetry.
//' @param q Dividend yield for the equity.
//' @param sigma_s Constant volatility of the stock.
//' @param sigma_x Constant volatility of the FX.
//' @param rho Correlation between the stock and the FX in `[-1, 1]`.
//' @param n ADR ratio (shares per ADR), usually >= 1.
//' @param type Either "call" or "put" (case-insensitive).
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' This is the standard Black–Scholes formula applied to ADR = (S·X)/n with
//' log-normal dynamics under the domestic measure. Let ADR_0 = (s_0 * x_0)/n and
//' sigma_adr as above. Then
//'   `d1 = [ln(ADR_0 / K) + (r_d - q + 0.5 * sigma_adr^2) * tau] / (sigma_adr * sqrt(tau))`,
//'   `d2 = d1 - sigma_adr * sqrt(tau)`.
//' The call is ADR_0 `e^{-q tau} N(d1) - K e^{-r_d tau} N(d2)`; the put uses put-call parity.
//'
//' @examples
//' european_option_cf_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01, 0.2, 0.1, 0.3, 1, "call")
//'
//' @export
// [[Rcpp::export]]
double european_option_cf_2d(double s_0,
                             double x_0,
                             double k,
                             double tau,
                             double r_d,
                             double r_f,
                             double q,
                             double sigma_s,
                             double sigma_x,
                             double rho,
                             int n,
                             String type) {

  // --- input validation ---
  if (n <= 0) stop("n must be >= 1");
  if (k <= 0.0) stop("k must be > 0");
  if (s_0 <= 0.0 || x_0 <= 0.0) stop("s_0 and x_0 must be > 0");
  if (tau < 0.0) stop("tau must be >= 0");
  if (rho < -1.0 || rho > 1.0) stop("rho must be within [-1, 1]");

  const double adr_0 = (s_0 * x_0) / static_cast<double>(n);
  const double sigma_adr2 = sigma_s * sigma_s
                          + sigma_x * sigma_x
                          + 2.0 * rho * sigma_s * sigma_x;
  const double sigma_adr = std::sqrt(std::max(0.0, sigma_adr2));

  // Edge case: tau == 0 -> payoff
  auto payoff_call = [&](double spot) { return std::max(spot - k, 0.0); };
  auto payoff_put  = [&](double spot) { return std::max(k - spot, 0.0); };

  const std::string t = type;
  const bool is_call = (t == "call" || t == "Call" || t == "CALL");

  if (tau == 0.0) {
    return is_call ? payoff_call(adr_0) : payoff_put(adr_0);
  }

  // If sigma_adr == 0, price is discounted intrinsic based on forward
  if (sigma_adr == 0.0) {
    const double F = adr_0 * std::exp((r_d - q) * tau);
    const double disc = std::exp(-r_d * tau);
    const double intrinsic_fwd = is_call ? std::max(F - k, 0.0)
                                         : std::max(k - F, 0.0);
    return disc * intrinsic_fwd;
  }

  const double sqrt_tau = std::sqrt(tau);
  const double d1 = (std::log(adr_0 / k) + (r_d - q + 0.5 * sigma_adr * sigma_adr) * tau)
                    / (sigma_adr * sqrt_tau);
  const double d2 = d1 - sigma_adr * sqrt_tau;

  const double Nd1 = R::pnorm5(d1, 0.0, 1.0, 1, 0);
  const double Nd2 = R::pnorm5(d2, 0.0, 1.0, 1, 0);
  const double Nm_d1 = R::pnorm5(-d1, 0.0, 1.0, 1, 0);
  const double Nm_d2 = R::pnorm5(-d2, 0.0, 1.0, 1, 0);

  const double disc_q  = std::exp(-q * tau);
  const double disc_rd = std::exp(-r_d * tau);

  if (is_call) {
    return adr_0 * disc_q * Nd1 - k * disc_rd * Nd2;
  } else {
    return k * disc_rd * Nm_d2 - adr_0 * disc_q * Nm_d1;
  }
}
