# LocalVolatility NEWS

## 2.0.0 (2025-10-24)

- 2D PDE drift and boundary alignment
  - Adopted domestic-measure S-drift consistently across 2D solvers:
    \(\mu_S = r_f - q - \rho\,\sigma_S\,\sigma_X\).
  - Updated far-field Dirichlet boundaries for product options to use
    \(\exp(-q\,\tau)\) on the \(S\cdot X\) term and \(\exp(-r_d\,\tau)\) on \(K\).
- American option solvers (2D) use penalty projection with split penalty per sub-step for stability.
- Added equivalence tests ensuring local-vol matrices filled with constants reproduce constant-vol prices (2D, European and American).
- Documentation overhaul: detailed roxygen blocks colocated with exported C++ functions; clarified grid sizes, drifts, boundaries, and examples.
- Minor typo fixes in parameter docs; standardized parameter names.

## 1.0 (2024-11-04)

- Initial release with pricing engines for 1D/2D options under local volatility.
- Utilities: Tavella–Randall nonuniform grid, Thomas algorithm, bilinear interpolation.
