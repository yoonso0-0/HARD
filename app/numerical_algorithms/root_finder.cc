
#include "root_finder.hh"

#include "../rad_function.hh"
#include <cmath>
#include <cstddef>

namespace hard::numerical_algorithms::root_finder {

// Options for root finders

// ------------------------------------------------------------------------------------------
double
Bisection::get_root(double c4,
  double c1,
  double c0,
  double upper_bracket_bound) {

  double an = 0.0;
  double bn = upper_bracket_bound;
  double cn;

  double fa, fb, fc;

  namespace energy_polynomial = hard::energy_polynomial;

  fa = energy_polynomial::func(an, c4, c1, c0);
  fb = energy_polynomial::func(bn, c4, c1, c0);

  if(fa * fb > 0.0) {
    flog_fatal("Check the initial bracket bounds for Bisection");
  }

  for(size_t i = 0; i < root_finder::max_iter; ++i) {

    cn = 0.5 * (an + bn);
    fc = energy_polynomial::func(cn, c4, c1, c0);

    double fa = energy_polynomial::func(an, c4, c1, c0);
    double fc = energy_polynomial::func(cn, c4, c1, c0);

    if(std::abs(fc) < root_finder::tol) {
      if(cn < 0.0) {
        flog_fatal("Bisection Method failed");
      }
      else {
        return cn;
      }
    }
    else if(fa * fc < 0.0) {
      // [a, c] is the new bracket
      bn = cn;
      fb = fc;
    }
    else if(fb * fc < 0.0) {
      // [c, b] is the new bracket
      an = cn;
      fa = fc;
    }
    else {
      flog_fatal("Bisection::Check the bounds");
    }
  }

  flog_fatal("Bisection method failed to converge");
  return -1.0;
}

// ------------------------------------------------------------------------------------------
double
NewtonRaphson::get_root(double c4,
  double c1,
  double c0,
  double upper_bracket_bound) {

  int i;
  double x = upper_bracket_bound * 0.5; // take this as an initial guess
  double dx{0.0};
  double f, fprime;

  for(i = 0; i < root_finder::max_iter; ++i) {

    f = energy_polynomial::func(x, c4, c1, c0);
    fprime = energy_polynomial::funcprime(x, c4, c1);

    dx = -f / fprime;

    if((std::abs(f) < tol) or (std::abs(dx) < tol * x)) {
      if(x < 0.0) {
        flog(error) << "Newton-Raphson method failed";
      }
      return x;
    }

    x = x + dx;
  }

  flog(error) << "Newton-Raphson method failed to converge with max_iter = "
              << root_finder::max_iter << "\n Last value : x = " << x
              << ", dx = " << dx << ", f = " << f << ", fprime = " << fprime;
  return -1.0; // No convergence
}

// ------------------------------------------------------------------------------------------
// double
// halleys_get_root(double c4, double c1, double c0, double upper_bracket_bound)
// {

//   double xn = upper_bracket_bound * 0.5; // Initial guess
//   double dx;

//   for(std::size_t i = 0; i < root_finder::max_iter; ++i) {
//     double fx = energy_polynomial::func(xn, c4, c1, c0);
//     double fprime = energy_polynomial::funcprime(xn, c4, c1);
//     double fprime2 = energy_polynomial::funcprime2(xn, c4);

//     dx = -(2.0 * fx * fprime) / (2.0 * fprime * fprime - fx * fprime2);

//     if((std::abs(fx) < tol) or (std::abs(dx) < tol * xn)) {
//       if(xn < 0.0) {
//         flog(error) << "Halley's method failed";
//       }
//       return xn;
//     }

//     xn = xn + dx;
//   }

//   flog_fatal("Halley's method failed to converge");
//   return -1.0; // No convergence
// }

} // namespace hard::numerical_algorithms::root_finder