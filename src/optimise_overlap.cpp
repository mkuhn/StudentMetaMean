#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double CBRT2 = std::cbrt(2);

// Find optimum point for product of likelihood functions for
// normal and Student's t distribution usign an exact numerical
// solution has been found using Wolfram Alpha
// [[Rcpp::export]]
NumericVector optimise_overlap(NumericVector meann, double sigman, double nu, double means, double sigmas) {
  const double k = sigman * sigman * (1+nu);
  const double b = means;
  const double bb = b*b;
  const double cc = sigmas*sigmas;
  const double d = nu;

  const double f1const = -bb+3*cc*d+3*k;
  const double f2const = 2*bb*b+18*b*cc*d-9*b*k;

  size_t n = meann.size();
  NumericVector out(n);

  for (size_t i = 0; i < n; ++i) {
    const double a = meann[i];
    const double aa = a*a;
    const double f1 = -aa+2*a*b+f1const;
    const double f2 = -2*aa*a+6*aa*b-6*a*bb-18*a*cc*d+9*a*k+f2const;
    const double f3 = f2+sqrt(4*f1*f1*f1+f2*f2);

    out[i] = -std::cbrt(f3)/(3*CBRT2) + (CBRT2*f1)/(3*std::cbrt(f3))+(a+2*b)/3;
  }

  return(out);
}


