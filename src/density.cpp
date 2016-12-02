#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double offset = (-log(2)-log(M_PI))/2;
const double df_limit_switch_to_normal = 10000;

// Density of normal or Student's t-distribution
inline double student_normal_pdf(const double& x, const double& sigma, const double& nu, const double& log_beta_precomp)
// // For testing:
// // [[Rcpp::export]]
// double student_normal_pdf(const double x, const double sigma, const double nu, const double log_beta_precomp)
{
  double result;

  if (nu > df_limit_switch_to_normal || nu == 0)
  {
    // normal distribution
    result = - x*x / 2 / (sigma*sigma) - log(sigma) + offset;
  }
  else
  {
    const double x_adj = x / sigma;
    double basem1 = x_adj * x_adj / nu;

    if (basem1 < 0.125)
    {
      result = -std::log1p(basem1);
    }
    else
    {
      result = -log(1 + basem1);
    }

    result *= (1+nu) / 2;
    result -= log(nu)/2 + log_beta_precomp + log(sigma);
  }
  return result;
}


/*** R

# tests for pdf

nu <- 10
x <- 0
sigma <- 2

print(dt(x/sigma, nu, log=T)-log(sigma))
print(student_normal_pdf(x, sigma, nu, log(beta(nu/2,0.5))))

nu <- 10
x <- 10
sigma <- 2

print(dt(x/sigma, nu, log=T)-log(sigma))
print(student_normal_pdf(x, sigma, nu, log(beta(nu/2,0.5))))

x <- 1
sigma <- 2

print(dnorm(x, 0, sigma, log=T))
print(student_normal_pdf(x, sigma, 0, NA))

*/



// Combined density of individual Student's t-distribution and overlapping distribution of means
// [[Rcpp::export]]
NumericVector dcombined(NumericVector mean1, NumericVector mean2, double sigma_normal, NumericVector nu,
                    NumericVector sigma, NumericVector log_beta_precomp) {


  const size_t n1 = mean1.size();
  const size_t n2 = mean2.size();

  const double log_sigma_normal = log(sigma_normal);
  const double sigma_normal2_2 = 2*sigma_normal*sigma_normal;
  const double N_offset = n2 * offset;

  NumericVector out(n1, 0);

  for (size_t i = 0; i < n1; i++) {
    const double m1 = mean1[i];

    double s = 0;
    double n = 0;

    for (size_t j = 0; j < n2; j++) {
      const double delta = mean2[j] - m1;
      s += student_normal_pdf(delta, sigma[j], nu[j], log_beta_precomp[j]);
      n -= delta*delta;
    }

    n = n/sigma_normal2_2 + N_offset - (n2-1) * log_sigma_normal;
    out[i] = exp(s)+exp(n);
  }

  return(out);
}


