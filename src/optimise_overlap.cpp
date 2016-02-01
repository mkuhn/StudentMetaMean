// [[Rcpp::depends(BH)]]

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

#include <boost/math/tools/minima.hpp>
#include <functional>

using namespace Rcpp;
using namespace std::placeholders;

double dStudentTimesNormal(double x, double meann, double sigman, double nu, double means, double sigmas) {
  double xm = x - meann;
  double xs = (x - means)/sigmas;
  return xm*xm + sigman * sigman * (1 + nu) * log(xs*xs + nu);
}

// [[Rcpp::export]]
double optimise_overlap(double meann, double sigman, double nu, double means, double sigmas) {
  double a = meann;
  double b = means;
  if (meann > means) {
    a = means;
    b = meann;
  }
  return(boost::math::tools::brent_find_minima(
      std::bind(dStudentTimesNormal, _1, meann, sigman, nu, means, sigmas),
      meann, means, 16).first);
}


