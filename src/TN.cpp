
#include <Rcpp.h>
using namespace Rcpp;


// the code for univariate truncated standard normal sampling is based on
// MatLab code by Z.I. Botev, see https://web.maths.unsw.edu.au/~zdravkobotev/
// and R package TruncatedNormal


static const double A = 0.4;
static const double minusA = -0.4;
static const double B = 2.05;


// Rayleigh rejection sampling
double nt(const double l, const double u) {
  double x;
  const double c = 0.5*l*l;
  const double f = std::expm1(c - 0.5*u*u);
  do {
    x = c - std::log1p(f * R::runif(0, 1));
  } while (x * std::pow(R::runif(0, 1), 2) > c);
  return std::sqrt(2*x);
}

// simple rejection sampling
double trnd(const double l, const double u) {
  double x;
  do {
    x = R::rnorm(0, 1);
  } while (x < l || x > u);
  return x;
}

//’ Generate a random value from a standardized univariate truncated normal distribution
//’
//’ @param l lower truncation bound.
//’ @param u upper truncation bound.
//’ @return A single draw from the standard univariate truncated normal distribution.
// [[Rcpp::export(rng=true)]]
double Crtuvn(const double l, const double u) {
  double out;
  if (l > A) {
    out = nt(l, u);
  } else if (u < minusA) {
    out = -nt(-u, -l);
  } else {
    if (std::abs(u - l) > B) {
      out = trnd(l, u);
    } else {
      out = R::qnorm(
        R::pnorm(l, 0, 1, true, false) + 
          (R::pnorm(u, 0, 1, true, false) - R::pnorm(l, 0, 1, true, false)) * R::runif(0, 1),
        0, 1, true, false
      );
    }
  }
  return out;
}

//’ Generate a Gibbs cycle for a standardized multivariate truncated normal distribution
//’
//’ @param v start value.
//’ @param Ut sparse matrix encoding the inequality constrained linear combinations.
//’ @param ustar vector encoding the (initial) lower bounds corresponding to Ut.
//’ @param eps small numerical value to stabilize the Gibbs sampler.
//’ @return A single draw from the standardized multivariate truncated normal distribution.
// [[Rcpp::export(rng=true)]]
NumericVector Crtmvn_Gibbs(const NumericVector & v, const SEXP Ut, const NumericVector & ustar, const double eps) {
  double a, b, vi, x, temp;
  if (!Rf_isS4(Ut) || !Rf_inherits(Ut, "dgCMatrix")) stop("Ut is not a dgCMatrix");
  const IntegerVector Utp(as<S4>(Ut).slot("p"));
  const IntegerVector Uti(as<S4>(Ut).slot("i"));
  const NumericVector Utx(as<S4>(Ut).slot("x"));
  NumericVector u(clone(ustar));
  int n = v.size();
  NumericVector out = no_init(n);
  // loop over variables
  for (int i = 0; i < n; i++) {
    a = R_NegInf;
    b = R_PosInf;
    vi = v[i];
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      u[Uti[j]] += Utx[j] * vi;
    }
    // loop over constraints that variable i is involved in
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      x = Utx[j];
      if (x > eps) {
        temp = u[Uti[j]]/x;
        if (temp > a) a = temp;
      } else if (x < -eps) {
        temp = u[Uti[j]]/x;
        if (temp < b) b = temp;
      }
    }
    if (a == R_NegInf && b == R_PosInf) {
      out[i] = R::rnorm(0, 1);
    } else if (a == b) {
      out[i] = a;
    } else if (a < b) {
      out[i] = Crtuvn(a, b);
    } else {
      // this seems a numerically stable way to deal with numerical inaccuracy:
      if (a < v[i]) {
        out[i] = a;
      } else if (b > v[i]) {
        out[i] = b;
      } else {
        out[i] = v[i];
      }
    }
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      u[Uti[j]] -= Utx[j] * out[i];
    }
  }
  return out;
}

//’ Generate a vector of truncated normal variates for use in probit binomial model
//’
//’ @param mu vector of normal means.
//’ @param y response vector.
//’ @return A vector of truncated normal variates.
// [[Rcpp::export(rng=true)]]
NumericVector CrTNprobit(const NumericVector & mu, const NumericVector & y) {
  double l, u;
  const int n = mu.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; i++) {
    l = y[i] == 1 ? -mu[i] : R_NegInf;
    u = y[i] == 0 ? -mu[i] : R_PosInf;
    out[i] = mu[i] + Crtuvn(l, u);
  }
  return out;
}
