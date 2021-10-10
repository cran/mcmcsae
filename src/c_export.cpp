
#include <Rcpp.h>


extern "C" SEXP CHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP LDL, SEXP super, SEXP Imult, SEXP m);
extern "C" SEXP CHMf_solve(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHMf_solve_matrix(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHMf_spsolve(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHM_update_inplace(SEXP object, SEXP parent, SEXP mult);
extern "C" void CHM_options();


// [[Rcpp::export]]
SEXP cCHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP LDL, SEXP super, SEXP Imult, SEXP m) {
  return CHM_dsC_Cholesky(a, perm, LDL, super, Imult, m);
}

// [[Rcpp::export]]
SEXP cCHMf_solve(SEXP a, SEXP b, SEXP system) {
  return CHMf_solve(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHMf_solve_matrix(SEXP a, SEXP b, SEXP system) {
  return CHMf_solve_matrix(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHMf_spsolve(SEXP a, SEXP b, SEXP system) {
  return CHMf_spsolve(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHM_update_inplace(SEXP object, SEXP parent, SEXP mult) {
  return CHM_update_inplace(object, parent, mult);
}

// [[Rcpp::export]]
void cCHM_options() {
  CHM_options();
}
