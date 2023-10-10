
// need to start cholmod, in addition to one started by Matrix
// --> probably two cholmod workspaces; how (memory-)inefficient is this?

#include <R_ext/RS.h>
#include <cholmod.h>
#include <Matrix_stubs.c>


extern cholmod_common c;  // see mcmcsae_init.c


// make sure mcmcsae_init.c includes this:
/*
#include "Matrix.h"

cholmod_common c;

// and at the end:
void R_init_mcmcsae(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

  M_R_cholmod_start(&c);
}

void R_unload_mcmcsae(DllInfo *dll) {
  M_cholmod_finish(&c);
}
*/


// m method, integer
void chm_set_ordering(const int m) {
  if (m == -1) {
    // natural ordering, i.e. no permutation
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_NATURAL; c.postorder = FALSE;
  } else if (m == 0) {
    c.default_nesdis = TRUE;
    c.nmethods = 0;  // the default, but without METIS since that does not seem to be available in Matrix
  } else if (m == 1) {
    // only AMD
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_AMD; c.postorder = TRUE;
  } else if (m == 2) {
    // natural ordering, but with postordering
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_NATURAL; c.postorder = TRUE;
  } else if (m == 3) {
    // most extensive search
    c.nmethods = 9;
  }
}

// Cholesky of dsCMatrix
// see Matrix package's internal_chm_factor and dsCMatrix_Cholesky in dsCMatrix.c
// added argument m: ordering method (integer)
SEXP CHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP super, SEXP Imult, SEXP m) {
  CHM_FR L;
  CHM_SP A = AS_CHM_SP__(a);
  double beta[2] = {0, 0};
  beta[0] = asReal(Imult);
  R_CheckStack();

  int iSuper = asLogical(super),
      iPerm  = asLogical(perm);
  int im     = asInteger(m);
  if ((im < -1) || (im > 3)) error("Cholesky ordering method must be an integer between -1 and 3");

  // NA --> let CHOLMOD choose
  if (iSuper == NA_LOGICAL)	iSuper = -1;

  c.final_ll = 1;
  c.supernodal = (iSuper > 0) ? CHOLMOD_SUPERNODAL :
    ((iSuper < 0) ? CHOLMOD_AUTO : CHOLMOD_SIMPLICIAL);

  if (iPerm) {
    chm_set_ordering(im);
  } else {  // no permutation, m ignored in this case
    chm_set_ordering(-1);
  }

  //printf("c.final_ll = %d", c.final_ll);
  L = M_cholmod_analyze(A, &c);
  if (!M_cholmod_factorize_p(A, beta, (int*)NULL, 0, L, &c))
    error("Cholesky factorization failed");

  return M_chm_factor_to_SEXP(L, 1 /* free */);
}

// there is no M_chm_dense_to_SEXP in Matrix_stubs so write it here
// based on Matrix package's chm_dense_to_matrix in chm_common.c
SEXP chm_dense_to_matrixSEXP(CHM_DN a) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocMatrix(REALSXP, a->nrow, a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  M_cholmod_free_dense(&a, &c);
  UNPROTECT(1);
  return ans;
}

// basically chm_dense_to_vector, see chm_common.c of Matrix package
SEXP chm_dense_to_vectorSEXP(CHM_DN a) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocVector(REALSXP, a->nrow * a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  M_cholmod_free_dense(&a, &c);
  UNPROTECT(1);
  return ans;
}

// dense vector solve
SEXP CHMf_solve(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  int n = LENGTH(b);
  CHM_DN B = N_AS_CHM_DN(REAL(b), n, 1);

  int sys = asInteger(system);
  R_CheckStack();

  if (!(sys--)) error("invalid system argument");

  SEXP ans = chm_dense_to_vectorSEXP(M_cholmod_solve(sys, L, B, &c));
  return ans;
}

// dense matrix solve
SEXP CHMf_solve_matrix(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  int* dim = INTEGER(getAttrib(b, R_DimSymbol));
  CHM_DN B = N_AS_CHM_DN(REAL(b), dim[0], dim[1]);

  int sys = asInteger(system);
  R_CheckStack();

  if (!(sys--)) error("invalid system argument");

  SEXP ans = chm_dense_to_matrixSEXP(M_cholmod_solve(sys, L, B, &c));
  return ans;
}

SEXP CHMf_spsolve(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  CHM_SP B = AS_CHM_SP__(b);
  int sys = asInteger(system);
  R_CheckStack();

  if (!(sys--)) error("invalid system argument");

  SEXP ans = M_chm_sparse_to_SEXP(
    M_cholmod_spsolve(sys, L, B, &c), 1, 0, 0, "", R_NilValue);
  return ans;
}

// destructive Cholesky
SEXP CHM_update_inplace(SEXP object, SEXP parent, SEXP mult) {
  CHM_FR L = AS_CHM_FR(object);
  CHM_SP A = AS_CHM_SP__(parent);
  R_CheckStack();

  M_chm_factor_update(L, A, asReal(mult));
  return R_NilValue;
}
