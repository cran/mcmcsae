
// need to start cholmod, in addition to one started by Matrix
// --> probably two cholmod workspaces; how (memory-)inefficient is this?

#include <Matrix.h>
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
  c.error_handler = M_R_cholmod_error;
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
SEXP CHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP LDL, SEXP super, SEXP Imult, SEXP m) {
  CHM_FR L;
  CHM_SP A = AS_CHM_SP__(a);
  R_CheckStack();

  int iSuper = asLogical(super),
      iPerm  = asLogical(perm),
      iLDL   = asLogical(LDL);
  int im     = asInteger(m);
  if ((im < -1) || (im > 3)) error("Cholesky ordering method must be an integer between -1 and 3");

  double dImult = asReal(Imult);

  // NA --> let CHOLMOD choose
  if (iSuper == NA_LOGICAL)	iSuper = -1;

  c.final_ll = (iLDL == 0) ? 1 : 0;
  c.supernodal = (iSuper > 0) ? CHOLMOD_SUPERNODAL :
    ((iSuper < 0) ? CHOLMOD_AUTO : CHOLMOD_SIMPLICIAL);

  if (iPerm) {
    chm_set_ordering(im);
  } else {  // no permutation, m ignored in this case
    chm_set_ordering(-1);
  }
  L = M_cholmod_analyze(A, &c);
  if (!M_cholmod_factorize_p(A, &dImult, (int*)NULL, 0, L, &c))
    error("Cholesky factorization failed");

  return M_chm_factor_to_SEXP(L, 1 /* free */);
}

// there is no M_chm_dense_to_SEXP in Matrix_stubs so write it here
// based on Matrix package's chm_dense_to_matrix in chm_common.c
SEXP chm_dense_to_matrixSEXP(CHM_DN a, int dofree) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocMatrix(REALSXP, a->nrow, a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  if (dofree > 0) {
    M_cholmod_free_dense(&a, &c);
  } else if (dofree < 0) Free(a);
  UNPROTECT(1);
  return ans;
}

// basically chm_dense_to_vector, see chm_common.c of Matrix package
SEXP chm_dense_to_vectorSEXP(CHM_DN a, int dofree) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocVector(REALSXP, a->nrow * a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  if (dofree > 0) {
    M_cholmod_free_dense(&a, &c);
  } else if (dofree < 0) Free(a);
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

  SEXP ans = chm_dense_to_vectorSEXP(
    M_cholmod_solve(sys, L, B, &c), 1);
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

  SEXP ans = chm_dense_to_matrixSEXP(
    M_cholmod_solve(sys, L, B, &c), 1);
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

void CHM_options() { int a=0;
  Rprintf("/* parameters for symbolic/numeric factorization and update/downdate */ \n");
  Rprintf("dbound = %lf \n", c.dbound);
  Rprintf("grow0 = %lf \n", c.grow0);
  Rprintf("grow1 = %lf \n", c.grow1);
  Rprintf("grow2 = %d \n", c.grow2);
  Rprintf("maxrank = %d \n", c.maxrank);
  Rprintf("supernodal_switch = %lf \n", c.supernodal_switch);
  Rprintf("supernodal = %d \n", c.supernodal);
  Rprintf("final_asis = %d \n", c.final_asis);
  Rprintf("final_super = %d \n", c.final_super);
  Rprintf("final_ll = %d \n", c.final_ll);
  Rprintf("final_pack = %d \n", c.final_pack);
  Rprintf("final_monotonic = %d \n", c.final_monotonic);
  Rprintf("final_resymbol = %d \n", c.final_resymbol);
  Rprintf("zrelax = (%lf, %lf, %lf) \n", c.zrelax[0], c.zrelax[1], c.zrelax[2]);
  Rprintf("nrelax = (%d, %d, %d) \n", c.nrelax[0], c.nrelax[1], c.nrelax[2]);
  Rprintf("prefer_zomplex = %d \n", c.prefer_zomplex);
  Rprintf("prefer_upper = %d \n", c.prefer_upper);
  Rprintf("quick_return_if_not_posdef = %d \n", c.quick_return_if_not_posdef);
  Rprintf("prefer_binary = %d \n", c.prefer_binary);
  Rprintf("/* printing and error handling options */ \n");
  Rprintf("print = %d \n", c.print);
  Rprintf("precise = %d \n", c.precise);
  Rprintf("try_catch = %d \n", c.try_catch);
  Rprintf("/* ordering options */ \n");
  Rprintf("nmethods = %d \n", c.nmethods);
  Rprintf("current = %d \n", c.current);
  Rprintf("selected = %d \n", c.selected);
  // TODO print options for all methods (struct)
  Rprintf("postorder = %d \n", c.postorder);
  Rprintf("default_nesdis = %d \n", c.default_nesdis);
  Rprintf("/* METIS workarounds */ \n");
  Rprintf("metis_memory = %lf \n", c.metis_memory);
  Rprintf("metis_dswitch = %lf \n", c.metis_dswitch);
  Rprintf("metis_nswitch = %d \n", c.metis_nswitch);
  Rprintf("/* workspace */ \n");
  Rprintf("nrow = %d \n", c.nrow);
  Rprintf("iworksize = %d \n", c.iworksize);
  Rprintf("xworksize = %d \n", c.xworksize);
  Rprintf("itype = %d \n", c.itype);
  Rprintf("dtype = %d \n", c.dtype);
  Rprintf("no_workspace_reallocate = %d \n", c.no_workspace_reallocate);
  Rprintf("/* statistics */ \n");
  Rprintf("status = %d \n", c.status);
  Rprintf("fl = %lf \n", c.fl);
  Rprintf("lnz = %lf \n", c.lnz);
  Rprintf("anz = %lf \n", c.anz);
  Rprintf("modfl = %lf \n", c.modfl);
  Rprintf("malloc_count = %d \n", c.malloc_count);
  Rprintf("memory_usage = %d \n", c.memory_usage);
  Rprintf("memory_inuse = %d \n", c.memory_inuse);
  Rprintf("nrealloc_col = %lf \n", c.nrealloc_col);
  Rprintf("nrealloc_factor = %lf \n", c.nrealloc_factor);
  Rprintf("ndbounds_hit = %lf \n", c.ndbounds_hit);
  Rprintf("rowfacfl = %lf \n", c.rowfacfl);
  Rprintf("aatfl = %lf \n", c.aatfl);
  Rprintf("called_nd = %d \n", c.called_nd);
  Rprintf("blas_ok = %d \n", c.blas_ok);
}
