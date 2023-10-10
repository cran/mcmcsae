#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "Matrix.h"

cholmod_common c;


/* .Call calls */
extern SEXP _mcmcsae_add_diagC(void *, void *);
extern SEXP _mcmcsae_add_vector(void *, void *);
extern SEXP _mcmcsae_Cbacksolve(void *, void *);
extern SEXP _mcmcsae_CbacksolveM(void *, void *);
extern SEXP _mcmcsae_cCHM_dsC_Cholesky(void *, void *, void *, void *, void *);
extern SEXP _mcmcsae_cCHM_update_inplace(void *, void *, void *);
extern SEXP _mcmcsae_cCHMf_solve(void *, void *, void *);
extern SEXP _mcmcsae_cCHMf_solve_matrix(void *, void *, void *);
extern SEXP _mcmcsae_cCHMf_spsolve(void *, void *, void *);
extern SEXP _mcmcsae_Ccholesky(void *);
extern SEXP _mcmcsae_Ccreate_sparse_crossprod_sym_template(void *, void *, void *, void *);
extern SEXP _mcmcsae_Cdense_crossprod_sym(void *, void *);
extern SEXP _mcmcsae_Cdense_crossprod_sym0(void *);
extern SEXP _mcmcsae_Cdense_crossprod_sym2(void *, void *);
extern SEXP _mcmcsae_Cdense_dense_crossprod(void *, void *);
extern SEXP _mcmcsae_Cdense_dense_prod(void *, void *);
extern SEXP _mcmcsae_Cdense_diag_crossprod(void *, void *);
extern SEXP _mcmcsae_Cdense_diag_prod(void *, void *);
extern SEXP _mcmcsae_Cdense_kron(void *, void *);
extern SEXP _mcmcsae_Cdense_numeric_crossprod(void *, void *);
extern SEXP _mcmcsae_Cdense_numeric_prod(void *, void *);
extern SEXP _mcmcsae_Cdense_sparse_crossprod(void *, void *);
extern SEXP _mcmcsae_Cdense_sparse_prod(void *, void *);
extern SEXP _mcmcsae_Cdense_sparse_tcrossprod(void *, void *);
extern SEXP _mcmcsae_Cdense_sparseS_prod(void *, void *);
extern SEXP _mcmcsae_Cdense_tab_tcrossprod(void *, void *);
extern SEXP _mcmcsae_Cdiag(void *);
extern SEXP _mcmcsae_Cdiag_sparse_prod(void *, void *);
extern SEXP _mcmcsae_CdiagU(void *);
extern SEXP _mcmcsae_Cforwardsolve(void *, void *);
extern SEXP _mcmcsae_CforwardsolveM(void *, void *);
extern SEXP _mcmcsae_Cnnz_per_col_scps_template(void *, void *, void *);
extern SEXP _mcmcsae_copy_vector(void *);
extern SEXP _mcmcsae_CrCRT(void *, void *, void *);
extern SEXP _mcmcsae_Crepgen(void *, void *, void *);
extern SEXP _mcmcsae_Crgig(void *, void *, void *, void *);
extern SEXP _mcmcsae_Crnorm(void *, void *, void *);
extern SEXP _mcmcsae_CrPGapprox(void *, void *, void *, void *);
extern SEXP _mcmcsae_Crtmvn_Gibbs_dense(void *, void *, void *, void *);
extern SEXP _mcmcsae_Crtmvn_Gibbs_sparse(void *, void *, void *, void *);
extern SEXP _mcmcsae_Crtmvn_slice_Gibbs_dense(void *, void *, void *, void *);
extern SEXP _mcmcsae_Crtmvn_slice_Gibbs_sparse(void *, void *, void *, void *);
extern SEXP _mcmcsae_CrTNprobit(void *, void *);
extern SEXP _mcmcsae_Crtuvn(void *, void *);
extern SEXP _mcmcsae_Cscale_dense(void *, void *);
extern SEXP _mcmcsae_Cscale_sparse(void *, void *);
extern SEXP _mcmcsae_Csparse_crossprod_sym(void *, void *);
extern SEXP _mcmcsae_Csparse_crossprod_sym2(void *, void *);
extern SEXP _mcmcsae_Csparse_dense_crossprod(void *, void *);
extern SEXP _mcmcsae_Csparse_dense_crossprod_sym(void *, void *);
extern SEXP _mcmcsae_Csparse_dense_prod(void *, void *);
extern SEXP _mcmcsae_Csparse_diag_crossprod_sym(void *, void *);
extern SEXP _mcmcsae_Csparse_numeric_crossprod(void *, void *);
extern SEXP _mcmcsae_Csparse_numeric_prod(void *, void *);
extern SEXP _mcmcsae_Csparse_sym_twist(void *, void *);
extern SEXP _mcmcsae_CsparseS_dense_prod(void *, void *);
extern SEXP _mcmcsae_CsparseS_numeric_prod(void *, void *);
extern SEXP _mcmcsae_Ctab(void *, void *, void *, void *, void *);
extern SEXP _mcmcsae_Ctab_dense_crossprod(void *, void *);
extern SEXP _mcmcsae_Ctab_dense_prod(void *, void *);
extern SEXP _mcmcsae_Ctab_numeric_crossprod(void *, void *);
extern SEXP _mcmcsae_Ctab_numeric_prod(void *, void *, void *);
extern SEXP _mcmcsae_Ctab_unary_crossprod(void *);
extern SEXP _mcmcsae_Ctab2dgC(void *);
extern SEXP _mcmcsae_Ctab2mat(void *);
extern SEXP _mcmcsae_diagC(void *);
extern SEXP _mcmcsae_dotprodC(void *, void *);
extern SEXP _mcmcsae_fast_aggrC(void *, void *, void *);
extern SEXP _mcmcsae_inverseSPD(void *);
extern SEXP _mcmcsae_log1pexpC(void *);
extern SEXP _mcmcsae_mv_update(void *, void *, void *, void *);
extern SEXP _mcmcsae_prec2se_cor(void *);
extern SEXP _mcmcsae_sparse_sum_x(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _mcmcsae_TMVN_HMC_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_mcmcsae_add_diagC",                             (DL_FUNC) &_mcmcsae_add_diagC,                              2},
    {"_mcmcsae_add_vector",                            (DL_FUNC) &_mcmcsae_add_vector,                             2},
    {"_mcmcsae_Cbacksolve",                            (DL_FUNC) &_mcmcsae_Cbacksolve,                             2},
    {"_mcmcsae_CbacksolveM",                           (DL_FUNC) &_mcmcsae_CbacksolveM,                            2},
    {"_mcmcsae_cCHM_dsC_Cholesky",                     (DL_FUNC) &_mcmcsae_cCHM_dsC_Cholesky,                      5},
    {"_mcmcsae_cCHM_update_inplace",                   (DL_FUNC) &_mcmcsae_cCHM_update_inplace,                    3},
    {"_mcmcsae_cCHMf_solve",                           (DL_FUNC) &_mcmcsae_cCHMf_solve,                            3},
    {"_mcmcsae_cCHMf_solve_matrix",                    (DL_FUNC) &_mcmcsae_cCHMf_solve_matrix,                     3},
    {"_mcmcsae_cCHMf_spsolve",                         (DL_FUNC) &_mcmcsae_cCHMf_spsolve,                          3},
    {"_mcmcsae_Ccholesky",                             (DL_FUNC) &_mcmcsae_Ccholesky,                              1},
    {"_mcmcsae_Ccreate_sparse_crossprod_sym_template", (DL_FUNC) &_mcmcsae_Ccreate_sparse_crossprod_sym_template,  4},
    {"_mcmcsae_Cdense_crossprod_sym",                  (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym,                   2},
    {"_mcmcsae_Cdense_crossprod_sym0",                 (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym0,                  1},
    {"_mcmcsae_Cdense_crossprod_sym2",                 (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym2,                  2},
    {"_mcmcsae_Cdense_dense_crossprod",                (DL_FUNC) &_mcmcsae_Cdense_dense_crossprod,                 2},
    {"_mcmcsae_Cdense_dense_prod",                     (DL_FUNC) &_mcmcsae_Cdense_dense_prod,                      2},
    {"_mcmcsae_Cdense_diag_crossprod",                 (DL_FUNC) &_mcmcsae_Cdense_diag_crossprod,                  2},
    {"_mcmcsae_Cdense_diag_prod",                      (DL_FUNC) &_mcmcsae_Cdense_diag_prod,                       2},
    {"_mcmcsae_Cdense_kron",                           (DL_FUNC) &_mcmcsae_Cdense_kron,                            2},
    {"_mcmcsae_Cdense_numeric_crossprod",              (DL_FUNC) &_mcmcsae_Cdense_numeric_crossprod,               2},
    {"_mcmcsae_Cdense_numeric_prod",                   (DL_FUNC) &_mcmcsae_Cdense_numeric_prod,                    2},
    {"_mcmcsae_Cdense_sparse_crossprod",               (DL_FUNC) &_mcmcsae_Cdense_sparse_crossprod,                2},
    {"_mcmcsae_Cdense_sparse_prod",                    (DL_FUNC) &_mcmcsae_Cdense_sparse_prod,                     2},
    {"_mcmcsae_Cdense_sparse_tcrossprod",              (DL_FUNC) &_mcmcsae_Cdense_sparse_tcrossprod,               2},
    {"_mcmcsae_Cdense_sparseS_prod",                   (DL_FUNC) &_mcmcsae_Cdense_sparseS_prod,                    2},
    {"_mcmcsae_Cdense_tab_tcrossprod",                 (DL_FUNC) &_mcmcsae_Cdense_tab_tcrossprod,                  2},
    {"_mcmcsae_Cdiag",                                 (DL_FUNC) &_mcmcsae_Cdiag,                                  1},
    {"_mcmcsae_Cdiag_sparse_prod",                     (DL_FUNC) &_mcmcsae_Cdiag_sparse_prod,                      2},
    {"_mcmcsae_CdiagU",                                (DL_FUNC) &_mcmcsae_CdiagU,                                 1},
    {"_mcmcsae_Cforwardsolve",                         (DL_FUNC) &_mcmcsae_Cforwardsolve,                          2},
    {"_mcmcsae_CforwardsolveM",                        (DL_FUNC) &_mcmcsae_CforwardsolveM,                         2},
    {"_mcmcsae_Cnnz_per_col_scps_template",            (DL_FUNC) &_mcmcsae_Cnnz_per_col_scps_template,             3},
    {"_mcmcsae_copy_vector",                           (DL_FUNC) &_mcmcsae_copy_vector,                            1},
    {"_mcmcsae_CrCRT",                                 (DL_FUNC) &_mcmcsae_CrCRT,                                  3},
    {"_mcmcsae_Crepgen",                               (DL_FUNC) &_mcmcsae_Crepgen,                                3},
    {"_mcmcsae_Crgig",                                 (DL_FUNC) &_mcmcsae_Crgig,                                  4},
    {"_mcmcsae_Crnorm",                                (DL_FUNC) &_mcmcsae_Crnorm,                                 3},
    {"_mcmcsae_CrPGapprox",                            (DL_FUNC) &_mcmcsae_CrPGapprox,                             4},
    {"_mcmcsae_Crtmvn_Gibbs_dense",                    (DL_FUNC) &_mcmcsae_Crtmvn_Gibbs_dense,                     4},
    {"_mcmcsae_Crtmvn_Gibbs_sparse",                   (DL_FUNC) &_mcmcsae_Crtmvn_Gibbs_sparse,                    4},
    {"_mcmcsae_Crtmvn_slice_Gibbs_dense",              (DL_FUNC) &_mcmcsae_Crtmvn_slice_Gibbs_dense,               4},
    {"_mcmcsae_Crtmvn_slice_Gibbs_sparse",             (DL_FUNC) &_mcmcsae_Crtmvn_slice_Gibbs_sparse,              4},
    {"_mcmcsae_CrTNprobit",                            (DL_FUNC) &_mcmcsae_CrTNprobit,                             2},
    {"_mcmcsae_Crtuvn",                                (DL_FUNC) &_mcmcsae_Crtuvn,                                 2},
    {"_mcmcsae_Cscale_dense",                          (DL_FUNC) &_mcmcsae_Cscale_dense,                           2},
    {"_mcmcsae_Cscale_sparse",                         (DL_FUNC) &_mcmcsae_Cscale_sparse,                          2},
    {"_mcmcsae_Csparse_crossprod_sym",                 (DL_FUNC) &_mcmcsae_Csparse_crossprod_sym,                  2},
    {"_mcmcsae_Csparse_crossprod_sym2",                (DL_FUNC) &_mcmcsae_Csparse_crossprod_sym2,                 2},
    {"_mcmcsae_Csparse_dense_crossprod",               (DL_FUNC) &_mcmcsae_Csparse_dense_crossprod,                2},
    {"_mcmcsae_Csparse_dense_crossprod_sym",           (DL_FUNC) &_mcmcsae_Csparse_dense_crossprod_sym,            2},
    {"_mcmcsae_Csparse_dense_prod",                    (DL_FUNC) &_mcmcsae_Csparse_dense_prod,                     2},
    {"_mcmcsae_Csparse_diag_crossprod_sym",            (DL_FUNC) &_mcmcsae_Csparse_diag_crossprod_sym,             2},
    {"_mcmcsae_Csparse_numeric_crossprod",             (DL_FUNC) &_mcmcsae_Csparse_numeric_crossprod,              2},
    {"_mcmcsae_Csparse_numeric_prod",                  (DL_FUNC) &_mcmcsae_Csparse_numeric_prod,                   2},
    {"_mcmcsae_Csparse_sym_twist",                     (DL_FUNC) &_mcmcsae_Csparse_sym_twist,                      2},
    {"_mcmcsae_CsparseS_dense_prod",                   (DL_FUNC) &_mcmcsae_CsparseS_dense_prod,                    2},
    {"_mcmcsae_CsparseS_numeric_prod",                 (DL_FUNC) &_mcmcsae_CsparseS_numeric_prod,                  2},
    {"_mcmcsae_Ctab",                                  (DL_FUNC) &_mcmcsae_Ctab,                                   5},
    {"_mcmcsae_Ctab_dense_crossprod",                  (DL_FUNC) &_mcmcsae_Ctab_dense_crossprod,                   2},
    {"_mcmcsae_Ctab_dense_prod",                       (DL_FUNC) &_mcmcsae_Ctab_dense_prod,                        2},
    {"_mcmcsae_Ctab_numeric_crossprod",                (DL_FUNC) &_mcmcsae_Ctab_numeric_crossprod,                 2},
    {"_mcmcsae_Ctab_numeric_prod",                     (DL_FUNC) &_mcmcsae_Ctab_numeric_prod,                      3},
    {"_mcmcsae_Ctab_unary_crossprod",                  (DL_FUNC) &_mcmcsae_Ctab_unary_crossprod,                   1},
    {"_mcmcsae_Ctab2dgC",                              (DL_FUNC) &_mcmcsae_Ctab2dgC,                               1},
    {"_mcmcsae_Ctab2mat",                              (DL_FUNC) &_mcmcsae_Ctab2mat,                               1},
    {"_mcmcsae_diagC",                                 (DL_FUNC) &_mcmcsae_diagC,                                  1},
    {"_mcmcsae_dotprodC",                              (DL_FUNC) &_mcmcsae_dotprodC,                               2},
    {"_mcmcsae_fast_aggrC",                            (DL_FUNC) &_mcmcsae_fast_aggrC,                             3},
    {"_mcmcsae_inverseSPD",                            (DL_FUNC) &_mcmcsae_inverseSPD,                             1},
    {"_mcmcsae_log1pexpC",                             (DL_FUNC) &_mcmcsae_log1pexpC,                              1},
    {"_mcmcsae_mv_update",                             (DL_FUNC) &_mcmcsae_mv_update,                              4},
    {"_mcmcsae_prec2se_cor",                           (DL_FUNC) &_mcmcsae_prec2se_cor,                            1},
    {"_mcmcsae_sparse_sum_x",                          (DL_FUNC) &_mcmcsae_sparse_sum_x,                           9},
    {"_mcmcsae_TMVN_HMC_C",                            (DL_FUNC) &_mcmcsae_TMVN_HMC_C,                            14},
    {NULL, NULL, 0}
};

void R_init_mcmcsae(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  
  M_R_cholmod_start(&c);
}

void R_unload_mcmcsae(DllInfo *dll) {
  M_cholmod_finish(&c);
}
