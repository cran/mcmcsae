#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "Matrix.h"

cholmod_common c;


/* .Call calls */
extern SEXP _mcmcsae_add_diagC(SEXP, SEXP);
extern SEXP _mcmcsae_Cbacksolve(SEXP, SEXP);
extern SEXP _mcmcsae_CbacksolveM(SEXP, SEXP);
extern SEXP _mcmcsae_cCHM_dsC_Cholesky(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_cCHM_options();
extern SEXP _mcmcsae_cCHM_update_inplace(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_cCHMf_solve(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_cCHMf_solve_matrix(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_cCHMf_spsolve(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_crossprod_sym(SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_crossprod_sym0(SEXP);
extern SEXP _mcmcsae_Cdense_crossprod_sym2(SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_diag_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_kron(SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_numeric_crossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Cdense_numeric_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Cdiag(SEXP);
extern SEXP _mcmcsae_Cdiag_sparse_prod(SEXP, SEXP);
extern SEXP _mcmcsae_CdiagU(SEXP);
extern SEXP _mcmcsae_Cforwardsolve(SEXP, SEXP);
extern SEXP _mcmcsae_CforwardsolveM(SEXP, SEXP);
extern SEXP _mcmcsae_Cmatrix_sparse_tcrossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Cmatrix_sparseS_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Cmatrix_tab_tcrossprod(SEXP, SEXP);
extern SEXP _mcmcsae_CrCRT(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Crepgen(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Crgig(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Crnorm(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_CrPGapprox(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Crtmvn_Gibbs(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_CrTNprobit(SEXP, SEXP);
extern SEXP _mcmcsae_Crtuvn(SEXP, SEXP);
extern SEXP _mcmcsae_Cscale_dense(SEXP, SEXP);
extern SEXP _mcmcsae_Cscale_sparse(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_crossprod_sym(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_crossprod_sym2(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_dense_crossprod_sym(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_diag_crossprod_sym(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_matrix_crossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_matrix_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_numeric_crossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_numeric_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Csparse_sym_twist(SEXP, SEXP);
extern SEXP _mcmcsae_CsparseS_matrix_prod(SEXP, SEXP);
extern SEXP _mcmcsae_CsparseS_numeric_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Ctab(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Ctab_matrix_crossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Ctab_matrix_prod(SEXP, SEXP);
extern SEXP _mcmcsae_Ctab_numeric_crossprod(SEXP, SEXP);
extern SEXP _mcmcsae_Ctab_numeric_prod(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_Ctab_unary_crossprod(SEXP);
extern SEXP _mcmcsae_Ctab2dgC(SEXP);
extern SEXP _mcmcsae_Ctab2mat(SEXP);
extern SEXP _mcmcsae_diagC(SEXP);
extern SEXP _mcmcsae_dotprodC(SEXP, SEXP);
extern SEXP _mcmcsae_fast_aggrC(SEXP, SEXP, SEXP);
extern SEXP _mcmcsae_get_col_dgC(SEXP, SEXP);
extern SEXP _mcmcsae_inverseSPD(SEXP);
extern SEXP _mcmcsae_log1pexpC(SEXP);
extern SEXP _mcmcsae_prec2se_cor(SEXP);
extern SEXP _mcmcsae_sparse_sum_x(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mcmcsae_add_diagC",                   (DL_FUNC) &_mcmcsae_add_diagC,                   2},
    {"_mcmcsae_Cbacksolve",                  (DL_FUNC) &_mcmcsae_Cbacksolve,                  2},
    {"_mcmcsae_CbacksolveM",                 (DL_FUNC) &_mcmcsae_CbacksolveM,                 2},
    {"_mcmcsae_cCHM_dsC_Cholesky",           (DL_FUNC) &_mcmcsae_cCHM_dsC_Cholesky,           6},
    {"_mcmcsae_cCHM_options",                (DL_FUNC) &_mcmcsae_cCHM_options,                0},
    {"_mcmcsae_cCHM_update_inplace",         (DL_FUNC) &_mcmcsae_cCHM_update_inplace,         3},
    {"_mcmcsae_cCHMf_solve",                 (DL_FUNC) &_mcmcsae_cCHMf_solve,                 3},
    {"_mcmcsae_cCHMf_solve_matrix",          (DL_FUNC) &_mcmcsae_cCHMf_solve_matrix,          3},
    {"_mcmcsae_cCHMf_spsolve",               (DL_FUNC) &_mcmcsae_cCHMf_spsolve,               3},
    {"_mcmcsae_Cdense_crossprod_sym",        (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym,        2},
    {"_mcmcsae_Cdense_crossprod_sym0",       (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym0,       1},
    {"_mcmcsae_Cdense_crossprod_sym2",       (DL_FUNC) &_mcmcsae_Cdense_crossprod_sym2,       2},
    {"_mcmcsae_Cdense_diag_prod",            (DL_FUNC) &_mcmcsae_Cdense_diag_prod,            2},
    {"_mcmcsae_Cdense_kron",                 (DL_FUNC) &_mcmcsae_Cdense_kron,                 2},
    {"_mcmcsae_Cdense_numeric_crossprod",    (DL_FUNC) &_mcmcsae_Cdense_numeric_crossprod,    2},
    {"_mcmcsae_Cdense_numeric_prod",         (DL_FUNC) &_mcmcsae_Cdense_numeric_prod,         2},
    {"_mcmcsae_Cdiag",                       (DL_FUNC) &_mcmcsae_Cdiag,                       1},
    {"_mcmcsae_Cdiag_sparse_prod",           (DL_FUNC) &_mcmcsae_Cdiag_sparse_prod,           2},
    {"_mcmcsae_CdiagU",                      (DL_FUNC) &_mcmcsae_CdiagU,                      1},
    {"_mcmcsae_Cforwardsolve",               (DL_FUNC) &_mcmcsae_Cforwardsolve,               2},
    {"_mcmcsae_CforwardsolveM",              (DL_FUNC) &_mcmcsae_CforwardsolveM,              2},
    {"_mcmcsae_Cmatrix_sparse_tcrossprod",   (DL_FUNC) &_mcmcsae_Cmatrix_sparse_tcrossprod,   2},
    {"_mcmcsae_Cmatrix_sparseS_prod",        (DL_FUNC) &_mcmcsae_Cmatrix_sparseS_prod,        2},
    {"_mcmcsae_Cmatrix_tab_tcrossprod",      (DL_FUNC) &_mcmcsae_Cmatrix_tab_tcrossprod,      2},
    {"_mcmcsae_CrCRT",                       (DL_FUNC) &_mcmcsae_CrCRT,                       3},
    {"_mcmcsae_Crepgen",                     (DL_FUNC) &_mcmcsae_Crepgen,                     3},
    {"_mcmcsae_Crgig",                       (DL_FUNC) &_mcmcsae_Crgig,                       4},
    {"_mcmcsae_Crnorm",                      (DL_FUNC) &_mcmcsae_Crnorm,                      3},
    {"_mcmcsae_CrPGapprox",                  (DL_FUNC) &_mcmcsae_CrPGapprox,                  4},
    {"_mcmcsae_Crtmvn_Gibbs",                (DL_FUNC) &_mcmcsae_Crtmvn_Gibbs,                4},
    {"_mcmcsae_CrTNprobit",                  (DL_FUNC) &_mcmcsae_CrTNprobit,                  2},
    {"_mcmcsae_Crtuvn",                      (DL_FUNC) &_mcmcsae_Crtuvn,                      2},
    {"_mcmcsae_Cscale_dense",                (DL_FUNC) &_mcmcsae_Cscale_dense,                2},
    {"_mcmcsae_Cscale_sparse",               (DL_FUNC) &_mcmcsae_Cscale_sparse,               2},
    {"_mcmcsae_Csparse_crossprod_sym",       (DL_FUNC) &_mcmcsae_Csparse_crossprod_sym,       2},
    {"_mcmcsae_Csparse_crossprod_sym2",      (DL_FUNC) &_mcmcsae_Csparse_crossprod_sym2,      2},
    {"_mcmcsae_Csparse_dense_crossprod_sym", (DL_FUNC) &_mcmcsae_Csparse_dense_crossprod_sym, 2},
    {"_mcmcsae_Csparse_diag_crossprod_sym",  (DL_FUNC) &_mcmcsae_Csparse_diag_crossprod_sym,  2},
    {"_mcmcsae_Csparse_matrix_crossprod",    (DL_FUNC) &_mcmcsae_Csparse_matrix_crossprod,    2},
    {"_mcmcsae_Csparse_matrix_prod",         (DL_FUNC) &_mcmcsae_Csparse_matrix_prod,         2},
    {"_mcmcsae_Csparse_numeric_crossprod",   (DL_FUNC) &_mcmcsae_Csparse_numeric_crossprod,   2},
    {"_mcmcsae_Csparse_numeric_prod",        (DL_FUNC) &_mcmcsae_Csparse_numeric_prod,        2},
    {"_mcmcsae_Csparse_sym_twist",           (DL_FUNC) &_mcmcsae_Csparse_sym_twist,           2},
    {"_mcmcsae_CsparseS_matrix_prod",        (DL_FUNC) &_mcmcsae_CsparseS_matrix_prod,        2},
    {"_mcmcsae_CsparseS_numeric_prod",       (DL_FUNC) &_mcmcsae_CsparseS_numeric_prod,       2},
    {"_mcmcsae_Ctab",                        (DL_FUNC) &_mcmcsae_Ctab,                        5},
    {"_mcmcsae_Ctab_matrix_crossprod",       (DL_FUNC) &_mcmcsae_Ctab_matrix_crossprod,       2},
    {"_mcmcsae_Ctab_matrix_prod",            (DL_FUNC) &_mcmcsae_Ctab_matrix_prod,            2},
    {"_mcmcsae_Ctab_numeric_crossprod",      (DL_FUNC) &_mcmcsae_Ctab_numeric_crossprod,      2},
    {"_mcmcsae_Ctab_numeric_prod",           (DL_FUNC) &_mcmcsae_Ctab_numeric_prod,           3},
    {"_mcmcsae_Ctab_unary_crossprod",        (DL_FUNC) &_mcmcsae_Ctab_unary_crossprod,        1},
    {"_mcmcsae_Ctab2dgC",                    (DL_FUNC) &_mcmcsae_Ctab2dgC,                    1},
    {"_mcmcsae_Ctab2mat",                    (DL_FUNC) &_mcmcsae_Ctab2mat,                    1},
    {"_mcmcsae_diagC",                       (DL_FUNC) &_mcmcsae_diagC,                       1},
    {"_mcmcsae_dotprodC",                    (DL_FUNC) &_mcmcsae_dotprodC,                    2},
    {"_mcmcsae_fast_aggrC",                  (DL_FUNC) &_mcmcsae_fast_aggrC,                  3},
    {"_mcmcsae_get_col_dgC",                 (DL_FUNC) &_mcmcsae_get_col_dgC,                 2},
    {"_mcmcsae_inverseSPD",                  (DL_FUNC) &_mcmcsae_inverseSPD,                  1},
    {"_mcmcsae_log1pexpC",                   (DL_FUNC) &_mcmcsae_log1pexpC,                   1},
    {"_mcmcsae_prec2se_cor",                 (DL_FUNC) &_mcmcsae_prec2se_cor,                 1},
    {"_mcmcsae_sparse_sum_x",                (DL_FUNC) &_mcmcsae_sparse_sum_x,                9},
    {NULL, NULL, 0}
};

void R_init_mcmcsae(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    
    M_R_cholmod_start(&c);
    // TODO maybe add error handler and set other cholmod options
}

void R_unload_mcmcsae(DllInfo *dll) {
    M_cholmod_finish(&c);
}
