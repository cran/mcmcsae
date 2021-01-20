
// RcppEigen.h also pulls in Rcpp.h
#include <RcppEigen.h>

using namespace Rcpp;


// For code that never uses random number generation,
// use option rng=false for a tiny speed-up


//’ Inverse of a symmetric positive definite (ordinary dense) matrix
//’ 
//’ @param M a symmetric positive definite matrix.
//’ @return The inverse of \code{M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd inverseSPD(const Eigen::Map<Eigen::MatrixXd> & M) {
  const int d = M.rows();
  if (M.cols() != d) stop("not a square matrix");
  return M.selfadjointView<Eigen::Upper>().ldlt().solve(Eigen::MatrixXd::Identity(d,d));
}

//’ Dense triangular solves, in particular solve of Lt system, with vector right hand side
//’
//’ @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//’ @param y a numeric vector.
//’ @return The solution of \code{Mx=y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cbacksolve(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & y) {
  if (M.cols() != y.size()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().solve(y);
}

//’ Dense triangular solves, in particular solve of L system, with vector right hand side.
//’
//’ @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//’ @param y a numeric vector.
//’ @return The solution of \code{M'x=y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cforwardsolve(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & y) {
  if (M.cols() != y.size()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().transpose().solve(y);
}

//’ Dense triangular solves, in particular solve of Lt system, with matrix right hand side
//’
//’ @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//’ @param y a numeric matrix.
//’ @return The solution of \code{Mx=y}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd CbacksolveM(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (M.cols() != y.rows()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().solve(y);
}

//’ Dense triangular solves, in particular solve of L system, with matrix right hand side
//’
//’ @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//’ @param y a numeric matrix.
//’ @return The solution of \code{M'x=y}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd CforwardsolveM(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (M.cols() != y.rows()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().transpose().solve(y);
}

//’ Inner product of two vectors
//’
//’ @param x a numeric vector.
//’ @param y a numeric vector.
//’ @return The inner product \code{x'y}.
// [[Rcpp::export(rng=false)]]
double dotprodC(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::Map<Eigen::VectorXd> & y) {
  if (x.size() != y.size()) stop("incompatible dimensions");
  return x.dot(y);
}

//’ Fast summation of a vector by group
//’
//’ @param x a numeric vector.
//’ @param group an integer vector of the same length as \code{x} defining the groups.
//’ @param n the number of groups, assumed to be labeled \code{1:n}.
//’ @return A vector with totals of \code{x} by group.
// [[Rcpp::export(rng=false)]]
NumericVector fast_aggrC(const NumericVector & x, const IntegerVector & group, int n) {
  const int xsize = x.size();
  if (xsize != group.size()) stop("incompatible dimensions");
  NumericVector out(n);
  for (int i = 0; i < xsize; i++) {
    out[group[i] - 1] += x[i];
  }
  return out;
}

//’ Matrix product of a sparse matrix with a vector
//’
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @param y a numeric vector.
//’ @return The matrix product \code{Ay}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Csparse_numeric_prod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> & y) {
  if (A.cols() != y.size()) stop("incompatible dimensions");
  return A * y;
}

//’ Matrix product of a symmetric sparse (dsC) matrix with a vector
//’
//’ @param A a numeric symmetric compressed, sparse, column-oriented matrix.
//’ @param y a numeric vector.
//’ @return The matrix product \code{Ay}.
// NB it is assumed that the symmetric matrix is stored in the upper triangular part
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd CsparseS_numeric_prod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> & y) {
  if (A.cols() != y.rows()) stop("incompatible dimensions");
  return A.selfadjointView<Eigen::Upper>() * y;
}

//’ Matrix product of a dense matrix with a vector
//’
//’ @param A a dense matrix.
//’ @param y a vector.
//’ @return The product \code{Ay}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cdense_numeric_prod(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & y) {
  if (A.cols() != y.size()) stop("incompatible dimensions");
  return A * y;
}

//’ Crossproduct of a dense matrix with a vector
//’
//’ @param A a dense matrix.
//’ @param y a vector.
//’ @return The product \code{A'y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cdense_numeric_crossprod(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & y) {
  if (A.rows() != y.size()) stop("incompatible dimensions");
  return A.adjoint() * y;
}

//’ Crossproduct of a sparse matrix with a vector
//’
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @param y a numeric vector.
//’ @return The crossproduct \code{A'y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Csparse_numeric_crossprod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> & y) {
  if (A.rows() != y.size()) stop("incompatible dimensions");
  return A.adjoint() * y;
}

//’ Matrix product of two dense matrices
//’
//’ @param A a numeric matrix.
//’ @param B a numeric matrix.
//’ @return The matrix product \code{AB}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cmatmat(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::MatrixXd> & B) {
  int n = A.cols();
  if (B.rows() != n) stop("incompatible matrices");
  return A * B;
}

//’ Matrix product of a sparse matrix with a matrix
//’
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @param y a matrix.
//’ @return The matrix product \code{Ay}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Csparse_matrix_prod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (A.cols() != y.rows()) stop("incompatible dimensions");
  return A * y;
}

//’ Matrix product of a matrix with a sparse matrix
//’
//’ @param A a matrix.
//’ @param B a numeric compressed, sparse, column-oriented matrix.
//’ @return The matrix product \code{AB}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cmatrix_sparse_prod(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::MappedSparseMatrix<double> & B) {
  if (A.cols() != B.rows()) stop("incompatible dimensions");
  return A * B;
}


//’ Matrix product of a symmetric sparse (dsC) matrix with a matrix
//’
//’ @param A a numeric symmetric compressed, sparse, column-oriented matrix.
//’ @param y a matrix.
//’ @return The matrix product \code{Ay}.
// NB it is assumed that the symmetric matrix is stored in the upper triangular part
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd CsparseS_matrix_prod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (A.cols() != y.rows()) stop("incompatible dimensions");
  return A.selfadjointView<Eigen::Upper>() * y;
}

//’ Matrix product of a matrix with a symmetric sparse (dsC) matrix
//’
//’ @param M a matrix.
//’ @param Q a numeric symmetric compressed, sparse, column-oriented matrix.
//’ @return The matrix product \code{MQ}.
// NB it is assumed that the symmetric matrix is stored in the upper triangular part
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cmatrix_sparseS_prod(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::MappedSparseMatrix<double> & Q) {
  if (M.cols() != Q.rows()) stop("incompatible dimensions");
  return M * Q.selfadjointView<Eigen::Upper>();
}

//’ Matrix product of a matrix with a diagonal matrix represented by a vector
//’
//’ @param M a dense matrix.
//’ @param d a numeric vector representing a diagonal matrix.
//’ @return The matrix product \code{M diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_diag_prod(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & d) {
  if (M.cols() != d.size()) stop("incompatible dimensions");
  return M * d.asDiagonal();
}

//’ Matrix product of the transpose of a sparse matrix with a matrix
//’
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @param y a matrix.
//’ @return The matrix product \code{A'y}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Csparse_matrix_crossprod(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (A.rows() != y.rows()) stop("incompatible dimensions");
  return A.transpose() * y;
}

//’ Product of a matrix with the transpose of a sparse matrix
//’
//’ @param y a (dense) matrix.
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @return The tcrossproduct \code{yA'}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cmatrix_sparse_tcrossprod(const Eigen::Map<Eigen::MatrixXd> & y, const Eigen::MappedSparseMatrix<double> & A) {
  if (y.cols() != A.cols()) stop("incompatible dimensions");
  return y * A.transpose();
}

//’ Product of a diagonal and sparse matrix
//’
//’ @param x a numeric vector, representing a diagonal matrix.
//’ @param A a numeric compressed, sparse, column-oriented matrix.
//’ @return The product \code{diag(x) A}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Cdiag_sparse_prod(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::MappedSparseMatrix<double> & A) {
  if (x.size() != A.rows()) stop("incompatible dimensions");
  return x.asDiagonal() * A;
}

// Used in closure for sparse matrix addition to compute x-slot.
// [[Rcpp::export(rng=false)]]
NumericVector sparse_sum_x(const int n,
      const IntegerVector & ind1, const IntegerVector & ind2,
      const NumericVector & M1x, const NumericVector & M2x,
      const bool UD1, const bool UD2,
      const double w1, const double w2) {
  NumericVector out(n);
  const int n1 = ind1.size();
  if (UD1) {
    for (int i = 0; i < n1; i++) {
      out[ind1[i]] = w1;
    }
  } else {
    for (int i = 0; i < n1; i++) {
      out[ind1[i]] = w1 * M1x[i];
    }
  }
  const int n2 = ind2.size();
  if (n2 > 0) {
    if (UD2) {
      for (int i = 0; i < n2; i++) {
        out[ind2[i]] += w2;
      }
    } else {
      for (int i = 0; i < n2; i++) {
        out[ind2[i]] += w2 * M2x[i];
      }
    }
  }
  return(out);
}

//’ Extract the diagonal of a dense matrix
//’
//’ @param A a numeric dense matrix.
//’ @return The diagonal of A as a vector.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd diagC(const Eigen::Map<Eigen::MatrixXd> & A) {
  return A.diagonal();
}

//’ Add a vector to the diagonal of a dense matrix
//’
//’ @param A a numeric dense matrix.
//’ @param d a numeric vector.
//’ @return Matrix A with x added to its diagonal.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd add_diagC(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & d) {
  if (d.size() != A.rows()) stop("incompatible dimensions");
  Eigen::MatrixXd out = A;
  out.diagonal() += d;
  return out;
}

//’ Compute the symmetric crossprod \code{M'M} for dense matrix \code{M}
//’
//’ @param M a numeric dense matrix.
//’ @return The product \code{M'M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_crossprod_sym0(const Eigen::Map<Eigen::MatrixXd> & M) {
  const int n = M.cols();
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(M.transpose());
}

//’ Compute the symmetric crossprod \code{M'diag(q)M} for dense matrix \code{M}
//’
//’ @param M a numeric dense matrix.
//’ @param q a numeric vector.
//’ @return The product \code{M'diag(q)M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_crossprod_sym(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & q) {
  if (M.rows() != q.size()) stop("incompatible input");
  const int n = M.cols();
  Eigen::VectorXd sq = q.array().sqrt();
  Eigen::MatrixXd sqM = sq.asDiagonal() * M;
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(sqM.transpose());
}

//’ Compute the symmetric crossprod \code{P'QP} for permutation matrix \code{P} and symmetric sparse Q
//’
//’ @param Q a symmetric sparse matrix, represented by its upper triangle.
//’ @param p an integer vector representing a permutation matrix P.
//’ @return The product \code{P'QP}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_sym_twist(const Eigen::MappedSparseMatrix<double> Q, const Eigen::Map<Eigen::VectorXi> & p) {
  if (Q.rows() != p.size()) stop("incompatible dimensions");
  Eigen::SparseMatrix<double> out;
  out = Q.selfadjointView<Eigen::Upper>().twistedBy( p.asPermutation() );
  return out.triangularView<Eigen::Upper>();
}

//’ Compute the crossprod of two dense matrices, known to result in a symmetric matrix
//’
//’ @param A a numeric dense matrix.
//’ @param B a numeric dense matrix, assumed to be of the form \code{QA} where \code{Q} is symmetric.
//’ @return The product \code{A'B}.
// [[Rcpp::export(rng=false)]]
NumericMatrix Cdense_crossprod_sym2(const NumericMatrix & A, const NumericMatrix & B) {
  const int n = A.ncol();
  if (B.ncol() != n) stop("incompatible dimensions");
  const int rank = A.nrow();
  if (B.nrow() != rank) stop("incompatible dimensions");
  NumericMatrix out = no_init(n, n);
  double z;
  int i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      z = 0;
      for (k = 0; k < rank - 4; k += 4) {
        z += (A(k,i) * B(k,j) + A(k+1,i) * B(k+1,j)) + (A(k+2,i) * B(k+2,j) + A(k+3,i) * B(k+3,j));
      }
      for (; k < rank; k++) {
        z += A(k,i) * B(k,j);
      }
      out(i,j) = z;
      out(j,i) = z;
    }
  }
  return out;
}

//’ Create a unit Diagonal Matrix
//’
//’ @param n an integer representing the dimension.
//’ @return The unit diagonal matrix of class \code{ddiMatrix} and dimension \code{n}.
// [[Rcpp::export(rng=false)]]
SEXP CdiagU(const int n) {
  S4 out("ddiMatrix");
  out.slot("Dim") = IntegerVector::create(n, n);
  out.slot("diag") = "U";
  return out;
}

//’ Create a Diagonal Matrix
//’
//’ @param x a numeric vector representing the matrix diagonal.
//’ @return The diagonal matrix of class \code{ddiMatrix} with diagonal \code{x}.
// [[Rcpp::export(rng=false)]]
SEXP Cdiag(const NumericVector x) {
  S4 out("ddiMatrix");
  const int n = x.size();
  out.slot("Dim") = IntegerVector::create(n, n);
  out.slot("x") = x;  // NB points to x, no copy
  return out;
}

//’ Compute scaled dense matrix
//’
//’ @param A a (dense) square matrix.
//’ @param d a vector of size the number of rows or columns of \code{A}, or a scalar.
//’ @return The sparse matrix \code{DAD} with \code{D=diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cscale_dense(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & d) {
  if (d.size() == 1) {
    return std::pow(d[0], 2) * A;
  } else {
    return d.asDiagonal() * A * d.asDiagonal();
  }
}

//’ Compute scaled sparse matrix
//’
//’ @param A a sparse matrix.
//’ @param d a vector of size the number of rows or columns of \code{A}, or a scalar.
//’ @return The sparse matrix \code{D'AD} with \code{D=diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Cscale_sparse(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> & d) {
  if (d.size() == 1) {
    return std::pow(d[0], 2) * A;
  } else {
    return d.asDiagonal() * A * d.asDiagonal();
  }
}

//’ Compute symmetric crossproduct of sparse and dense matrices
//’
//’ @param A a sparse matrix.
//’ @param Q a dense matrix, assumed to be symmetric.
//’ @return The (dense) matrix \code{A'QA}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Csparse_dense_crossprod_sym(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::MatrixXd> & Q) {
  if (A.rows() != Q.rows()) stop("incompatible dimensions");
  return (A.adjoint() * Q * A).selfadjointView<Eigen::Upper>();
}

//’ Compute symmetric crossproduct of sparse matrices
//’
//’ @param A a sparse matrix.
//’ @param Q a sparse matrix, assumed to be symmetric.
//’ @return The sparse matrix \code{A'QA}, in upper triangular view.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_crossprod_sym(const Eigen::MappedSparseMatrix<double> & A, const Eigen::MappedSparseMatrix<double> & Q) {
  if (A.rows() != Q.rows()) stop("incompatible dimensions");
  return (A.adjoint() * Q.selfadjointView<Eigen::Upper>() * A).triangularView<Eigen::Upper>();
}

//’ Compute symmetric 'weighted' crossproduct of sparse matrices
//’
//’ @param A a sparse matrix.
//’ @param Q a vector, representing a diagonal matrix.
//’ @return The sparse matrix \code{A'diag(Q)A}, in upper triangular view.
// NB here it seems faster to first assign QA to a temporary sparse matrix
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_diag_crossprod_sym(const Eigen::MappedSparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> & Q) {
  if (Q.size() != A.rows()) stop("incompatible dimensions");
  Eigen::SparseMatrix<double> temp = Q.asDiagonal() * A;
  return (A.adjoint() * temp).triangularView<Eigen::Upper>();
}

//’ Compute crossproduct of two sparse matrices, known to yield a symmetric (sparse) matrix
//’
//’ @param A a sparse matrix.
//’ @param B a sparse matrix, assumed to be the product of a symmetric matrix with \code{A}.
//’ @return The sparse matrix \code{A'B}, in upper triangular view.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_crossprod_sym2(const Eigen::MappedSparseMatrix<double> & A, const Eigen::MappedSparseMatrix<double> & B) {
  if (A.rows() != B.rows()) stop("incompatible dimensions");
  return (A.adjoint() * B).triangularView<Eigen::Upper>();
}

// [[Rcpp::export(rng=false)]]
Rcpp::List prec2se_cor(const Eigen::Map<Eigen::MatrixXd> & Q) {
  Eigen::MatrixXd V = inverseSPD(Q);
  Eigen::VectorXd se = V.diagonal().array().sqrt();
  Eigen::VectorXd se_inv = se.array().inverse();
  V = se_inv.asDiagonal() * V * se_inv.asDiagonal(); // overwrite V by the correlation matrix
  const int d = se.size();
  Eigen::VectorXd cor(d*(d-1)/2);
  int ind = 0;
  for (int i = 1; i < d; i++) {
    cor.segment(ind, i) = V.col(i).head(i);
    ind += i;
  }
  return Rcpp::List::create(Rcpp::Named("se") = se,
                            Rcpp::Named("cor") = cor);
}

//’ Extract a column from a sparse matrix
//’
//’ @param A a sparse matrix (\code{dgCMatrix}).
//’ @param j the (1-based) column index.
//’ @return A (dense) vector with the \code{j}th column of \code{A}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd get_col_dgC(const Eigen::MappedSparseMatrix<double> & A, const int j) {
  return A.col(j - 1);
}

//’ Compute log(1 + exp(x)) in a numerically robust way
//’
//’ @param x input vector.
//’ @return The computed vector log(1 + exp(x)).
//’ @references
//’   M. Maechler (2012)
//’     Accurately Computing log(1 − exp(− |a|)) Assessed by the Rmpfr package.
//’     URL: \url{https://cran.r-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
// [[Rcpp::export(rng=false)]]
NumericVector log1pexpC(const NumericVector & x) {
  const int n = x.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; i++) {
    if (x[i] <= -37) {
      out[i] = std::exp(x[i]);
    } else if (x[i] <= 18) {
      out[i] = std::log1p(std::exp(x[i]));
    } else if (x[i] <= 33.3) {
      out[i] = x[i] + std::exp(-x[i]);
    } else {
      out[i] = x[i];
    }
  }
  return out;
}

//’ Compute Kronecker product of two dense matrices
//’
//’ @param M1 numeric dense matrix.
//’ @param M2 numeric dense matrix.
//’ @return The Kronecker product of M1 and M2 as a dense matrix.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_kron(const Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::MatrixXd> & M2) {
  return kroneckerProduct(M1, M2);  //.selfadjointView<Eigen::Lower>() for symmetric matrices --> slower
}

//’ Generalisation of rep function useful for sparse kronecker product templates
//’
//’ @param v numeric vector.
//’ @param n integer vector.
//’ @param M2 numeric vector.
//’ @return A vector with replicated multiples of v's elements.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Crepgen(const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::VectorXi> & n, const Eigen::Map<Eigen::VectorXd> & M2) {
  const int q = M2.size();
  Eigen::VectorXd out(q * v.size());
  Eigen::VectorXd x;
  int ind = 0;
  int ind_out = 0;
  int ni;
  for (int i = 0; i < n.size(); i++) {
    ni = n.coeff(i);
    x = v.segment(ind, ni);
    for (int j = 0; j < q; j++) {
      out.segment(ind_out, ni) = M2.coeff(j) * x;
      ind_out += ni;
    }
    ind += ni;
  }
  return out;
}
