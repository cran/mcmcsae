
check_ny <- function(ny, data) {
  if (is.null(ny)) return(1L)  # Bernoulli or categorical
  if (is.character(ny)) {
    if (length(ny) != 1L) stop("wrong input for 'ny'")
    ny <- eval_in(ny, data)
  } else {
    if (!is.numeric(ny)) stop("wrong input for 'ny'")
    if (is.integer(data) && length(data) == 1L) {
      # data interpreted as sample size
      n <- data
    } else {
      n <- nrow(data)
    }
    if (all(length(ny) != c(1L, n))) stop("wrong length for 'ny'")
    if (anyNA(ny)) stop("missings in 'ny' not allowed")
    if (any(ny < 0)) stop("'ny' cannot be negative")
  }
  if (any(ny - as.integer(ny) > sqrt(.Machine$double.eps))) warn("non-integral values in 'ny' have been rounded")
  as.integer(round(ny))
  # or should we allow non-integral ny? allowed for model fitting but not for prediction by rbinom
  # note that y in create_sampler is currently not rounded, but y <= ny is checked
  # maybe add argument round and round only in case of prediction
}

check_ry <- function(ry, data) {
  if (is.null(ry)) return(1L)  # default value
  if (is.character(ry)) {
    if (length(ry) != 1L) stop("wrong input for 'ry'")
    ry <- eval_in(ry, data)
  } else {
    if (!is.numeric(ry)) stop("wrong input for 'ry'")
    if (is.integer(data) && length(data) == 1L) {
      # data interpreted as sample size
      n <- data
    } else {
      n <- nrow(data)
    }
    if (all(length(ry) != c(1L, n))) stop("wrong length for 'ry'")
    if (anyNA(ry)) stop("missings in 'ry' not allowed")
    if (any(ry <= 0)) stop("'ry' must be positive")
  }
  ry
}
