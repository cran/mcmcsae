
has_response <- function(formula) length(formula) == 3L

get_response <- function(formula, data=NULL) {
  tf <- terms(formula, data=data)
  ind <- attr(tf, "response")
  if (ind > 0L) {
    vars <- as.list(attr(tf, "variables"))[-1L]
    eval(vars[[ind]], data, environment(formula))
  } else {
    NULL
  }
}

get_offset <- function(formula, data=NULL) {
  if (is.integer(data) && length(data) == 1L) data <- NULL
  tf <- terms(formula, data=data)
  ind <- attr(tf, "offset")
  if (is.null(ind)) return(NULL)
  if (length(ind) > 1L) stop("only a single offset allowed in 'formula'")
  vars <- as.list(attr(tf, "variables"))[-1L]
  eval(vars[[ind]], data, environment(formula))
}

get_vars <- function(formula, rhs.only=TRUE) {
  tf <- terms(formula)
  vars <- as.list(attr(tf, "variables"))[-1L]
  if (rhs.only) {
    if (attr(tf, "response") > 0L || !is.null(attr(tf, "offset"))) {
      if (attr(tf, "response") > 0L) {
        vars <- vars[-c(attr(tf, "response"), attr(tf, "offset"))]
      } else {
        vars <- vars[-attr(tf, "offset")]
      }
    }
  }
  vars
}

get_types <- function(formula, specials=c("gen", "mec", "reg", "bart")) {
  vars <- get_vars(formula)
  if (!length(vars)) return(NULL)
  sapply(vars, function(x) match.arg(as.character(x[[1L]]), specials))
}

get_coefnames <- function(formula, specials=c("gen", "mec", "reg", "bart"), data=NULL) {
  vars <- get_vars(formula)
  if (!length(vars)) return(NULL)
  parnames <- sapply(vars, function(x) if (is.null(x$name)) NA_character_ else x$name)
  types <- sapply(vars, function(x) match.arg(as.character(x[[1L]]), specials))
  parnames[is.na(parnames)] <- paste0(types[is.na(parnames)], which(is.na(parnames)))
  check_mod_names(parnames)
  parnames
}

has_explicit_intercept <- function(formula) {
  fstr <- tail(as.character(formula), 1L)
  grepl("(^|\\+)\\s*\\(*\\s*1\\s*\\)*\\s*(\\+|$)", fstr)
}

standardize_formula <- function(formula, specials=c("reg", "mec", "gen", "bart"), default="reg", data=NULL) {
  # interpret everything not in special terms as a default component
  tf <- terms(formula, keep.order=TRUE, specials=specials, data=data)
  # NB ~ . - var does not warn if var is not in data
  idx <- unlist(attr(tf, "specials"), use.names=FALSE)  # variable indices of special terms
  if (length(idx)) {
    fac <- attr(tf, "factors")
    for (i in seq_along(idx)) {
      term.idx <- which(fac[idx[i], ] > 0)  # translate to term indices
      if (length(term.idx) != 1L) stop("cannot parse formula")
      idx[i] <- term.idx
    }
    remainder <- attr(tf, "term.labels")[-idx]
  } else {
    remainder <- attr(tf, "term.labels")
  }
  e <- environment(formula)
  if (length(remainder) || has_explicit_intercept(formula) || (!length(idx) && attr(tf, "intercept") == 1L)) {
    if (length(remainder)) {
      if (attr(tf, "intercept") == 0L) {
        formula <- paste0(default, "( ~ 0 +", paste(remainder, collapse=" + "), ")")
      } else {
        formula <- paste0(default, "( ~", paste(remainder, collapse=" + "), ")")
      }
    } else {
      formula <- paste0(default, "( ~ 1)")
    }
    if (length(idx)) {
      funpart <- paste(attr(tf, "term.labels")[idx], collapse=" + ")
      formula <- paste(formula, funpart, sep=" + ")
    }
    if (!is.null(attr(tf, "offset"))) {
      if (length(attr(tf, "offset")) > 1L) stop("only one offset allowed")
      formula <- paste(rownames(attr(tf, "factors"))[[attr(tf, "offset")]], "+", formula)
    }
    if (attr(tf, "response") > 0L) {
      formula <- as.formula(paste(deparse(attr(tf, "variables")[[attr(tf, "response") + 1L]]), "~", formula), env=e)
    } else {
      formula <- as.formula(paste("~", formula), env=e)
    }
  }
  formula
}

to_mclist <- function(formula) {
  vars <- get_vars(formula)
  parnames <- get_coefnames(formula)
  mod <- list()
  for (m in seq_along(vars)) {
    mc <- vars[[m]]
    if (is.null(names(mc)) || names(mc)[2L] == "") names(mc)[2L] <- "formula"  # first argument is formula
    mod[[parnames[m]]] <- mc
  }
  mod
}

to_Vmclist <- function(formula) {
  vars <- get_vars(formula)
  parnames <- get_coefnames(formula, c("vreg", "vfac"))
  mod <- list()
  for (m in seq_along(vars)) {
    mc <- vars[[m]]
    #if (is.null(mc$factor)) mc$factor <- "local_"  # by default apply a local factor
    mod[[parnames[m]]] <- mc
  }
  mod
}
