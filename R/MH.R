
# functions for Metropolis-Hastings updates


###############################
# degrees of freedom parameter

# RW proposal on log(df) with scale tau
draw_df_MH_RW <- function(r, df.current, Q.current, df.mod) {
  # draw new degrees of freedom value from proposal
  df.star <- exp(rnorm(1L, sd=df.mod$tau)) * df.current
  # compute log-acceptance-ratio
  log_ar <- r * (  0.5 * df.star * log(0.5 * df.star) - 0.5 * df.current * log(0.5 * df.current)
                   + lgamma(0.5 * df.current) - lgamma(0.5 * df.star) ) +
    df.mod$alpha0 * log(df.star / df.current) +
    (df.current - df.star) * (df.mod$beta0 + 0.5 * sum(Q.current - log(Q.current)))
  if (log(runif(1L)) < log_ar) df.star else df.current
}

# random walk with drift: Metropolis Adjusted Langevin Algorithm
draw_df_MH_mala <- function(r, df.current, Q.current, df.mod) {
  # compute mala drift
  temp <- df.mod$beta0 + 0.5 * sum(Q.current - log(Q.current))
  drift <- df.mod$alpha0 + 0.5 * r * df.current * (1 - log(0.5 * df.current) - digamma(0.5 * df.current)) - df.current * temp
  # draw new degrees of freedom value from mala proposal
  df.star <- exp(df.mod$tau * rnorm(1L, mean=0.5 * drift)) * df.current
  # compute log-acceptance-ratio
  log_ar <- r * ( 0.5 * df.star * log(0.5 * df.star) - 0.5 * df.current * log(0.5 * df.current)
                  + lgamma(0.5 * df.current) - lgamma(0.5 * df.star) ) +
    + (df.current - df.star) * temp + df.mod$alpha0 * log(df.star/df.current)
  if (log(runif(1L)) < log_ar) df.star else df.current
}


###################
# Leroux parameter

# draw Leroux parameter and also update
# p: current state, mc: model component, coef_raw: raw coefficients (vector if q0=1, otherwise matrix)
draw_Leroux <- function(p, mc, coef_raw) {
  L <- p[[mc$name_Leroux]]  # current value
  # draw a candidate value (proposal)
  #   and initialize log acceptance probability
  if (mc$Leroux_type == "default") {
    L.star <- runif(1L)
    log_ar <- 0
  } else {
    L.star <- rbeta(1L, mc$Leroux$a.star, mc$Leroux$b.star)
    log_ar <- (mc$Leroux$a - mc$Leroux$a.star) * log(L.star / L) + (mc$Leroux$b - mc$Leroux$b.star) * log((1 - L.star) / (1 - L))
  }
  # compute log full conditional posterior in L.star
  Ldiff <- L.star - L
  QA_diff <- mc$mat_sum_Leroux(mc$QA, mc$idL, Ldiff, -Ldiff)
  tr_diff <- switch(mc$var,
    unstructured = sum(crossprod_sym(coef_raw, QA_diff) * p[[mc$name_Qraw]]),
    diagonal = sum(.colSums(coef_raw * (QA_diff %*% coef_raw), mc$l, mc$q0) * p[[mc$name_Qraw]]),
    scalar =
      if (mc$q0 == 1L)
        dotprodC(coef_raw, QA_diff %m*v% coef_raw) * p[[mc$name_Qraw]]
      else
        sum(coef_raw * (QA_diff %*% coef_raw)) * p[[mc$name_Qraw]]
  )
  # for now, assume a uniform prior on Leroux parameter
  detQA.star <- mc$det(2*L.star, 1 - 2*L.star)
  log_ar <- log_ar + 0.5 * ( mc$q0 * (detQA.star - p[[mc$name_detQA]]) - tr_diff )
  if (log(runif(1L)) < log_ar) {  # accept
    p[[mc$name_Leroux]] <- L.star
    p[[mc$name_detQA]] <- detQA.star
  }
  p
}
