################################################################################
# EDITS unm_glm FUNCTION OF R PACKAGE unmconf VERSION 1.0.0
# CALLS THE NEW FUNCTION MyUnm_glm

# EDIT 1: When one unmeasured confounder is specified U=(U_1 ), the function incorrectly omits U1 from the substantive analysis of interest and fits Y|X,C. We modified the function so that it correctly fits Y|X,C,U1.
# LINE 35: placed ! in the correct place
# EDIT 2: The function inserts a line-break in the middle of the list (after 50L) of covariates for glm. It is a problem when the covariate names are long and/or there are many confounders. We modified the function so that it is inserted after 500L (maximum allowed).
# LINES 37, 59, 116: increased width.cutoff
# EDIT 3: We modified the function to allow us to fix seed numbers so that our results can be reproduced.
# LINE 183: added ... to allow additional arguments
################################################################################

MyUnm_glm <- function (form1, form2 = NULL, form3 = NULL, family1 = binomial(), 
                       family2 = NULL, family3 = NULL, data, n.iter = 2000, n.adapt = 1000, 
                       thin = 1, n.chains = 4, filename = tempfile(fileext = ".jags"), 
                       quiet = getOption("unm_quiet"), progress.bar = getOption("unm_progress.bar"), 
                       code_only = FALSE, priors, response_nuisance_priors, response_params_to_track, 
                       confounder1_nuisance_priors, confounder1_params_to_track, 
                       confounder2_nuisance_priors, confounder2_params_to_track, 
                       ...) 
{
  g <- glue::glue
  if (!inherits(family2, "family")) 
    family2 <- list(family = "none")
  if (!inherits(family3, "family")) 
    family3 <- list(family = "none")
  y <- deparse(form1[[2]])
  u1 <- if (inherits(form2, "formula")) 
    deparse(form2[[2]])
  else NULL
  u2 <- if (inherits(form3, "formula")) 
    deparse(form3[[2]])
  else NULL
  if ((!is.null(u1) || !is.null(u2)) && grepl(paste(g("\\b{u1}\\b"), # edited
                                                    g("\\b{u2}\\b"), sep = "|"),
                                              deparse(form1[[3]], width.cutoff=500L))) { # edited
    conf_piece <- "+ inprod(U[i,], lambda)"
  }
  else {
    conf_piece <- ""
  }
  response_model_code <- switch(family1$family, gaussian = g("{y}[i] ~ dnorm(mu_{y}[i], tau_{y})\n{'    '}mu_{y}[i] <- inprod(X[i,], beta) {conf_piece}"), 
                                binomial = g("{y}[i] ~ dbern(p_{y}[i])\n{'    '}logit(p_{y}[i]) <- inprod(X[i,], beta) {conf_piece}"), 
                                Gamma = g("{y}[i] ~ dgamma(alpha_{y}, d[i])\n{'    '}d[i] <- alpha_{y} / mu_{y}[i]\n{'    '}log(mu_{y}[i]) <- inprod(X[i,], beta) {conf_piece}"), 
                                poisson = g("{y}[i] ~ dpois(mu_{y}[i])\n{'    '}log(mu_{y}[i]) <- log_e[i] + inprod(X[i,], beta) {conf_piece}"))
  if (missing(response_nuisance_priors)) {
    response_nuisance_priors <- switch(family1$family, gaussian = g("tau_{y} ~ dgamma(0.001, 0.001)"), 
                                       binomial = " ", Gamma = g("alpha_{y} ~ dgamma(.1, .1)"), 
                                       poisson = " ")
    response_params_to_track <- switch(family1$family, gaussian = g("tau_{y}"), 
                                       binomial = character(0), Gamma = g("alpha_{y}"), 
                                       poisson = character(0))
  }
  else {
    response_nuisance_priors <- g(response_nuisance_priors)
    response_params_to_track <- g(response_params_to_track)
  }
  if (!is.null(u2) && grepl(paste(g("\\b{u2}\\b")), deparse(form2[[3]], width.cutoff=500L))) { # edited
    conf2_piece <- "+ inprod(U[i, 2], zeta)"
  }
  else {
    conf2_piece <- ""
  }
  confounder1_model_code <- switch(family2$family, gaussian = g("U[i,1] ~ dnorm(mu_{u1}[i], tau_{u1})\n                   {'    '}mu_{u1}[i] <- inprod(W[i,], gamma) {conf2_piece}"), 
                                   binomial = g("U[i,1] ~ dbern(p_{u1}[i])\n                   {'    '}logit(p_{u1}[i]) <- inprod(W[i,], gamma) {conf2_piece}"), 
                                   none = "")
  if (missing(confounder1_nuisance_priors)) {
    confounder1_nuisance_priors <- switch(family2$family, 
                                          gaussian = g("tau_{u1} ~ dgamma(0.001, 0.001)"), 
                                          binomial = " ", none = "")
    confounder1_params_to_track <- switch(family2$family, 
                                          gaussian = g("tau_{u1}"), binomial = character(0), 
                                          none = character(0))
  }
  else {
    confounder1_nuisance_priors <- g(confounder1_nuisance_priors)
    confounder1_params_to_track <- g(confounder1_params_to_track)
  }
  confounder2_model_code <- switch(family3$family, gaussian = g("U[i,2] ~ dnorm(mu_{u2}[i], tau_{u2})\n                 {'    '}mu_{u2}[i] <- inprod(V[i,], delta)"), 
                                   binomial = g("U[i,2] ~ dbern(p_{u2}[i])\n                 {'    '}logit(p_{u2}[i]) <- inprod(V[i,], delta)"), 
                                   none = "")
  if (missing(confounder2_nuisance_priors)) {
    confounder2_nuisance_priors <- switch(family3$family, 
                                          gaussian = g("tau_{u2} ~ dgamma(0.001, 0.001)"), 
                                          binomial = " ", none = "")
    confounder2_params_to_track <- switch(family3$family, 
                                          gaussian = g("tau_{u2}"), binomial = character(0), 
                                          none = character(0))
  }
  else {
    confounder2_nuisance_priors <- g(confounder2_nuisance_priors)
    confounder2_params_to_track <- g(confounder2_params_to_track)
  }
  current_na_action <- getOption("na.action")
  on.exit(options(na.action = current_na_action), add = TRUE)
  options(na.action = "na.pass")
  (X <- model.matrix(form1, data = data))
  if (!is.null(u1)) 
    (W <- model.matrix(form2, data = data))
  else (W <- NULL)
  if (!is.null(u2)) 
    (V <- model.matrix(form3, data = data))
  else (V <- NULL)
  if (is.null(u1)) {
    (U <- NULL)
  }
  else if (is.null(u2)) {
    (U <- X[, c(u1), drop = FALSE])
  }
  else {
    (U <- X[, c(u1, u2), drop = FALSE])
  }
  (X <- X[, setdiff(colnames(X), colnames(U)), drop = FALSE])
  if (!is.null(u2) && !is.null(W) && grepl(paste(g("\\b{u2}\\b")), 
                                           deparse(form2[[3]], width.cutoff=500L))) { # edited
    (U2 <- W[, c(u2), drop = FALSE])
  }
  else {
    U2 <- NULL
  }
  (W <- W[, setdiff(colnames(W), colnames(U2)), drop = FALSE])
  p_be <- ncol(X)
  p_la <- ncol(U)
  p_ga <- ncol(W)
  p_ze <- ncol(U2)
  p_de <- ncol(V)
  X_vars <- colnames(X)
  U_vars <- colnames(U)
  W_vars <- colnames(W)
  U2_vars <- colnames(U2)
  V_vars <- colnames(V)
  pretty_X_vars <- gsub("\\(Intercept\\)", "1", X_vars)
  pretty_W_vars <- gsub("\\(Intercept\\)", "1", W_vars)
  pretty_V_vars <- gsub("\\(Intercept\\)", "1", V_vars)
  jags_coefs <- c(g("beta[{1:p_be}]"), if (!is.null(p_la)) g("lambda[{1:p_la}]"), 
                  if (!is.null(p_ga)) g("gamma[{1:p_ga}]"), if (!is.null(p_ze)) g("zeta[{1:p_ze}]"), 
                  if (!is.null(p_de)) g("delta[{1:p_de}]"))
  real_coefs <- c(g("beta[{X_vars}]"), g("lambda[{U_vars}]"), 
                  g("gamma[{W_vars}]"), g("zeta[{U2_vars}]"), g("delta[{V_vars}]"))
  real_coefs <- gsub("\\(Intercept\\)", "1", real_coefs)
  real_to_jags_coefs <- function(x) {
    dict <- structure(jags_coefs, names = real_coefs)
    x[x %in% real_coefs] <- unname(dict[x])[x %in% real_coefs]
    x
  }
  jags_to_real_coefs <- function(x) {
    dict <- structure(real_coefs, names = jags_coefs)
    x[x %in% jags_coefs] <- unname(dict[x])[x %in% jags_coefs]
    x
  }
  default_prior1 <- switch(family1$family, gaussian = "dnorm(0, .001)", 
                           binomial = "dnorm(0, .1)", Gamma = "dnorm(0, .001)", 
                           poisson = "dnorm(0, .1)", none = "")
  default_prior2 <- switch(family2$family, gaussian = "dnorm(0, .001)", 
                           binomial = "dnorm(0, .1)", none = "")
  default_prior3 <- switch(family3$family, gaussian = "dnorm(0, .001)", 
                           binomial = "dnorm(0, .1)", none = "")
  prior_struct1 <- rep(default_prior1, length(jags_coefs[grepl("(beta|lambda)", 
                                                               jags_coefs)]))
  prior_struct2 <- rep(default_prior2, length(jags_coefs[grepl("(gamma|zeta)", 
                                                               jags_coefs)]))
  prior_struct3 <- rep(default_prior3, length(jags_coefs[grepl("(delta)", 
                                                               jags_coefs)]))
  jags_priors <- c(prior_struct1, prior_struct2, prior_struct3)
  names(jags_priors) <- jags_coefs
  if (!missing(priors)) 
    jags_priors[real_to_jags_coefs(names(priors))] <- priors
  jags_priors <- g("{names(jags_priors)} ~ {jags_priors} \t\t# = {real_coefs}")
  jd <- list(n = nrow(X), X = X, U = U, W = W, V = V, log_e = model.offset(model.frame(form1, 
                                                                                       data = data)))
  jd[[y]] <- data[[y]]
  jd <- drop_nulls(jd)
  code <- g("\n  model {\n    # models\n    for (i in 1:n) {\n      {{response_model_code}}\n\n      {{confounder1_model_code}}\n      {{confounder2_model_code}}\n    }\n\n    # priors\n    {{paste(jags_priors, collapse = '\n    ')}}\n    {{response_nuisance_priors}}\n    {{confounder1_nuisance_priors}}\n    {{confounder2_nuisance_priors}}\n  }\n  ", 
            .open = "{{", .close = "}}")
  if (code_only) {
    cat(code)
    cat("\n")
    return(invisible())
  }
  cat(code, file = filename)
  jm <- jags.model(filename, data = jd, n.adapt = n.adapt, 
                   n.chains = n.chains, quiet = quiet, ...) # edited
  params_of_interest <- unique(sub("\\[.+\\]", "", real_coefs))
  samps <- coda.samples(jm, c(params_of_interest, response_params_to_track, 
                              confounder1_params_to_track, confounder2_params_to_track), 
                        n.iter = n.iter, thin = thin, progress.bar = progress.bar)
  names <- dimnames(samps[[1]])[[2]]
  names[names == "beta"] <- "beta[1]"
  names[names == "lambda"] <- "lambda[1]"
  names[names == "gamma"] <- "gamma[1]"
  names[names == "delta"] <- "delta[1]"
  names[names == "zeta"] <- "zeta[1]"
  names <- jags_to_real_coefs(names)
  names <- gsub("(\\w+)_(\\w+)", "\\1[\\2]", names)
  samps <- lapply(samps, function(mcmc) {
    dimnames(mcmc) <- list(NULL, names)
    mcmc
  })
  class(samps) <- c("unm_int", "unm_mod", "mcmc.list")
  attr(samps, "file") <- filename
  attr(samps, "code") <- code
  attr(samps, "form1") <- form1
  attr(samps, "family1") <- family1
  attr(samps, "form2") <- form2
  attr(samps, "family2") <- family2
  attr(samps, "form3") <- form3
  attr(samps, "family3") <- family3
  attr(samps, "call") <- match.call()
  attr(samps, "jm") <- jm
  attr(samps, "n.iter") <- n.iter
  attr(samps, "n.adapt") <- n.adapt
  attr(samps, "thin") <- thin
  attr(samps, "n.chains") <- n.chains
  samps
}
