#' #' Split weights across geographies
#'
#' \code{geoweight} calculates state weights for each household in a microdata
#' file that add up to the household total weight, such that weighted state
#' totals for selected characteristics hit or come close to desired targets
#'
#' \code{geoweight} uses the solver \code{\link[nleqslv]{nleqslv}} or the solver
#' \code{\link[minpack.lm]{nls.lm}} depending on user choice.
#'
#' The default method, LM, uses \code{nls.lm} as it appears to be the most robust
#' of the methods, rarely failing and often producing a better optimum than
#' Broyden or Newton. However, in some circumstances one of the latter may work
#' better. It is hard to define guidelines for when a particular method will be
#' better. The Broyden method can be faster or more robust than the Newton method
#' but generally requires many more iterations than the Newton method, although
#' iterations will be faster.
#'
#'
#' @param wh Household weights, 1 per household, numeric vector length h. Each
#'   household's geography weights must sum to its household weight.
#' @param xmat Data for households. Matrix with 1 row per household and 1 column
#'   per characteristic (h x k matrix). Columns can be named.
#' @param targets Targeted values. Matrix with 1 row per geographic area and 1
#'   column per characteristic. If columns are named, names must match column
#'   names of `xmat`. Rownames can be used to identify geographic areas. If
#'   unnamed, rows will be named geo1, geo2, ..., geo_s
#' @param dweights Difference weights: weights to be applied to Weighting
#'   factors for targets (h x k matrix).
#' @param betavec optional vector of initial guess at parameters, length s * k;
#'   default is zero for all
#' @param method optional parameter for approach to use; must be one of c('LM',
#'   'Broyden', 'Newton'); default is 'LM'
#' @param maxiter maximum number of iterations; integer; defaults vary by method:
#'   LM (default):        200
#'   Broyden:            2000
#'   Newton:              200
#' @param optlist list of options that will update nelsqv or nls.lm options
#'   respectively
#' @param quiet c(TRUE, FALSE) FALSE is default; TRUE provides newlsqv or nls.lm
#'   output
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#'   \item{solver_message}{message from the solver that was used}
#'   \item{etime}{elapsed time}
#'   \item{beta_opt_mat}{s x k matrix of optimal parameters}
#'   \item{whs}{h x s matrix of state weights for each household, computed
#'     using the optimal parameters}
#'   \item{wh}{the input vector of household total weights, length h}
#'   \item{xmat}{matrix of data for households, h x k}
#'   \item{dweights}{optional vector of weighting factors for targets, length s * k}
#'   \item{output}{list of output from the solver that was used}
#' }
#' @examples
#' # Example 1: Determine state weights for a simple problem with random data
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#'
#' res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw)
#'
#' res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, method = 'Newton')
#'
#' res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, method = 'Broyden')
#'
#' res1
#' res2
#' res3
#' c(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
#'
#' # verify that the state weights produce the desired targets
#' t(res1$whs) %*% p$xmat
#' p$targets
#'
#' # Example 2: Determine state weights for a problem using ACS data
#' # In this case, we actually know a set of weights that work (because they
#' # are in the ACS) but we want to see if we can construct state weights
#' # as if we did not know the true weights.
#' library(tidyverse)
#' data(acs)
#' glimpse(acs)
#'
#' # create a subset that has just one income group (incgroup ==5), with
#' # variables we want to target, including indicator variables where we want
#' # to target the weighted number of records for which a particular variable
#' data_df <- acs %>%
#'   filter(incgroup == 5) %>%
#'   select(stabbr, pwgtp, pincp, wagp, ssip) %>%
#'   # create the indicator variables
#'   mutate(nrecs = 1, # indicator used for number of records
#'          wagp_nnz = (wagp != 0) * 1.,
#'          ssip_nnz = (ssip != 0) * 1.)
#' data_df # 1,000 records
#'
#' wh <- data_df$pwgtp
#'
#' targets_df <- data_df %>%
#'   pivot_longer(-c(pwgtp, stabbr)) %>%
#'   mutate(wtd_value = value * pwgtp) %>%
#'   group_by(stabbr, name) %>%
#'   summarise(wtd_value = sum(wtd_value), .groups = "drop") %>%
#'   pivot_wider(values_from = wtd_value)
#' targets_df
#'
#' targets <- targets_df %>%
#'   select(-stabbr) %>%
#'   as.matrix
#' rownames(targets) <- targets_df$stabbr
#' targets
#'
#' xmat <- data_df %>%
#'   select(all_of(colnames(targets))) %>%
#'   as.matrix
#' xmat
#'
#' opts_LM <- list(ptol = 1e-8, ftol = 1e-8)
#' resx2 <- geoweight(wh = wh, xmat = xmat, targets = targets, optlist = opts_LM, quiet = TRUE)
#' names(resx2)
#' resx2$solver_message
#' resx2$h; resx2$s; resx2$k
#' resx2$method
#' resx2$etime
#' resx2$sse_unweighted
#' resx2$sse_weighted
#' resx2$targets
#' resx2$targets_calc
#' resx2$targets_diff %>% round(1)
#' resx2$targets_pctdiff %>% round(1)
#'
#' # let's compare the estimated weights to the actual weights, which we know in this case
#' # (but would not know in a real-world application)
#' df <- tibble(wh = resx2$wh) %>%
#'   cbind(resx2$whs) %>%
#'   mutate(person = row_number()) %>%
#'   pivot_longer(cols = any_of(state.abb),
#'                names_to = "stabbr_solved",
#'                values_to = "weight_solved") %>%
#'  left_join(data_df %>%
#'              select(stabbr_true = stabbr, pwgtp) %>%
#'              mutate(person = row_number()),
#'              by = "person")
#'
#' # Example 3: Repeat the above, but add some random noise to the targets to
#' # make them harder to hit. Nothing else is changed.
#'
#' set.seed(1234)
#' targetsx3 <- targets * (1 + rnorm(length(targets), mean=0, sd=.01))
#' targetsx3
#'
#' resx3 <- geoweight(wh = wh, xmat = xmat, targets = targetsx3, optlist = opts_LM, quiet = TRUE)
#' resx3$solver_message
#' resx3$h; resx3$s; resx3$k
#' resx3$method
#' resx3$etime
#' resx3$sse_unweighted
#' resx3$sse_weighted
#' resx3$targets
#' resx3$targets_calc
#' resx3$targets_diff %>% round(1)
#' resx3$targets_pctdiff %>% round(1)
#'
#' # Note that this was harder to solve and we do not hit the targets exactly.
#'
#' @export
geoweight <- function(wh,
                      xmat,
                      targets,
                      dweights = get_dweights(targets),
                      betavec = rep(0, length(targets)),
                      method = 'LM',
                      maxiter = NULL,
                      optlist = NULL,
                      quiet = TRUE){

  args <- as.list(environment()) # gets explicit and default arguments
  stopifnot(method %in% c('LM', 'Broyden', 'Newton'))

  # dweights will be set so that targets are normalized to 100, unless alternative
  #   difference weights are supplied
  # betavec will set to 0 unless supplied
  # Global choices under the Newton and Broyden methods:
  # c("dbldog", "pwldog", "cline", "qline", "gline", "hook", "none")

  if(!is.null(maxiter) & !(is.null(optlist))) {
    print("CAUTION: maxiter and optlist both supplied. maxiter will override any iteration limit included in optlist")
  }

  opts_LM <- minpack.lm::nls.lm.control(maxiter = 200,
                                        nprint = 1,
                                        factor = 100,
                                        ptol = .001)

  opts_Broyden <- list(maxit = 2000,
                       trace = 1,
                       allowSingular = TRUE,
                       ftol = 1e-04)

  opts_Newton <- list(maxit = 200,
                      trace = 1,
                      allowSingular = TRUE,
                      ftol = 1e-04)

  if(quiet){
    opts_LM$nprint <- 0
    opts_Broyden$trace <- 0
    opts_Newton$trace <- 0
  }

  h <- nrow(xmat)
  s <- nrow(targets)
  k <- ncol(targets)

  snames <- rownames(targets)
  knames <- colnames(targets)

  t1 <- proc.time()
  if (method == "LM") {
    opts_LM <- purrr::list_modify(opts_LM, !!!optlist) # splice lists
    if(!is.null(maxiter)) opts_LM$maxiter <- maxiter
    output <- minpack.lm::nls.lm(par = betavec,
                                 fn = diff_vec,
                                 control = opts_LM,
                                 wh = wh, xmat = xmat, targets = targets,
                                 dweights = dweights)
    beta_opt <- output$par
    opts_used <- opts_LM
  } else if (method == "Broyden") {
    opts_Broyden <- purrr::list_modify(opts_Broyden, !!!optlist) # splice lists
    if(!is.null(maxiter)) opts_Broyden$maxit <- maxiter
    output <- nleqslv::nleqslv(x = betavec, fn = diff_vec, jac = NULL,
                               wh = wh, xmat = xmat, targets = targets,
                               dweights = dweights,
                               method = c("Broyden"),
                               global = c("dbldog"), xscalm = c("auto"),
                               jacobian = FALSE,
                               control = opts_Broyden)
    beta_opt <- output$x
    opts_used <- opts_Broyden

    } else if (method == "Newton") {
      opts_Newton <- purrr::list_modify(opts_Newton, !!!optlist) # splice lists
      if(!is.null(maxiter)) opts_Newton$maxit <- maxiter
      output <- nleqslv::nleqslv(x = betavec, fn = diff_vec, jac = NULL,
                                 wh = wh, xmat = xmat, targets = targets,
                                 dweights = dweights,
                                 method = c("Newton"),
                                 global = c("dbldog"), xscalm = c("auto"), # dbldog pwldog cline
                                 jacobian = FALSE,
                                 control = opts_Newton)
      beta_opt <- output$x
      opts_used <- opts_Newton
    }
  t2 <- proc.time()

  # define additional values to be returned
  solver_message <- output$message
  etime <- t2 - t1
  sse_unweighted <- sse(beta_opt, wh, xmat, targets)
  sse_weighted <- sse(beta_opt, wh, xmat, targets, dweights)
  beta_opt_mat <- vtom(beta_opt, s)
  rownames(beta_opt_mat) <- snames
  colnames(beta_opt_mat) <- knames

  whs <- get_weights(beta_opt_mat,
                     get_delta(wh, beta_opt_mat, xmat),
                     xmat)
  colnames(whs) <- snames

  targets_calc <- t(whs) %*% xmat
  targets_diff <- targets_calc - targets
  targets_pctdiff <- targets_diff / targets * 100

  keepnames <- c("h", "s", "k", "method", "opts_used", "solver_message", "etime",
                 "sse_unweighted", "sse_weighted",
                 "beta_opt_mat",
                 "targets", "targets_calc", "targets_diff", "targets_pctdiff",
                 "whs", "wh", "xmat", "dweights", "output")
  result <- list()
  for(var in keepnames) result[[var]] <- get(var)

  result
}
