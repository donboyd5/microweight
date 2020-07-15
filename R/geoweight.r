#' Split weights across geographies
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
#' @param wh vector of household total weights, length h (see h, s, k definitions below)
#' @param xmat h x k matrix of data for households
#' @param targets s x k matrix of desired target values
#' @param dweights optional vector of weighting factors for targets, length s * k
#' @param method optional parameter for approach to use; must be one of
#' c('LM', 'Broyden', 'Newton'); default is 'LM'
#' @param betavec optional vector of initial guess at parameters, length s * k;
#'   default is zero for all
#' @param maxiter integer; defaults: Broyden (2000), Newton (200), LM (200)
#' @param opts list of options that will update nelsqv or nls.lm options respectively
#' @param quiet c(TRUE, FALSE) FALSE is default; TRUE provides newlsqv or nls.lm output
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
#'   \item{output}{list of output from the solver that was used}
#' }
#' @examples
#' # Example 1: Determine state weights for a simple problem with random data
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#'
#' res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, quiet=TRUE)
#'
#' res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, method = 'Newton', quiet=TRUE)
#'
#' res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, method = 'LM', quiet=TRUE)
#'
#' res1
#' res2
#' res3
#' c(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
#'
#' # verify that the state weights produce the desired targets
#' t(res2$whs) %*% p$xmat
#' p$targets
#' @export
geoweight <- function(wh, xmat, targets, dweights = get_dweights(targets),
                      method = 'LM', betavec = rep(0, length(targets)),
                      maxiter = NULL, opts = NULL, quiet = FALSE){

  # dweights will be set so that targets are normalized to 100, unless alternative
  #   difference weights are supplied
  # betavec will set to 0 unless supplied
  # Global choices under the Newton and Broyden methods:
  # c("dbldog", "pwldog", "cline", "qline", "gline", "hook", "none")

  if(!is.null(maxiter) & !(is.null(opts))) {
    print("CAUTION: maxiter and opts both supplied. maxiter will override any iteration limit included in opts.")
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
    opts_LM <- purrr::list_modify(opts_LM, !!!opts) # splice lists
    if(!is.null(maxiter)) opts_LM$maxiter <- maxiter
    output <- minpack.lm::nls.lm(par = betavec,
                                 fn = diff_vec,
                                 control = opts_LM,
                                 wh = wh, xmat = xmat, targets = targets,
                                 dweights = dweights)
    beta_opt <- output$par
    opts_used <- opts_LM
  } else if (method == "Broyden") {
    opts_Broyden <- purrr::list_modify(opts_Broyden, !!!opts) # splice lists
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
      opts_Newton <- purrr::list_modify(opts_Newton, !!!opts) # splice lists
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
                 "whs", "wh", "output")
  result <- list()
  for(var in keepnames) result[[var]] <- get(var)

  result
}
