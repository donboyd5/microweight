#' Reweight a microdata file.
#'
# This is a comment. It will not show in the documentation.
#'
#' Calculate new weights for each household in a microdata file
#' so that (1) selected variables, weighted with the new weights and summed, hit
#' or come close to desired targets, and (2) a measure of distortion based on
#' how much the new weights differ from an initial set of weights is minimized.
#'
#'
#' @param iweights Numeric vector of initial household weights, length h.
#' @param targets Numeric vector of desired target values, length k.
#' @param target_names Character vector of names for targets, length k.
#' @param tol Numeric vector of additive tolerances for targets, length k. The
#'   solver will seek to hit the targets, plus or minus tol.
#' @param xmat Matrix of data for households, dimension h x k.
#' @param xlb Numeric vector of lower bounds for the ratio of new weights to
#'   initial weights. Default is 0.
#' @param xub Numeric vector of upper bounds for the ration of new weights to
#'   initial weights. Default is 50.
#' @param maxiter Integer value for maximum number of iterations. Default is 50.
#' @param optlist Named list of allowable
#'   \href{https://coin-or.github.io/Ipopt/OPTIONS.html}{IPOPT options}.
#'
#' @returns A list with the following elements:
#'
#' \describe{
#' \item{solver_message}{The message produced by IPOPT. See
#' \href{https://coin-or.github.io/Ipopt/OUTPUT.html}{IPOPT output}.}
#' \item{etime}{Elapsed time.} \item{objective}{The objective function value at
#' the solution.} \item{weights}{Numeric vector of new weights.}
#' \item{targets_df}{Data frame with target names and values, values at the
#' starting point, values at the solution, tolerances, differences, and percent
#' differences. The suffix diff indicates the difference from a target and the
#' suffix pdiff indicates the percent difference from a target.}
#' \item{result}{List with output from the solver that was used.}
#' }
#'
#' @details
#' \code{reweight} uses the \href{https://coin-or.github.io/Ipopt/}{IPOPT
#' solver} in the package \code{\link[ipoptr]{ipoptr}}. The problem is set up as
#' a nonlinear program with constraints. The constraints are the desired
#' targets. The user can set tolerances around these targets in which case they
#' are inequality constraints. By default the distortion measure to be minimized
#' is the sum of squared differences between the ratio of each new weight to the
#' corresponding initial weight and 1. The user can provide alternative
#' distortion measures.
#'
#' @examples
#' # Example 1: Determine new weights for a simple problem with random data
#' p <- make_problem(h=30, s=1, k=4)
#' # we'll have 30 households (30 weights) and 4 targets
#'
#' # break out needed components and prepare them for reweight
#' iweights <- p$wh
#'
#' targets <- as.vector(p$targets)
#' target_names <- paste0("targ", 1:length(targets))
#' # the targets created by make_problem are hit exactly when using the initial
#' # weights so perturb them slightly so that we have to find new weights
#' set.seed(1234)
#' noise <- rnorm(n=length(targets), mean = 0, sd = .02)
#' targets <- targets * (1 + noise)
#'
#' tol <- .005 * abs(targets)
#' xmat <- p$xmat
#' colnames(xmat) <- target_names
#'
#' res <- reweight(iweights = iweights, targets = targets,
#'                 target_names = target_names, tol = tol,
#'                 xmat = xmat)
#'
#' res
#'
#' # Example 2: Determine new weights for a small problem using ACS data
#' library(tidyverse)
#' data(acsbig)
#' data(acs_targets)

#' # let's focus on income group 5 and create and then try to hit targets for:
#' #    number of records (nrecs -- to be created based on the weight, pwgtp)
#' #    personal income (pincp)
#' #    wages (wagp)
#' #    number of people with wages (wagp_nnz -- to be created)
#' #    supplemental security income (ssip)
#' #    number of people with supplemental security income (ssip_nnz -- to
#'        be created)
#' # we also need to get pwgtp - the person weight for each record, which
#'        will be our initial weight
#' # for each "number of" variable we need to create an indicator variable that
#' # defines whether it is true for that record
#'
#' # get the data and prepare it
#' data_df <- acsbig %>%
#'   filter(incgroup == 5) %>%
#'   select(pwgtp, pincp, wagp, ssip) %>%
#'   # create the indicator variables
#'   mutate(nrecs = 1, # indicator used for number of records
#'          wagp_nnz = (wagp != 0) * 1.,
#'          ssip_nnz = (ssip != 0) * 1.)
#' data_df # 10,000 records
#'
#' # prepare targets: in practice we would get them from an external source but
#' # in this case we'll get actual sums on the file and perturb them randomly
#' # so that targets differ from initial sums
#' set.seed(1234)
#' targets_df <- data_df %>%
#'   pivot_longer(-pwgtp) %>%
#'   mutate(wtd_value = value * pwgtp) %>%
#'   group_by(name) %>%
#'   summarise(wtd_value = sum(wtd_value), .groups = "drop") %>%
#'   mutate(target = wtd_value * (1 + rnorm(length(.), mean=0, sd=.02)))
#' targets_df # in practice we'd make sure that targets make sense (e.g.,
#'   not negative)
#'
#' iweights <- data_df$pwgtp
#' targets <- targets_df$target
#' target_names <- targets_df$name
#'
#' tol <- .005 * abs(targets) # use 0.5% as our tolerance
#' xmat <- data_df %>%
#'   # important that columns be in the same order as the targets
#'   select(all_of(target_names)) %>%
#'   as.matrix
#'
#' res <- reweight(iweights = iweights, targets = targets,
#'                 target_names = target_names, tol = tol,
#'                 xmat = xmat)
#' names(res)
#' res$solver_message
#' res$etime
#' res$objective
#' res$targets_df
#'
#' @export
reweight <- function(iweights,
                     targets,
                     target_names,
                     tol,
                     xmat,
                     xlb=0,
                     xub=50,
                     maxiter = 50,
                     optlist = NULL,
                     method="auglag"){

  # stopifnot(all(!is.na(x$eval_jac_g(x$x0))))
  # stopifnot(is.function(x$eval_f))
  run_input_checks(method)
  if(!is.null(optlist)){ # check options list
    if(optlist$file_print_level > 0) {
      stopifnot("valid output_file name needed if file_print_level > 0" = !is.null(optlist$output_file))
    }
  }


  t1 <- proc.time()

  # define ipopt options
  if(!is.null(maxiter) & !(is.null(optlist))) {
    print("CAUTION: maxiter and opts both supplied. maxiter will override any iteration limit included in optlist.")
  }

  opts <- list(print_level = 0,
               file_print_level = 0, # integer
               max_iter= maxiter,
               linear_solver = "ma27", # mumps pardiso ma27 ma57 ma77 ma86 ma97
               jac_c_constant = "yes", # does not improve on moderate problems equality constraints
               jac_d_constant = "yes", # does not improve on  moderate problems inequality constraints
               hessian_constant = "yes") # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
               # output_file = NULL)

  # update default list just defined with any options specified in
  opts <- purrr::list_modify(opts, !!!optlist)
  if(!is.null(maxiter)) opts$max_iter <- maxiter # give priority to directly passed maxiter

  cc_sparse <- get_cc_sparse(xmat, target_names, iweights)

  inputs <- get_inputs(iweights,
                       targets, target_names, tol,
                       cc_sparse,
                       xlb, xub)

  result <- ipoptr::ipoptr(x0 = inputs$x0,
                   lb = inputs$xlb,
                   ub = inputs$xub,
                   eval_f = eval_f_xm1sq, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq
                   eval_grad_f = eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
                   eval_g = eval_g, # constraints LHS - a vector of values
                   eval_jac_g = eval_jac_g,
                   eval_jac_g_structure = inputs$eval_jac_g_structure,
                   eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
                   eval_h_structure = inputs$eval_h_structure,
                   constraint_lb = inputs$clb,
                   constraint_ub = inputs$cub,
                   opts = opts,
                   inputs = inputs)
  t2 <- proc.time()


  # define additional values to be returned
  solver_message <- result$message
  etime <- t2 - t1
  objective <- result$objective
  weights <- iweights * result$solution

  targets_df <- tibble::tibble(targnum = 1:length(targets),
                       targname = target_names,
                       target = targets) %>%
    dplyr::mutate(targinit = eval_g(inputs$x0, inputs),
           targcalc = result$constraints,
           targinit_diff = targinit - target,
           targtol = tol,
           targcalc_diff = targcalc - target,
           targinit_pdiff = targinit_diff / target * 100,
           targtol_pdiff = abs(targtol / target) * 100,
           targcalc_pdiff = targcalc_diff / target * 100)

  keepnames <- c("solver_message",
                 "etime",
                 "objective",
                 "weights",
                 "targets_df",
                 "result")
  output <- list()
  for(var in keepnames) output[[var]] <- get(var)

  output
}


run_input_checks <- function(method){
  stopifnot(method %in% c("auglag", "ipopt"))
  if(method == "auglag"){
    print("all good auglag")
  } else if(method == "ipopt") {
    print("all good ipopt")
  }
}



# Create sparse constraints coefficients data frame
# # @keyword internal
#' @export
get_cc_sparse <- function(xmat, target_names, iweights) {
  #.. dense matrix of constraint coefficients for the nation
  cc_sparse <- tidyr::as_tibble(xmat) %>%
    dplyr::select(dplyr::all_of(target_names)) %>%
    dplyr::mutate(iweight=iweights,
           j=dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(target_names),
                 names_to="cname",
                 values_to = "value") %>%
    dplyr::mutate(nzcc = iweight * value) %>%
    dplyr::filter(nzcc != 0) %>%
    dplyr::mutate(i=match(cname, target_names)) %>%
    dplyr::select(i, j, cname, nzcc, iweight, value) %>%
    dplyr::arrange(i, j) # this ordering is crucial for the Jacobian

  return(cc_sparse)
}

# Get inputs
# # @keyword internal
#' @export
get_inputs <- function(iweights,
                       targets,
                       target_names,
                       tol,
                       cc_sparse,
                       xlb, xub){

  inputs <- list()
  inputs$iweight <- iweights
  inputs$cc_sparse <- cc_sparse
  inputs$constraints <- targets
  inputs$constraint_names <- target_names
  inputs$n_variables <- length(inputs$iweight)
  inputs$n_constraints <- length(inputs$constraints)
  inputs$clb <- targets - tol
  inputs$cub <- targets + tol
  inputs$i_heq <- which(tol==0)
  inputs$i_hin <- which(tol!=0)

  # finally, add xlb, xub, x0, and the relevant structures
  inputs$xlb <- rep(xlb, inputs$n_variables)
  inputs$xub <- rep(xub, inputs$n_variables)
  inputs$x0 <- rep(1, inputs$n_variables)

  # for now:
  inputs$objscale <- 1 # divisor to scale objective function

  inputs$eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
  inputs$eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian

  inputs
}

