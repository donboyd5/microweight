#' Reweight a microdata file
#'
#' \code{reweight} calculates new weights for each household in a microdata
#' file that hit or come close to desired targets, while minimizing a measure of
#' distortion
#'
#' \code{reweight} uses the IPOPT solver in the package \code{\link[ipoptr]{ipoptr}}
#'
#' The problem is set up as a nonlinear program with constraints. The constraints are
#' the desired targets. The user can set tolerances around these targets in which case
#' they are inequality constraints. By default the distortion measure to be minimized
#' is the sum of squared differences between 1 and the ratio of new weights to the initial weight.
#' The user can provide alternative distortion measures.
#'
#'
#' @param iweights vector of initial household weights
#' @param xmat matrix of data for households, dimension h x k
#' @param targets named vector of desired target values, length k
#' @param maxiter integer; default 50
#' @param opts list of options that will update ipopt options
#' @param quiet c(TRUE, FALSE) FALSE is default; TRUE provides newlsqv or nls.lm output
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{nweights}{vector of new weights}
#'   \item{iweights}{vector of original weights}
#'   \item{result}{list of output from the solver that was used}
#' }
#' @examples
#' # Example 1: Determine state weights for a simple problem with random data
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#'
#' res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
#'   dweights = dw, quiet=TRUE)
#' @export
#'
#'

reweight <- function(iweights,
                     targets, target_names, tol,
                     xmat,
                     xlb=0, xub=50,
                     maxiter = 50,
                     optlist = NULL){
  # stopifnot(all(!is.na(x$eval_jac_g(x$x0))))
  # stopifnot(is.function(x$eval_f))
  if(optlist$file_print_level > 0) {
    stopifnot("valid output_file name needed if file_print_level > 0" = !is.null(optlist$output_file))
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

  result <- ipoptr(x0 = inputs$x0,
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
    mutate(targinit = eval_g(inputs$x0, inputs),
           targcalc = result$constraints,
           targinit_diff = targinit - target,
           targtol = tol,
           targcalc_diff = targcalc - target,
           targinit_pdiff = targinit_diff / target * 100,
           targtol_pdiff = abs(targtol / target) * 100,
           targcalc_pdiff = targcalc_diff / target * 100)

  keepnames <- c("solver_message", "etime", "objective", "weights",
                 "targets_df",
                 "result")
  output <- list()
  for(var in keepnames) output[[var]] <- get(var)

  output
}


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
    arrange(i, j) # this ordering is crucial for the Jacobian

  return(cc_sparse)
}

get_inputs <- function(iweights,
                       targets, target_names, tol,
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
