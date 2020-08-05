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
#' @param method "auglag" (default) or "ipopt" (requires installation of ipoptr)
#' @param quiet TRUE (default) or FALSE
#'
#' @returns A list with the following elements:
#'
#' \describe{
#' \item{solver_message}{The message produced by IPOPT. See
#' \href{https://coin-or.github.io/Ipopt/OUTPUT.html}{IPOPT output}.}
#' \item{etime}{Elapsed time.}
#'  \item{objective}{The objective function value at
#' the solution.}
#' \item{weights}{Numeric vector of new weights.}
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
#' res <- reweight(iweights = iweights,
#'                 targets = targets,
#'                 target_names = target_names,
#'                 tol = tol,
#'                 xmat = xmat,
#'                 method="ipopt")
#' res$etime
#' res$objective_unscaled
#' res$targets_df
#'
#' res2 <- reweight(iweights = iweights,
#'                  targets = targets,
#'                  target_names = target_names,
#'                  tol = tol,
#'                  xmat = xmat,
#'                  method="auglag")
#' res2$etime
#' res2$objective_unscaled
#' res2$targets_df
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
#' #       be created)
#' # we also need to get pwgtp - the person weight for each record, which
#' #       will be our initial weight
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
#' # in practice we'd make sure that targets make sense (e.g., not negative)
#' targets_df
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
#'                 xmat = xmat,
#'                 method = "ipopt")
#' names(res)
#' res$solver_message
#' res$etime
#' res$objective
#' res$targets_df
#'
#' res2 <- reweight(iweights = iweights, targets = targets,
#'                 target_names = target_names, tol = tol,
#'                 xmat = xmat,
#'                 method="auglag")
#' names(res2)
#' res2$solver_message
#' res2$etime
#' res2$objective
#' res2$targets_df
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
                     method = "auglag",
                     quiet = TRUE){

  # document these:
  # 'define_jac_g_structure_sparse' 'eval_f_xm1sq' 'eval_g'
  # 'eval_grad_f_xm1sq' 'eval_h_xm1sq' 'eval_jac_g' 'get_cc_sparse'
  # 'get_inputs'

  stopifnot(method %in% c("auglag", "ipopt"))

  check_args(method)

  t1 <- proc.time()
  cc_sparse <- get_cc_sparse(xmat, target_names, iweights)

  inputs <- get_inputs(iweights,
                       targets, target_names, tol,
                       cc_sparse,
                       xlb, xub)

  if(method == "auglag"){
    output <- call_auglag(iweights,
                          targets,
                          target_names,
                          tol,
                          xmat,
                          xlb,
                          xub,
                          maxiter,
                          optlist,
                          quiet,
                          inputs)
    } else if(method == "ipopt") {
      output <- call_ipopt(iweights,
                           targets,
                           target_names,
                           tol,
                           xmat,
                           xlb,
                           xub,
                           maxiter,
                           optlist,
                           quiet,
                           inputs)
  }
  t2 <- proc.time()

  # define additional values to be returned
  # solver_message <- result$message
  etime <- t2 - t1

  # use if statements below, even though not needed, in case result names change
  if(method == "auglag"){
    solver_message <- output$result$message
    objective_unscaled <- output$result$objective * output$inputs$objscale
    weights <- iweights * output$result$solution
    xfinal <- output$result$solution

  } else if(method == "ipopt"){
    solver_message <- output$result$message
    objective_unscaled <- output$result$objective * output$inputs$objscale
    weights <- iweights * output$result$solution
    xfinal <- output$result$solution
  }

  targets_df <- tibble::tibble(targnum = 1:length(targets),
                       targname = target_names,
                       target = targets) %>%
    dplyr::mutate(targinit = eval_g(inputs$x0, inputs),
           targcalc = eval_g(xfinal, inputs),
           targinit_diff = .data$targinit - .data$target,
           targtol = tol,
           targcalc_diff = .data$targcalc - .data$target,
           targinit_pdiff = .data$targinit_diff / .data$target * 100,
           targtol_pdiff = abs(.data$targtol / .data$target) * 100,
           targcalc_pdiff = .data$targcalc_diff / .data$target * 100)

  keepnames <- c("solver_message",
                 "etime",
                 "objective_unscaled",
                 "weights",
                 "xfinal",
                 "targets_df")
  # output <- list()
  for(var in keepnames) output[[var]] <- get(var)
  # reorder the list
  ordered_names <- c(keepnames, "result", "opts", "inputs")
  output <- output[ordered_names]

  output
}


check_args <- function(method){
  if(method == "auglag"){
    print("auglag check arguments")
  } else if(method == "ipopt"){
  #   if(!is.null(maxiter) & !(is.null(optlist))) {
  #     print("CAUTION: maxiter and opts both supplied. maxiter will override any iteration limit included in optlist.")
  #   }
  #   if(!is.null(optlist)){ # check options list
  #     if(optlist$file_print_level > 0) {
  #       stopifnot("valid output_file name needed if file_print_level > 0" = !is.null(optlist$output_file))
  #     }
  #   }
  # }
    print("ipopt check")
  }
}


call_ipopt <- function(iweights,
                       targets,
                       target_names,
                       tol,
                       xmat,
                       xlb,
                       xub,
                       maxiter,
                       optlist,
                       quiet,
                       inputs){
  # check inputs

  # define ipopt options
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

  output <- list()
  output$inputs <- inputs
  output$opts <- opts
  output$result <- result
  output
}


call_auglag <- function(iweights,
                        targets,
                        target_names,
                        tol,
                        xmat,
                        xlb,
                        xub,
                        maxiter,
                        optlist,
                        quiet,
                        inputs){
  # check inputs

  # add inputs needed for auglag
  inputs$objscale <- length(iweights)
  inputs$gscale <- targets / 10
  # inputs$objscale <- 1
  # inputs$gscale <- rep(1, length(targets))

  jac <- matrix(0, nrow=inputs$n_constraints, ncol=inputs$n_variables)
  f <- function(i){
    rep(i, length(inputs$eval_jac_g_structure[[i]]))
  }
  iidx <- lapply(1:length(inputs$eval_jac_g_structure), f) %>% unlist
  jidx <- unlist(inputs$eval_jac_g_structure)
  indexes <- cbind(iidx, jidx)
  jac[indexes] <- inputs$cc_sparse$nzcc

  jac_scaled <- sweep(jac, MARGIN=1, FUN="/", STATS=inputs$gscale)
  jac <- jac_scaled

  inputs$jac_heq <- jac[inputs$i_heq, ]
  inputs$jac_hin <- rbind(-jac[inputs$i_hin, ], jac[inputs$i_hin, ])

  local_opts <- list()
  local_opts$algorithm <- "NLOPT_LD_LBFGS"
  local_opts$xtol_rel <- 1.0e-4

  opts <- list()
  opts$algorithm <- "NLOPT_LD_AUGLAG"
  opts$ftol_rel <- 1.0e-4
  opts$maxeval <-  5000
  opts$print_level <- 0
  opts$local_opts <- local_opts

  result <- nloptr::nloptr(x0=inputs$x0,
                eval_f=eval_f_xm1sq,
                eval_grad_f = eval_grad_f_xm1sq,
                lb = inputs$xlb, ub = inputs$xub,
                eval_g_ineq = hin_fn,
                eval_jac_g_ineq = hin.jac_fn,
                eval_g_eq = heq_fn,
                eval_jac_g_eq = heq.jac_fn,
                opts = opts,
                inputs = inputs)

  output <- list()
  output$inputs <- inputs
  output$opts <- opts
  output$result <- result
  output
}


#' Create Sparse Constraints Coefficients Data Frame.
#'
#' Construct a data frame that represents a sparse matrix with nonzero
#' constraint coefficients. A constraint coefficient is the amount by which a
#' constraint changes if an element of the vector `x` changes by one unit. This
#' stores only those coefficients that are nonzero, in triplet form, where `i`
#' is an index identifying which constraint is represented, `j` is an index
#' identifying which element of the vector `x` is represented, and `nzcc` is the
#' value of the nonzero constraint coefficient.
#'
#' The data value in `xmat` multiplied by the initial weight determines how much
#' a constraint will change with a unit change in `x`, the ratio of the new
#' weight to the initial weight. For example, suppose household number 7 has
#' income of $10,000 then. Suppose that this household is in row 7 of `xmat` and
#' that column 3 of `xmat` corresponds to income. If this household has an
#' initial weight of 20, then the constraint coefficient for this cell - how
#' much changing the `x` value for this household will change the constraint for
#' total income - is 200,000 (20 x 100,000). The constraints coefficient data
#' frame will have a row where `i` (the index for the constraint, not the household)
#' is 3, `j` (the index for the household) is 7, and `nzcc` is 200000.
#'
#' @param xmat Numeric matrix of unweighted data values in dense
#'   form where each row is a household, each column corresponds to a
#'   constraint, and each cell is a data value.
#' @param target_names Character vector with names for each constraint. One element per
#' column in `xmat`.
#' @param iweights Numeric vector of initial weights in the microdata file, one per row in `xmat`.
#'
#' @returns A data frame representing a sparse matrix of constraint
#'   coefficients, with columns:
#'   \describe{
#'   \item{i}{the constraint number}
#'   \item{j}{an index into x}
#'   \item{cname}{constraint name, from `target_names`}
#'   \item{nzcc}{the nonzero constraint coefficient - the initial weight
#'   multiplied by value}
#'   \item{iweight}{initial weight for this cell}
#'   \item{value}{data value}
#'   }
#'
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
    dplyr::mutate(nzcc = .data$iweight * .data$value) %>%
    dplyr::filter(.data$nzcc != 0) %>%
    dplyr::mutate(i=match(.data$cname, target_names)) %>%
    dplyr::select(.data$i, .data$j, .data$cname, .data$nzcc,
                  .data$iweight, .data$value) %>%
    dplyr::arrange(.data$i, .data$j) # this ordering is crucial for the Jacobian

  return(cc_sparse)
}

#' Create List of Inputs for Optimization Functions.
#'
#' Most users will never use this function. Can be useful for debugging.
#'
#' @param iweights Numeric vector of initial weights in the microdata file.
#' @param targets Numeric vector of targets.
#' @param target_names Character vector of names for targets.
#' @param tol Numeric vector of tolerances, one per target.
#' @param cc_sparse Data frame with nonzero constraint coefficients in sparse format.
#' @param xlb Numeric vector of lower bounds for the `x` vector.
#' @param xub Numeric vector of upper bounds for the `x` vector.
#'
#' @export
get_inputs <- function(iweights,
                       targets,
                       target_names,
                       tol,
                       cc_sparse,
                       xlb,
                       xub){

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


# functions for nloptr
evg_scaled <- function(x, inputs){
  eval_g(x, inputs) / inputs$gscale
}


heq_fn <- function(x, inputs){
  evg_scaled(x, inputs)[inputs$i_heq] - inputs$targets[inputs$i_heq] / inputs$gscale[inputs$i_heq]
}

hin_fn <- function(x, inputs){
  # NLopt always expects constraints to be of the form myconstraint (x) ≤ 0,
  g <- evg_scaled(x, inputs)[inputs$i_hin]

  c(inputs$clb[inputs$i_hin]  / inputs$gscale[inputs$i_hin] - g,
    g - inputs$cub[inputs$i_hin]  / inputs$gscale[inputs$i_hin])
}


heq.jac_fn <- function(x, inputs){
  inputs$jac_heq
}


hin.jac_fn <- function(x, inputs){
  inputs$jac_hin
}


# Note that in the main package file, I have: @import dplyr so that we can access dplyr functions
# Throughout this file, the use of dplyr will generate a "no visible binding for global variable"
# note. However, it is not an error - it is ok.
# use .data$ before the variables that offend
# must  @importFrom rlang .data
