#' Calculate household-specific constants
#'
#' \code{get_delta} calculates household-specific constants for the poisson model
#'
#' See (Khitatrakun, Mermin, Francis, 2016, p.5)
#'
#' @param wh vector of household total weights, length h (see h, s, k definitions below)
#' @param beta s x k matrix of parameters to poisson model
#' @param xmat h x k matrix of data for households
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @returns A vector of household-specific constants of length h
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' beta <- matrix(0, nrow=p$s, ncol=p$k)
#' get_delta(p$wh, beta, p$xmat)
#' @export
get_delta <- function(wh, beta, xmat) {
    # calculate the household-specific constants (1 per household) given: wh: vector of constants for each household beta: s x k matrix of poisson
    # model coefficients that apply to all households xmat: h x k matrix of characteristics for each household

    # See (Khitatrakun, Mermin, Francis, 2016, p.5)

    # Note: we cannot let beta %*% xmat get too large!! or exp will be Inf and problem will bomb. It will get large when a beta element times an
    # xmat element is large, so either beta or xmat can be the problem.  not sure if this will need to be exported

    beta_x <- exp(beta %*% t(xmat))
    log(wh / colSums(beta_x))  # denominator is sum for each person
}


#' Calculate state weights for each household
#'
#' \code{get_weights} calculates state weights for each household
#'
#' See (Khitatrakun, Mermin, Francis, 2016, p.4)
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @param beta s x k matrix of parameters to poisson model
#' @param delta vector of household-specific constants, length h
#' @param xmat h x k matrix of data for households
#'
#' @returns h x s matrix of state weights for each household
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' beta <- matrix(0, nrow=p$s, ncol=p$k)
#' delta <- get_delta(p$wh, beta, p$xmat)
#' get_weights(beta, delta, p$xmat)
#' @export
get_weights <- function(beta, delta, xmat) {
    # calculate state-specific weights for each household that will add to their return an h x s matrix with weight for each state s for each
    # household h for each household the state weights will add to the total weight

    # beta: s x k matrix of poisson model coefficients that apply to all households delta: h-length vector of household-specific constants xmat: h
    # x k matrix of characteristics for each household

    # See (Khitatrakun, Mermin, Francis, 2016, p.4)

    # calc s x h matrix: each row has the sum over k of beta[s_i, k] * x[h_j, k] for each household where s_i is the state in row i each column is
    # a specific household
    beta_x <- beta %*% t(xmat)  #
    # add the delta vector of household constants to every row of beta_x and transpose
    beta_xd <- apply(beta_x, 1, function(mat) mat + delta)
    exp(beta_xd)  # exponentiate to calculate the weights
}


#' Calculate targets matrix
#'
#' \code{targets_mat} calculates s x k matrix of targets
#'
#' See (Khitatrakun, Mermin, Francis, 2016)
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#' @param betavec vector parameters to poisson model, length s * k (see h, s, k definitions below)
#' @param wh vector of household total weights, length h
#' @param xmat h x k matrix of data for households
#' @param s number of states (or other geographies or subgroups)
#'
#' @returns s x k matrix of targets
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' betavec <- rep(0, p$s * p$k)
#' targets_mat(betavec, p$wh, p$xmat, p$s)
#' @export
targets_mat <- function(betavec, wh, xmat, s) {
    # return an s x k matrix of calculated targets, given a beta vector, household weights, and x matrix
    beta <- vtom(betavec, s)
    delta <- get_delta(wh, beta, xmat)
    whs <- get_weights(beta, delta, xmat)
    t(whs) %*% xmat
}


#' Calculate targets vector
#'
#' \code{targets_vec} calculates vector of targets of length s * k
#'
#' See (Khitatrakun, Mermin, Francis, 2016)
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @param betavec vector parameters to poisson model, length s * k (see h, s, k definitions below)
#' @param wh vector of household total weights, length h
#' @param xmat h x k matrix of data for households
#' @param s number of states (or other geographies or subgroups)
#'
#' @returns vector of targets of length s * k
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' betavec <- rep(0, p$s * p$k)
#' targets_vec(betavec, p$wh, p$xmat, p$s)
#' @export
targets_vec <- function(betavec, wh, xmat, s) {
    # return a vector of calculated targets and corresponding values calculated given a beta vector, household weights, and x matrix
    calc_targets <- targets_mat(betavec, wh, xmat, s)
    as.vector(calc_targets)
}


#' Calculate vector of differences between desired targets and calculated targets
#'
#' \code{targets_mat} calculates vector of differences of length s * k
#'
#' See (Khitatrakun, Mermin, Francis, 2016)
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @param betavec vector parameters to poisson model, length s * k (see h, s, k definitions below)
#' @param wh vector of household total weights, length h
#' @param xmat h x k matrix of data for households
#' @param targets s x k matrix of desired target values
#' @param dweights optional length s * k vector of weights to multiply the vector of targets
#'
#' @returns vector of differences of length s * k
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#' betavec <- rep(0, p$s * p$k)
#' diff_vec(betavec, p$wh, p$xmat, p$targets, dw)
#' @export
diff_vec <- function(betavec, wh, xmat, targets, dweights = rep(1, length(targets))) {
    # return a vector of length s*k of differences between desired targets and calculated targets given a beta vector, household weights, and x
    # matrix of k characteristics
    calc_targets <- targets_mat(betavec, wh, xmat, nrow(targets))
    d <- targets - calc_targets  # these are both s x k matrices
    as.vector(d) * as.vector(dweights)
}


#' Sum of squared differences between desired targets and calculated targets
#'
#' \code{sse} calculates sum of squared differences
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @param betavec vector parameters to poisson model, length s * k (see h, s, k definitions below)
#' @param wh vector of household total weights, length h
#' @param xmat h x k matrix of data for households
#' @param targets s x k matrix of desired target values
#' @param dweights optional length s * k vector of weights to multiply the vector of targets
#'
#' @returns scalar: sum of squared differences
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#' betavec <- rep(0, p$s * p$k)
#' sse(betavec, p$wh, p$xmat, p$targets, dw)
#' @export
sse <- function(betavec, wh, xmat, targets, dweights = NULL) {
    # return a single value - sse (sum of squared errors)
    sum(diff_vec(betavec, wh, xmat, targets, dweights)^2)
}


#' Calculate weights to be multiplied by each target so that it hits a goal
#'
#' \code{get_dweights} calculates difference weights
#'
#' \describe{
#'   \item{h}{number of households (or individuals, records, tax returns, etc.)}
#'   \item{s}{number of states (or other geographies or subgroups)}
#'   \item{k}{number of characteristics each household has}
#' }
#'
#' @param targets s x k matrix of desired target values
#' @param goal optional scalar -- weights multiplied by targets will equal this
#'
#' @returns scalar: sum of squared differences
#' @examples
#' p <- make_problem(h=10, s=3, k=2)
#' dw <- get_dweights(p$targets)
#' dw * p$targets
#' @export
get_dweights <- function(targets, goal = 100) {
    # difference weights - a weight to be applied to each target in the difference function so that it hits its goal
    dw <- ifelse(targets!=0, goal / targets, 1)
    as.vector(dw)
}


#' Construct a starting point from a sample of the data
#'
#' \code{get_starting_point} returns a list with a starting point constructed
#'  from a sample and also returns the optimization results from that sample.
#'
#'
#' @param p a problem
#'
#' @returns list
#' @examples
#' p <- make_problem(h=100, s=3, k=2)
#' spoint <- get_starting_point(p)
#' @export
#' @importFrom rlang .data
get_starting_point <- function(p){
    # @importFrom  rlang .data needed in package so that we can use .data so that
    # wh does not trigger warning see below use of .data
    # see https://www.r-bloggers.com/no-visible-binding-for-global-variable/
    prows <- max(round(.1 * p$h), 5 * p$s * p$k) # n of rows we want
    if(is.null(colnames(p$xmat))) colnames(p$xmat) <- paste0('k', 1:ncol(p$xmat))

    if(prows >= p$h){print("Problem too small"); return(NULL)}

    print("Building small problem for starting point...")

    set.seed(1234)
    whx_df <-  tibble::as_tibble(p$xmat) %>%
        dplyr::mutate(wh=p$wh) %>%
        dplyr::sample_n(prows)

    p2 <- p
    p2$wh <- whx_df$wh
    p2$xmat <- whx_df %>%
        dplyr::select(-.data$wh) %>% # .data avoids R CMD check warning
        as.matrix()
    p2$targets <- p$targets * nrow(p2$xmat) / nrow(p$xmat)
    p2$h <- nrow(p2$xmat)

    print("Solving small problem...")
    res <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets,
                     method = 'Newton', maxiter = 50, quiet=TRUE)

    splist <- list()
    splist$spoint <- as.vector(res$beta_opt_mat)
    splist$result <- res

    return(splist)
}

