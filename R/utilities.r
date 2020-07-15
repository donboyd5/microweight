
# create problems ------------------------------------------------------------

#' Create Problem of Desired Size
#'
#' Will have households, states, and household characteristics.
#'
#' @param h number of households to create.
#' @param s number of states to create.
#' @param k number of characteristics each household will have.
#' @name make_problem
#' @returns A list with elements of a problem.
#' @examples
#' make_problem(4, 3, 2)
#' @export
make_problem <- function(h, s, k) {
    # create a problem of a chosen size h: # of households s: # of states k: # of characteristics per household

    # returns a list with items created below

    # example call: make_problem(8, 3, 2)

    set.seed(1234)
    xmat <- matrix(stats::runif(h * k), nrow = h, byrow = TRUE)

    set.seed(1234)
    whs <- matrix(stats::runif(h * s, 10, 20), nrow = h, byrow = TRUE)

    wh = rowSums(whs)
    ws = colSums(whs)

    targets <- t(whs) %*% xmat  # s x k

    keepnames <- c("h", "s", "k", "xmat", "wh", "ws", "whs", "targets")
    problem <- list()
    for (var in keepnames) problem[[var]] <- get(var)
    problem
}


# utilities for vectors and matrices -------------------------------------------

#' Convert Vector of States and Characteristics to Matrix With States as Rows
#'
#' Matrix will have same format as data matrix.
#'
#' @param v vector to be converted to matrix.
#' @param s number of states.
#' @name vtom
#' @returns A matrix with states as rows and characteristics as columns.
#' @examples
#' n_states <- 3
#' n_chars <- 2
#' set.seed(1)
#' value_vec <- runif(n_states * n_chars)
#' vtom(value_vec, n_states)
#' @export
vtom <- function(v, s) {
    # vector to matrix in the same ordering as a beta matrix
    matrix(v, nrow = s, byrow = FALSE)
}

