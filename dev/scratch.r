
# Determine new weights for a small problem using ACS data
library(microweight)
library(tidyverse)
data(acs)

data_df <- acs %>%
  filter(incgroup == 5) %>%
  select(stabbr, pwgtp, pincp, wagp, ssip) %>%
  # create the indicator variables
  mutate(nrecs = 1, # indicator used for number of records
         wagp_nnz = (wagp != 0) * 1.,
         ssip_nnz = (ssip != 0) * 1.)
data_df # 1,000 records

wh <- data_df$pwgtp

set.seed(1234)
targets_df <- data_df %>%
  pivot_longer(-c(pwgtp, stabbr)) %>%
  mutate(wtd_value = value * pwgtp) %>%
  group_by(stabbr, name) %>%
  summarise(wtd_value = sum(wtd_value), .groups = "drop") %>%
  mutate(wtd_value = wtd_value * (1 + rnorm(length(.), mean=0, sd=.01))) %>%
  pivot_wider(values_from = wtd_value)
targets_df

targets <- targets_df %>%
  select(-stabbr) %>%
  as.matrix
rownames(targets) <- targets_df$stabbr
targets

xmat <- data_df %>%
  select(all_of(colnames(targets))) %>%
  as.matrix
xmat

opts_LM <- list(ptol = 1e-16, ftol = 1e-16)
resx <- geoweight(wh = wh, xmat = xmat, targets = targets, optlist = opts_LM, quiet = FALSE)
names(resx)
resx$solver_message
resx$h; resx$s; resx$k
resx$method
resx$etime
resx$sse_unweighted
resx$sse_weighted
resx$targets
resx$targets_calc
resx$targets_diff %>% round(1)
resx$targets_pctdiff %>% round(1)






opts <- list(output_file = here::here("test.out"),
             file_print_level = 5,
             linear_solver = "ma57")

res <- reweight(iweights = iweights,
                targets = targets,
                target_names = target_names,
                tol = tol,
                xmat = xmat,
                method = "ipopt",
                maxiter = 10,
                optlist = opts)

res$solver_message
res$objective_unscaled
res$etime
names(res$result)


library(tidyverse)
library(microweight)
library(nloptr)
# library(alabama)

# can we use alabama

p <- make_problem(h=50000, s=1, k=50)
# we'll have 30 households (30 weights) and 4 targets

# break out needed components and prepare them for reweight
iweights <- p$wh

targets <- as.vector(p$targets)
target_names <- paste0("targ", 1:length(targets))
# the targets created by make_problem are hit exactly when using the initial
# weights so perturb them slightly so that we have to find new weights
set.seed(1234)
noise <- rnorm(n=length(targets), mean = 0, sd = .05)
targets <- targets * (1 + noise)

tol <- .005 * abs(targets)
tol[c(1, 3)] <- 0 # force some equality constraints
tol
# tol <- rep(0, length(targets))
xmat <- p$xmat
colnames(xmat) <- target_names

# eval_f_xm1sq
# eval_grad_f_xm1sq
# eval_h_xm1sq
# eval_g
# eval_jac_g

cc_sparse <- get_cc_sparse(xmat, target_names, iweights)

inputs <- get_inputs(iweights,
                     targets, target_names, tol,
                     cc_sparse,
                     xlb=0, xub=50)
inputs$objscale <- length(iweights)
inputs$gscale <- targets / 100

inputs$objscale <- 1
inputs$gscale <- rep(1, length(targets))

jac <- matrix(0, nrow=inputs$n_constraints, ncol=inputs$n_variables)
f <- function(i){
  rep(i, length(inputs$eval_jac_g_structure[[i]]))
}
iidx <- lapply(1:length(inputs$eval_jac_g_structure), f) %>% unlist
jidx <- unlist(inputs$eval_jac_g_structure)
indexes <- cbind(iidx, jidx)
jac[indexes] <- inputs$cc_sparse$nzcc
jac_bak <- jac
# jac

jac_scaled <- sweep(jac, MARGIN=1, FUN="/", STATS=inputs$gscale)
# jac[1:5, 1:8]
# inputs$gscale[1:5]
# jac_scaled[1:5, 1:8]
jac <- jac_scaled

jac_heq <- jac[inputs$i_heq, ]
# jac_hin <- rbind(jac[inputs$i_hin, ], -jac[inputs$i_hin, ])
jac_hin <- rbind(-jac[inputs$i_hin, ], jac[inputs$i_hin, ])


# scale evalg
evg_scaled <- function(x, inputs){
  eval_g(x, inputs) / inputs$gscale
}


heq_fn <- function(x, inputs){
  evg_scaled(x, inputs)[inputs$i_heq] - targets[inputs$i_heq] / inputs$gscale[inputs$i_heq]
}
heq_fn(inputs$x0, inputs)


hin_fn <- function(x, inputs){
  # NLopt always expects constraints to be of the form myconstraint (x) ≤ 0,
  g <- evg_scaled(x, inputs)[inputs$i_hin]

  c(inputs$clb[inputs$i_hin]  / inputs$gscale[inputs$i_hin] - g,
    g - inputs$cub[inputs$i_hin]  / inputs$gscale[inputs$i_hin])
}
hin_fn(inputs$x0, inputs)



heq.jac_fn <- function(x, inputs){
  jac_heq
}

heq.jac_fn <- function(x, inputs){
  jac_heq
}

hin.jac_fn <- function(x, inputs){
  jac_hin
}

inputs$clb[inputs$i_hin]; inputs$cub[inputs$i_hin]
eval_g(inputs$x0, inputs)[inputs$i_hin]

res <- reweight(iweights = iweights, targets = targets,
                target_names = target_names, tol = tol,
                xmat = xmat, method="ipopt")
res$solver_message
res$etime
res$objective

a <- proc.time()
tmp <- nloptr(x0=tmp2$solution, # x0=inputs$x0, # x0=tmp$solution,
               eval_f=eval_f_xm1sq,
               eval_grad_f = eval_grad_f_xm1sq,
               lb = inputs$xlb, ub = inputs$xub,
               eval_g_ineq = hin_fn, # hin_fn_nlopt
               eval_jac_g_ineq = hin.jac_fn, # hin.jac_fn_nlopt
               eval_g_eq = heq_fn,
               eval_jac_g_eq = heq.jac_fn,
               opts = list(algorithm="NLOPT_LD_AUGLAG", # NLOPT_LD_AUGLAG_EQ NLOPT_LD_AUGLAG
                           # xtol_rel=1.0e-6,
                           ftol_rel=1.0e-4,
                           maxeval = 5000,
                           print_level = 1,
                           # nlopt_set_vector_storage=1e6,
                           local_opts = list(algorithm="NLOPT_LD_LBFGS", # NLOPT_LD_LBFGS NLOPT_LD_LBFGS NLOPT_LD_SLSQP (SLOW)
                                             xtol_rel  = 1.0e-4)),
               inputs = inputs)
b <- proc.time()
b - a

# nloptr.print.options()
# https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-gradient-based-optimization
#  ?'nloptr-package'
tmpbak <- tmp
tmp2 <- tmp

# names(tmp2)
tmp2$message
tmp2$objective
tmp2$iterations
tmp2$num_constraints_eq
tmp2$num_constraints_ineq
# outer: NLOPT_LD_AUGLAG, NLOPT_LD_AUGLAG_EQ, NLOPT_LD_SLSQP
# NLOPT_LD_SLSQP does not require local optimizer; TOO LONG


eval_f_xm1sq(res$result$solution, inputs)
# eval_f_xm1sq(tmp$par, inputs)
eval_f_xm1sq(tmp2$solution, inputs)

# targets
inputs$clb; inputs$cub
eval_g(res$result$solution, inputs)
# eval_g(tmp$par, inputs)
eval_g(tmp2$solution, inputs)

(eval_g(tmp2$solution, inputs) / eval_g(res$result$solution, inputs) * 100 - 100) %>% round(2)
(eval_g(tmp2$solution, inputs) / targets * 100 - 100) %>% round(5)

# cor(tmp$par, res$result$solution)
cor(tmp2$solution, res$result$solution)

# (tmp$par - res$result$solution) %>% round(2)
(tmp2$solution - res$result$solution) %>% round(2)


# d <- check.derivatives(inputs$x0,
#                        func=eval_f_xm1sq,
#                        func_grad=eval_grad_f_xm1sq,
#                        check_derivatives_tol = 1e-04,
#                        check_derivatives_print = "errors",
#                        func_grad_name = "grad_f",
#                        inputs=inputs)

# library(numDeriv)
# j <- jacobian(func=evg_scaled, x=inputs$x0, inputs=inputs)
# dim(j)
# ivals <- 1:5
# jvals <- 1:8
#
# j[ivals, jvals]
# jac_scaled[ivals, jvals]
# sum((j - jac_scaled)^2)

# j <- jacobian(func=hin_fn, x=inputs$x0, inputs=inputs)
# dim(j)
# ivals <- 1:6
# jvals <- 1:8
#
# j[ivals, jvals]
# jac_hin[ivals, jvals]
# sum((j - jac_hin)^2)



# a <- proc.time()
# tmp <- alabama::auglag(par=inputs$x0,
#                        fn=eval_f_xm1sq,
#                        gr = eval_grad_f_xm1sq,
#                        heq = heq_fn,
#                        heq.jac = heq.jac_fn,
#                        hin = hin_fn,
#                        hin.jac = hin.jac_fn,
#                        control.outer = list(itmax = 100,
#                                             trace=TRUE,
#                                             kkt2.check=FALSE,
#                                             # mu0=.9,
#                                             # sig0=.9, # .9 worked well
#                                             # eps = 1e-8,
#                                             method="nlminb"), # "nlminb L-BFGS-B
#                        control.optim = list(trace=0, maxit=20), # , factr=1e-10, pgtol=1e-10 # inner loop, optim options
#                        inputs=inputs)
# b <- proc.time()
# b - a







data(acsbig)
data(acs_targets)

# let's focus on income group 5 and create and then try to hit targets for:
#    number of records (nrecs -- to be created based on the weight, pwgtp)
#    personal income (pincp)
#    wages (wagp)
#    number of people with wages (wagp_nnz -- to be created)
#    supplemental security income (ssip)
#    number of people with supplemental security income (ssip_nnz -- to be created)
# we also need to get pwgtp - the person weight for each record, which will be
#   our initial weight
# for each "number of" variable we need to create an indicator variable that
# defines whether it is true for that record

# get the data and prepare it
data_df <- acsbig %>%
  filter(incgroup == 5) %>%
  select(pwgtp, pincp, wagp, ssip) %>%
  # create the indicator variables
  mutate(nrecs = 1, # indicator used for number of records
         wagp_nnz = (wagp != 0) * 1.,
         ssip_nnz = (ssip != 0) * 1.)
data_df # 10,000 records

# prepare targets: in practice we would get them from an external source but
# in this case we'll get actual sums on the file and perturb them randomly
# so that targets differ from initial sums
set.seed(1234)
targets_df <- data_df %>%
  pivot_longer(-pwgtp) %>%
  mutate(wtd_value = value * pwgtp) %>%
  group_by(name) %>%
  summarise(wtd_value = sum(wtd_value), .groups = "drop") %>%
  mutate(target = wtd_value * (1 + rnorm(length(.), mean=0, sd=.02)))
targets_df # in practice we'd make sure that targets make sense (e.g., not negative)


iweights <- data_df$pwgtp
targets <- targets_df$target
target_names <- targets_df$name

tol <- .005 * abs(targets) # use 0.5% as our tolerance
xmat <- data_df %>%
  # important that they be in the same order as the targets
  select(all_of(target_names)) %>%
  as.matrix

res <- reweight(iweights = iweights, targets = targets,
                target_names = target_names, tol = tol,
                xmat = xmat, xlb = 0.1)
names(res)
res$solver_message
res$etime
res$objective
res$targets_df

sum(iweights)
sum(res$weights)

sum(iweights==0)
sum(res$weights==0)

quantile(iweights)
quantile(res$weights)
quantile(res$result$solution)

sort(res$result$solution, decreasing = TRUE)
sum(res$result$solution)


glimpse(acs)
glimpse(acs_targets)

glimpse(acsbig)



library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr

options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(microweight)

library(optimx)

library(mize)

p <- make_problem(h=10, s=3, k=2)
p <- make_problem(h=100, s=10, k=5)
p <- make_problem(h=5000, s=20, k=8)


dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))
x <- rep(1, length(p$targets))

library(numDeriv)

system.time(hmat <- hessian(sse, betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
hmat
hmat2

(hmat - hmat2) %>% round(3)



# gradient of the sse function
system.time(a <- numDeriv::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
system.time(b <- pracma::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)) # much faster
a
b
sum(a)
sum(b)
sum((a-b)^2)

library(mize)
ff <- function(x, wh, xmat, targets, dweights=NULL){
  sse(betavec=x, wh=wh, xmat=xmat, targets=targets, dweights = dweights)
}

gf <- function(x, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=x, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}

ff(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)
gf(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)

fgm <- list(fn = ff,
           gr = gf)
mize(par=betavec, fg=fgm, method = "L-BFGS", wh = p$wh, xmat = p$xmat, targets = p$targets,
     dweights = dw)

res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw)

res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw, method = 'Newton')

res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
  dweights = dw, method = 'LM')


# try something new
library(nloptr)
# allmeth <- c("BFGS", "CG", "Nelder-Mead", "L-BFGS-B", "nlm", "nlminb",
#              "lbfgsb3", "Rcgmin", "Rtnmin", "Rvmmin", "spg", "ucminf",
#              "newuoa", "bobyqa", "nmkb", "hjkb", "hjn", "lbfgs", "subplex")
library(lbfgsb3c)
library(subplex)
library(ucminf)
# newuoa slow; spg slow; bobyqa slow; nmkb no progress; ucminf slow; nlm slow; Rcgmin [needs gradient]

betavec <- rep(0, length(p$targets))
a <- proc.time()
res4 <- optimx::optimr(par = betavec,
                       fn = sse,
                       method = "L-BFGS-B", # L-BFGS-B
                       control = list(trace=1,
                                      maxit=2000,
                                      kkt=FALSE,
                                      starttests=FALSE),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a
str(res4)


gf <- function(betavec, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=betavec, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}
gf(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)

a <- proc.time()
res5 <- optimx::optimr(par = betavec,
                       fn = sse,
                       method = "L-BFGS-B", # L-BFGS-B
                       gr=gf, # "grcentral", #gf,
                       control = list(trace=1,
                                      maxit=2000,
                                      kkt=FALSE,
                                      starttests=FALSE),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a
str(res5)

betavec <- rep(0, length(p$targets))
a <- proc.time()
res5 <- subplex(par = betavec,
                       fn = sse,
                       control = list(maxit=200),
                       wh = p$wh, xmat = p$xmat, targets = p$targets,
                       dweights = dw)
b <- proc.time()
b - a

res4$par
as.vector(res3$beta_opt_mat)




library(dfoptim)
p
system.time(resrv <- hjk(rep(0, length(p$targets)), sse, control = list(info=TRUE), wh=p$wh, xmat=p$xmat, targets=p$targets))
system.time(resrv2 <- nmk(rep(0, length(p$targets)), sse, control = list(), wh=p$wh, xmat=p$xmat, targets=p$targets))
resrv$par
resrv2$par

library(rootSolve)
A <- matrix(nrow = 500, ncol = 500, runif(500*500))
B <- runif(500)
jfun <- function (x) A
fun <- function(x) A %*%x - B
fun <- function(x) (A %*%x - B) %>% as.vector
fun2 <- function(x, z) (A %*%x - B) %>% as.vector
system.time(X <- multiroot(start = 1:500, f = fun,
                           jactype = "fullusr", jacfunc = jfun))
system.time(X <- multiroot(start = 1:500, f = fun))

system.time(X <- multiroot(start = 1:500, f = fun, parms=list(z=10)))


assert_that(is.character(x))
assert_that(is.matrix(A))

f <- function(x, z){
  c(x[1]+10, x[2]-5)^2 * z
}
f(c(2, 3), 10)
multiroot(start = c(2, 3), f, z=10)

system.time(X <- multiroot(start = 1:500, f = fun))

sum(fun(X$root)^2)
fun(1:500)

p2 <- make_problem(h=4000, s=30, k=15)

system.time(v1 <- multiroot(start=rep(0, length(p2$targets)), f=diff_vec,
                            atol=1e-4,
                            wh=p2$wh, xmat=p2$xmat, targets=p2$targets))
max(abs(v1$f.root))

v2 <- geoweight(wh = p2$wh, xmat = p2$xmat, targets = p2$targets,
          method = 'Newton')
max(abs(v2$output$fvec))
v2$etime

p2$targets
targets_mat(as.vector(v2$beta_opt_mat), p2$wh, p2$xmat, p2$s)
targets_mat(as.vector(v1$root), p2$wh, p2$xmat, p2$s)



v2 <- geoweight(wh = p2$wh, xmat = p2$xmat, targets = p2$targets,
                dweights = rep(1, length(p2$targets)), method = 'Newton')

v3 <- multiroot(diff_vec, rep(0, length(p$targets)), maxiter = 100,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))

v4 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                dweights = rep(1, length(p$targets)), method = 'Newton')

v4 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets, method = 'Newton')


p2$targets
targets_mat(as.vector(v2$beta_opt_mat), p2$wh, p2$xmat, p2$s)
targets_mat(as.vector(v1$root), p2$wh, p2$xmat, p2$s)


multiroot(diff_vec, rep(0, length(p$targets)), maxiter = 100,
          rtol = 1e-6, atol = 1e-8, ctol = 1e-8,
          useFortran = TRUE, positive = FALSE,
          jacfunc = NULL, jactype = "fullint",
          verbose = FALSE, bandup = 1, banddown = 1,
          parms = NULL, wh=p$wh, xmat=p$xmat, targets=p$targets)

rosenbrock <- function(x){
  n <- length(x)
  sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
}
par0 <- rep(0, 10)
hjk(par0, rosenbrock)

hjkb(c(0, 0, 0), rosenbrock, upper = 0.5)

library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

# library(readxl)
# library(grDevices)
# library(knitr)
# library(kableExtra)
# library(janitor)
# library(ipoptr)


glimpse(acs)
count(acs, stabbr)
count(acs, incgroup)

acsdf <- acs
acsdf <- acsbig
glimpse(acsdf)

acs_prep <- acsdf %>%
  filter(incgroup == 2) %>%
  mutate(pop = 1,
         young = agep < 21,
         workage = agep %in% 22:64,
         married = mar==1,
         female = sex==1,
         ssip_nnz = ssip != 0) %>%
  # convert to indicator variables
  mutate(across(c(young, workage, married, female, ssip_nnz), as.integer)) %>%
  select(serialno, stabbr, pwgtp,
         pop, pincp, wagp, ssip, ssip_nnz, young, workage, married, female)
glimpse(acs_prep)

xmat <- acs_prep %>% select(-serialno, -stabbr, -pwgtp) %>%
  as.matrix()

wh <- acs_prep$pwgtp

targs_df <- acs_prep %>%
  group_by(stabbr) %>%
  summarise(across(pop:female, ~ sum(. * pwgtp)), .groups = "drop")
targets <- as.matrix(targs_df[ , -1])
rownames(targets) <- targs_df$stabbr

identical(colnames(xmat), colnames(targets))

colnames(xmat) <- NULL
colnames(targets) <- NULL

p <- list()
p$h <- nrow(xmat)
p$s <- nrow(targets)
p$k <- ncol(targets)
p$wh <- acs_prep$pwgtp
p$xmat <- xmat

p$h; p$s; p$k

set.seed(1234)
noise <- rnorm(n=length(targets), mean=0, sd=.01)
# noise <- 0
p$targets <- targets * (1 + noise)

res1 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Broyden')
res2 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton')
res3 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM')

rbind(res1$solver_message, res2$solver_message, res3$solver_message)
cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)

res1$targets_diff %>% round(1)
res2$targets_diff %>% round(1)
res3$targets_diff %>% round(1)

res1$targets_pctdiff %>% round(1)
res2$targets_pctdiff %>% round(1)
res3$targets_pctdiff %>% round(1)
c(max(abs(res1$targets_pctdiff)),
  max(abs(res2$targets_pctdiff)),
  max(abs(res3$targets_pctdiff))) %>% round(3)

cbind(res1$etime, res2$etime, res3$etime)

quantile(res3$beta_opt_mat)


spoint <- get_starting_point(p)
spoint$result$solver_message
spoint$result$etime

res1s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'Broyden')
res1s$sse_weighted
res1s$sse_unweighted
res1s$etime
res1$beta_opt_mat
res1s$beta_opt_mat

res2s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'Newton')
res2s$sse_weighted
res2s$sse_unweighted
res2s$etime
res2$beta_opt_mat
res2s$beta_opt_mat

res3s <- geoweight(wh=p2$wh, xmat=p2$xmat, targets=p2$targets, method = 'LM')
res3s$sse_weighted
res3s$sse_unweighted
res3s$etime
res3$beta_opt_mat
res3s$beta_opt_mat

res2_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton', betavec=spoint$spoint)
res2_2step$etime

res3_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM', betavec=spoint$spoint)
res3_2step$etime

rbind(res1$solver_message, res2$solver_message, res3$solver_message,
      res2_2step$solver_message, res3_2step$solver_message)

cbind(res1$etime, res2$etime, res3$etime, res2_2step$etime, res3_2step$etime)

cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted,
      res2_2step$sse_unweighted, res3_2step$sse_unweighted)

cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted,
      res2_2step$sse_weighted, res3_2step$sse_weighted)


res1$targets_diff %>% round(1)
res2$targets_diff %>% round(1)
res3$targets_diff %>% round(1)
res2_2step$targets_diff %>% round(1)
res3_2step$targets_diff %>% round(1)

res1$targets_pctdiff %>% round(1)
res2$targets_pctdiff %>% round(1)
res3$targets_pctdiff %>% round(1)
res3_2step$targets_pctdiff %>% round(1)

c(max(abs(res1$targets_pctdiff)),
  max(abs(res2$targets_pctdiff)),
  max(abs(res3$targets_pctdiff)),
  max(abs(res2_2step$targets_pctdiff)),
  max(abs(res3_2step$targets_pctdiff))
  ) %>% round(3)


res3$beta_opt_mat %>% round(4)
res2_2step$beta_opt_mat %>% round(4)
res3_2step$beta_opt_mat %>% round(4)

summary(res3$whs)
summary(res2_2step$whs)

# now do tax data ----



# could we use ipopt?? how fast is gradient?? ----
dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))
# gradient of the sse function
# system.time(a <- numDeriv::grad(sse, x, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw))
system.time(b <- pracma::grad(sse, betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)) # much faster

gf <- function(x, wh, xmat, targets, dweights=NULL){
  pracma::grad(sse, x0=x, wh=wh, xmat=xmat, targets=targets, dweights=dweights)
}

inputs <- list()
inputs$wh <- p$wh
inputs$xmat <- p$xmat
inputs$targets <- p$targets
inputs$dweights <- dw

library(ipoptr)

opts_ipopt <- list("print_level" = 0,
                   "file_print_level" = 5, # integer
                   "print_user_options" = "yes",
                   "max_iter"= 50,
                   "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
                   # "mehrotra_algorithm" = "yes", # default no
                   "hessian_approximation" = "exact", # default exact; limited-memory
                   # "limited_memory_update_type" = "sr1", # default bfgs; sr1 docs say "not working well" but this really helps
                   "obj_scaling_factor" = 3, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
                   # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
                   # "accept_every_trial_step" = "yes",
                   # "derivative_test" = "first-order",
                   "output_file" = here::here("temp3.out"))

set.seed(1)
xran <- runif(length(betavec), -10, 10)

# maybe fast hessian?

ip1 <- ipoptr(x0 = betavec,
              eval_f = sse, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq eval_f_wfs
              eval_grad_f = gf, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
              lb=rep(-100, length(betavec)),
              ub=rep(100, length(betavec)),
              # eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
              # eval_h_structure = inputs$eval_h_structure,
              opts = opts_ipopt,
              wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw)


x2 <- ip1$solution
# res3$beta_opt_mat

iptarg <-  targets_mat(x2, p$wh, p$xmat, p$s)
iptarg


fn <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/data/acs_10krecs_5states.rds"
df <- readRDS(fn)
acs <- df %>%
  mutate(incgroup=ntile(pincp, 10)) %>%
  select(-pop, -nrecs)
usethis::use_data(acs)
quantile(df$pincp)


library(minpack.lm) # nls.lm
library(nleqslv) # nleqslv

p <- make_problem(h=10, s=50, k=10)

p$targets
dw <- get_dweights(p$targets)
betavec <- rep(0, length(p$targets))

v1 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="Newton")

v2 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="Broyden")

v3 <- geoweight(betavec, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=dw,
        method="levmarq")

v1$etime
v2$etime
v3$etime

(diff_vec(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)
(diff_vec(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)
(diff_vec(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) / as.vector(p$targets) * 100) %>% round(3)

diff_vec(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)
diff_vec(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)
diff_vec(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets) %>% round(3)

sse(v1$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)
sse(v2$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)
sse(v3$betavec, wh=p$wh, xmat=p$xmat, targets=p$targets)

dt <- t4 - t3
dt

cbind(v1$x, v2$x, v3$x)

first <- list(a="1", b="2")
second <- list(b="3", d="4")

Map(c, first, second)


p <- make_problem(h=1000, s=20, k=10)
dw <- get_dweights(p$targets)

res1 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                  dweights = get_dweights(p$targets), maxiter=20)

res2 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
                  dweights = get_dweights(p$targets), method = 'Newton', maxiter=10)

res3 <- geoweight(wh = p$wh, xmat = p$xmat, targets = p$targets,
   dweights = get_dweights(p$targets), method = 'LM', maxiter=20, opts=list(factor=10))


res1
res2
res3
c(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)
c(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)

if (!assertthat::is.number(n)) { # make sure n is just one value and numeric
  stop("n should be a number of length 1!")
}
assert_that(is.character(x))

start <- rep(0, length(p$targets))
a <- proc.time()
res4 <- multiroot(start=start,
                  f=diff_vec, verbose=TRUE,
                  wh=p$wh, xmat=p$xmat, targets=p$targets)
b <- proc.time()
b - a
res4

library(BB)
start <- rep(0, length(p$targets))
t1 <- proc.time()
res5 <- BBoptim(par=start, fn=sse,
                method=3,
                wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))
t2 <- proc.time()
t2 - t1

res6 <- spg(start, fn=sse, control=list(), quiet=FALSE, alertConvergence=TRUE,
            wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))

res1$beta_opt_mat

res1$targets_calc
res2$targets_calc
res3$targets_calc
targets_mat(res5$par, p$wh, p$xmat, p$s)
targets_mat(res6$par, p$wh, p$xmat, p$s)
targets_mat(as.vector(res3$beta_opt_mat), p$wh, p$xmat, p$s)
p$targets

res1$etime; res2$etime; res3$etime; (b - a)
ssew <- sse(betavec=res4$root, wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=get_dweights(p$targets))


res4 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=rep(1, length(targets)))
res5 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "Newton", dweights=rep(1, length(targets)))
res6 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "LM", dweights=rep(1, length(targets)))

cbind(res1$sse_weighted, res2$sse_weighted, res3$sse_weighted)
cbind(res1$sse_unweighted, res2$sse_unweighted, res3$sse_unweighted)

res1$etime; res2$etime; res3$etime



# res2d <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=rep(1, length(targets)), method = "Newton")

out <- res1$targets_calc
out <- res2$targets_calc
# out <- res2d$targets_calc
out <- res3$targets_calc

out <- res4$targets_calc
out <- res5$targets_calc
out <- res6$targets_calc


out
p$targets

(out - p$targets) %>% round(4)

(out / p$targets * 100 - 100) %>% round(4)


res1$etime
res2$etime
res3$etime

res2$whs %>% round(2) %>% head(10)
whs2 %>% round(2) %>% head(10)

st <- "CA"
st <- "FL"
cor(whs2[, st], res2$whs[, st])


library(purrr)
x <- list(x = 1:10, y = 4, z = list(a = 1, b = 2))
str(x)
# Update values
str(list_modify(x, a = 1))
x
str(list_modify(x, z = 5))

a <- list(x = 1:10, y = 4)
str(a)
b <- list(y=7, m=1:3)
str(b)
str(list_modify(a, !!!b))

str(list_merge(a, b))
str(update_list(a, b))


df <- tibble(
  x = 1:3,
  y = c("a", "d,e,f", "g,h")
)
df %>%
  transform(y = strsplit(y, ",")) %>%
  unnest(y)


f <- function(a, b, c){
  args <- as.list(match.call())[-1]
  print(str(args))
  print(names(args))
  print(args)
}

f(6, "cat", 1.3)

