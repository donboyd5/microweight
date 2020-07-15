
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

library(microweight)

# get and save a small ACS sample ----
fn <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/data/acs_10krecs_5states.rds"
fn2 <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/data/acs_100krecs_20states.rds"

df <- readRDS(fn)
acs <- df %>%
  mutate(incgroup=ntile(pincp, 10)) %>%
  select(-pop, -nrecs)
usethis::use_data(acs)
quantile(df$pincp)

df <- readRDS(fn2)
acsbig <- df %>%
  mutate(incgroup=ntile(pincp, 10)) %>%
  select(-pop, -nrecs)
usethis::use_data(acsbig, overwrite = TRUE)
quantile(df$pincp)


# create an ACS problem ----
glimpse(acs)
count(acs, stabbr)
count(acs, incgroup)

acs_df <- acs %>%
  # filter(incgroup == 2) %>%
  mutate(pop = 1,
         young = agep < 21,
         workage = agep %in% 22:64,
         married = mar==1,
         female = sex==1,
         ssip_nnz = ssip != 0) %>%
  # convert to indicator variables
  mutate(across(c(young, workage, married, female, ssip_nnz), as.integer)) %>%
  select(serialno, incgroup, stabbr, pwgtp,
         pop, pincp, wagp, ssip, ssip_nnz, young, workage, married, female)
glimpse(acs_df)

acs_targets <- acs_df %>%
  group_by(incgroup, stabbr) %>%
  summarise(across(pop:female, ~ sum(. * pwgtp)), .groups = "drop")
usethis::use_data(acs_targets)

# create a big ACS problem ----
glimpse(acsbig)
count(acsbig, stabbr)
count(acsbig, incgroup)

acs_df2 <- acsbig %>%
  # filter(incgroup == 2) %>%
  mutate(pop = 1,
         young = agep < 21,
         workage = agep %in% 22:64,
         married = mar==1,
         female = sex==1,
         ssip_nnz = ssip != 0) %>%
  # convert to indicator variables
  mutate(across(c(young, workage, married, female, ssip_nnz), as.integer)) %>%
  select(serialno, incgroup, stabbr, pwgtp,
         pop, pincp, wagp, ssip, ssip_nnz, young, workage, married, female)
glimpse(acs_df2)

acsbig_targets <- acs_df2 %>%
  group_by(incgroup, stabbr) %>%
  summarise(across(pop:female, ~ sum(. * pwgtp)), .groups = "drop")
usethis::use_data(acsbig_targets)


#.. test out the problem ---
ig <- 2

acs_prep <- acsbig %>%
  dplyr::filter(incgroup == ig) %>%
  dplyr::mutate(pop = 1,
         young = agep < 21,
         workage = agep %in% 22:64,
         married = mar==1,
         female = sex==1,
         ssip_nnz = ssip != 0) %>%
  # convert to indicator variables
  dplyr::mutate(dplyr::across(c(young, workage, married, female, ssip_nnz), as.integer)) %>%
  dplyr::select(stabbr, pwgtp, pop, pincp, wagp, ssip, ssip_nnz, young, workage, married, female)


targs_df <- acsbig_targets %>%
  dplyr::filter(incgroup==ig)

# tdf2 <- targs_df %>%
#   dplyr::mutate(dplyr::across(c(pincp, wagp, ssip, ssip_nnz, young, workage, married, female),
#                               ~ifelse(.==0, pop * mean(. / pop), .)))

whs2 <- acs_prep %>%
  select(stabbr, pwgtp) %>%
  mutate(rn=row_number()) %>%
  arrange(stabbr) %>%
  pivot_wider(names_from = stabbr, values_from = pwgtp, values_fill=0) %>%
  arrange(rn) %>%
  select(-rn) %>%
  as.matrix

targets <- as.matrix(targs_df[ , -c(1:2)])
# targets <- as.matrix(tdf2[ , -c(1:2)])

rownames(targets) <- targs_df$stabbr

xmat <- acs_prep %>%
  dplyr::select(pop:female) %>%
  as.matrix()

wh <- acs_prep$pwgtp

identical(colnames(xmat), colnames(targets))

targets
targets_mat(rep(0, length(targets)), wh, xmat, nrow(targets))
# t(whs2) %*% xmat
# (t(whs2) %*% xmat) - targets

p <- list()
p$h <- nrow(xmat)
p$s <- nrow(targets)
p$k <- ncol(targets)
p$wh <- acs_prep$pwgtp
p$xmat <- xmat
p$targets <- targets

res1 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, maxiter=10)
res2 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "Newton", maxiter=10)
res3 <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = "LM", maxiter=10, opts=list(factor=15))

res2d <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, dweights=rep(1, length(p$targets)), method = "Newton")

p$targets %>% round(2)
res1$targets_calc %>% round(2)
res2$targets_calc %>% round(2)
res3$targets_calc %>% round(2)
res2d$targets_calc %>% round(2)


res3$targets_diff %>% round(2)
res3$targets_pctdiff %>% round(2)
