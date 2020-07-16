
# TODO ----
# step 7 create new state weights
# step 8 test results


# Problem targets ----
# N09600_nnz and especially A09600_sum bad for most states in AGI_STUB 1
# also bad in other stubs, especially small states


# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r


# libraries ----
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(purrr)
# library(repurrrsive) ??
library(multidplyr)

library(knitr)
library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")
library(ipoptr)

library(microweight) # devtools::install_github("donboyd5/microweight")

devtools::session_info()

# functions ----
source(here::here("dev", "functions_puf_tpc.r"))
source(here::here("dev", "functions_puf_ipoptr.r"))
source(here::here("dev", "functions_ipopt.r"))
source(here::here("dev", "functions_targets.r"))


# globals ----
dbdir <- "C:/Users/donbo/Dropbox (Personal)/"

# useful cut points for quantiles
probs1_99 <- c(0, .01, .05, .10, .25, .5, .75, .9, .95, .99, 1)

# here is Peter's list of desired targets
desired_target_vars <- c("N1_nnz", "MARS1_nnz", "MARS2_nnz",
                         "A00100_sum", "posAGI_nnz",
                         "A00200_sum", "N00200_nnz",
                         "A01000_sum", "N01000_nnz",
                         "A04470_sum", "N04470_nnz",
                         "A17000_sum", "N17000_nnz",
                         "A04800_sum", "N04800_nnz",
                         "A05800_sum", "N05800_nnz",
                         "A09600_sum", "N09600_nnz",
                         "A00700_sum", "N00700_nnz")


# 1. get targets, examine ----
puf_targ <- read_csv(file = paste0(dbdir, 'puf_2017_targets.csv')) %>%
  select(-X1)
ht(puf_targ)
(available_target_vars <- setdiff(names(puf_targ), c("STATE", "AGI_STUB")))

setdiff(desired_target_vars, available_target_vars)
setdiff(available_target_vars, desired_target_vars)

#.. 1-a. check the targets ----
# get US totals for each income range and in total, for each target
count(puf_targ, STATE)
ustargtot <- puf_targ %>%
  select(AGI_STUB, all_of(available_target_vars)) %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(available_target_vars), function(x) sum(x) / 1e3), .groups = "drop") %>%
  janitor::adorn_totals()

# useful for comparison vs puf values
ustargtot[, 1:7] %>%
  kable(digits=0, format.args = list(big.mark = ','))

# get each state as % of all states
shares <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(available_target_vars)) %>%
  group_by(AGI_STUB) %>%
  mutate(across(all_of(available_target_vars), ~ . / sum(.) * 100))
shares %>% filter(AGI_STUB==2)

# which have a state share of total that is far from their pop share?
shareslong <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(available_target_vars)) %>%
  pivot_longer(-c(AGI_STUB, STATE), names_to="targname", values_to="target") %>%
  group_by(AGI_STUB, targname) %>%
  mutate(target_pct=target / sum(target, na.rm=TRUE) * 100) %>%
  group_by(AGI_STUB, STATE) %>%
  mutate(pop_pct=target_pct[targname=="N1_nnz"],
         diff=target_pct - pop_pct,
         pdiff=diff / pop_pct * 100,
         apdiff=abs(pdiff)) %>%
  ungroup

shareslong %>%
  arrange(AGI_STUB, desc(apdiff)) %>%
  group_by(AGI_STUB) %>%
  filter(row_number() < 5)

# bad stub 10: NH, MS, RI, AR, IA WY
shareslong %>%
  filter(AGI_STUB==10, STATE=="WY")

# SD 2 (9600) WV 3 (9600)
shareslong %>%
  filter(AGI_STUB==2, STATE=="MI")


# worst 20 -- all are in incgroup 1
shareslong %>%
  filter(AGI_STUB == 1) %>%
  filter(STATE == "NY") %>%
  arrange(AGI_STUB, desc(apdiff)) %>%
  filter(row_number() < 20)

# check n returns
mars <- puf_targ %>%
  select(AGI_STUB, STATE, N1_nnz, MARS1_nnz, MARS2_nnz) %>%
  mutate(resid=N1_nnz - MARS1_nnz - MARS2_nnz,
         pct=resid / N1_nnz * 100)

mars %>% arrange(pct) # DC group 1 is negative

quantile(mars$pct)

mars %>% arrange(-pct)


# 2. get puf ----
puf_2017 <- read_csv(file = paste0(dbdir, 'puf_2017_filter.csv'))
glimpse(puf_2017)
ht(puf_2017)

# put an id on the file that allows us to identify original records, then sort
# and put a sortid on the file
puf_2017b <- puf_2017 %>%
  mutate(recid=row_number()) %>%
  arrange(AGI_STUB, recid) %>%
  mutate(sortid=row_number())
glimpse(puf_2017b)


# 3. prepare puf ----
available_target_vars

# get vector of variables for each kind of target, and remove suffix from name
nnz_vars <-remove_suffix("_nnz", available_target_vars) # variables for which we want the weighted number of nonzero values
sum_vars <- remove_suffix("_sum", available_target_vars) # variables for which we want the weighted sum
sumneg_vars <- remove_suffix("_sumneg", available_target_vars) # variables for which we want the weighted sum

pufprep <- puf_2017b %>%
  mutate_at(nnz_vars,
            list(nnz = ~ 1 * (. != 0))) %>%
  mutate_at(sum_vars,
            list(sum = ~ . * (. != 0))) %>%
  mutate_at(sumneg_vars,
            list(sumneg = ~ . * (. < 0)))
pufprep

# prepare a stripped-down version - an important file
pufstrip <- pufprep %>%
  select(sortid, AGI_STUB, s006, contains("_nnz"), contains("_sum"), contains("_sumneg"))
glimpse(pufstrip)

(puf_target_vars <- setdiff(names(pufstrip), c("sortid", "AGI_STUB", "s006"))) # 21 target vars

# get puf totals and compare them to selected target totals, created earlier
puftots <- pufstrip %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(puf_target_vars), function(x) sum(x * s006) / 1e3), .groups="drop") %>%
  janitor::adorn_totals()

names(ustargtot)
cols <- 1:8
cols <- c(1, 9:14)
cols <- c(1, 15:21)
(vars <- names(ustargtot)[cols])
puftots[, vars] %>%
  kable(digits=0, format.args = list(big.mark = ','))

ustargtot[, vars] %>%
  kable(digits=0, format.args = list(big.mark = ','))
# good, they are the close or identical
# note that there are some NA values in the targets and the corresponding
# puf values are not always zero -- esp. A04470_sum and A17000_sum


# 4. define targets ----
# pick a subset of possible_target_vars
(possible_target_vars <- intersect(puf_target_vars, available_target_vars)) # 21 targets

# These are the possible targets until I go in and get more
# "N1_nnz"
# "MARS1_nnz", "MARS2_nnz"
# "posAGI_nnz",
# "A00100_sum"
# "N00200_nnz", "A00200_sum" # Wages
# "N00700_nnz", "A00700_sum" # State income tax refunds
# "N01000_nnz", "A01000_sum" # Net capital gain or loss
# "N04470_nnz", "A04470_sum" # Total deductions (standard or itemized)
# "N04800_nnz", "A04800_sum" # Taxable income
# "N05800_nnz", "A05800_sum" # Income tax before credits
# "N09600_nnz", "A09600_sum" # Alternative minimum tax
# "N17000_nnz", "A17000_sum" # Total medical and dental expense deduction

# ivars <- c(1, 4, 6)
# ivars <- c(1, 3, 5, 13, 14)
# ivars <- c(1:4, 6:21) # 21 is top
# (target_vars <- possible_target_vars[ivars])
# target_vars <- c("N1_nnz", "MARS2_nnz", "posAGI_nnz", "N00200_nnz",
#                  "A00100_sum", "A00200_sum")
# target_vars <- setdiff(possible_target_vars, c("MARS1_nnz", "MARS2_nnz", "A09600_sum", "N09600_nnz")) # A09600_sum N09600_nnz
# target_vars <- setdiff(possible_target_vars, c("MARS1_nnz", "MARS2_nnz")) # A09600_sum N09600_nnz
# A09600_sum N09600_nnz have zeroes -- must adjust

# define target_vars_df which will have different target vars for different
# STATE-AGI_STUB combinations
# we need to define fewer targets for income group 1, I think
possible_target_vars

# minimalist
# default_targets <- c("N1_nnz", "MARS2_nnz",
#                "posAGI_nnz", "A00100_sum",
#                "N00200_nnz", "A00200_sum")
# incgroup1_targets <- c("N1_nnz", "MARS2_nnz",
#                "A00100_sum",
#                "N00200_nnz", "A00200_sum",
#                "N09600_nnz", "A09600_sum")

# best attempt
default_targets <- c("N1_nnz", "MARS2_nnz",
               "posAGI_nnz", "A00100_sum",
               "N00200_nnz", "A00200_sum",
               "N00700_nnz", "A00700_sum",
               "N01000_nnz", "A01000_sum",
               "N04470_nnz", "A04470_sum",
               "N04800_nnz", "A04800_sum",
               "N05800_nnz", "A05800_sum",
               "N09600_nnz", "A09600_sum",
               "N17000_nnz", "A17000_sum")

# basic targets for income group 1 and for OA
basic_targets <- c("N1_nnz", "MARS2_nnz",
                       "A00100_sum",
                       "N00200_nnz", "A00200_sum",
                       "N05800_nnz", "A05800_sum",
                       "N09600_nnz", "A09600_sum")

# create lists of each kind of target
# targs_default <- list(default_targets)
# targs_ig1 <- list(incgroup1_targets)

# data frame, each AGI_STUB has a vector of target vars
# get bad targets
targets_bad <- shareslong %>%
  filter(is.na(target) | target==0) %>%
  select(AGI_STUB, STATE, targname) %>%
  group_by(AGI_STUB, STATE) %>%
  summarise(badtargets=list(as.vector(targname)), .groups="drop")

remove_bad <- function(target_vec, badtargets) {
  diff <- setdiff(unlist(target_vec), unlist(badtargets))
  list(diff)
  }
target_vars_df <- expand_grid(STATE=unique(puf_targ$STATE),
                              AGI_STUB=unique(puf_targ$AGI_STUB)) %>%
  left_join(targets_bad, by = c("STATE", "AGI_STUB")) %>%
  rowwise() %>%
  mutate(target_vec=case_when(AGI_STUB == 1 ~ list(basic_targets),
                              STATE == "OA" ~ list(basic_targets),
                              TRUE ~ list(default_targets)),
         target_vec=remove_bad(target_vec, badtargets))
target_vars_df %>% filter(STATE=="OA")

df <- target_vars_df %>% filter(STATE=="OA", AGI_STUB==1)
basic_targets
df$badtargets[[1]]
df$target_vec[[1]]


# 5. prepare initial state weights one at a time ----
# I found in prior investigation that LM is the best method
# notes on parallel:
# https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html
# https://cfss.uchicago.edu/notes/split-apply-combine/

#.. 5-a Loop through a set of states ----

#.... test the one_state function ----
# OA 1 RI 1 DC 1 OA 5 AL 2 NY 1 WY 2
st <- "WY"; ig <- 2 # NH MS RI
st <- "SD"; ig <- 2 #
st <- "WV"; ig <- 3 #
st <- "OA"; ig <- 5 #
st <- "NY"; ig <- 1 #

df <- target_vars_df %>% filter(STATE==st, AGI_STUB==ig)
basic_targets
df$badtargets[[1]]
df$target_vec[[1]]

# # SD 2 (9600) WV 3 (9600)
tvdf2 <- target_vars_df %>%
  mutate(target_vec=ifelse(STATE==st & AGI_STUB==ig, targs_ig1, list(target_vec)))
tvdf2 %>% filter(STATE==st & AGI_STUB==ig)
tvdf2 <- target_vars_df

res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, quiet=FALSE)

# res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, method='Newton', quiet=FALSE)
# res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, method='Broyden', quiet=FALSE)

res$h; res$s; res$k
res$etime
res$solver_message
res$sse_unweighted
res$sse_weighted
# res$targets %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_calc %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_diff %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_pctdiff %>% round(2)
# cbind(res$sortid, res$wh)

comp <- bind_rows(as_tibble(res$targets, rownames = "stabbr") %>%
                    mutate(type="target"),
                  as_tibble(res$targets_calc, rownames = "stabbr") %>%
                    mutate(type="calc")) %>%
  pivot_longer(-c(stabbr, type)) %>%
  pivot_wider(names_from = type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100)
comp %>% filter(stabbr != 'ZZ') %>%
  kable(digits=c(0, 0, 0, 0, 0, 1), format.args=list(big.mark = ','))


#.... loop ----
cluster <- new_cluster(6)

# CAUTION: unfortunately no progress reporting when run in parallel
parallel <- TRUE
# parallel <- FALSE

if(parallel){
  # set the latest versions of functions, etc. up for the run
  cluster_copy(cluster, c('one_state', 'puf_targ',
                          'pufstrip', 'target_vars_df')) # functions and data not in a library
  cluster_library(cluster, c("dplyr", "microweight", "tidyr", "purrr"))
}

a <- proc.time()
res_df <- puf_targ %>%
  select(AGI_STUB, STATE) %>%
  # filter(AGI_STUB %in% 2:3, STATE %in% c("AL", "CA", "IL")) %>%
  # filter(AGI_STUB == 1) %>%
  # filter(STATE == "OA") %>%
  group_by(AGI_STUB, STATE) %>%
  nest()  %>%
  {if (parallel) partition(., cluster) else .} %>%
  mutate(res=map(STATE, one_state, incgroup=AGI_STUB, target_vars_df=target_vars_df, quiet=TRUE)) %>%
  {if (parallel) collect(.) else .} %>%
  ungroup %>%
  arrange(AGI_STUB, STATE) # must sort if parallel
b <- proc.time()
b - a # seconds
(b - a) / 60 # minutes 69 minutes

# save results as the above can take a long time to run
# system.time(saveRDS(res_df, here::here("ignore", "res_df.rds"))) # 33 secs
# res_df <- readRDS(here::here("ignore", "res_df.rds"))

names(res_df$res[[1]])

res_df %>%
  unnest_wider(res)

# examine solver messages and sse_weighted
tmp <- res_df %>%
  ungroup %>%
  hoist(res, "h", "solver_message", "sse_weighted")
tmp
count(tmp, h)
count(tmp, solver_message)
count(tmp, AGI_STUB, solver_message)
tmp %>%
  group_by(AGI_STUB) %>%
  summarise(sum(is.na(sse_weighted)))
quantile(tmp$sse_weighted, na.rm=TRUE)

tmp %>% filter(AGI_STUB==10, STATE=="RI")

tmp %>%
  arrange(-sse_weighted)

rec <- res_df %>%
  filter(AGI_STUB==2, STATE=="OA") %>%
  .$res
rec[[1]]$targets_pctdiff %>% round(2)


#.. 5-b retrieve state weights and examine puf ----
# res_df %>%
#   unnest_wider(res) %>%
#   mutate(bvl=map(beta_opt_mat, as.vector)) %>%
#   unnest_longer(bvl) %>%
#   select(AGI_STUB, STATE, bvl)

# iglook <- 10

#.. retrieve original national weights, one set per AGI_STUB
wts_wh <- res_df %>%
  select(AGI_STUB, res) %>%
  group_by(AGI_STUB) %>%
  slice(n=1) %>%
  ungroup %>% # ungroup is needed but I don't know why
  hoist(res, "sortid", "wh") %>%
  select(AGI_STUB, sortid, wh) %>%
  unnest(c(sortid, wh))

wts_wh %>%
  group_by(AGI_STUB) %>%
  slice_head(n=4)

#.. retrieve initial state weights
# res_df$res[[1]]$whs[1:10, ] # examine the whs matrix; AGI_STUB 1 is bad
wts_whs <- res_df %>%
  unnest_wider(res) %>%
  select(AGI_STUB, STATE, sortid, whs) %>%
  mutate(whs=map(whs, function(mat) mat[, 1])) %>% # get whs vector for each state
  unnest(c(sortid, whs)) %>%
  pivot_wider(names_from = STATE, values_from=whs)

#.. combine weights and calculate new national weight as sum of state weights
wts_all <- wts_wh %>%
  left_join(wts_whs, by = c("AGI_STUB", "sortid")) %>%
  mutate(wh_sum=rowSums(across(-c(AGI_STUB, sortid, wh)))) %>%
  select(AGI_STUB, sortid, wh, wh_sum, everything())
glimpse(wts_all)

# examine the new national weight in comparison to original
wts_all %>%
  select(AGI_STUB, sortid, wh, wh_sum) %>%
  mutate(ratio = wh_sum / wh) %>%
  group_by(AGI_STUB) %>%
  slice_head(n=4)

wts_all %>%
  mutate(ratio=wh_sum / wh) %>%
  group_by(AGI_STUB) %>%
  do(qtiledf(.$ratio, probs=probs1_99))

# compare calculated national targets with new state-sum weights to targets ----
calctargs <- pufstrip %>%
  right_join(wts_all %>% select(AGI_STUB, sortid, wh_sum),
             by = c("sortid", "AGI_STUB")) %>%
  pivot_longer(-c(sortid, AGI_STUB, s006, wh_sum)) %>%
  mutate(target=s006 * value, calc=wh_sum * value) %>%
  group_by(AGI_STUB, name) %>%
  summarise(target=sum(target, na.rm=TRUE),
            calc=sum(calc, na.rm=TRUE),
            .groups="drop") %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100)
calctargs

calctargs %>%
  arrange(-abs(pdiff))

# compare calculated targets to target values, by state and AGI stub ----
calc_targ <- calc_targets(wtsdf=wts_all, pufdf=pufstrip)

target_comp <- bind_rows(puf_targ %>% mutate(type="target"),
                         calc_targ %>% mutate(type="calc")) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE, type),
               names_to = "targname") %>%
  pivot_wider(names_from=type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  left_join(target_vars_df %>%
              select(AGI_STUB, STATE, targname=target_vec) %>%
              unnest(targname) %>%
              mutate(targeted=TRUE),
            by = c("AGI_STUB", "STATE", "targname")) %>%
  mutate(targeted=ifelse(is.na(targeted), FALSE, TRUE))

target_comp %>%
  filter(STATE=="CT") %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="NY", AGI_STUB==10) %>%
  arrange(-abs(pdiff))


# 6. reweight the national file to hit national targets ----
# get the national targets
ustargs <- puf_targ %>%
  group_by(AGI_STUB) %>%
  select(-STATE) %>%
  summarise(across(.fns=sum), .groups="drop")
ustargs

#.. test the function ----
d <- stub_opt(1)
d$solver_message
quantile(d$result$solution)
d$weights_df


#.. loop through the stubs ----
# cluster <- new_cluster(6)

# CAUTION: unfortunately no progress reporting when run in parallel
# parallel <- TRUE
parallel <- FALSE

if(parallel){
  # set the latest versions of functions, etc. up for the run
  cluster_copy(cluster,
               c('stub_opt', 'ustargs', 'wts_all',
                 'pufstrip', 'get_cc_national', 'get_inputs_national',
                 'get_conbounds', 'scale_inputs', 'define_jac_g_structure_sparse',
                 'eval_f_xm1sq', 'eval_grad_f_xm1sq',
                 'eval_g', 'eval_jac_g',
                 'eval_h_xm1sq')) # functions and data not in a library
  cluster_library(cluster, c("dplyr", "ipoptr", "plyr", "tidyr", "purrr"))
}

a <- proc.time()
opt_df <- tibble(AGI_STUB=1:10) %>%
  group_by(AGI_STUB) %>%
  nest()  %>%
  {if (parallel) partition(., cluster) else .} %>%
  mutate(output=map(AGI_STUB, stub_opt)) %>%
  {if (parallel) collect(.) else .} %>%
  ungroup %>%
  arrange(AGI_STUB) # must sort if parallel
b <- proc.time()
b - a # seconds
(b - a) / 60 # minutes 69 minutes

# saveRDS(opt_df, here::here("ignore", "opt_df.rds"))
# opt_df <- readRDS(here::here("ignore", "opt_df.rds"))


names(opt_df$output[[1]])

tmp <- opt_df %>%
  hoist(output, "solver_message")
tmp

# rlang::last_error()
# rlang::last_trace()


# 7. calculate adjusted state weights that hit the adjusted national total ----

# mutate(whs=map(whs, function(mat) mat[, 1])) %>% # get whs vector for each state
# wts_whs <- res_df %>%
#   unnest_wider(res) %>%
#   select(AGI_STUB, STATE, sortid, whs) %>%
#   mutate(whs=map(whs, function(mat) mat[, 1])) %>% # get whs vector for each state
#   unnest(c(sortid, whs)) %>%
#   pivot_wider(names_from = STATE, values_from=whs)

names(opt_df$output[[1]])
# opt_df$output[[1]]$weights_df # AGI_STUB sortid iweight     x whs_adj
wts_adj <- opt_df %>%
  hoist(output, "weights_df") %>%
  select(AGI_STUB, weights_df) %>%
  group_by(AGI_STUB) %>%
  mutate(weights_df=map(weights_df, function(x) x %>% select(-AGI_STUB))) %>%
  unnest(col=weights_df) %>%
  ungroup %>%
  left_join(wts_all, by = c("AGI_STUB", "sortid")) %>%
  mutate(across(-c(AGI_STUB, sortid, iweight, x, wh_adj, wh, wh_sum),
                function(wt) wt * x)) %>%
  mutate(wh_adjsum=rowSums(across(-c(AGI_STUB, sortid, iweight, x, wh_adj, wh, wh_sum)))) %>%
  select(AGI_STUB, sortid, iweight, x, starts_with("wh"), everything())


calc_targ_adj <- calc_targets(wtsdf=wts_adj, pufdf=pufstrip)


target_comp_adj <- bind_rows(puf_targ %>% mutate(type="target"),
                             calc_targ_adj %>% mutate(type="calc")) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE, type),
               names_to = "targname") %>%
  pivot_wider(names_from=type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  left_join(target_vars_df %>%
              select(AGI_STUB, STATE, targname=target_vec) %>%
              unnest(targname) %>%
              mutate(targeted=TRUE),
            by = c("AGI_STUB", "STATE", "targname")) %>%
  mutate(targeted=ifelse(is.na(targeted), FALSE, TRUE))

target_comp_adj %>%
  filter(STATE=="CT") %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="CT") %>%
  arrange(-abs(pdiff))

target_comp_adj %>%
  filter(STATE=="NY", AGI_STUB==10) %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="NY", AGI_STUB==10) %>%
  arrange(-abs(pdiff))

wts_all_natladj <- wts_all %>%
  left_join(iweights2 %>% select(sortid, wh_sumadj), by = "sortid") %>%
  select(AGI_STUB, sortid, wh, wh_sum, wh_sumadj, everything())


# 8. construct final state weights ----
#.. 7-a use poisson method ----
# targets
# grab the targets for this STATE-AGI_STUB combination
ig <- 1

wh <- wts_all_natladj %>%
  filter(AGI_STUB==ig) %>%
  .$wh_sumadj # the new national weight -- use either wh_sum or wh_sumadj

# target_vars <- default_targets
target_vars <- target_vars_df %>%
  filter(STATE=="AL", AGI_STUB==ig) %>%
  .$target_vec %>%
  unlist

targets <- puf_targ %>%
  filter(AGI_STUB==ig) %>%
  select(AGI_STUB, STATE, all_of(target_vars))

targmat <- targets %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(targmat) <- targets$STATE
targmat

xmat <- pufstrip %>%
  filter(AGI_STUB==ig) %>%
  select(all_of(target_vars)) %>%
  as.matrix

# res1 <- geoweight(wh = wh, xmat = xmat, targets = targmat, method = 'Broyden')
res2 <- geoweight(wh = wh, xmat = xmat, targets = targmat, method = 'LM') # best
# res3 <- geoweight(wh = wh, xmat = xmat, targets = targmat, method = 'Newton')

res1$targets_pctdiff %>% round(2)
res2$targets_pctdiff %>% round(2)
res3$targets_pctdiff %>% round(2)

res2$h; res2$s; res2$k
res2$opts_used
res2$solver_message
res2$etime
res2$sse_unweighted

#.. 7-b use reweight method ----
# get the initial state weights from step 5
# rescale so that they hit the new national weights from step 6
# use ipopt with state constraints and adding-up constraints
# let's keep states in alpha order throughout, even OA

# prepare initial weights
glimpse(wts_all2)
iweights3 <- wts_all_natladj %>%
  filter(AGI_STUB==2) %>%
  pivot_longer(-c(AGI_STUB, sortid, wh, wh_sum, wh_sumadj),
               names_to = "STATE",
               values_to = "whs_raw") %>%
  group_by(AGI_STUB, sortid) %>%
  mutate(iweight=whs_raw / sum(whs_raw) * wh_sumadj) %>%
  ungroup %>%
  arrange(AGI_STUB, sortid, STATE) %>%
  mutate(j = row_number())
ht(iweights3)


# prepare targets
target_vars <- default_targets

targagg <- puf_targ %>%
  filter(AGI_STUB==2) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) %>%
  arrange(AGI_STUB, STATE) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE), names_to = "var_calc", values_to = "value") %>%
  mutate(cname=paste0(STATE, "_", var_calc),
         targtype="aggregate") %>%
  separate(var_calc, c('var', 'calc')) %>%
  select(cname, value, AGI_STUB, STATE, var, calc, targtype)
ht(targagg)

targadd <- wts_all_natladj %>%
  filter(AGI_STUB == 2) %>%
  arrange(sortid) %>%
  mutate(cname=paste0("p", str_pad(sortid, width=8, side="left", pad="0")),
         targtype = "addup",
         var=as.character(sortid),
         STATE="ZZ") %>%
  select(cname, value=wh_sumadj, AGI_STUB, STATE, var, targtype)
ht(targadd)

targdf <- bind_rows(targagg, targadd) %>%
  arrange(desc(targtype), cname) %>%
  mutate(i = row_number())
ht(targdf)


# prepare income-group data
glimpse(pufstrip)
igdata <- pufstrip %>%
  filter(AGI_STUB==2)


cc_sparse <- get_cc_states(.incgroup_data = igdata,
                             .target_vars = target_vars,
                             .iweights = iweights3,
                             .targets_df = targdf)

inputs <- get_inputs(.targets_df = targdf,
                     .iweights = iweights3,
                     .cc_sparse = cc_sparse,
                     .targtol=.15, .xub=20, .conscaling=FALSE, scale_goal=1)
names(inputs)
inputs$n_variables; inputs$n_constraints; inputs$n_targets

condf <- targdf %>%
  mutate(tol=ifelse(targtype == "addup", 0, .005),
         tol=ifelse(STATE == "OA", .06, tol),
         tol=ifelse(cname == "OA_A00200_sum", .25, tol)) %>%
  mutate(clb = value - abs(tol * value),
         cub = value + abs(tol * value))
ht(condf)
count(condf, tol)
condf %>% filter(tol > .01)
# A00200_sum

inputs$clb <- condf$clb
inputs$cub <- condf$cub

# For LARGE problems use linear_solver=ma77, obj_scaling=1, and mehrotra_algorithm=yes
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 50,
             "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             #"linear_system_scaling" = "mc19",
             #"linear_scaling_on_demand" = "no", # default is yes -- no means ALWAYS scale
             # "ma57_automatic_scaling" = "yes", # if using ma57
             # "ma57_pre_alloc" = 3, # 1.05 is default; even changed, cannot allocate enough memory, however
             # "ma77_order" = "amd",  # metis; amd -- not clear which is faster
             "mehrotra_algorithm" = "yes",
             #"obj_scaling_factor" = 1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             #"nlp_scaling_method" = "none", # NO - use default gradient_based, none, equilibration-based
             # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
             "jac_c_constant" = "yes", # does not improve on moderate problems equality constraints
             "jac_d_constant" = "yes", # does not improve on  moderate problems inequality constraints
             "hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "limited_memory_update_type" = "bfgs", # default bfgs; sr1 docs say "not working well" but this really helps
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = here::here("test7.out"))


# setwd(here::here("temp1"))
# getwd()
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
names(result)
quantile(result$solution, probs=c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))

ht(inputs$cc_sparse)



# 8. check results ----
# [TO COME]

