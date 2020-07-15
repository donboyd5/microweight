
remove_suffix <- function(suffix, target_vars) {
  # return a named vector with the base names of target variables
  var_backend <- str_extract(target_vars, "_.*$") # gets the ending part of each variable name
  base_vars <- target_vars[which(var_backend==suffix)] %>% str_remove(suffix)
  names(base_vars) <- base_vars
  base_vars
}


one_state <- function(st, incgroup, target_vars_df, quiet=TRUE, method='LM'){
  print(sprintf("State: %s, Income group: %i", st, incgroup))

  # grab the targets for this STATE-AGI_STUB combination
  target_vars <- target_vars_df %>%
    filter(STATE==st, AGI_STUB==incgroup) %>%
    .$target_vec %>%
    unlist

  targs <- puf_targ %>%
    filter(AGI_STUB==incgroup) %>%
    select(AGI_STUB, STATE, all_of(target_vars)) %>%
    mutate(STATE=ifelse(STATE==st, STATE, "ZZ")) %>%
    group_by(AGI_STUB, STATE) %>%
    summarise(across(all_of(target_vars), sum), .groups="drop")

  targmat <- targs %>%
    select(all_of(target_vars)) %>%
    as.matrix
  rownames(targmat) <- targs$STATE
  targmat

  xmat <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    select(all_of(target_vars)) %>%
    as.matrix

  sortid <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    .$sortid

  wh <- pufstrip %>%
    filter(AGI_STUB==incgroup) %>%
    .$s006

  p <- list()
  p$s <- nrow(targmat)
  p$k <- ncol(targmat)
  p$h <- length(wh)
  p$wh <- wh
  p$targets <- targmat
  p$xmat <- xmat

  result <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = method, quiet=quiet)
  result$sortid <- sortid

  return(result)
}

