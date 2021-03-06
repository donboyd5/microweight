#' Targets for Extract from the American Community Survey (ACS).
#'
#' weighted values for `acsbig` data frame. Can be used as targets.
#'
#' @format A data frame with 1 record for each combination of 10 income deciles and 20
#' states (AK, AL, AR, AZ, CA, CO, CT, DE, FL, GA, HI, IA, ID, IL, IN, KS, KY, LA, MD, ME).
#' \url{https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/PUMS_Data_Dictionary_2013-2017.pdf}:
#' \describe{
#'   \item{incgroup}{income decile among individuals in the data}
#'   \item{stabbr}{state postal abbreviation}
#'   \item{pop}{sum of person weights}
#'   \item{pincp}{personal income}
#'   \item{wagp}{wage income}
#'   \item{ssip}{Supplemental Security Income}
#'   \item{ssip_nnz}{weighted number of records for which Supplemental Security Income is nonzero}
#'   \item{young}{weighted number of records for which age < 21}
#'   \item{workage}{weighted number of records for which age is in 21 to 64 years}
#'   \item{married}{weighted number of married records (mar==1)}
#'   \item{female}{weighted number of female records (sex==1)}
#' }
#' @source \url{calculated from data at https://www.census.gov/programs-surveys/acs/}
"acsbig_targets"
