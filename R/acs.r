#' Extract from the American Community Survey (ACS).
#'
#' tibble with 10,000 person records from 5 states (CA, FL, IL, NY, TX)
#' from the 2013-2017 5-year ACS.
#'
#' @format A data frame with 10,000 rows and 17 variables. The most important
#' variables are described below. For details see ACS documentation at
#' \url{https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/PUMS_Data_Dictionary_2013-2017.pdf}:
#' \describe{
#'   \item{serialno}{person id}
#'   \item{sporder}{person order within household}
#'   \item{stabbr}{state postal abbreviation}
#'   \item{pwgtp}{person weight}
#'   \item{adjinc}{income adjustment factor}
#'   \item{mar}{marital status 1=married, ..., 5=never married or under 15 years old}
#'   \item{sex}{1=male, 2=female}
#'   \item{agep}{age}
#'   \item{pincp}{personal income}
#'   \item{wagp}{wage income}
#'   \item{intp}{interest income}
#'   \item{pap}{public assistance income}
#'   \item{retp}{retirement income}
#'   \item{ssip}{Supplemental Security Income}
#'   \item{ssp}{Social Security income}
#'   \item{otherincp}{other income}
#'   \item{incgroup}{income decile among individuals in the data}
#'   ... variables ending in p, except otherincp, are income items per ACS documentation
#'   (otherincp is calculated as the difference between pincp and other income
#'   items).
#' }
#' @source \url{https://www.census.gov/programs-surveys/acs/}
"acs"
