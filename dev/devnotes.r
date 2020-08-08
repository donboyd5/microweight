
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

devtools::session_info()

# resources ----
# http://r-pkgs.had.co.nz/
# Developing packages: https://support.rstudio.com/hc/en-us/articles/200486488?version=1.3.1039&mode=desktop
# Build, test: https://support.rstudio.com/hc/en-us/articles/200486508-Building-Testing-and-Distributing-Packages
# Documenting: https://support.rstudio.com/hc/en-us/articles/200532317-Writing-Package-Documentation
# Documenting: https://cran.r-project.org/doc/manuals/R-exts.html#Writing-R-documentation-files
# CRAN release: https://cran.r-project.org/web/packages/policies.html

# Good intro: https://johnmuschelli.com/smi_2019/index.html
# Easier to use: https://github.com/muschellij2/smi_2019/blob/master/index.pdf

# https://cran.r-project.org/doc/contrib/Leisch-CreatingPackages.pdf
# http://portal.stats.ox.ac.uk/userdata/ruth/APTS2012/Rcourse10.pdf
# https://cran.r-project.org/doc/manuals/R-exts.html



# examples to learn from ----
# data package:  https://github.com/hadley/fueleconomy
# code package:  https://github.com/hadley/cubelyr


# general advice ----



# workflow ----
# Ctrl+Shift+D Update documentation
# Ctrl+Shift+B Build and reload


# documentation ----
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Sectioning
# https://support.rstudio.com/hc/en-us/articles/200532317-Writing-Package-Documentation
# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html


# package development tools ----
library(usethis)
library(devtools)
library(assertthat)


# package setup ----
#.. create a  readme shell ----
usethis::use_readme_rmd()

# usethis::use_gpl3_license()
# usethis::use_gpl3_license(name = find_name())


#.. add packages etc. to DESCRIPTION ----
usethis::use_package("dplyr") # Defaults to imports
usethis::use_package("ipoptr", "Suggests")
usethis::use_package("minpack.lm") # Defaults to imports
usethis::use_package("nleqslv") # Defaults to imports
usethis::use_package("nloptr")
usethis::use_package("purrr")
usethis::use_package("rlang") # so that we can use .data I think
usethis::use_package("stats")
usethis::use_package("tibble")
usethis::use_package("tidyr")

usethis::use_roxygen_md()
usethis::use_pipe(export = TRUE)

# usethis::use_package("tidyverse") # Defaults to imports ERROR

# examples of rbuildignore (DO NOT RUN) ----
# use_build_ignore("./data/sourceData/", escape = TRUE, pkg = ".")
# use_build_ignore("./R/DataConversionPrograms/", escape = TRUE, pkg = ".")


# package coding ----
#.. code formatting ----
# install.packages("formatR")
formatR::tidy_dir("R")

# install.packages("lintr")
lintr::lint_package()


# vignettes ----
# http://r-pkgs.had.co.nz/vignettes.html
# devtools::use_vignette("my-vignette")

# This will:
#
#   Create a vignettes/ directory.
#   Add the necessary dependencies to DESCRIPTION (i.e. it adds knitr to the
#     Suggests and VignetteBuilder fields).
#   Draft a vignette, vignettes/my-vignette.Rmd.
#
# The draft vignette has been designed to remind you of the important parts of
# an R Markdown file. It serves as a useful reference when you’re creating a
# new vignette.
#
# Once you have this file, the workflow is straightforward:
#
#   Modify the vignette.
#   Press Ctrl/Cmd + Shift + K (or click ) to knit the vignette and preview the output.
#
# There are three important components to an R Markdown vignette:
#
#   The initial metadata block.
#   Markdown for formatting text.
#   Knitr for intermingling text, code and results.

# You can build all vignettes from the console with devtools::build_vignettes(),
# but this is rarely useful. Instead use devtools::build() to create a package
# bundle with the vignettes included. RStudio’s “Build & reload” does not build
# vignettes to save time. Similarly, devtools::install_github() (and friends) will
# not build vignettes by default because they’re time consuming and may require additional
# packages. You can force building with devtools::install_github(build_vignettes = TRUE).
# This will also install all suggested packages.


# Do NOT RERUN the lines below as they will write over old vignette files - use
# ONLY to create initial vignette
# usethis::use_vignette("reweighting")
# usethis::use_vignette("geoweighting")


# devtools::build() will create package here:
# "C:/RPrograms PC/Packages/microweight_0.1.0.tar.gz"

# install with:
#  install.packages("C:/RPrograms PC/Packages/microweight_0.1.0.tar.gz", repos = NULL, type="source", dependencies = TRUE)
browseVignettes(package = "microweight")

# manual ----
devtools::build_manual(pkg = ".", path = NULL)


# CRAN policies ----
# https://cran.r-project.org/web/packages/policies.html
# CRAN submission
# https://cran.r-project.org/web/packages/submission_checklist.html


