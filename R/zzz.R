##' @useDynLib mmsu, .registration = TRUE
##' @importFrom odin odin
##' @importFrom utils modifyList
##' @importFrom ICDMM equilibrium_init_create model_param_list_create
NULL

.onLoad <- function(libname, pkgname) {
  # Check and install ICDMM package
  if (!requireNamespace("ICDMM", quietly = TRUE)) {
    devtools::install_github("mrc-ide/deterministic-malaria-model")
  }
}
