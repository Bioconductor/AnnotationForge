.onLoad <- function(libname, pkgname) {
  ## load the utils::data
  where = asNamespace(pkgname)
  utils::data(list = pkgname, package = pkgname, envir = where)
}
