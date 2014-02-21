.onLoad <- function(libname, pkgname) {
  ## load the data
  where = asNamespace(pkgname)
  data(list = pkgname, package = pkgname, envir = where)
}
