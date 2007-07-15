.First.lib <- function(lib, pkg)
  library.dynam("sjevor", pkg, lib)

.Last.lib <- function (libpath) {
  ## Run when the package is being unloaded.  This allows us to test
  ## packages within same session when dynlib is updated.
  library.dynam.unload("sjevor", libpath)
}
