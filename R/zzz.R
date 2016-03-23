.onLoad <- function(libname, pkgname) {
  # Do a messed up hack to fix the poLCA package.
  fixPoLCA = ".Csave = .C; .C = function(...) { a = list(...); a$PACKAGE='poLCA'; do.call(.Csave, a)};"
  trace("poLCA.probHat.C", tracer = parse(text = fixPoLCA), where = poLCA::poLCA, at = 1, print = F)
  trace("poLCA.postClass.C", tracer = parse(text = fixPoLCA), where = poLCA::poLCA, at = 1, print = F)
  trace("poLCA.ylik.C", tracer = parse(text = fixPoLCA), where = poLCA::poLCA, at = 1, print = F)
}