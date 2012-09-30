
#Cannot work with NAMESPACE, because environment-variable has to be set first...
.onLoad <- function(lib, pkg, ...) {
  Sys.setenv(RCHEMCPPPATH= system.file("chemcpp", package=pkg) )
  library.dynam(chname="Rchemcpp",package=pkg,lib.loc=lib)
}

#loadRcppModules(direct=TRUE)


#loadModule("Ratom", TRUE)
#loadModule("Rbond", TRUE)
loadModule("Rmolecule", TRUE)
loadModule("Rmoleculeset", TRUE)
#loadModule("Rnode", TRUE)
#loadModule("Rring", TRUE)

loadModule("Relements", TRUE)

loadModule("subtreehelper",TRUE)
loadModule("spectrumhelper",TRUE)
loadModule("spectrum3Dhelper",TRUE)

