THISPKG <- "Bureau"
.onAttach <- function(libname, pkgname) {
	version <- packageDescription("Bureau", fields="Version")
	packageStartupMessage(paste("
Welcome to Bureau version ", version, "\n", sep = "" ) )
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
