
#include "Rmolecule.h"



RCPP_MODULE(Rmolecule){
	Rcpp::class_<Rmolecule>( "Rmolecule")
	.constructor()
	.method( "addAtom", &Rmolecule::addAtom, "Add an atom" )
	.method( "linkAtoms", &Rmolecule::linkAtoms, "link two atoms (zero-indexed)" )
	.method( "writeSD", &Rmolecule::writeSD, "write an SD File")
	;
}









