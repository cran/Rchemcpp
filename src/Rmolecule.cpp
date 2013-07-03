
#include "Rmolecule.h"



RCPP_MODULE(Rmolecule){
	Rcpp::class_<Rmolecule>( "Rmolecule")
	.constructor()
	.method( "addAtom", &Rmolecule::addAtom, "Add an atom" )
	.method( "linkAtoms", &Rmolecule::linkAtoms, "link two atoms (zero-indexed)" )
	.method( "writeSD", &Rmolecule::writeSD, "write an SD file")
	
	.method( "listStringDescriptors", &Rmolecule::listStringDescriptors, "list the string-descriptors")
	//uncommon for Molecules ... therefore commented out...	
	//.method( "listFloatDescriptors", &Rmolecule::listFloatDescriptors, "list the float-descriptors")
	//.method( "listIntDescriptors", &Rmolecule::listIntDescriptors, "list the int-descriptors")
	
	.method( "getStringDescriptorValue", &Rmolecule::getStringDescriptorValue, "get the value of a string-descriptor")
	.method( "getStringDescriptorUnit", &Rmolecule::getStringDescriptorUnit, "get the unit of a string-descriptor")
	.method( "getStringDescriptorComment", &Rmolecule::getStringDescriptorComment, "get the comment of a string-descriptor")

	//uncommon for Molecules ... therefore commented out...
	//.method( "getFloatDescriptorValue", &Rmolecule::getFloatDescriptorValue, "get the value of a float-descriptor")
	//.method( "getFloatDescriptorUnit", &Rmolecule::getFloatDescriptorUnit, "get the unit of a float-descriptor")
	//.method( "getFloatDescriptorComment", &Rmolecule::getFloatDescriptorComment, "get the comment of a float-descriptor")
	//.method( "getIntDescriptorValue", &Rmolecule::getIntDescriptorValue, "get the value of a int-descriptor")
	//.method( "getIntDescriptorUnit", &Rmolecule::getIntDescriptorUnit, "get the unit of a int-descriptor")
	//.method( "getIntDescriptorComment", &Rmolecule::getIntDescriptorComment, "get the comment of a int-descriptor")

	.method( "setStringDescriptor", &Rmolecule::setStringDescriptor, "set the value of a string-descriptor")
	.method( "deleteStringDescriptor", &Rmolecule::deleteStringDescriptor, "delete a string-descriptor")

	;
}









