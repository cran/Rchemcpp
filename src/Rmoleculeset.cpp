
#include "Rmoleculeset.h"

#include<string>
#include<iostream>
#include<sstream>




RCPP_MODULE(Rmoleculeset){
	Rcpp::class_<Rmoleculeset>( "Rmoleculeset")
	.constructor()
//	.constructor<int,int>()


	.method( "addSD", &Rmoleculeset::addSD, "read molecule" )
	.method( "addSD2", &Rmoleculeset::addSD2, "read molecule" )
//Rcpp::List::create( Rcpp::_["aFileName"] Rcpp::_["genericAtomTypeFlag"] = 0, Rcpp::_["beginMolecule"] = -1, Rcpp::_["endMolecule"] = -1  ),		

	.method( "addKCF", &Rmoleculeset::addKCF, "read molecule" )
	.method( "addKCF2", &Rmoleculeset::addKCF2, "read molecule" )

	.method( "numMolecules", &Rmoleculeset::numMolecules, "return the number of molecules in the set")

	.method( "hideHydrogens", &Rmoleculeset::hideHydrogens, "hide hydrogens")

	.method( "setMorganLabels", &Rmoleculeset::setMorganLabels, "set Morgan labels")

	.method( "setKashimaKernelParam", &Rmoleculeset::setKashimaKernelParam, "set Kashima Kernel parameter")	

	.method( "initializeGram", &Rmoleculeset::initializeGram, "initialize the Gram matrix")

	.method( "initializeSelfKernel", &Rmoleculeset::initializeSelfKernel, "initialize the Self Kernel")

	.method( "setComparisonSetSelf", &Rmoleculeset::setComparisonSetSelf, "set comparison set to be the moleculeset itself;") 

	.method( "setComparisonSetCopy", &Rmoleculeset::setComparisonSetCopy, "set comparison set; argument-object is cloned! (and disposed of at destruction)")	

	.method( "getComparisonSet", &Rmoleculeset::getComparisonSet, "Pointer to the current comparison set; managed by chemcpp; Object becomes invalid, when the moleculeset runs ouf of scope.") //new method	

	.method( "normalizeTanimoto", &Rmoleculeset::normalizeTanimoto, "normalize Tanimoto")

	.method( "normalizeTanimotoMinMax", &Rmoleculeset::normalizeTanimotoMinMax, "normalize TanimotoMinMax")

	.method( "normalizeGram", &Rmoleculeset::normalizeGram, "normalize Gram")

	.method( "noTottersTransform", &Rmoleculeset::noTottersTransform, "transform to avoid totters")

	.method( "readPartialCharges", &Rmoleculeset::readPartialCharges, "add partial charges from file")

	.method( "setMorganChargesLabels", &Rmoleculeset::setMorganChargesLabels, "set Morgan charges threshold")

	.method( "gramCompute3D", &Rmoleculeset::gramCompute3D, "compute 3d gram")

	.method( "gramCompute", &Rmoleculeset::gramCompute, "compute gram")

	.method( "addMoleculeCopy", &Rmoleculeset::addMoleculeCopy, "Clone Molecule and add it")

	.method( "getGram", &Rmoleculeset::getGram, "Clone Molecule and add it")

	.method( "getGramNormal", &Rmoleculeset::getGramNormal, "Clone Molecule and add it")

	//currently not expicitly necessary
    	.method( "atomsLabelsListing", &Rmoleculeset::atomsLabelsListing, "get atom label listings")

    	.method( "bondsListing", &Rmoleculeset::bondsListing, "get bond listings")

	.method( "normalizeTanimoto_raw", &Rmoleculeset::normalizeTanimoto_raw, "normalize Tanimoto")

	.method( "normalizeGram_raw", &Rmoleculeset::normalizeGram_raw, "normalize Gram")

	.method( "writeGramMatrix", &Rmoleculeset::writeGramMatrix, "write Gram matrix to file")

	.method( "writeSelfKernelList", &Rmoleculeset::writeSelfKernelList, "write self kernel list")

	//Extension 1.0.7
	.method( "getMolByIndex", &Rmoleculeset::getMolByIndex, "get Molecule by Index (zero-based)")
	;
}

















