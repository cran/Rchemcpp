
#ifndef RMOLECULE_H
#define RMOLECULE_H

#include <Rcpp.h>
#include <RcppClassic.h>

#include <molecule.h>

class Rmolecule : public Molecule {
public:
	
	Rmolecule(){};
	
	void addAtom(string aSymbol)
	{Molecule::addAtom(aSymbol);}

	void linkAtoms( int firstAtom, int secondAtom, int aBondLabel)
	{Molecule::linkAtoms(firstAtom, secondAtom, aBondLabel);}

	void writeSD( string aFileName )
	{Molecule::writeSD(aFileName);}

	//...


//private:

};

#endif

