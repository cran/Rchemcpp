#ifndef CONSTANT_H
#define CONSTANT_H

#define MDL2000 2000
#define MDL2000JLP 2001
#define MDL3000 3000

#define NAVALUE -9999

#define INTEGER 0
#define FLOAT 1
#define STRING 3
#define DEFAULT 4

#define MULTITHREAD 1
#define MAXTHREADS 100

#ifndef CHEMCPPPATH
#endif

#define ELEMENTSFILENAME CHEMCPPPATH "/data/elements.csv"	// semicolumn (;) separated file with  each line describing an element

#define KEGGATOMS CHEMCPPPATH "/data/keggatoms.csv"

#define NUMELEMENTS 109	// maximum number of elements
#define GRAMATOMSIDENTITYFILENAME CHEMCPPPATH "/data/gramAtoms.binary.csv"           // KASHIMA ORIGINAL
#define GRAMATOMSDOTPRODUCT1FILENAME CHEMCPPPATH "/data/gramAtoms.dotproduct.csv"    // WARNING WRONG

#define GRAMATOMSCRENDOTPRODUCTFILENAME CHEMCPPPATH "/data/gramAtoms.cren.dotproduct.csv"
#define GRAMATOMSCRENPOLYNOMIALFILENAME CHEMCPPPATH "/data/gramAtoms.cren.polynomial.csv"
#define GRAMATOMSCRENGAUSSIANFILENAME CHEMCPPPATH "/data/gramAtoms.cren.gaussian.csv"
#define GRAMKEGGATOMSIDENTITYFILENAME CHEMCPPPATH "/data/gramKeggAtoms.kashima.csv"

#define PS "ps"  // name used to designate the Kashima Start Probability
#define PT "pt"  // name used to designate the Kashima Transition Probability
#define PQ "pq"  // name used to designate the Kashima Stop Probability

// ATOM PROPERTIES

#define FULLNAME "Name"		// name used to designate the Full name
#define SYMBOLNAME "Symbol"	// name used to designate the symbol letter in the elements file
#define AN "An"			// name used to designate the atomic number
#define AM "Am" 		// name used to designate the atomic mass
#define VWR "Vwr"		// name used to designate the van der Waals radius
#define CR "Cr"			// Covalent radius
#define AV "Av"             	// Atomic volume
#define MP "Mp"             	// Melting point
#define BP "Bp"             	// Boiling point
#define VALENCE "Valence"   	// Valence
#define EA1 "Ea1"           	// First electron affinity
#define IE1 "Ie1"           	// First ionization energy
#define IE2 "Ie2"           	// second ionization energy
#define IE3 "Ie3"           	// third ionization energy
#define IE4 "Ie4"           	// fourth ionization energy
#define IE5 "Ie5"           	// fifth ionization energy
#define IE6 "Ie6"           	// sixth ionization energy
#define IE7 "Ie7"           	// seventh ionization energy
#define IE8 "Ie8"           	// eigth ionization energy
#define EN "En"             	// Electronegativity

// MOLECULE PROPERTIES

#define ACTIVITY "activity" // stores the biological activity

// BOND LABELS

#define NOBOND 0
#define SINGLEBOND 1
#define DOUBLEBOND 2
#define TRIPLEBOND 3
#define AROMATICBOND 4
#define SINGLECYCLEBOND 5
#define DOUBLECYCLEBOND 6
#define TRIPLECYCLEBOND 7


// CALCULATIONS

#define GRAMDECIMALPRECISION 9		// number of decimals to print in the gram matrix

#define MAXSIMILARITYWARNING 0.999999	// number between 0 and 1. If the normalised kernel
					// exceeds this value between two molecules a warning
					// is printed

#define MINSIMILARITYWARNING 1-MAXSIMILARITYWARNING


#define DUPLICATEMWIDENTITY 0.0001      // in MoleculeSet::removeDuplicates, if mw1-mw2 < DUPLICATEMWIDENTITY then
                                        // two compounds are considered to have identical MW

#define COMPONENTINDEX "componentIndex" // name of the descriptor used for marking fragments in a
					// molecule.

#endif

#include <string>
#include <stdlib.h>

std::string GETChemcpppath();
#define CHEMCPPPATH GETChemcpppath() +	

