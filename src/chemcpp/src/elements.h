/****************************************************************************************
					  elements.h 
					-------------
    copyright            : (C) 2006 Jean-Luc Perret - Pierre Mahé
    email                : jean-luc.perret@unine.ch - pierre.mahe@ensmp.fr
 ***************************************************************************************/

/****************************************************************************************
 *                                                                         		*
 *	This program is free software; you can redistribute it and/or			*
 * 	modify it under the terms of the GNU Lesser General Public			*
 * 	License as published by the Free Software Foundation; either			*
 * 	version 2.1 of the License, or (at your option) any later version.		*
 *											*
 *	This program is distributed in the hope that it will be useful,			*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of			*
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU		*
 *	Lesser General Public License for more details.					*
 *											*
 *	You should have received a copy of the GNU Lesser General Public		*
 * 	License along with this library; if not, write to the Free Software		*
 * 	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA	*
 *											*
 ****************************************************************************************/


#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>

//using namespace std;

#include <datacontainer.h>
#include <atom.h>
#include <stringutils.h>

/** Instanciates chemical elements
  @author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
	@version 0.1
	@date 8 Jan 2004

	CLASS NAME: 	Elements

	FOR:					SNSF SPONSORED PROJECT

	PURPOSE:  		Instanciates the chemical elements (each of class Atom) and
								load their physico-chemical properties


	Physico-chemical properties are loaded from file ELEMENTSFILENAME in the chemcpp
	distribution.

	\warning Elements currently only instanciates the following elements: H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca with the full set of physico-chemical properties, and Br, I with only full name, atomic number and symbol.

	An Atom pointer to an Element can be obtained by using the elements global instance
	of class Elements, using the following syntax:

	 \code Atom* myPointer = elements["H"]; \endcode

	\warning note that to create a new Atom the copy operator should be used as described in "atom_example.cpp".

	The properties of the elements can be retrieve as pointers to a Descriptor by using
	the getIntDescriptor(), getFloatDescriptor() or getStringDescriptor() functions of the DataContainer class from which the class Atom is derived.
	To designate the physico-chemical properties use the macros defined in constants.h
	Here is a list of these macros:

	- FULLNAME: "Name" name used to designate the Full name
	- SYMBOLNAME: "Symbol" name used to designate the symbol letter in the elements file
	- AN: name used to designate the atomic number
	- AM: name used to designate the atomic mass
	- VWR: name used to designate the van der Waals radius
	- CR: Covalent radius
	- AV: Atomic volume
	- MP: Melting point
	- BP: Boiling point
	- VALENCE: Valence
	- EA1: First electron affinity
	- IE1: First ionization energy
	- IE2: second ionization energy
	- IE3: third ionization energy
	- IE4: fourth ionization energy
	- IE5: fifth ionization energy
	- IE6: sixth ionization energy
	- IE7: seventh ionization energy
	- IE8: eigth ionization energy
	- EN: Electronegativity

	for example to get the covalent radius value of an atom use:

	\code float myCR = myPointer->getFloatDescriptor(CR)->getValue(); \endcode


	\todo Implement deletion of KindDescriptors in destructor of class Elements

*/

class Elements : public DataContainer  {
public:

///@name General functions for building and manipulating the set of elements
//@{
	/** loads the elements from a file.
	*/
	Elements( string aFileName, string aGramAtomFileName );

	/** elements destructor.
			at the moment only a default desctructor, but could eventually delete all
			chemical elements this class created.
	*/
	~Elements();

	/** returns an atom pointer to an element designed using its chemical symbol.
	*/
	Atom* operator[](string aSymbol);

	/** returns a pointer to the element designed by aSymbol.
	*/
	Atom* getElement(string aSymbol);

	/** function to load a gram matrix to compare atoms.
			this matrix will then be used by the kashima kernel.
			by default without calling this function, the identity matrix
			is loaded by the elements constructor.
	*/
	void loadGramAtoms(string filename);

	/** returns the filename of the atom kernel matrix loaded.
		returns "default" if no kernel matrix was loaded.
	*/
	string getAtomKernelName();

	/** loads the definition of elements.
	*/
	void loadDefinition( string aFilename );

	/** returns the number of elements in periodicTable.
	*/
	uint numElements(){ return( periodicTable.size() ); }

//@}




// ****************************** //
// **** DEPRECATED FUNCTIONS **** //
// ****************************** //
	/** write a description of the elements to cout 
	*/
	//void describe();
	//Elements();



// SECOND TIMES PUBLIC ??????
//public:


	/** gram matrix used to compare atoms. can be set with loadGramAtom().
			this matrix is initialized to the identity matrix by the constructor.
	*/
	float gramAtom[NUMELEMENTS][NUMELEMENTS];

	/** filename of the atom kernel matrix loaded.
 	*/
	string gramAtomName;



private :

	/** hash of atoms indexed by their symbolic names (ex: H, Na,...).
	*/
	map<string, Atom*> periodicTable;



	
};



#endif

