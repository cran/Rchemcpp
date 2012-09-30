/****************************************************************************************
					  atom.h 
					----------
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



#ifndef ATOM_H
#define ATOM_H

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>

#include <vector>

#include <node.h>
#include <bond.h>
#include <ring.h>

using namespace std;

/** Atom class
		@author Dr Jean-Luc Perret ("luc@kuicr.kyoto-u.ac.jp"), Kyoto University, Japan
		@version 0.3
		@date 17 Jan 2004


		CLASS NAME: 	Atom

		FOR:					SNSF SPONSORED PROJECT

		PURPOSE:

		This class implements the notion of atom. Atoms have
			- a map of bonds
			- properties, 2 kinds:
				- descriptors
				- kind descriptors
					    
		- Some descriptors are unique to an atom instance and are called descriptors.
		- Some descriptors are common to a chemical element and are called kind descriptors.
	  

		These descriptors are implemented using the DataContainer class which takes care of
		memory allocation for these descriptors. Therefore descriptors can be added and removed
		at runtime using the DataContainer functions.

		Some atom descriptors need fast access and are thus directly implemented in the
		atom class using native data types (int or float). These is the case for id and should
		be done for ps and pq.

	 	Atoms should always be created using the copy operator using the atoms stored in the
		global elements vector as template (\#include \<elements.h\>).

		for example instanciating a new hydrogen atom should be done in this way:

			Atom* myatom = new Atom( elements["H"] );

*/


class Atom : public Node {

friend class Elements;
friend class Bond;
friend class Molecule;

/**
		\example atom_example.cpp
*/

public :

///@name Atom construction functions
//@{

	/** class constructor.
	*/
	Atom();

	/** class constructor. 
	*/
	Atom( const Atom& anAtom );

	/** class constructor.
	*/
	Atom( string aLabel );
	
	/** class destructor.
	*/
	~Atom();

	/** assignment operator.
	*/
	 Atom& operator=( const Atom& anAtom );

	/** adds a bond to the bond hash.
	*/
	Bond* addBond( Bond *aBond ) throw ( CError );

//@}

///@name Accessor functions
//@{
	
	/** returns the atom unique Id.
	*/
	int getId() const { return( id ); }

	/** returns the atom unique Id as a string.
	*/
	string getIdString();

	/** returns the atom Id in the molecule.
	*/
	int getIdInMolecule() const { return( idInMolecule ); }

	/** sets the Id of this atom in the molecule. 
	*/
	void setIdInMolecule( int anId );
	
	/** returns the element type.
	*/
	int getType() const { return( type ); }

	/** returns the value of the name descriptor in case atoms where names
			(used for reading the Mutag dataset).
	*/
	string getName(){ return( getStringDescriptor( "name" )->getValue() ); }

	/** returns the value of the symbol descriptor in case atoms where names
			(used for reading the Mutag dataset).
	*/
	string getSymbol(){ return( getStringDescriptor( "Symbol" )->getValue() ); }

	/** returns the element symbol. In case of mol files it is the same as getSymbol()
	but in case of KCF files, getSymbol() returns the Kegg atom type while
	getElementSymbol() returns the element atom type. 
	*/
	string getElementSymbol(){ return( getStringDescriptor( "ElementSymbol" )->getValue() ); }
	
	/** sets the value of the ElementSymbol descriptor. 
	*/
	void setElementSymbol( string anElement );

	/** returns the atomic number of the atom (ex: hydrogens should return 1). 
	*/
	int getAN() const { return ( an ); }

	/** returns true if the atom is 'generic'. 
	*/
	bool isGeneric(){ return genericAtomType; }

	/** returns true if the atom is a non terminal Carbon
		(a terminal carbon is a carbon atom attached to only one
		non hydrogen atom).
	*/
	bool isCSkeleton();

	/** returns true if the atom has valid coordinates. 
	*/
	bool hasCoordinates();

	/** returns the x coordinate of the atom.
	*/
	float getX(){ return( x ); }

	/** returns the y coordinate of the atom. 
	*/
	float getY(){ return( y ); }

	/** returns the z coordinate of the atom. 
	*/
	float getZ(){ return( z ); }

	/** sets the x, y, and z coordinates of the atom.
	*/
	void setCoordinates( float aX, float aY, float aZ );

	/** returns the number of Bonds.
	*/
	int numBonds(){ return ( bonds.size() ); }

	/** returns the number of hidden Bonds.
	*/
	int numHiddenBonds(){ return ( hiddenBonds.size() ); }

	/** returns the number of aromatic Bonds for this atom.
	*/
	int getNumAromaticBonds();

	/** returns the number of Bonds (synonymous to numBonds).
	*/
	int degree(){ return ( bonds.size() ); }

	/** returns true if the node was already visited. 
	*/
	bool wasVisited(){ return visited;}

	/** sets the node as visited. 
	*/
	void setVisited(){ visited = true;}

	/** sets the node as unvisited. 
	*/
	void unsetVisited(){ visited = false;}

	/** returns the Morgan index of order anOrder. Emits an error if it was not calculated before.
	*/
	int getMorganIndex( int order );

	/** returns the value of the Morgan index for which the diversity of Morgan indices is maximum
			in the current molecule. set by setUniqueMorganIndex(int).
	*/
	int getUniqueMorganIndex( bool silentError = false ) throw( CError );

	/** set the Morgan label of the atom.
	*/
	void setMorganLabel(string aLabel);

	/** sets the Morgan label (used in the Molecule::morganKernel function) to the
	concatenation of the atomSymbol and the Morgan index of iteration anOrder.
	*/
	void setMorganLabel( int anOrder );

	/** returns the value of the Morgan label (set by setMorganLabel()).
	*/
	string getMorganLabel( bool silentError = false ) throw( CError );

	/** sets the Perret label of the atom (if a carbon has more than 2 aromatic bonds
	    it is renamed CJ).
	*/
	void setPerretLabel() throw( CError );

	/** returns the value of the perret label (set by setPerretLabel(), called
		everytime a molecule is modified).
	*/
	string getPerretLabel( bool silentError = false ) throw( CError );

	/** adds a ring to the set of rings the atom is member of.
	*/
	void addRing( Ring* aRing ){ if( !hasRing( aRing ) ){ rings.push_back( aRing ); } }

	/** returns true if the atom is member of the given ring.
	*/
	bool hasRing( Ring* aRing );

	/** returns true if the atom is member of any ring.
	*/
	bool hasRing(){ if( numRings() > 0 ){ return( true ); }else{ return( false );};}
	
	/** returns the number of rings the atom is member of.
	*/
	int numRings(){ return( rings.size() ); }

	/** returns the transition probability to anAtom.
	returns 0 if no atoms exist with anAtom.
	*/
	double getKashimaPT( Atom* anAtom );

	/** returns the BFS vector.
	*/
	vector<Atom*>* getBFSVector(){ return( &BFSVector ); }
	
	/** returns the BFS Bond vector.
	*/
	vector<Bond*>* getBFSBondVector(){ return( &BFSBondVector ); }

	/** returns the BFSVector size.
	*/
	int getBFSVectorSize(){ return( BFSVector.size() ); }
	
	/** resets the BFSVector.
	*/
	void resetBFSVector(){ BFSVector.clear(); BFSBondVector.clear(); }



//@}


///@name Atom manipulation functions
//@{
	/** erases the atoms coordinates (set flagHasCoordinates to false). 
	*/
	void eraseCoordinates();

	/** hides a bond (but not reverse bond).
	*/
	void hideBond( map< Atom*, Bond* >::iterator aBondI );

	/** hides a bond (but not reverse bond).
	*/
	void hideBond( Bond* aBond );
	
	/** hides a bond (but not reverse bond).
	*/
	Bond* hideBond( Atom* aTarget );

	/** hides all bonds whose target is an hydrogen.
	*/
	int hideHydrogenBonds();

	/** hides all bonds to and from this Atom.
	*/
	void hideAllToFromBonds();

	/** hides all bonds to and from to aTarget Atom
		and returns a poiter to the 'to Bond' of aTarget atom.
	*/
	Bond* hideToFromBonds( Atom* aTarget );

	/** hides the first bond in bonds (both to and from bonds)
		and returns a pointer to the 'to Bond'.
	*/
	Bond* hideToFromFirstBond();

	/** restores all hidden bonds (for example, bonds to hydrogens).
	*/
	int restoreHiddenBonds();

	/** restores hidden bond aBond.
	*/
	void restoreHiddenBond( Bond* aBond ) throw( CError );

	/** restores hidden bond to atom aTarget. Throws a CError if no hidden bond to aTarget is found.
	*/
	void restoreHiddenBond( Atom* aTarget ) throw( CError );

	/** deletes all bonds of that atom (called from the destructor of molecule).
	*/
	void deleteBonds();

	/** deletes all hidden bonds of that atom (called from the destructor of molecule).
	*/
	void deleteHiddenBonds();

	/** returns true if a bond (of any type) exists with the target atom.
	*/
	bool bondExists( Bond* aBond );

	/** returns the sum of the bond types.
	*/
	long bondSum();

	/** unsets all bond flags.
	*/
	void unsetBondFlags();
	
	/** unsets all bond flags original.
	*/
	void unsetBondFlagsOriginal();

	/** sets the unique Morgan index of this atom to the value of anOrder iteration
			of the Morgan index calculation.
	*/
	void setUniqueMorganIndex( int anOrder );

	/** returns the sum of the Morgan indices of order anOrder of all neighboors.
	*/
	int getSumOfNeighboorMorganIndex( int anOrder );

	/** resets the Morgan index map. call this function whenever the molecule is modified!
	*/
	void resetMorganIndex();

	/** set the value of the partial charge of the instance to aValue
	 */
	void setPartialCharge(double aValue);

	/** include the sign of the partial charge in the Morgan label of the instance.
	If partial charge > threshold: sign = + ; otherwise, sign = -
 	*/
	void setMorganChargeLabel(double threshold);


	/** returns the label of the bond connecting this atom with otherAtom. returns NOBOND if there are no bonds connecting both atoms. 
	*/
	Bond* getBondWithTarget( Atom* otherAtom ) throw( CError );

	/** returns a pointer to the next unvisited Node.
		returns NULL pointer if all atoms were visited.
	*/
	Atom* nextUnvisitedAtom() throw( CError );

	/** BFS procedure for smallest ring detection. returns a pointer
		to the smallest ring containing the atom or NULL if the atom is not member
		of any ring. note you should check ring membership with hasRing().
	*/
	Ring* getRingBFS( vector<Atom*>* toVisit, vector<Bond*>* toVisitBond ) throw( CError );

	/** NO DOCUMENTATION
	*/
	void pushBFSVector( vector<Atom*>* aPath, vector<Bond*>* aBondPath );
//@}


///@name Iterators
//@{

	/** returns a pointer to the bonds map.
	*/
	map<Atom*, Bond*>& getBonds(){ return( bonds ); }
	/** returns an iterator to the first Bond.
	*/
	map<Atom*, Bond*>::iterator beginBond(){ return( bonds.begin() ); }
	/** returns a pointer to the last Bond.
	*/
	map<Atom*, Bond*>::iterator endBond(){ return( bonds.end() ); }

	/** returns an iterator to the first ring the atom is member of.
	*/
	vector<Ring*>::iterator beginRing(){ return( rings.begin() ); }

	/** returns an iterator to the last ring the atom is member of.
	*/
	vector<Ring*>::iterator endRing(){ return( rings.end() ); }

	/** NOT DOCUMENTED
	*/
	map<Atom*, Bond*>::iterator getBondIteratorWithTarget( Atom* otherAtom );



//@}


///@name Output functions
//@{
	/** returns a string describing the atom (atomic symbol + unique Id + memory location + number of bonds).
	*/
	string toString();

	/** returns a short string describing the atom (atomic symbol + unique Id).
	*/
	string toStringShort();



	/** prints a description of the atom to cout (describeShort + descriptors, but not the kind descriptors).
	*/
	void describe() throw( CError );
	
	/** prints a short description of the atom to cout (atomic symbol + unique Id).
	*/
	void describeShort();
	
	/** NO DOCUMENTATION
	*/
	string toStringBFSVector();

//@}



	/** stream operator for Atom.
	*/
	friend ostream& operator<<(ostream& os, const Atom& anAtom);
	
	/** NO DOCUMENTATION.
	*/
	static void getVectorIntersect( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result);




// ******************************* //
// **** DEPRECATED FUNCTIONS **** //
// ****************************** //

	//Atom(const Atom& anAtom);

	//Atom(const Atom* anAtom) : DataContainer( (DataContainer&) anAtom );



	/** returns an iterator to the first Hiden Bond
	*/
	//map<Atom*, Bond*>::iterator beginHiddenBond(){ return( hiddenBonds.begin() ); }
	/** returns an iterator to the last Hiden Bond
	*/
	//map<Atom*, Bond*>::iterator endHiddenBond(){ return( hiddenBonds.end() )

	/** same as getBondWithTarget but among hidden bonds
	*/
	//Bond* getHiddenBondWithTarget( Atom* otherAtom );
	/** same as getBondWithTarget but among saved bonds
	*/
	//Bond* getSavedBondWithTarget( Atom* otherAtom );



	/** returns the transition probability to anAtom, but looks among hidden bonds
	  	Returns 0 if no atoms exist with anAtom
	*/
	//float getHiddenKashimaPT( Atom* anAtom );
	/** returns the transition probability to anAtom, but looks among saved bonds
	  	Returns 0 if no atoms exist with anAtom
	*/
	//float getSavedKashimaPT( Atom* anAtom );

	//static void getVectorUnion( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result, vector<Atom*>* intersect );


private :
	/** atom unique Id.
	*/
	int id;

	/** atom unique Id in the molecule it belongs.
	*/
	int idInMolecule;

	/** atom type: start from Hydrogen at 0, incremented at the instanciation of elements.
		used by MoleculeUtils::atomKernelExternalMatrix().
	*/
	int type;

	/** atomic number.
			to enhance retrival speed this property is hard coded
			rather than using the DataContainer mechanism.
	*/
	int an;


	/** atom x coodinate.
	*/
	float x;

	/** atom y coodinate.
	*/
	float y;

	/** atom z coodinate.
	*/
	float z;

	/** flag indicating if the atom has valid coordinates.
	*/
	bool flagHasCoordinates;

	/** counter of the number of atoms which have been instanciated.
	*/
	static int counter;

	/** hash of bonds indexed by atom.
	*/
	map<Atom*, Bond*> bonds;

	/** vector of hidden bonds indexed by atom.
	*/
	map<Atom*, Bond*> hiddenBonds;

	/** vector of saved bonds indexed by atom.
			filled by calling saveAllBonds().
	*/
	//map<Atom*, Bond*> savedBonds;


	/** smallest rings of which this atom is member
		(set by Molecule::detectSSSR()).
	*/
	vector<Ring*> rings;


	/** flag to indicate if the atom was visited. Used in DFS graph algorithm.
	*/
	bool visited;

	/** indicates if this atom is part of the elements or if it is generic (created with constructor: Atom::Atom(string)).
	*/
	bool genericAtomType;


protected : // Protected methods

	/** this function is used when instanciating the atom for setting the atomic number.
	*/
	void setAN( int a ) { an = a; }

	/** resets the atom counter. Only the elements should call this function.
	*/
	static void resetCounter() { counter = 0; }

	/** sets the atom type.
	*/
	void setType( int aType ) { type = aType; }

	/** clear hidden bonds.
	*/
	//void clearHiddenBonds();

	/** hide all bonds leaving this atom.
	*/
	//void hiddeAllBonds();

	/** save all bonds leaving this atom.
	*/
	//void saveAllBonds();

	/** map containing the morgan indices.
	*/
	map< int, int > morganIndex;

	/** Morgan label of the atom set by the Atom::setMorganLabel() function.
	*/
	string morganLabel;

	/** Perret label of the atom.set by the Atom::setPerretLabel() function.
	*/
	string perretLabel;

	/** Morgan index value for the iteration with the biggest diversity for a molecule
	(set by Molecule::setUniqueMorganIndices()).
	*/
	int uniqueMorganIndex;


	/** this vector is used to store paths in a BFS search
		(stores atoms).
	*/
	vector<Atom*> BFSVector;

	/** this vector is used to store paths in a BFS search
		(stores bonds).
	*/
	vector<Bond*> BFSBondVector;

	/** partial charge of the atom.
	*/
	double partialCharge;


};



#endif
