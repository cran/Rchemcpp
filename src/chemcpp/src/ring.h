/****************************************************************************************
					  ring.h 
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


#ifndef RING_H
#define RING_H

#include <vector>

#include <cerror.h>


using namespace std;

class Atom;
class Bond;


/**
This class implements the notion of ring. It contains a list of atoms.

@author Jean-Luc Perret
*/
class Ring : public std::vector<Atom*>{
public:

///@name Ring construction functions
//@{
	/** class constructor.
	*/
	Ring();

	/** class constructor.
	*/
	Ring( vector<Atom*>*, vector<Bond*>* );

	/** class desctructor.
	*/
	~Ring();

	/** adds a bond to the ring.
		this function should only be called during ring creation from a file. Not for ring detection.
	*/
	void addBond( Bond* aBond, bool silentError = false ) throw ( CError );

	/** adds an atom to the ring.
	*/
	void addAtom( Atom* anAtom, bool silentError = false ) throw ( CError );

//@}

///@name Accessor functions
//@{

	/** returns the id of the ring.
	*/
	int getID(){ return( id ); }

	/** sest the id of the ring.
	*/
	void setID( int a ){ id = a; }

	/** returns the vector of bonds the ring is made of.
	*/
	vector<Bond*>* getBonds(){ return( &bonds ); }

//@}



///@name Ring manipulation functions
//@{
	/** returns true if the ring contains anAtom.
	*/
	bool hasAtom( Atom* anAtom );

	/** returns true if the ring contains aBond.
	*/
	bool hasBond( Bond* aBond );

	/** comparison operator.
	*/
	bool equals( Ring* anotherRing );

//@}


///@name Output functions
//@{

	/** returns a string shortly describing of the ring.
	*/
	string toStringShort();

	/** returns a string describing of the ring.
	*/
	string toString();

	/** gives a short description of the ring.
	*/
	void describeShort();

	/** gives a description of the ring.
	*/
	void describe();

//@}


// ****************************** //
// **** DEPRECATED FUNCTIONS **** //
// ****************************** //
	//Ring( vector<Atom*>* );
	//vector<Bond*>::iterator beginBond(){ return( bonds.begin() ) ; }
	//vector<Bond*>::iterator endBond(){ return( bonds.end() ) ; }

protected:

	/** vector of bonds the ring is made of.
	*/
	vector<Bond*> bonds;

	/** stores a unique id of the ring in the molecule.
	*/
	int id;

};

#endif
