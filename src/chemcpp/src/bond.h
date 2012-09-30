/****************************************************************************************
					  bond.h 
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



#ifndef BOND_H
#define BOND_H

//#include<sstream>

#include "constant.h"
#include "datacontainer.h"
#include "atom.h"

#include <ring.h>


class Atom;

/** Bond class
    @author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
		@version 0.3
		@date 17 Jan 2004

		CLASS NAME: 	Bond

		FOR:					SNSF SPONSORED PROJECT

		PURPOSE:

	  this class implements the notion of bond between two atoms
    A bond is directed: it has a source atom and a target atom,
    so a chemical bond should be represented by two instances of class Bond,
    one for each direction.
  
    A bond as implemented here contains a source and a target atom, a label
		describing the bond type (see constants.h).
		Since bond is a daughter of DataContainer other Descriptors may be added
		at runtime.

		A bond is derived from DataContainer, so it has all the associated
		descriptors, but in addition bonds inplement an atomDescriptors map
		allowing to store atoms as descriptors to the bond.
	  (in addition to source ant target atom).

		For faster access, some descriptors are directly implemented in the class bond
		using native data type (double and int). This is the case for pt and label.

  	*/



class Bond { //: public DataContainer  {  // bonds do not need to be a DataContainer

/**
		\example bond_example.cpp
*/


public:

///@name Bond construction functions
//@{
	/** class constructor.
	*/
	Bond();

	/** class constructor.
	*/
	Bond(Atom* aSource, Atom* aTarget, int aLabel, int aPerretLabel = NAVALUE, int aBondStereo = 0,  int aBondNotUsed = 0, int aBondTopology = 0, int aBondReactionCenter = 0 );

	/** class destructor.
	*/
	~Bond();

//@}



///@name Accessor functions
//@{

	/** returns a pointer to the target atom.
	*/
	Atom* getTarget() throw( CError );

	/** returns true if target atom is set.
	*/
	bool hasTarget(){ if( target == NULL ){ return( false ); }else{ return( true ); }}

	/** returns a pointer to the source atom.
	*/
	Atom* getSource() throw( CError );

	/** returns true if source atom is set.
	*/
	bool hasSource(){ if( source == NULL ){ return( false ); }else{ return( true ); }}

	/** returns the bond label.
	*/
	int getLabel(){ return(label);}

	/** returns the value of the Perret label (set by setPerretLabel, called
		everytime a molecule is modified).
	*/
	int getPerretLabel(){
		return( perretLabel );
	}

	/** sets the Perret label of the bond.
	*/
	void setPerretLabel();


	/** sets the transition probability used in the Kashima kernel.
		note that for performance reasons this function does not use the descriptors
		of DataContainer to store its values but rather uses a native double variable
		in the bond class.
	*/
	void setKashimaPT(double aPt){ pt = aPt; }

	/** returns the transition probability used in the Kashima kernel.
	*/
	double getKashimaPT(){ return( pt ); }

	/** returns true if the bond has been marked with setFlag().
		flag can be removed with unsetFlag().
	*/
	bool hasFlag(){ return flag; }

	/** sets 'flag' to true
	*/
	void setFlag(){ flag = true; }

	/** sets 'flag' to false
	*/
	void unsetFlag(){ flag = false; }

	/** returns true if the bond has been marked with setFlagOriginal().
		flag can be removed with unsetFlagOriginal().
	*/
	bool hasFlagOriginal(){ return flagOriginal; }
	
	/** sets 'flagOriginal' to true
	*/
	void setFlagOriginal(){ flagOriginal = true; }

	/** sets 'flagOriginal' to false
	*/
	void unsetFlagOriginal(){ flagOriginal = false; }


	/** adds a ring to the set of rings the bond is member of.
	*/
	void addRing( Ring* aRing ){ if( !hasRing( aRing ) ){ rings.push_back( aRing ); } }

	/** returns true if the bond is member of the given ring.
	*/
	bool hasRing( Ring* aRing );

	/** returns true if the bond is member of any ring.
	*/
	bool hasRing(){ if( numRings() > 0 ){ return( true ); }else{ return( false );};}
	
	/** returns the number of rings the bond is member of.
	*/
	int numRings(){ return( rings.size() ); }

	/** sets the stereo value.
	*/
	void setStereo( int a ){ stereo = a; }

	/** returns the stereo value.
	*/
	int getStereo(){ return( stereo ); }
	
	/** sets the topology value.
	*/
	void setTopology( int a ){ topology = a; }

	/** returns the topology value.
	*/
	int getTopology(){ return( topology ); }

	/** sets the reaction center value.
	*/
	void setReactionCenter( int a ){ reactionCenter = a; }

	/** returns the reaction center value.
	*/
	int getReactionCenter(){ return( reactionCenter ); }
	
	/** sets the 'not used' value.
	*/
	void setNotUsed( int a ){ notUsed = a; }

	/** returns the 'not used' value.
	*/
	int getNotUsed(){ return( notUsed ); }

//@}


///@name Bond manipulation functions
//@{
	/** reverses the bond direction.
	*/
	void reverse();

	/** since bonds are represented directional, returns the reciprocal bond.
	*/
	Bond* getReverse() throw( CError );

	/** NO DOCUMENTATION
	*/
	void hideToFrom();
	
	/** NO DOCUMENTATION
	*/
	void restoreToFrom();

//@}



///@name Output functions
//@{

	/** returns a string description of the bond 
		( = string description of source and target atoms and bond label).
	*/
	string toString();

	/** returns a short string description of the bond
		( = short string description of source and target atoms and bond label).
	*/
	string toStringShort();

	/** writes a string description of the bond to cout.
	*/
	void describe();

	/** writes a short string description of the bond to cout.
	*/
	void describeShort();

//@}



///@name Iterators
//@{


	/** iterator on the first element of the set of rings the bond is member of.
	*/
	vector<Ring*>::iterator beginRing(){ return( rings.begin() );}
	/** iterator on the last element of the set of rings the bond is member of.
	*/
	vector<Ring*>::iterator endRing(){ return( rings.end() );}

//@}


	// ****************************** //
	// **** DEPRECATED FUNCTIONS **** //
	// ****************************** //

	//string toStringBFSVector();
	//Ring* getRingBFS( vector<Bond*>* toVisit, Bond* init = NULL );


private: // Private attributes

	/** Source atom of this bond.
	*/
	Atom *source;

	/** Target atom of this bond.
	*/
	Atom *target;

	/** Label of the bond. Possible values are defined in constants.h
	*/
	int label;

	/** Perret labeling.
	*/
	int perretLabel;

	/** stereo information. one of:
		if single bond: STEREONOT, STEREOUP, STEREOEITHER, STERODOWN
		if double bond: STEROCISTRANSNOT, STEROCISTRANS
	*/
	int stereo;

	/** not used field (xxx in bond mol files).
	*/
	int notUsed;

	/** topology (rrr in bond mol files).
		DO NOT USE, just for storing informations contained in mol files.
	*/
	int topology;

	/** reaction center (ccc in mol files).
	*/
	int reactionCenter;

	/** transition probability (used in Kashima Kernel).
	*/
	double pt;

	/** flag which can be set using setFlag() and unset using unsetFlag().
	*/
	bool flag;

	/** flag which can be set using setFlagOriginal() and unset using unsetFlagOriginal().
	this flag is set if the bond is present in the original molecule file (backward bonds don't).
	*/
	bool flagOriginal;	
				
	/** atom Descriptors. Intended to store source and target atom of the bond.
			may include other descriptors via the addAtomDescriptor() function.
	*/
	//map< string, Descriptor< Atom >* > atomDescriptors;

	/** this vector is used to store paths in a BFS search
	*/
	//vector<Bond*> BFSVector;

	/** vector of rings the bond is member of.
	*/
	vector<Ring*> rings;


};

#endif
