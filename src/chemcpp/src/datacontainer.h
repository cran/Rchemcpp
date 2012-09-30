
/****************************************************************************************
					  datacontainer.h 
					-------------------
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

#ifndef DATACONTAINER_H
#define DATACONTAINER_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

using std::map;
using std::ostream;


#include "descriptor.h"
#include "cerror.h"

/** Container for scientific values of type int, float and string
	@author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
	@version 0.1
	@date 8 Jan 2004

	CLASS NAME: 	DataContainer

	FOR:					SNSF SPONSORED PROJECT

	PURPOSE:  		Base Class implementing a collection of Descriptors and a
								collection of kind Descriptors ( of types int, float, and string )

	DataContainer is a base class to be used to derive classs representing concepts
	which need to be described with real values (properties)

	We distinguish two kind of properties: unique descriptors and kind descriptors.
	Unique descriptors are properties whose values are unique to an instance of
	the concept.
	Kind descriptors are properties whose values are shared among a group of concept
	instances.

	For example the class Atom is derived from DataContainer. Atoms have unique
	descriptors like a name, a start and stop probability. On the other hand atoms
	have many kind descriptors which are shared among atoms belonging to the
	same chemical element. The atomic number of all hydrogens is the same, but
	is different for carbons. Kind descriptors allow for this.

	Kind desriptors are stored only once in memory and are then referenced by all
	class instances of the same kind can therefore be added at any time but should
	be removed with caution since many class instances may point to it. In chemcpp
	generic instances of Atoms are instanciated globally in the elements.h file by
	creation of an Elements instance and stored in the elements map. Kind descriptors
	are created at that time and will only be deleted with the deletion of the unique
	instance of Elements.

	\todo deletion of kindDescriptors by elements is not yet implemented

	Memory allocation for unique descriptors on the other hand can be managed by the
	DataContainer class itself since only one instance of the derived class points to
	them. So deletion of the DataContainer causes the deletion of the descriptors.

	\warning unique descriptors should not be referenced by something else than this
	container class.

  */


class DataContainer {
public:

///@name DataContainer construction functions
//@{ 
	/** class constructor.
	*/
	DataContainer();

	/** class constructor.
	*/
	DataContainer( DataContainer& aDataContainer );
  	//DataContainer& operator= (const DataContainer& aDataContainer);
	
	/** (virtual) class desctructor.
	*/
	virtual ~DataContainer();
//@}


///@name Accessor functions
//@{

	/** returns true if the datacontainer has a string descriptor with label aLabel.
	*/
	bool hasStringDescriptor( string aLabel );

	/** returns true if the datacontainer has a float descriptor with label aLabel. 
	*/
	bool hasFloatDescriptor( string aLabel );
	
	/** returns true if the datacontainer has an int descriptor with label aLabel. 
	*/
	bool hasIntDescriptor( string aLabel );

//@}


///@name DataContainer manipulation functions
//@{

  	/** adds an integer kind descriptor with a label, value, unit and comment to the int kind descriptor container.
			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< int >* addKindIntDescriptor( string aLabel, int aValue, string aUnit, string aComment );

	/** adds an existing int descriptor to the int kind descriptor container.
			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< int >* addKindIntDescriptor(Descriptor< int>* aDescriptor);

  	/** adds an float kind descriptor with a label, value, unit and comment to the float kind descriptor container.
			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< float >* addKindFloatDescriptor( string aLabel, float aValue, string aUnit, string aComment );

	/** adds an existing float descriptor to the float kind descriptor container.
			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< float >* addKindFloatDescriptor(Descriptor< float>* aDescriptor);

	/** adds a string kind descriptor with a label, value, unit and comment to the string kind descriptor container.
 			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< string >* addKindStringDescriptor( string aLabel, string aValue, string aUnit, string aComment );

	/** adds an existing string descriptor to the string kind descriptor container.
			WARNING this descriptor will be deleted with the destruction of the
			Elements class.
	*/
	Descriptor< string >* addKindStringDescriptor(Descriptor< string>* aDescriptor);

	/** adds an integer unique descriptor with a label, value, unit and comment.
			this descriptor will be deleted with the destruction of this class.
	*/
	Descriptor< int >*  addIntDescriptor( string aLabel, int aValue, string aUnit, string aComment );

  	/** adds a float unique descriptor with a label, value, unit and comment.
			this descriptor will be deleted with the destruction of this class.
	*/
	Descriptor< float >*  addFloatDescriptor( string aLabel, float aValue, string aUnit, string aComment );

	/** adds a string unique descriptor with a label, value, unit and comment.
			this descriptor will be deleted with the destruction of this class.
	*/
	Descriptor< string >*  addStringDescriptor( string aLabel, string aValue, string aUnit, string aComment );

	/** sets the value of an int label (unique or kind descriptor).
	*/
	Descriptor< int >* setIntDescriptor(string aLabel, int aValue, string aUnit, string aComment, bool addIfMissing, bool silentError );

	/** sets the value of float label (unique or kind descriptor).
	*/
	Descriptor< float >* setFloatDescriptor(string aLabel, float aValue, string aUnit, string aComment, bool addIfMissing, bool silentError );

	/** sets the value of a string label (unique or kind descriptor).
	*/
	Descriptor< string >* setStringDescriptor(string aLabel, string aValue, string aUnit, string aComment, bool addIfMissing, bool silentError);


	/** returns a int property (unique or kind descriptor).
	*/
	Descriptor< int >* getIntDescriptor( string aLabel, bool silentError = true ) throw( CError );


	/** fills a vector with the possible values of a int descriptor and returns the number of possible values.
	*/
	long getPossibleValuesInIntDescriptor( string aDescriptorName, vector< int >* );

	/** returns a float property (unique or kind descriptor).
	*/
	Descriptor< float >* getFloatDescriptor( string aLabel, bool silentError = true  ) throw( CError );

	/** returns a string property (unique or kind descriptor).
	*/
	Descriptor< string >* getStringDescriptor( string aLabel, bool silentError = true  ) throw( CError );

	/** deletes the unique descriptor aString from Descriptors.
			deletes the object and removes pointer from the Hash, does not remove kind descriptors.
	*/
	virtual bool deleteDescriptor( string aString, bool found = false );

	/** deletes all unique descriptors from ***Descriptors.
			deletes the objects and removes pointer from the Hash, does not remove kind descriptors.
	*/
	void deleteAllDescriptors();

	  /** adds a descriptor to the datacontainer. 
		If aName terminates with .string, a string descriptor will be added. 
		If the aName terminates with a .integer, an integer descriptor will be added. 
		If aName terminates with .float, a float descriptor will be added. 
	*/
	void addUnknownTypeDescriptor( string aName, string aValue );

//@}



///@name Iterators
//@{

	/** start iterator for int descriptors.
	*/
	map< string, Descriptor< int >* >::iterator beginIntDescriptor(){ return( intDescriptors.begin() ); }
	/** end iterator for int descriptors.
	*/
	map< string, Descriptor< int >* >::iterator endIntDescriptor(){ return( intDescriptors.end() ); }

	/** start iterator for float descriptors.
	*/
	map< string, Descriptor< float >* >::iterator beginFloatDescriptor(){ return( floatDescriptors.begin() ); }
	/** end iterator for float descriptors.
	*/
	map< string, Descriptor< float >* >::iterator endFloatDescriptor(){ return( floatDescriptors.end() ); }

	/** start iterator for string descriptors.
	*/
	map< string, Descriptor< string >* >::iterator beginStringDescriptor(){ return( stringDescriptors.begin() ); }
	/** end iterator for string descriptors.
	*/
	map< string, Descriptor< string >* >::iterator endStringDescriptor(){ return( stringDescriptors.end() ); }

//@}




///@name Output functions
//@{
	/** writes a description of all labels to cout (unique and kind descriptors).
	*/
	void describe() throw( CError );

	/** writes a description of unique descriptors only to cout (no kind descriptors).
	*/
	void describeShort() throw( CError );
//@}

protected:

	/** false if this datacontainer was constructed by copy of an existingdatacontainer. Useful for destruction.
	*/
	bool flagElement;

	/** called by the destructor if flagElement == true.
	*/
	void deleteAllKindDescriptors();


	// Data container provides two different Descriptors:
	//   - kind descriptors which are shared by a group of instances of DataContainers
	//   - descriptors which are unique to an instance of this class.

	//   Kind descriptors can be created on the heap by this class, but are never deleted.
	// 	 on the other hand descriptors are fully managed by this class. WARNING:
	//	 Descriptors are deleted by the destructor of this class so they
	//   should not be used by reference from elsewhere.


	// DESCRIPTORS
  	// Hash of pointers on Descriptor<> to store properties of an instance of
	// this class
	// Objects in the hash are created by the add***Descriptor function
	// Objects are deleted in the deleteDescriptor or the
	// deleteAllDescriptors functions (also called by the destructor)

	/** hash of pointers on Integer descriptors.
	*/       
	map< string, Descriptor< int >* > intDescriptors;

	/** hash of pointers on Float descriptors.
	*/     
	map< string, Descriptor< float >* > floatDescriptors;
       
	/** hash of pointers on String descriptors.
	*/ 
	map< string, Descriptor< string >* > stringDescriptors;


	// KIND DESCRIPTORS
	// Hashes of pointers on Descriptor<> to store properties shared by a group of
	// instances of this class. The values of these descriptors are shared among the group
	// of instances. Example: the kind descriptor symbol with value H
	// is shared among all hydrogen instances of the derived class atom.
  	//
	// Objects in these maps can be created by the function
	// addKind***Descriptor( string aLabel, float aValue, string aUnit, string aComment )
	// which returns a pointer to the created descriptor.
	// An already created descriptor can also be included in the map by the function
	// addKind***Descriptor( descriptor& )
	
	/** hash of pointers on Integer kind descriptors.
	*/  	
	map< string, Descriptor< int >* >* kindIntDescriptors;

	/** hash of pointers on Float kind descriptors.
	*/  
	map< string, Descriptor< float >* >* kindFloatDescriptors;

	/** hash of pointers on String kind descriptors.
	*/  
	map< string, Descriptor< string >* >* kindStringDescriptors; // Hash of Desriptor<string>

  	// Hashes of pointers on Descriptor<string>
	// use when datacontainer Class instances should only contain a reference to a Descriptor
  	// map< string, Descriptor< int >* > intPDescriptors;
  	// map< string, Descriptor< float >* > floatPDescriptors;
  	// map< string, Descriptor< string >* > stringPDescriptors;

};
#endif
