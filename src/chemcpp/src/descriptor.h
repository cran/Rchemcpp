/****************************************************************************************
					  descriptor.h 
					----------------
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



#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include<iostream>
#include<sstream>
using std::string;
using std::endl;
using std::stringstream;

#include <constant.h>
#include <cerror.h>

/** Template class for storing scientific values
  @author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
	@version 0.3
	@date 17 Jan 2004

	CLASS NAME: 	Descriptor

	FOR:					SNSF SPONSORED PROJECT

	PURPOSE:  		Storage type for a general scientific descriptor with label,
								value, unit and comment

	This template class implements the storage of a general scientific measurment
	It not only stores a value with its label, but also the units used, and a comment.
	
  the label, unit and comment are all strings. The value however can be of any type since
	its type is a tamplate T. Creating a new descriptor can be done in that way:

	(in this example we are interested in a descriptor for body temperature and we initialize
	its value to 37.5)

	Descriptor<float> myDescriptor = new Descriptor<float>("temperature", 37.5, "degrees C", "body temperature");

	If we don't want the descriptor to contain any value we can give it any value and call
	the setEmpty() function. After that any call to getValue() will throw an exception, until
	a value is assigned with setValue(T)

	Descriptors of type int, float and string are used by the DataContainer class.

  */

template <class T> class Descriptor {
public:
///@name Descriptor construction functions
//@{
	/** class constructor.
	*/
	Descriptor();

	/** class constructor.
	*/
	Descriptor(T avalue);

	/** class constructor.
	*/
  	Descriptor(string aLabel, T aValue, string aUnit, string aComment);

	//~Descriptor(){}
//@}


///@name Accessor functions
//@{
	/** sets the value of the descriptor.
	*/
	void setValue( T );
	
	/** sets the label of the descriptor.
	*/
	void setLabel( string aLabel ){ label = aLabel; }
	
	/** sets the unit of the descriptor.
	*/
	void setUnit( string aUnit ){ unit = aUnit; }
	
	/** sets the comment of the descriptor.
	*/
	void setComment( string aComment){ comment = aComment; }

	/** empties the value of the Descriptor. isEmpty() function will then return true.
	*/
  	void setEmpty(){ empty = true; }

	/** returns the value. Throws a CError exception if the descriptor has been setEmpty().
	*/
	T getValue( bool silentError = false ) throw( CError );
	
	/** returns the label string.
	*/
	string getLabel(){ return( label ); }
	
	/** returns the unit string.
	*/
	string getUnit(){ return( unit ); }
	
	/** returns the comment string.
	*/
	string getComment(){ return( comment ); }

	/** returns true is the value was not set, false otherwise.
	*/
	bool isEmpty(){ return(empty); }
//@}

///@name Output functions
//@{
	/** returns a string description of the descriptor.
	*/
	string toString();
	
	/** returns a string description of the descriptor omiting the comment.
	*/
	string toStringShort();

	/** prints a description of the descriptor to cout.
	*/
	void describe();

	/** prints a description of the descriptor to cout omiting the comment.
	*/
	void describeShort();
  
//@}

private:

	/** string representing the descriptor name, example: temperature.
	*/
	string label;

	/** variable representing the value of the descriptor, exemple: <float> 37.5.
	*/
	T value;

	/** string representing the descriptor unit, example: degree C.
	*/
	string unit;

	/** string representing a comment for the descriptor, example: body temperature.
	*/
	string comment;

	/** flag telling if the Descriptor is empty.
	*/
	bool empty;
 
};





template <class T>
Descriptor<T>::Descriptor(){
	setLabel( "no name" );
	setUnit( "no unit" );
	setComment( "no comment" );
	empty = true;
}

template <class T>
Descriptor<T>::Descriptor( T avalue ){
	//Descriptor();
	setLabel( "no name" );
	setUnit( "no unit" );
	setComment( "no comment" );
	empty = true;

	setValue( avalue );
}

template <class T>
Descriptor<T>::Descriptor( string aLabel, T aValue, string aUnit, string aComment ){
	setLabel( aLabel );
	setUnit( aUnit );
	setComment( aComment );
	setValue( aValue );
}


template <class T>
void Descriptor<T>::setValue( T avalue ){
	value = avalue;
	empty = false;
}

template <class T>
T Descriptor<T>::getValue( bool silentError ) throw( CError ){
	if( !isEmpty() ){
		return value;
	}else{
		// value is emtpy, so ERROR_NA exception is raised
		CError e( ERRORNA, getLabel() + " is empty" );
		if( silentError == false){
			e.describe();
		}
		throw( e );
	}
}

template <class T>
void Descriptor<T>::describe(){
	cout << toString() << endl;
}

template <class T>
void Descriptor<T>::describeShort(){
	cout << toStringShort() << endl;
}

template <class T>
string Descriptor<T>::toString(){
	stringstream out;
	if( empty == false ){
		out << getComment() << ": " << getLabel() << " = " << getValue() << " (" << getUnit() << ") ";
	}else{
		out << getComment() << ": " << getLabel() << " = " << "NA" << " (" << getUnit() << ") ";
	}
	return( out.str() );
}

template <class T>
string Descriptor<T>::toStringShort(){
	stringstream out;
	if( empty == false ){
		out << getLabel() << " = " << getValue() << " (" << getUnit() << ") ";
	}else{
		out << getLabel() << " = " << "NA" << " (" << getUnit() << ") ";
	}
	return( out.str() );
}

#endif
