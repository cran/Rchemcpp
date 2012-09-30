/****************************************************************************************
					  node.h 
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


#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <datacontainer.h>


/**
@author Jean-Luc Perret
*/
class Node : public DataContainer{
public:

///@name Node construction functions
//@{

	/** class constructor.
	*/
	Node();
	
	/** class constructor.
	*/
	Node( const DataContainer& aDataContainer );
	
	/** class destructor.
	*/
	~Node();

//@}

///@name Accessor functions
//@{

	/** returns the Label of the node.
	*/
	string getLabel(){ return label; }

	/** sets the Label of the node.
	*/
	void setLabel( string aString );

	/** sets the Start probability of the atom to aPs (used in the Kashima Kernel).
	*/
	void setKashimaPS(double aPs){ ps = aPs; flagHasPs = true; }

	/** returns the Start probability of the atom (used in the Kashima Kernel).
	*/
	double getKashimaPS( bool silentError = false );
	//{ return ( getFloatDescriptor(PS)->getValue( silentError ) ); }

	/** sets the Stop probability of the atom to aPq (used in the Kashima Kernel).
	*/
	void setKashimaPQ( double aPq ){ pq = aPq; flagHasPq = true; }

	/** returns the Stop probability of the atom ( used in the Kashima Kernel ).
	*/
	double getKashimaPQ( bool silentError = false );

	/** helper function for the calculation of the Kashima Kernel.
	*/
	int getRPosition(){ return( rPosition ); }

	/** helper function for the calculation of the Kashima Kernel.
	*/
	void setRPosition( int aR ){ rPosition = aR; }

//@}



///@name Output functions
//@{
	/** gives a short (string) description of the node.
	*/
	string toStringShort();
//@}



// ****************************** //
// **** DEPRECATED FUNCTIONS **** //
// ****************************** //
	//void resetFlagHasPQ();
	//void setFlagHasPQ( bool );
	//bool hasPQ(){ return flagHasPq; }


protected:

	/** Label of the node.
	*/
	string label;

	/** stop probability for random walks (graph kernel).
	*/
	double pq;

	/** flag indicating wether the stop probability variable is initialized or not.
	*/
	bool flagHasPq;

	/** start probability for random walks (graph kernel).
	*/
	double ps;

	/** flag indicating wether the start probability variable is initialized or not.
	*/
	bool flagHasPs;

	/** helper variable for the calculation of the Kashima Kernel.
	*/
	int rPosition;

private:

};

#endif
