/****************************************************************************************
					  node.cpp 
					------------
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


#include "node.h"



Node::Node() : DataContainer()
{
}

Node::Node( const DataContainer& aDataContainer ) : DataContainer( (DataContainer&) aDataContainer )
{
}

Node::~Node()
{
}

void Node::setLabel( string aString ){
	label = aString;
}


double Node::getKashimaPQ( bool silentError ){
	if( !flagHasPq ){
		stringstream out;
		out << "No value was defined for pq in atom " << toStringShort() << " ";
		CError e( ERRORNA, out.str() );
		if( !silentError ){
			e.describe();
		}
		throw(e);
	}
	return( pq );
}

double Node::getKashimaPS( bool silentError ){
	if( !flagHasPs ){
		stringstream out;
		out << "No value was defined for ps in atom " << toStringShort() << " ";
		CError e( ERRORNA, out.str() );
		if( !silentError ){
			e.describe();
		}
		throw(e);
	}
	return( ps );
}



string Node::toStringShort(){
	stringstream out;
	out << "Node " << getLabel() << endl;
	return( out.str() );
}



/*!
    \fn Node::resetFlagHasPQ()
    the node has no pq value
 */
//void Node::resetFlagHasPQ()
//{
//	flagHasPq = false;
//}


/*!
    \fn Node::setFlagHasPQ( bool )
 */
//void Node::setFlagHasPQ( bool aBool )
//{
//	flagHasPQ = aBool
//}


