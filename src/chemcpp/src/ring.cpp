/****************************************************************************************
					  ring.cpp 
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


#include "ring.h"

#include<string>
#include<iostream>
#include<sstream>

#include <constant.h>

Ring::Ring( vector<Atom*>* atomsList, vector<Bond*>* bondsList ){
	vector<Atom*>::iterator ai;
	for( ai = atomsList->begin(); ai != atomsList->end(); ai++ ){
		push_back( *ai );
	}

	vector<Bond*>::iterator bi;
	for( bi = bondsList->begin(); bi != bondsList->end(); bi++ ){
		bonds.push_back( *bi );
	}
	#ifdef DEBUG
	cout << "a7h     created ring with " << size() << " atoms and " << bonds.size() << " bonds" << endl;
	#endif
}

Ring::Ring(){
}

Ring::~Ring()
{
}

bool Ring::equals( Ring* anotherRing ){
	// two rings are equal if they are the same size and contain the same atoms
	if( size() != anotherRing->size() ){
		return( false );
	}
	// they are the same size so let's check if they contain the same atoms
	vector<Atom*>::iterator ai;
	for( ai = begin(); ai != end(); ai++ ){
		if( anotherRing->hasAtom( (*ai) ) ){
		}else{
			return( false );
		}
	}
	return( true );
}

bool Ring::hasAtom( Atom* anAtom ){
	vector< Atom* >::iterator ai;
	for( ai = begin(); ai != end(); ai++ ){
		if( (*ai) == anAtom ){
			return( true );
		}
	}
	return( false );
}

bool Ring::hasBond( Bond* aBond ){
	vector< Bond* >::iterator ai;
	for( ai = bonds.begin(); ai != bonds.end(); ai++ ){
		if( (*ai) == aBond ){
			return( true );
		}
	}
	return( false );
}


string Ring::toStringShort(){
	stringstream out;
	out << "Ring of size " << size();
	return( out.str() );
}

string Ring::toString(){
	stringstream out;
	out << toStringShort() << " whith " << bonds.size() << " bonds";
	return( out.str() );
}

void Ring::describeShort(){
	cout << toStringShort() << endl;
}

void Ring::describe(){
	cout << toString() << endl;
}

void Ring::addBond( Bond* aBond, bool silentError ) throw( CError ){
	if( !hasBond( aBond ) ){
	// add the bond to the vector of bonds
		bonds.push_back( aBond );
	}else{
		if( !silentError ){
			CError e = CError( BONDALREADYEXISTS, "Ring::addBond: cannot add bond because bond already exist in ring" );
			throw( e );
		}
	}
}

void Ring::addAtom( Atom* anAtom, bool silentError ) throw( CError ){
	if( !hasAtom( anAtom ) ){
	// add the bond to the vector of bonds
		push_back( anAtom );
	}else{
		if( !silentError ){
			CError e = CError( ATOMALREADYEXISTS, "cannot add atom Ring::addAtom:  because atom already exist in ring" );
			throw( e );
		}
	}
}
