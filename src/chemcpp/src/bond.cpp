/****************************************************************************************
					  bond.cpp 
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


#include "bond.h"

Bond::Bond(){
	pt = 0;
	label = NAVALUE;
	perretLabel = NAVALUE;
	stereo = 0;
	notUsed = 0;
	topology = 0;
	reactionCenter = 0;
}


Bond::Bond(Atom* aSource, Atom* aTarget, int aLabel, int aPerretLabel, int aBondStereo,  int aBondNotUsed, int aBondTopology, int aBondReactionCenter ){
	source = aSource;
	target = aTarget;
	label = aLabel;
	perretLabel = aPerretLabel;
	pt = 0;

	stereo = aBondStereo;
	notUsed = aBondNotUsed;
	topology = aBondTopology;
	reactionCenter = aBondReactionCenter;


	//addFloatDescriptor("pt" ,0.0, "probability", "transition probability for calculation of the Kashima Kernel");
}


Bond::~Bond(){
	//cout << "DELETING BOND" << endl;
	// delete atom descriptors
	/*map<const string, Descriptor<Atom>* >::iterator iti;
	for( iti = atomDescriptors.begin(); iti != atomDescriptors.end(); iti++ ){
			#ifdef DEBUG
				cout << "deleting atom descriptor [ " << (*iti).second << "]" << endl;
			#endif
			delete (*iti).second;
  		(atomDescriptors).erase(iti);
	}*/
}


/** reverses the bond direction */
void Bond::reverse(){
	Atom* toto;
	toto = target;
	target = source;
	source = toto;
}


/** returns a string description of the bond */
string Bond::toString(){
	stringstream out;
  out << getSource()->toStringShort() << " -" << getLabel() << "/" << getPerretLabel() << "- " << getTarget()->toStringShort();

  	if( rings.size() > 0 ){
		out << " in " << rings.size() << " ring of size";
		vector<Ring*>::iterator ri;
		for( ri = rings.begin(); ri != rings.end(); ri++ ){
			out << " " << (*ri)->size();
		}
	}else{
		out << " in no ring";
	}
	return( out.str() );
}



string Bond::toStringShort(){
	stringstream out;
  out << getSource()->toStringShort() << " -" << getLabel() << "- " << getTarget()->toStringShort();
	return( out.str() );
}



void Bond::describe(){
	cout << "    " << this->toString() << endl;
	//DataContainer::describeShort();
}



void Bond::describeShort(){
	cout << "    " << this->toStringShort() << endl;
}



Bond* Bond::getReverse() throw( CError ){
	try{
		Bond* result = getTarget()->getBondWithTarget( getSource() );
		return( result );
	}catch( CError e ){
		throw( e );
	}
}



bool Bond::hasRing( Ring* aRing ){
	vector<Ring*>::iterator ri;
	for( ri = rings.begin(); ri != rings.end(); ri++ ){
		if( (*ri) == aRing ){
			return( true );
		}
	}
	return( false );
}



void Bond::setPerretLabel(){
	perretLabel = label;
	if( hasRing() ){
		// only labels for bonds member of a ring a changed
		if( label != AROMATICBOND ){
			// only non aromatic cycles need renaming
			switch( label ){
				case SINGLEBOND:
					perretLabel = SINGLECYCLEBOND;
					break;
				case DOUBLEBOND:
					perretLabel = DOUBLECYCLEBOND;
					break;
				case TRIPLEBOND:
					perretLabel = TRIPLECYCLEBOND;
					break;
			}
		}
	}
}



Atom* Bond::getTarget() throw( CError ){
	if( !hasTarget() ){
		stringstream out;
		out << "Bond::getTarget: Bond " << toStringShort() << " has no target";
		CError e( ATOMNOTFOUND , out.str() );
		e.describe();
		throw(e);
	}else{
		return( target );
	}
}



Atom* Bond::getSource() throw( CError ){
	if( !hasSource() ){
		stringstream out;
		out << "Bond::getSource: Bond " << toStringShort() << " has no source";
		CError e( ATOMNOTFOUND , out.str() );
		e.describe();
		throw(e);
	}else{
		return( source );
	}
}



void Bond::hideToFrom(){
	getSource()->hideBond( getTarget() );
	getTarget()->hideBond( getSource() );
}



void Bond::restoreToFrom(){
	getSource()->restoreHiddenBond( getTarget() );
	getTarget()->restoreHiddenBond( getSource() );
}





/** Adds an atom descriptor intended to store stoms, for example source ant target are automatically added to a bond */
//void Bond::addAtomDescriptor(string aLabel, Atom* aValue, string aUnit, string aComment){
//	atomDescriptors[aLabel] = new Descriptor<Atom>( aLabel, *aValue, aUnit, aComment );
//}


/*Descriptor< Atom >* Bond::getAtomDescriptor(string aLabel){

	//bool found;

	// look among atom descriptors
  map<const string, Descriptor<Atom>* >::iterator its;
	for( its = atomDescriptors.begin(); its != atomDescriptors.end(); its++ ){
   	if( (*its).first == aLabel ){
			//found = true;
			return ( (*its).second );
		}
	}

  // if no corresponding descriptors found, throw exception
	CError e( MISSINGDESCRIPTOR, "no Atom descriptor " + aLabel );
	e.describe();
	throw(e);
}*/


/*bool Bond::deleteDescriptor( string aLabel, bool found ){
	// first check for the descriptor among the atom descriptors
	map<const string, Descriptor<Atom>* >::iterator iti;
	for( iti = (atomDescriptors).begin(); iti != (atomDescriptors).end(); iti++ ){
    if( (*iti).first == aLabel ){
			#ifdef DEBUG
				cout << "deleting Atom descriptor ";
			#endif
			(*iti).second->describe();
			delete (*iti).second;
   		(atomDescriptors).erase(iti);
			found = true;
		}
	}

	// then look in the remaining descriptors
  found = DataContainer::deleteDescriptor( aLabel, found );

	

	return(found);	
}*/

/** returns 1 if two bonds have the same label, 0 otherwise */
/*float Bond::binKernel(Bond* aBond){
	if( aBond->getLabel() == this->getLabel() ) {
		return 1;
	}else{
		return 0;
	}
}  */


/*string Bond::toStringBFSVector(){
	stringstream out;

	vector<Bond*>::iterator i;
	for( i = BFSVector.begin(); i!= BFSVector.end(); i++){
		out << (*i)->toStringShort() << ", ";
	}
	return( out.str() );
}*/


/*Ring* Bond::getRingBFS( vector<Bond*>* toVisit, Bond* init ){
	// set BFSVector for this atom
	//if( (*toVisit->begin()) == this ){
	toVisit->erase( toVisit->begin() );
	//}

	// update BFSPath
	//
	vector<Bond*>* ppath = getSource()->getBFSVector();
	vector<Bond*>::iterator ai;
	for( ai = ppath->begin(); ai != ppath->end(); ai++ ){
		BFSVector.push_back( (*ai) );
	}
	BFSVector.push_back( this );

	cout << "a4h    IN GETRINGBFS for atom " << toStringShort() << endl;
	cout << "a4h    with BFSVector of size " << getBFSVectorSize() << ": " << toStringBFSVector() << endl;

	map<Atom*, Bond*>::iterator i;
	cout << "a4h    checking target neighbour bonds" << endl;
	for( i = getTarget()->bonds.begin(); i != getTarget()->bonds.end(); i++ ){
	// for all neighbours except where we just come from check if they were already visited
		cout << "a4h       checking " << (*i).second->toStringShort() << endl;

		if( (*i).first != getSource() ){   // if the bond is not a bond returning were we just come from
			cout << "a4h       is not a return bond" << endl;
			if( (*i).second->getBFSVectorSize() > 0){ // if the bond was  already visited
				cout << "a4h    " << (*i).first->toStringShort() << " was already visited" << endl;
				// check if the ring is valid
				vector<Bond*> vectorIntersect;
				getVectorIntersect( getBFSVector(), (*i).second->getBFSVector(), &vectorIntersect );
				cout << "a4h    ring has intersect of size " << vectorIntersect.size() << endl;
				if( vectorIntersect.size() == 1 ){
					// if yes create the ring and return a pointer to it
					vector<Bond*> vectorUnion;
					getVectorUnion( getBFSVector(), (*i).second->getBFSVector(), &vectorUnion, &vectorIntersect );
					return( new Ring( &vectorUnion ) );

				}else{
					// if no abandon this branch
				}
			}else{
				// it the bond was not yet visited, add it to the list of bonds to visit and set BFSVector of target atom
				cout << "a4h    " << (*i).first->toStringShort() << " not yet visited, will continue this branch" << endl;

				for( ai = BFSVector.begin(); ai != BFSVector.end(); ai++ ){
					getTarget()->BFSVector.push_back( (*ai) );
				}

				cout << "a4h    target now has BFSVector" << endl;
				getTarget()->describe();

				// if no continue with this branch
				toVisit->push_back((*i).second);
			}
		}else{
			cout << "a4h       this bond returns to where we come from skipping" << endl;

		}
	}

	// call getRingBFS for the next bond in the queue
	if( toVisit->size() > 0 ){
		vector<Bond*>::iterator nextI = toVisit->begin();
		Bond* next = *nextI;
		toVisit->erase( nextI );
		Ring* ring = next->getRingBFS( toVisit, this );
	}else{
		// did not find a ring
		return( NULL );
	}

	// if no ring was found return NULL
	return( NULL );
}*/





