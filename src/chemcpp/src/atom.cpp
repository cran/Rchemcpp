/****************************************************************************************
					  atom.cpp 
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




#include "atom.h"

#include <moleculeutils.h>

#include <queue>

//#define DEBUG 1
//#define DEBUGPERRETLABEL 1

int Atom::counter = 0;

/** atom class constructor */
Atom::Atom() : Node() {
	#ifdef DEBUG
		cout << "Atom::Atom()" << endl;
	#endif
	Atom::counter++;
	id = Atom::counter;

	x = 0;
	y = 0;
	z = 0;
	flagHasCoordinates = false;

	an = 0;

	pq = 0;
	flagHasPq = false;
	ps = 0;
	flagHasPs = false;

	genericAtomType = false;

	#ifdef DEBUG
		cout << "Atom::Atom() DONE" << endl;
	#endif
}


/* the following function should be called to create an atom outside the elements
	it can be used to define a particular atom type. The only property therefore available
	from it is its label and its unique id.
	NOTE: atoms created this way are deleted when the molecule they are part of is deleted
*/
Atom::Atom( string aLabel ) : Node() {

	#ifdef DEBUG
		cout << "Atom::Atom( string aLabel )" << endl;
	#endif

	Atom::counter++;
	id = Atom::counter;

	x = 0;
	y = 0;
	z = 0;
	flagHasCoordinates = false;

	an = -1;

	pq = 0;
	flagHasPq = false;
	ps = 0;
	flagHasPs = false;

	label = aLabel;
	genericAtomType = true;

	#ifdef DEBUG
		cout << "Atom::Atom() DONE" << endl;
	#endif

}

Atom::Atom(const Atom& anAtom) : Node( (DataContainer&) anAtom ) { //: DataContainer( (DataContainer&) anAtom ){
//Atom::Atom(const Atom& anAtom) : DataContainer( anAtom ){

	#ifdef DEBUG
		cout << "Atom::Atom(const Atom* anAtom)" << endl;
	#endif

	Atom::counter++;
	id = Atom::counter;

	x = 0;
	y = 0;
	z = 0;
	flagHasCoordinates = false;

	this->setAN( anAtom.getAN() );
	this->setType( anAtom.getType() );
	this->setLabel( anAtom.label );   // MODIFIED JLP BUG

	#ifdef DEBUG
		cout << "created atom " << getSymbol() << getId() << endl;
	#endif

	ps = anAtom.ps;
	flagHasPs = anAtom.flagHasPs;
	pq = anAtom.pq;
	flagHasPq = anAtom.flagHasPq;

	//copy bonds
	//map<Atom*, Bond*>::iterator bit;
	//for ( bit = anAtom.bonds.begin(); bit != anAtom.bonds.end(); bit++ ){
		//Bond* newBond = new Bond( *(*bit) );
		//Atom* newAtom = new Atom( *elements["H"] );
		//atoms.push_back( newAtom );
		//this->addBond( newBond );
	//}

	genericAtomType = false;

	resetMorganIndex();
	#ifdef DEBUG
		//cout << "created atom " << toStringShort() << endl;
		cout << "Atom::Atom(const Atom* anAtom) DONE" << endl;
	#endif

}

Atom& Atom::operator=(const Atom& anAtom) {

  #ifdef DEBUG
		cout << "COPYING ATOM" << endl;
	#endif

	if(this != &anAtom){

		/*this->kindStringDescriptors = anAtom.kindStringDescriptors;
		this->kindIntDescriptors = anAtom.kindIntDescriptors;
		this->kindFloatDescriptors = anAtom.kindFloatDescriptors;*/

		// Copy all non kind descriptors
		// first copy int descriptors

		//cout << "copying An: " << anAtom.getAN() << endl;

		this->setAN( anAtom.getAN() );
		this->setType( anAtom.getType() );

		x = anAtom.x;
		y = anAtom.y;
		z = anAtom.z;
		flagHasCoordinates = anAtom.flagHasCoordinates;

		ps = anAtom.ps;
		flagHasPs = anAtom.flagHasPs;

		pq = anAtom.pq;
  		flagHasPq = anAtom.flagHasPq;

		// increment the atom id and set the id for this atom
		Atom::counter++;
		id = Atom::counter;
	}

	resetMorganIndex();

	return( *this );
}

void Atom::describe() throw( CError ) {
	cout << getLabel() << endl;
	cout << getMorganLabel() << endl;
	cout << getPerretLabel() << endl;
	cout << numBonds() << " bonds " << endl;
	cout << numHiddenBonds() << " hidden bonds" << endl;

	DataContainer::describe();
	cout << toString() << endl;
	cout << "coordinates: " << x << ", " << y << ", " << z << endl;
	cout << "An = " << getAN() << endl;
	cout << "Type = " << getType() << endl;
}
void Atom::describeShort() {
	cout << this->toString() << endl;
}

Bond* Atom::addBond(Bond *aBond) throw ( CError ){
	#ifdef DEBUG
		cout << "adding bond [ " << aBond->toStringShort() << "]" << endl;
	#endif
	if( bondExists( aBond ) ){
		//cout << "bond exists" << endl;

		stringstream out;
		out << "Bond " << aBond->toString() << " already exists ";
		CError e( BONDALREADYEXISTS, out.str() );
		e.describe();
		throw(e);

		return( bonds[aBond->getTarget()] );
	}else{

		bonds[aBond->getTarget()] = aBond;
		return( aBond );
	}
}

Atom::~Atom(){
	// delete all bonds
	#ifdef DEBUG
		cout << "Atom::~Atom()" << endl;
	#endif

	deleteBonds();

	#ifdef DEBUG
		cout << "Atom::~Atom() DONE" << endl;
	#endif
}

/** returns a string describing the atom */
string Atom::toString(){
	stringstream out;

  	out << this->getStringDescriptor(SYMBOLNAME)->getValue() << this->getIdString() << " ";
	out << "PerretLabel: " << getPerretLabel() << " ";
	out << "Morgan: " << getMorganIndex(1) << " " << getMorganIndex(2) << " " << getMorganIndex(3) << " ";
	if( hasCoordinates() ){
		out << " (" << x << ", " << y << ", " << z << ") ";
	}else{
		//out << this->getStringDescriptor(SYMBOLNAME)->getValue() << this->getIdString() << " Morgan: " << getMorganIndex(1) << " " << getMorganIndex(2) << " " << getMorganIndex(3) << " " << " at " << this << " (" << numBonds() << " bonds)";
	}

	out << " at " << this << " (" << numBonds() << " bonds)";

	try{
		if( hasRing() ){
			out << ", in " << numRings() << " ring of size:";
			vector<Ring*>::iterator ri;
			for( ri = beginRing(); ri != endRing(); ri++ ){
				out << " " << (*ri)->size();
			}
			out << " ";
		}else{
			out << ", no rings ";
		}
	}catch( CError e ){
		out << " rings not detected " << endl;
	}


	return( out.str() );
}

/** returns a short string describing the atom */
string Atom::toStringShort(){
	stringstream out;
	out << this->getStringDescriptor(SYMBOLNAME)->getValue() << this->getIdString();
	return( out.str() );
}




/** convenience function to retrieve the id as a string */
string Atom::getIdString(){
	stringstream out;
  out << getId();
	return( out.str() );
}


ostream& operator << (ostream& os, const Atom& anAtom) {
	os << "ATOM " << anAtom.getId();
  return os;
}


/** hide all bonds whose target is an hydrogen */
int Atom::hideHydrogenBonds(){

  //deleteHiddenBonds();

  //cout << "entering ATOM::hideHydrogenBonds" << endl;
	int i = 0;
	map<Atom*, Bond*>::iterator bond;

	for( bond = beginBond(); bond != endBond(); ){
	  //cout << "checking " << (*bond).second->toString() << endl;
		if( (*bond).first->getAN() == 1 ){
		  //cout << "ERASING" << endl;
			//hiddenBonds.push_back( (*bond) );
			hiddenBonds[(*bond).first] = (*bond).second;
			i++;
			map<Atom*, Bond*>::iterator aBond = bond;
			bond++;
			bonds.erase((*aBond).first);
			if( bond != beginBond() ){
				bond--;
			}
			//cout << "DONE" << endl;
		}else{
		  bond++;
		}
	}

	//cout << "leaving Atom::hideHydrogenBonds" << endl;

	return(i);
}
/** restore all hidden bonds (for example to hydrogens) */
int Atom::restoreHiddenBonds(){

	#ifdef DEBUG
		cout << "tttty Atom::restoreHiddenBonds()" << endl;
	#endif

	int i = 0;
	map<Atom*, Bond*>::iterator bond;

	for( bond = hiddenBonds.begin(); bond != hiddenBonds.end(); bond++ ){
		// restore bonds
		bonds[(*bond).first] = (*bond).second;
		i++;
	}
	hiddenBonds.clear();
	//hiddenBonds.erase( hiddenBonds.begin(),  hiddenBonds.end() );

	//cout << "  there are now " << hiddenBonds.size() << " hidden bonds" << endl;
	//cout << "  there are now " << bonds.size() << " bonds" << endl;

	return(i);
}

void Atom::restoreHiddenBond( Bond* aBond ) throw( CError ){
	map<Atom*, Bond*>::iterator bi;
	for( bi = hiddenBonds.begin(); bi != hiddenBonds.end(); bi++ ){
		if( (*bi).second == aBond ){
			bonds[(*bi).first] = (*bi).second;
			hiddenBonds.erase( bi );
			return;
		}
	}
	stringstream out;
	out << "Bond " << aBond->toString() << " not found among hidden bounds";
	CError e( BONDNOTFOUND, out.str() );
	e.describe();
	throw(e);

}

void Atom::restoreHiddenBond( Atom* aTarget ) throw( CError ){

	Bond* aBond = hiddenBonds[aTarget];
	if( hiddenBonds.find( aTarget ) == hiddenBonds.end() ){
		stringstream out;
		out << "Bond " << aBond->toString() << " not found among hidden bounds";
		CError e( BONDNOTFOUND, out.str() );
		e.describe();
		throw(e);
	}else{
		#ifdef DEBUG
			cout << "Atom::restoreHiddenBond( Atom* ): " << aBond->toStringShort() << endl;
		#endif
		bonds[aBond->getTarget()] = aBond;
	}

}



/** sets the x, y, and z coordinate of the atom */
void Atom::setCoordinates( float aX, float aY, float aZ ){

	//cout << " Atom::setCoordinates: setting coordinates to : " << aX << ", " << aY << ", " << aZ << endl;

	x = aX;
	y = aY;
	z = aZ;

	flagHasCoordinates = true;
}

/** erases the atoms coordinates (set flagHasCoordinates to false) */
void Atom::eraseCoordinates(){
	flagHasCoordinates = false;
}
/** returns true if the atom has valid coordinates */
bool Atom::hasCoordinates(){
	return( flagHasCoordinates );
}


/** unset all bond flags */
void Atom::unsetBondFlags(){
	map<Atom*, Bond*>::iterator bond;
	for( bond = beginBond(); bond != endBond(); bond++ ){
		(*bond).second->unsetFlag();
	}
}
/** unset all bond flags */
void Atom::unsetBondFlagsOriginal(){
	map<Atom*, Bond*>::iterator bond;
	for( bond = beginBond(); bond != endBond(); bond++ ){
		(*bond).second->unsetFlagOriginal();
	}
}

bool Atom::bondExists(Bond* aBond){
	map< Atom*, Bond* >::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		if( (*i).first == aBond->getTarget() ){
			return( true );
		}
	}
	return( false );
}


/** sets the id of this atom in the molecule */
void Atom::setIdInMolecule( int anId ){
	idInMolecule = anId;
}
/** sets the value of the ElementSymbol descriptor */
void Atom::setElementSymbol( string anElement ){
	setStringDescriptor("ElementSymbol", anElement, "", "", true, true);
}

/** returns a pointer to the next unvisited Node
		throws an exception if all atoms were visited
*/
Atom* Atom::nextUnvisitedAtom() throw( CError ){
	map<Atom*, Bond*>::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		//if(!(*i).second->getSource()->wasVisited()){
		//	return( (*i).second->getSource() );
		//}
		//cout << " was " << (*i).first->toStringShort() << " visited ? ";
		if( !(*i).first->wasVisited() ){
			//cout << "no" << endl;
			return( (*i).first );
		}//else{
			//cout << "yes" << endl;
		//}
	}

	CError e( ALLNEIGHBOURSVISITED, "all neighbours to atom " + toStringShort() + " were already visited" );
	throw(e);
}

/** hides a bond (but not reverse bond)
*/
void Atom::hideBond( map< Atom*, Bond* >::iterator aBondI ){
	hiddenBonds[(*aBondI).first] = (*aBondI).second;
	bonds.erase( aBondI );
}

void Atom::hideBond( Bond* aBond ){
	hiddenBonds[ aBond->getTarget() ] = aBond;
	bonds.erase( aBond->getTarget() );
}

Bond* Atom::hideBond( Atom* aTarget ){
	Bond* result = bonds[ aTarget ];

	#ifdef DEBUG
		cout << "Atom::hideBond( Atom* ): hiding bond " << result->toStringShort() << endl;
	#endif


	hiddenBonds[ aTarget ] = result;
	bonds.erase( aTarget );
	return( result );
}

/** hides all bonds to and from this Atom
*/
void Atom::hideAllToFromBonds(){
  #ifdef DEBUG
  cout << "Atom::hideAllToFromBonds: hiding bonds from " << toStringShort() << " to " << bonds.size() << " atoms "<< endl;
  #endif

  int i = 0;
  map<Atom*, Bond*>::iterator bond;

	for( bond = beginBond(); bond != endBond(); ){
			i++;
			map<Atom*, Bond*>::iterator aBond = bond;
			bond++;
			hideToFromBonds( (*aBond).first );
			if( bond != beginBond() ){
				bond--;
			}
	}

  #ifdef DEBUG
	cout << "Atom::hideAllToFromBonds done" << endl;
  #endif
}

// hides the first bond in bonds (both to and from bonds)
//
Bond* Atom::hideToFromFirstBond(){
	return( hideToFromBonds( ( *bonds.begin() ).first ) );
}

/** hides all bonds to and from to aTarget Atom
*/
Bond* Atom::hideToFromBonds( Atom* aTarget ){
	// hide from bond
        #ifdef DEBUG
  		cout << "hiding bond from " << aTarget->toStringShort() << " to this" << endl;
        #endif
	aTarget->hideBond( this );

        // hide to bond
	#ifdef DEBUG
	cout << "hiding bond from this to " << aTarget->toStringShort() << endl;
	#endif

	Bond* result = hideBond( aTarget );
	#ifdef DEBUG
	cout << "hiding bond to done" << endl;
	#endif
	return( result );
}


/** returns the morgan index of order anOrder. Emits an error if it was not calculated before
*/
int Atom::getMorganIndex( int order ){

	//order = order - 1;

	if( order < 1 ){
		stringstream out;
	  out << "Atom::getMorganIndex Morgan index cannot be calculated for an order < 1 (here : " << order << ")";
		CError e(VALUENOTALLOWED, out.str() );
		e.describe();
		throw(e);
	}

	if( morganIndex.count( order ) == 1 ){
		// the morgan index was already calculated
		return( morganIndex[ order ] );

	}else{
		// morgan as not been computed yet

		int res = 0;
		if( order > 1 ){
			// the morgan index is the sum of the connectivity values of the neighboors

			res = getSumOfNeighboorMorganIndex( order - 1 );
			morganIndex[ order ] = res;
			return( res );

		}else{
			// the first order morgan index is just the number of neighboors

			res = numBonds();
			morganIndex[ 1 ] = res;
			return( res );

		}
	}
}

/** returns the sum of the morgan indices of order anOrder of all neighboors.
*/
int Atom::getSumOfNeighboorMorganIndex( int anOrder ){
	// for each neighboor add the morgan index of order anOrder
	int sum = 0;
	map<Atom*, Bond*>::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		sum += (*i).first->getMorganIndex( anOrder );
	}
	return( sum );
}

long Atom::bondSum(){
	// for each neighboor add the morgan index of order anOrder
	long sum = 0;
	map<Atom*, Bond*>::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		sum += (*i).second->getLabel();
	}
	return( sum );
}


/** reset the morgan index map. Call this function whenever the molecule is modified!
*/
void Atom::resetMorganIndex(){
	morganLabel = "";
	uniqueMorganIndex = -1;
	morganIndex.clear();
}

/** sets the morgan label (used in the Molecule::morganKernel function) to the concatenation of the atomSymbol and the morgan index of iteration anOrder
*/
void Atom::setMorganLabel( int anOrder ){

	stringstream out;
	if( anOrder > 0 ){
		out << getLabel() << getMorganIndex( anOrder );
	}else{
		out << getLabel();
	}
	morganLabel = out.str();

}

void Atom::setPerretLabel() throw( CError ){
	stringstream out;
	out << getLabel();

	try{
		int numR = numRings();
		#ifdef DEBUGPERRETLABEL
			cout << "Atom::setPerretLabel for atom: " << toStringShort() << endl;
			cout << "   rings: " << numR << endl;
		#endif
		if( numR > 1 ){
			if( numR > 2 ){
				out << "K";
			}else{
				out << "J";
			}
		}
		perretLabel = out.str();

		#ifdef DEBUGPERRETLABEL
			cout << "   -> " << perretLabel << endl;
		#endif

		//set the perretLabel for all bonds in this atom
		map<Atom*,Bond*>::iterator bi;
		for( bi = bonds.begin(); bi != bonds.end(); bi++ ){
			(*bi).second->setPerretLabel();
		}
	}catch( CError e ){
		CError e(NOTCALCULATED, "Atom::setPerretLabel: ring membership was not calculated, please use Molecule::setectSSSR before calling Molecule::setPerretLabel");
		e.describe();
		throw(e);
	}
}

int Atom::getNumAromaticBonds(){
	map<Atom*, Bond*>::iterator i;
	int numAromaticBonds = 0;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		if( (*i).second->getLabel() == AROMATICBOND ){
			numAromaticBonds++;
		}
	}
	return( numAromaticBonds );
}


/** sets the unique morgan index of this atom to the value of anOrder iteration of the Morgan index calculation */
void Atom::setUniqueMorganIndex( int anOrder ){
	uniqueMorganIndex = getMorganIndex( anOrder );
}

string Atom::getMorganLabel( bool silentError ) throw( CError ){
	/*if(morganLabel == ""){
		// morgan label was not set, emit error
		CError e(NOTCALCULATED, "Atom::getMorganLabel: morganLabel was not calculated, use Molecule::setMorganLabels( int )");
		if( silentError == false){
			e.describe();
		}
		throw(e);
	} */
	return( morganLabel );
}

string Atom::getPerretLabel( bool silentError ) throw( CError ){
	return( perretLabel );
}

int Atom::getUniqueMorganIndex( bool silentError ) throw( CError ){
	if(uniqueMorganIndex == -1){
		// unique morgan index was not set, emit error
		CError e(NOTCALCULATED, "Atom::getUniqueMorganIndex: uniqueMorganIndex was not calculated, use Molecule::setUniqueMorganIndices()");
		if( silentError == false){
			e.describe();
		}
		throw(e);
	}
	return( uniqueMorganIndex );	
}

/** returns the label of the bond connecting this atom with otherAtom. Returns NOBOND if there are no bonds connecting both atoms */
Bond* Atom::getBondWithTarget( Atom* otherAtom ) throw( CError ){
	if( bonds.find( otherAtom ) == bonds.end() ){
		CError e( BONDNOTFOUND, "No bond with target " + otherAtom->toStringShort() + " in atom " + toStringShort() );
		throw( e );
	}
	return( bonds[ otherAtom ] );
}
map<Atom*, Bond*>::iterator Atom::getBondIteratorWithTarget( Atom* otherAtom ){
	return( bonds.find( otherAtom ) );
}



/** returns the transition probability to anAtom.
		Returns 0 if no atoms exist with anAtom */
double Atom::getKashimaPT( Atom* anAtom ){
	try{
		Bond* ab = getBondWithTarget( anAtom );
		return( ab->getKashimaPT() );
	}catch( CError ){
		return( 0 );
	}
}

/** delete all bonds in of that atom (called from the destructor of molecule) */
void Atom::deleteBonds(){
	//cout << "Atom::deleteBonds - deleting bonds for atom " << toStringShort() << endl;
	map<Atom*, Bond*>::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		#ifdef DEBUG
			cout << "deleting bond " << (*i).second->toStringShort() << endl;
		#endif
	delete (*i).second; // delete bonds
	}
	bonds.clear();
}

void Atom::deleteHiddenBonds(){
	//cout << "Atom::deleteBonds - deleting bonds for atom " << toStringShort() << endl;
	map<Atom*, Bond*>::iterator i;
	for( i = hiddenBonds.begin(); i != hiddenBonds.end(); i++ ){
		#ifdef DEBUG
			cout << "deleting bond " << (*i).second->toStringShort() << endl;
		#endif
	delete (*i).second; // delete bonds
	}
	hiddenBonds.clear();
}


bool Atom::isCSkeleton()
{
	// check if this atom is a carbon
	int numNonH = 0;

	if( getSymbol() == "C" ){
		//cout << " got a carbon " << endl;
		// if true, check the neighbours if there are more than one non H than this is a
		// non terminal carbon
		map<Atom*, Bond*>::iterator i;
		for( i = bonds.begin(); i != bonds.end(); i++ ){
			//cout << "    linked to " << (*i).first->getSymbol() << endl;
			if( (*i).first->getSymbol() != "H" ){
				//cout << ".";
				numNonH++;
				//if( numNonH > 1){
				//	cout << "    NOT A SKELETON" << endl;
				//	return( true );
				//}
			}
		}
		if( numNonH > 1 ){
			//cout << "    SKELETON" << endl;
			return( true );
		}else{
			//cout << "    NOT SKELETON" << endl;
			return( false );

		}

	}else{
		return( false );
	}
}

bool Atom::hasRing( Ring* aRing ){
	vector<Ring*>::iterator ri;
	for( ri = rings.begin(); ri != rings.end(); ri++ ){
		if( (*ri) == aRing ){
			return( true );
		}
	}
	return( false );
}


Ring* Atom::getRingBFS( vector<Atom*>* toVisit, vector<Bond*>* toVisitBond ) throw( CError ){
	#ifdef DEBUG
	cout << "a5h     ATOM::GETRINGBFS now at " << toStringShort() << endl;
	#endif

	if( toVisit->size() == 0 ){
		#ifdef DEBUG
		cout << "a5h     ROOT" << endl;
		#endif
		toVisit = new vector<Atom*>;
		BFSVector.push_back( this );  // update the path for this
	}else{
		BFSVector.push_back( (*toVisit->begin()) );  // update the path for this atom
		toVisit->erase( toVisit->begin() ); // erase the bond we just used from toVisit
	}

	if( toVisitBond->size() == 0 ){
		toVisitBond = new vector<Bond*>;
	}else{
		BFSBondVector.push_back( (*toVisitBond->begin()) );  // update the path for this atom
		toVisitBond->erase( toVisitBond->begin() ); // erase the bond we just used from toVisit
	}

	#ifdef DEBUG
	cout << "a5h     with path " << toStringBFSVector() << endl;
	cout << "        (length = " << BFSVector.size() << ")" << endl;
	cout << "        (bond length = " << BFSBondVector.size() << ")" << endl;
	cout << "a5h     ---->" << endl;
	#endif

	// add neighbour bonds with unvisited atoms to toVisit
	map<Atom*, Bond*>::iterator i;
	for( i = bonds.begin(); i != bonds.end(); i++ ){
		#ifdef DEBUG
		cout << "a5h     checking " << (*i).second->toStringShort() << endl;
		cout << "        with target " << (*i).first->toStringShort() << endl;
		#endif
		//cout << "        cannot go to " << (*(BFSVector.end()-2))->toStringShort() << endl;

		if( BFSVector.size() > 1 && (*i).first == (*(BFSVector.end()-2)) ){
			#ifdef DEBUG
			cout << "        cannot go to " << (*(BFSVector.end()-2))->toStringShort() << endl;
			cout << "a5h      return bond, skipping " << endl;
			#endif
		}else{
			//if( (*i).first->degree() > 1 ){
				//cout << "a5h     has degree >1 " << endl;
				// check if target atom was visited
				if( (*i).first->getBFSVectorSize() == 0 ){

					// target atom was not visited, add the bond to toVisit and set the path of the target
					#ifdef DEBUG
					cout << "a5h    adding Bond " << (*i).second->toStringShort() << endl;
					#endif
					toVisit->push_back( (*i).first ); // add this atom  to toVisit
					toVisitBond->push_back( (*i).second ); // add this bond to the list of bonds to follow.

					BFSBondVector.push_back( (*i).second );


					(*i).second->getTarget()->pushBFSVector( &BFSVector, &BFSBondVector ); // set the target's path
				}else{
					#ifdef DEBUG
					cout << "a5h     TARGET WAS VISITED BEFORE" << endl;
					#endif
											// set the target's path

					// target atom was visited. Check if it is a valid ring
					// (if the paths of this and target atom
					// form an intersection with only one atom)

					vector<Atom*> vectorIntersect;
					getVectorIntersect( getBFSVector(),
							(*i).first->getBFSVector(),
							 &vectorIntersect
						);
					#ifdef DEBUG
					cout << "a5h    ring has intersect of size " << vectorIntersect.size() << endl;
					#endif
					if( vectorIntersect.size() == 1 ){
						// if yes create the ring and return a pointer to it
																			BFSVector.push_back( (*i).second->getTarget() );
						#ifdef DEBUG
						cout << "a5h    ADDING LAST BOND: " << (*i).second->toStringShort() << endl;
						#endif
						BFSBondVector.push_back( (*i).second );
						#ifdef DEBUG
						cout << "a5h    NOW CREATING RING" << endl;
						cout << "a5h    merging atoms: " << getBFSVector()->size() << " " << (*i).first->getBFSVector()->size() << endl;
						#endif

						vector<Atom*> vectorUnion;
						//getVectorUnion( getBFSVector(), (*i).first->getBFSVector(), &vectorUnion, &vectorIntersect );
						MoleculeUtils::mergeSet( getBFSVector(), (*i).first->getBFSVector(), &vectorUnion );
						#ifdef DEBUG
						cout << "a5h    union vector has size " << vectorUnion.size() << endl;

						cout << "a5h    merging bonds: " << getBFSBondVector()->size() << " " << (*i).first->getBFSBondVector()->size() << endl;
						#endif

						vector<Bond*> bondVectorUnion;
						MoleculeUtils::mergeBondSet( getBFSBondVector(), (*i).first->getBFSBondVector(), &bondVectorUnion );

						// check that all bond source and target atoms are member of the ring

#ifdef DEBUG
						cout << "a5h    bond union vector has size " << bondVectorUnion.size() << endl;
						#endif


						vector<Bond*> bondVectorMembers;
						MoleculeUtils::selectRingMemberBonds( &bondVectorUnion, &vectorUnion, &bondVectorMembers );


						#ifdef DEBUG
						cout << "a5h    among which " << bondVectorMembers.size() << " are ring members " << endl;
						#endif

						// add reverse bonds to the set of rings
						vector<Bond*>::iterator bi;
						vector<Bond*> reverseBonds;

						for( bi = bondVectorMembers.begin(); bi != bondVectorMembers.end(); bi++ ){
							try{
								Bond* reverse = (*bi)->getReverse();
																				#ifdef DEBUG
									cout << "a5h    added reverse bond " << reverse->toStringShort() << endl;
								#endif
								reverseBonds.push_back( reverse );

							}catch( CError ){
							}
						}

						for( bi = reverseBonds.begin(); bi != reverseBonds.end(); bi++ ){
							bondVectorMembers.push_back( (*bi) );
						}

						delete toVisit;
						delete toVisitBond;

						return( new Ring( &vectorUnion, &bondVectorMembers ) );
					}else{
						// it is not a valid ring, skip this bond

					}
				}

			//}else{
			//	cout << "a5h    skipping Bond " << (*i).second->toStringShort() << " because " << (*i).second->getTarget()->toStringShort() << " has degree 1, cannot be member of a cycle" << endl;
			//}
		}
	}
	#ifdef DEBUG
	cout << "a5h " << toVisit->size() << " in toVisit" << endl;
	#endif
	if( toVisit->size() > 0 ){
		// proceed to the nex toVisit Bond target atom
		return( (*toVisit->begin())->getRingBFS( toVisit, toVisitBond ) );
	}else{
		#ifdef DEBUG
		cout << "a5h NO RING THROWING ATOMNOTINRING" << endl;
		#endif
		delete toVisit;
		delete toVisitBond;

		CError e = CError( ATOMNOTINRING, "Atom " + toStringShort() + " is not a ring member " );
		throw( e );
	}
}

void Atom::pushBFSVector( vector<Atom*>* aPath, vector<Bond*>* aBondPath ){

	#ifdef DEBUG
	cout << "b1h    atom " << toStringShort() << " now has bfsVector:" << endl;
	#endif

	vector<Atom*>::iterator ai;
	for( ai = aPath->begin(); ai != aPath->end(); ai++ ){
		BFSVector.push_back(*ai);
		#ifdef DEBUG
		(*ai)->describeShort();
		#endif
	}

	#ifdef DEBUG
	cout << " and bfsBondVector " << endl;
	#endif
	vector<Bond*>::iterator bi;
	for( bi = aBondPath->begin(); bi != aBondPath->end(); bi++ ){
		BFSBondVector.push_back(*bi);
		#ifdef DEBUG
		(*bi)->describeShort();
		#endif
	}

}



void Atom::getVectorIntersect( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result){
	vector<Atom*>::iterator ai;
	for ( ai = v1->begin(); ai != v1->end(); ai++ ){
		vector<Atom*>::iterator aj;
		for ( aj = v2->begin(); aj != v2->end(); aj++ ){
			if( (*aj) == (*ai) ){
				result->push_back( (*ai) );
			}
		}
	}
}

/*void Atom::getVectorUnion( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result, vector<Atom*>* intersect ){

	cout << "a9h    UNION" << endl;

	vector<Atom*> all;

	vector<Atom*>::iterator ai;
	for ( ai = v1->begin(); ai != v1->end(); ai++ ){
		all.push_back( *ai );
		cout << (*ai)->toStringShort() << ",";
	}
	for ( ai = v2->begin(); ai != v2->end(); ai++ ){
		all.push_back( *ai );
		cout << (*ai)->toStringShort() << ",";
	}
	cout << endl;

	MoleculeUtils::mergeSet( &all, intersect, result );

}*/

string Atom::toStringBFSVector(){
	stringstream out;

	vector<Atom*>::iterator i;
	for( i = BFSVector.begin(); i!= BFSVector.end(); i++){
		out << (*i)->toStringShort() << ", ";
	}
	return( out.str() );
}





/** set the Morgan label of the instance to aLabel
 */
void Atom::setMorganLabel( string aLabel ){
  morganLabel = aLabel;
}


/** set the value of the partial charge of the instance to aValue
 */
void Atom::setPartialCharge(double aValue){
  partialCharge = aValue;
}

/** include the sign of the partial charge in the Morgan label of the instance.
If partial charge > threshold: sign = + ; otherwise, sign = -
 */
void Atom::setMorganChargeLabel(double threshold){

    if(partialCharge > threshold){
	setMorganLabel(getMorganLabel() + "+");
    }
    if(partialCharge < -1.0*threshold){
        setMorganLabel(getMorganLabel() + "-");
    }

}


