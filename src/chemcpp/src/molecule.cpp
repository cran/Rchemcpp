/****************************************************************************************
					  molecule.cpp 
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



#include "molecule.h"
#include "moleculeutils.h"

//#define DEBUGK 1
//#define DEBUG 1
//#define DEBUGPERRETLABEL 1
//#define DEBUGSSSR 1



int Molecule::counter = 0;

Molecule::Molecule() : DataContainer(){

	fastPT = new map<Atom*, map<Atom*, double>* >;
	fastPTSave = new map<Atom*, map<Atom*, double>* >;
	fastPTNext = NULL;

	Molecule::counter++;
	id = Molecule::counter;

	selfKernelCalculated = false;

	addStringDescriptor( "name" ,"", "", "molecule name" );
	addStringDescriptor( "comment" ,"", "", "comment" );
	addStringDescriptor( "comment2" ,"", "", "comment2" );
	addStringDescriptor( "comment3" ,"", "", "comment3" );

	//setIntDescriptor( ACTIVITY ,98, "", "biological activity" );

	selectedFlag = false;

	chiral=false;
	location = "";

	flagHasSSSRDetected = false;

	moleculeChanged(true,true);

	activity = 0;
	flagActivity = false;

	originalFormat = NAVALUE;
}

Molecule::Molecule( Molecule& aMolecule , bool bool_resetMorganIndex ) : DataContainer( (DataContainer&) aMolecule ){

	fastPT = new map<Atom*, map<Atom*, double>* >;
	fastPTSave = new map<Atom*, map<Atom*, double>* >;
	fastPTNext = NULL;

	#ifdef DEBUG
		cout << "Molecule COPY operator" << endl;
	#endif

	//setStringDescriptor( "name" ,"", "", "molecule name" );
	//setStringDescriptor( "comment" ,"", "", "comment" );
	addStringDescriptor( "name" ,"", "", "molecule name" );
	addStringDescriptor( "comment" ,"", "", "comment" );
	addStringDescriptor( "comment2" ,"", "", "comment2" );
	addStringDescriptor( "comment3" ,"", "", "comment3" );

	this -> selectedFlag = aMolecule.selectedFlag;
	this -> chiral = aMolecule.chiral;
	this -> maxMorganIteration = aMolecule.maxMorganIteration;
	this -> selfKernel = aMolecule.selfKernel;
	this -> selfKernelCalculated = aMolecule.selfKernelCalculated;
	this -> sortDescriptorName = aMolecule.sortDescriptorName;
	this -> sortDescriptorType = aMolecule.sortDescriptorType;

	this ->	activity = aMolecule.activity;
	this -> flagActivity = aMolecule.flagActivity;

	this->location = aMolecule.location;
	this->originalFormat = aMolecule.originalFormat;

	this-> flagHasSSSRDetected = aMolecule.flagHasSSSRDetected;


	// increment the molecule id and set the id for this molecule
	Molecule::counter++;
	id = Molecule::counter;

	selfKernelCalculated = false;

	// add atoms and make a map of corresponding atoms in original molecule
	map<Atom*, Atom*> match; // old and new atom
	vector<Atom*>::iterator iit;
	for( iit = aMolecule.beginAtom(); iit != aMolecule.endAtom(); iit++ ){
		Atom* newAtom = new Atom( *(*iit) );
		atoms.push_back( newAtom );
		match[(*iit)] = newAtom;

		newAtom->setCoordinates( (*iit)->getX(), (*iit)->getY(), (*iit)->getZ() );
		newAtom->setAN( (*iit)->getAN() );
		newAtom->setType( (*iit)->getType() );
		newAtom->setMorganLabel( (*iit)->getMorganLabel() );

	}

	// add bonds
	// THIS SHOULD BE IMPLEMENTED YET
	map< Atom*, Bond*>::iterator bi;
	for( iit = aMolecule.beginAtom(); iit != aMolecule.endAtom(); iit++ ){
		for( bi = (*iit)->beginBond(); bi != (*iit)->endBond(); bi++ ){
			Bond* newBond = new Bond(
							match[(*iit)],
							match[(*bi).first],
							(*bi).second->getLabel()
						);
			if( (*bi).second->hasFlagOriginal() ){
				newBond->setFlagOriginal();
			}else{
				newBond->unsetFlagOriginal();
			}

			match[(*iit)]->addBond( newBond );

		}
		//match[(*iit)]->describeShort();

	}

	// copy int descriptors
	map<const string, Descriptor<int>* >::iterator dii;
	for( dii = aMolecule.beginIntDescriptor(); dii != aMolecule.endIntDescriptor(); dii++ ){
		if( (*dii).second->isEmpty() ){
			Descriptor<int>* d = addIntDescriptor(
					(*dii).second->getLabel(),
					-1,
					(*dii).second->getUnit(),
					(*dii).second->getComment() );
			d->setEmpty();

		}else{
			addIntDescriptor(
					(*dii).second->getLabel(),
					(*dii).second->getValue(),
					(*dii).second->getUnit(),
					(*dii).second->getComment() );
		}
	}

	// copy float descriptors
	map<const string, Descriptor<float>* >::iterator dif;
	for( dif = aMolecule.beginFloatDescriptor(); dif != aMolecule.endFloatDescriptor(); dif++ ){
		if( (*dif).second->isEmpty() ){
			Descriptor<float>* d = addFloatDescriptor(
					(*dif).second->getLabel(),
					-1.0,
					(*dif).second->getUnit(),
					(*dif).second->getComment());
			d->setEmpty();

		}else{

			addFloatDescriptor(
					(*dif).second->getLabel(),
					(*dif).second->getValue(),
					(*dif).second->getUnit(),
					(*dif).second->getComment());
		}
	}

	// copy string descriptors
	map<const string, Descriptor<string>* >::iterator dis;
	for( dis = aMolecule.beginStringDescriptor(); dis != aMolecule.endStringDescriptor(); dis++ ){
		if( (*dis).second->isEmpty() ){
			Descriptor<string>* d = addStringDescriptor(
					(*dis).second->getLabel(),
					"",
					(*dis).second->getUnit(),
					(*dis).second->getComment() );
			d->setEmpty();

		}else{

			addStringDescriptor(
					(*dis).second->getLabel(),
					(*dis).second->getValue(),
					(*dis).second->getUnit(),
					(*dis).second->getComment() );
		}
	}

	//cout << "WARNING: the molecule copy operator does not copy fastPQ, fastPS and fastPT!!!" << endl;

	if(bool_resetMorganIndex == true)
	  {
	    resetMorganIndex();
	  }




}


Molecule::Molecule(
			Molecule& m1,
			Molecule& m2,
			double (*pt2AtomKernel)( Atom*, Atom* ),
			double (*pt2BondKernel)( Bond*, Bond* )
			) : DataContainer(){

	atoms.clear();

	fastPT = new map<Atom*, map<Atom*, double>* >;
	fastPTSave = new map<Atom*, map<Atom*, double>* >;
	fastPTNext = NULL;

	Molecule::counter++;
	id = Molecule::counter;

	selfKernelCalculated = false;

	addStringDescriptor( "name" ,"", "", "molecule name" );
	addStringDescriptor( "comment" ,"", "", "comment" );
	addStringDescriptor( "comment2" ,"", "", "comment2" );
	addStringDescriptor( "comment3" ,"", "", "comment3" );
	//setIntDescriptor( ACTIVITY ,98, "", "biological activity" );

	selectedFlag = false;
	chiral=false;

	activity = 0;
	flagActivity = false;

	flagHasSSSRDetected = false;

	vector<Atom*>::iterator iit;
	vector<Atom*>::iterator jit;

	vector< Atom*> source1;
	vector< Atom*> source2;
	vector< Atom*> newAtoms;

	// get a list of all matching atoms in both graphs
	// and add all atoms matching in both graphs
	#ifdef DEBUG
		cout << "adding atoms" << endl;
	#endif
	int ac = 0;
	for ( iit = m1.beginAtom(); iit != m1.endAtom(); iit++ ){
		for ( jit = m2.beginAtom(); jit != m2.endAtom(); jit++ ){
			if ( pt2AtomKernel( *iit, *jit ) == 1 ){
				#ifdef DEBUG
					cout << "  " << ac << endl;
					cout << "    " << (*iit)->toStringShort() << endl;
				#endif

				Atom* newAtom = new Atom( *(*iit) );
				this->addAtom( newAtom , true, true);

				newAtoms.push_back( newAtom );
				ac++;
				source1.push_back( *iit );
				source2.push_back( *jit );



				try{
					double apq = (*iit)->getKashimaPQ( true ) * (*jit)->getKashimaPQ( true );
					newAtom->setKashimaPQ( apq );
					fastPQ[newAtom] = apq;
				} catch( CError e ){
					cout << "WARNING: please set PQ before computing fused graph" << endl;
				}

				try{
					double aps = (*iit)->getKashimaPS( true ) * (*jit)->getKashimaPS( true );
					newAtom->setKashimaPS( aps );
					fastPS[newAtom] = aps;
				} catch( CError e ){
					cout << "WARNING: please set PS before computing fused graph" << endl;
				}

			}
		}
  	}

	#ifdef DEBUG
		cout << "adding atoms DONE" << endl;
		cout << "  " << m1.toStringShort() << " - " << m2.toStringShort() << ": product graph has " << ac << " atoms" << endl;
	#endif



	Bond* bond1 = NULL;
	Bond* bond2 = NULL;
	Bond* bond1r = NULL;
	Bond* bond2r = NULL;

	for( int i = 0; i < numAtoms()-1; i++ ){
		//cout << " " << i << endl;
		for( int j = i+1; j < numAtoms(); j++ ){
			//cout << "  " << j << endl;

			//cout << source1[i]->toStringShort() << " " << source2[i]->toStringShort() << endl;
			//cout << source1[j]->toStringShort() << " " << source2[j]->toStringShort() << endl;

			try{
				bond1 = source1[i]->getBondWithTarget( source1[j] );
				bond2 = source2[i]->getBondWithTarget( source2[j] );

				bond1r = source1[j]->getBondWithTarget( source1[i] );
				bond2r = source2[j]->getBondWithTarget( source2[i] );

				//cout << "checking " << source1[i]->toStringShort() << " with " << source1[j]->toStringShort() << " -> " << bond1->getLabel() << endl;
				//cout << "checking " << source2[i]->toStringShort() << " with " << source2[j]->toStringShort() << " -> " << bond2->getLabel() << endl;

				//if( bondType1 == bondType2 ){
				if( pt2BondKernel( bond1, bond2 ) ){
					//cout << "Linking " << endl;

					Atom* aSource = newAtoms[i];
					Atom* aTarget = newAtoms[j];

					Bond* forwardBond = new Bond( aSource, aTarget, 1 );
					Bond* backwardBond = new Bond( aTarget, aSource, 1 );

					double fpt = bond1->getKashimaPT() * bond2->getKashimaPT();
					double bpt = bond1r->getKashimaPT() * bond2r->getKashimaPT();

					forwardBond->setKashimaPT( fpt );
					backwardBond->setKashimaPT( bpt );

					backwardBond->unsetFlagOriginal();
					forwardBond->setFlagOriginal();

					aSource->addBond( forwardBond );
					aTarget->addBond( backwardBond );

					map<Atom*, double>* aMap = (*fastPT)[aSource];
					if( aMap == NULL ){
						aMap = new map<Atom*, double>;
						(*aMap)[aTarget] = fpt;
						(*fastPT)[aSource] = aMap;
					}else{
						(*aMap)[aTarget] = fpt;
					}

					aMap = (*fastPT)[aTarget];
        				if( aMap == NULL ){
						aMap = new map<Atom*, double>;
						(*aMap)[aSource] = bpt;
						(*fastPT)[aTarget] = aMap;
					}else{
						(*aMap)[aSource] = bpt;
					}


					aMap = (*fastPTSave)[aSource];
        				if( aMap == NULL ){
						aMap = new map<Atom*, double>;
						(*aMap)[aTarget] = fpt;
						(*fastPTSave)[aSource] = aMap;
					}else{
						(*aMap)[aTarget] = fpt;
					}

					aMap = (*fastPTSave)[aTarget];
					if( aMap == NULL ){
						aMap = new map<Atom*, double>;
						(*aMap)[aSource] = bpt;
						(*fastPTSave)[aTarget] = aMap;
					}else{
						(*aMap)[aSource] = bpt;
					}

					//cout << aSource->toStringShort() << " - " << aTarget->toStringShort() << ": " << fpt << " " << bpt << endl;

					//cout << (*fastPTSave[aSource])[aTarget] << " " <<(*fastPTSave[aTarget])[aSource]<< endl;
					//cout << (*fastPT[aSource])[aTarget] << " " <<(*fastPT[aTarget])[aSource]<< endl;

					//linkAtoms( newAtoms[i], newAtoms[j], bondType1 );

				}
			}catch( CError ){
			}
			//cout << "####" << endl;
		}


	}


	source1.clear();
	source2.clear();
	newAtoms.clear();

	location = "";
	originalFormat = NAVALUE;

	//resetMorganIndex();
	//cout << "DONE" << endl;

}




// 3D product graph constructor  //
//------------------------------//
Molecule::Molecule(
			Molecule& m1,
			Molecule& m2,
			double (*pt2AtomKernel)( Atom*, Atom*),
			double (*pt2BondKernel)( float, float, float ),
			float edgeKernelParameter
			) : DataContainer(){

  // step 0 : initialization
  // -----------------------
  atoms.clear();
  Molecule::counter++;
  id = Molecule::counter;
  
  addStringDescriptor( "name" ,"", "", "molecule name" );
  addStringDescriptor( "comment" ,"", "", "comment" );
  addStringDescriptor( "comment2" ,"", "", "comment2" );
  addStringDescriptor( "comment3" ,"", "", "comment3" );
  //setIntDescriptor( ACTIVITY ,98, "", "biological activity" );
  
  selfKernelCalculated = false;
  selectedFlag = false;
  chiral=false;
  activity = 0;
  flagActivity = false;
  flagHasSSSRDetected = false;
  

        // NOT USED YET, but need to be declared when deleting molecule
	fastPT = new map<Atom*, map<Atom*, double>* >;
	fastPTSave = new map<Atom*, map<Atom*, double>* >;
	fastPTNext = NULL;


  
  // step 1 : build atoms
  // --------------------
  vector<Atom*>::iterator iit;
  vector<Atom*>::iterator jit;
  
  vector< Atom*> source1;
  vector< Atom*> source2;
  vector< Atom*> newAtoms;
  
  // get a list of all matching atoms in both graphs
  // and add all atoms matching in both graphs
  int ac = 0;
  double k_at1_at2 = 0.0;
  
  for ( iit = m1.beginAtom(); iit != m1.endAtom(); iit++ ){
    for ( jit = m2.beginAtom(); jit != m2.endAtom(); jit++ ){
      k_at1_at2 = pt2AtomKernel( *iit, *jit);
      if ( k_at1_at2 != 0 ){
	Atom* newAtom = new Atom( *(*iit) );
	this->addAtom( newAtom , true, true);
	
	newAtoms.push_back( newAtom );
	ac++;
	source1.push_back( *iit );
	source2.push_back( *jit );
      }
    }
  }
  
  // step 2 : build edges
  // --------------------
  // --> compute bonds + adjacency matrix
  double k_e1_e2 = 0;
  double k_n1_n2 = 0;
  double k_n2_n1 = 0;
  float d1 = 0;
  float d2 = 0;

  //init adjacency matrix + walks matrix
  adjacency = new vector< vector<double> >;
  walks = new vector< vector<double> >;

  for(int i = 0 ; i < numAtoms() ; i++){
    adjacency->push_back(vector<double>() );
    walks->push_back(vector<double>() );
    for(int j = 0 ; j < numAtoms() ; j++){
      (*adjacency)[i].push_back(0.0);
      (*walks)[i].push_back(0.0);
    }
  }
  

  for( int i = 0; i < numAtoms()-1; i++ ){
    for( int j = i+1; j < numAtoms(); j++ ){
      
      d1 = atomicDistance(source1[i],source1[j]);
      d2 = atomicDistance(source2[i],source2[j]);

      k_e1_e2 = pt2BondKernel(d1,d2,edgeKernelParameter);

	// WARNING : new test to be sure to introduce a correct edge
	if( source1[i] == source1[j] || source2[i] == source2[j] ){
		k_e1_e2 = 0;
	} 
	//else{
	//	cout << "distance 1 : " << d1 << endl;
	//	cout << "distance 2 : " << d2 << endl; 
	//}


      k_n1_n2 = pt2AtomKernel(source1[i],source2[i]);
      k_n2_n1 = pt2AtomKernel(source1[j],source2[j]);

      //cout << "#### atom kernel value : " << k_n2_n1 << endl;
      //cout << "#### edge kernel value : " << k_e1_e2 << endl;


      // 1 - set adjacency matrix + walks matrix
      setAdjacency(i, j, k_n1_n2 * k_e1_e2);
      setAdjacency(j, i, k_n2_n1 * k_e1_e2);

      setWalks(i, j, k_n1_n2 * k_e1_e2);
      setWalks(j, i, k_n2_n1 * k_e1_e2);

      // 2 - introdce atomic bonds : a pair of bonds : forward (i->j) + backward (j->i)
      if(k_e1_e2 != 0){
	Atom* aSource = newAtoms[i];
	Atom* aTarget = newAtoms[j];
	
	Bond* forwardBond = new Bond( aSource, aTarget, 1 );
	Bond* backwardBond = new Bond( aTarget, aSource, 1 );
	
	forwardBond->setKashimaPT( k_n1_n2 * k_e1_e2 );
	backwardBond->setKashimaPT( k_n2_n1 * k_e1_e2 );
	
	//backwardBond->unsetFlagOriginal();
	//forwardBond->setFlagOriginal();
	
	aSource->addBond( forwardBond );
	aTarget->addBond( backwardBond );
      }	
      
    }
  }


}



// function to compute the Euclidian distance between a pair of atoms
float Molecule::atomicDistance(Atom* atom1, Atom* atom2){

  float distance;

  distance = 
    (atom1->getX() - atom2->getX())*(atom1->getX() - atom2->getX()) 
    +  (atom1->getY() - atom2->getY())*(atom1->getY() - atom2->getY()) 
    +  (atom1->getZ() - atom2->getZ())*(atom1->getZ() - atom2->getZ());

  return( sqrt(distance) );

}



/** this function should be called everytime a molecule is modified
*/
void Molecule::compute(){
	#ifdef DEBUG
	cout << "NOW IN MOLECULE::COMPUTE " << endl;
	#endif
	if( !hasSSSRDetected() ){
	//	detectSSSR();
	}
	setMorganLabels(0);
	try{
	  setPerretLabels();
	}catch( CError e ){
		cout << "Warning: SSSR was not computed in molecule " << toStringShort() << " so cannot set Perret Labels, skippink!" << endl;
	}
}



void Molecule::moleculeChanged( bool resetSSSR, bool bool_resetMorganIndex ){
    if( bool_resetMorganIndex ){
        resetMorganIndex();
     }
    if( resetSSSR ){
        flagHasSSSRDetected = false;
    }
}




Molecule& Molecule::operator=( const Molecule& aMolecule ) {

	#ifdef DEBUG
		cout << "Molecule = operator" << endl;
	#endif

	if( this != &aMolecule ){

		this->kindStringDescriptors = aMolecule.kindStringDescriptors;
		this->kindIntDescriptors = aMolecule.kindIntDescriptors;
		this->kindFloatDescriptors = aMolecule.kindFloatDescriptors;

		this -> selectedFlag = aMolecule.selectedFlag;

		this -> location = aMolecule.location;

		// Copy all non kind descriptors
		// first copy int descriptors

		map<const string, Descriptor<int>* >::iterator iti;
		for( iti = (intDescriptors).begin(); iti != (intDescriptors).end(); iti++ ){
				addIntDescriptor( (*iti).second->getLabel() , (*iti).second->getValue(), (*iti).second->getUnit(), (*iti).second->getComment() );
		}
		map<const string, Descriptor<float>* >::iterator itf;
		for( itf = (floatDescriptors).begin(); itf != (floatDescriptors).end(); itf++ ){
				addFloatDescriptor( (*itf).second->getLabel() , (*itf).second->getValue(), (*itf).second->getUnit(), (*itf).second->getComment() );
		}
		map<const string, Descriptor<string>* >::iterator its;
		for( its = (stringDescriptors).begin(); its != (stringDescriptors).end(); its++ ){
				addStringDescriptor( (*its).second->getLabel() , (*its).second->getValue(), (*its).second->getUnit(), (*its).second->getComment() );
		}

		// increment the atom id and set the id for this atom
		Molecule::counter++;
		id = Molecule::counter;

//		selfKernelCalculated = false;
	}

	//resetMorganIndex(); ***???

	return( *this );
}


Molecule::~Molecule(){

	// delete all atoms
	erase();

	delete fastPT;
	delete fastPTSave;
	delete fastPTNext;
}


Atom* Molecule::addAtom( string aSymbol, bool resetSSSR ) throw( CError ){

	/*ToUpper      up(std::locale::classic());
 	ToLower      down(std::locale::classic());

	std::transform (aSymbol.begin(), aSymbol.begin()+1, aSymbol.begin(), up);
	if( aSymbol.length() > 1 ){
		std::transform (aSymbol.begin()+1, aSymbol.end(), aSymbol.begin()+1, down);
	} */

	string h = aSymbol.substr(0, 1);
	string r = "";
	if( aSymbol.length() > 1 ){
		r = aSymbol.substr( 1, aSymbol.length()-1 );
	}

	h = StringUtils::toUpper( h );
	r = StringUtils::toLower( r );
	aSymbol = h + r;


	Atom* anAtom = NULL;

	try{
		anAtom = new Atom( *elements[aSymbol] );
	}catch( CError e ){
		//cerr << "WARNING could not add element " << aSymbol << endl;
		throw( e );
	}


	atoms.push_back(anAtom);
	moleculeChanged( resetSSSR , true);

	#ifdef DEBUG
		//cout << "created atom " << anAtom->toString() << endl;
	#endif
	return( anAtom );

}


Atom* Molecule::addAtom( Atom* anAtom, bool resetSSSR, bool resetMorganIndex ) throw( CError ){
    // if atom does not already exist, then add the new atom
    if( atomExists( anAtom ) ){
        //cout << "ERROR ATOM EXISTS!!!!" << endl;

        stringstream out;
        out << "Atom " << anAtom->toString() << " already exists ";
        CError e( ATOMALREADYEXISTS, out.str() );
        e.describe();
        throw(e);
        //return( anAtom );
    }else{
        atoms.push_back( anAtom );
        //cout << "adding new atom " << anAtom->toString() << endl;
        moleculeChanged( resetSSSR, resetMorganIndex );
        return( anAtom );
    }
}




/** check whether Atom exists in the atom vector */
bool Molecule::atomExists(Atom* anAtom){
	vector< Atom* >::iterator i;
	for( i = atoms.begin(); i != atoms.end(); i++ ){
		if( *i == anAtom ){
			return( true );
		}
	}
	return( false );
}

Atom* Molecule::getAtom(int anId) throw( CError ){

	vector<Atom*>::iterator atom;
	// for each atom
	for( atom = beginAtom(); atom != endAtom(); atom++ ){
		if( (*atom)->getId() == anId ){
			return( *atom );
		}
	}
	stringstream out;
	out << "No atom with id " << anId << " in molecule " << getId();
	CError e(ERRORATOMNOTFOUND, out.str() );
	e.describe();
	describeLong();
	throw(e);
}


Atom* Molecule::getAtomByIndex(int ind){
	return atoms[ind];
}




Bond* Molecule::linkAtoms( Atom* aSource, Atom* aTarget, int aBondLabel, int aBondStereo,  int aBondNotUsed, int aBondTopology, int aBondReactionCenter, bool resetSSSR ){
	if( !atomExists( aSource ) ){
		addAtom( aSource , true, true);
	}
	if( !atomExists( aTarget ) ){
		addAtom( aTarget , true , true );
	}

	Bond* forwardBond = new Bond(aSource, aTarget, aBondLabel, NAVALUE, aBondStereo, aBondNotUsed, aBondTopology, aBondReactionCenter );
	Bond* backwardBond = new Bond(aTarget, aSource, aBondLabel, NAVALUE, aBondStereo, aBondNotUsed, aBondTopology, aBondReactionCenter);

	backwardBond->unsetFlagOriginal();
	forwardBond->setFlagOriginal();

	aSource->addBond(forwardBond);
	aTarget->addBond(backwardBond);

	moleculeChanged( resetSSSR, true );
	return( forwardBond );
}

Bond* Molecule::linkAtomsNoReturn( Atom* aSource, Atom* aTarget, int aBondLabel, int aBondStereo,  int aBondNotUsed, int aBondTopology, int aBondReactionCenter ){
	//addAtom( aSource );
	//addAtom( aTarget );

	Bond* forwardBond = new Bond( aSource, aTarget, aBondLabel, NAVALUE, aBondStereo, aBondNotUsed, aBondTopology, aBondReactionCenter );
	//forwardBond->setFlagOriginal();

	forwardBond = aSource->addBond( forwardBond );

	//resetMorganIndex();
	// WARNING !!!! 
	// --> modified by PM ??
	//   --> PM needs this modification in function Molecule::threeDtransform 
	//       in order not to loose the MI related to te 2D structure when switching to 3D structure (i.e., completely connected graph)

	return( forwardBond );
}


Bond* Molecule::linkAtoms( int aSource, int aTarget, int aBondLabel, int aBondStereo,  int aBondNotUsed, int aBondTopology, int aBondReactionCenter, bool resetSSSR ) throw( CError ){
	// atoms must exists otherwise cannot link

	if( max(aSource, aTarget) > numAtoms()-1 ){
		stringstream out;
		out << "molecule " << toString() << " has only " << numAtoms() << " atoms, so cannot link atom " << aSource << " to " << aTarget ;
		CError e = CError(NOTENOUGHATOMSINMOLECULE, out.str() );
		e.describe();
		throw(e);
	}else{
		try{
			Bond* forwardBond = linkAtoms( atoms[aSource], atoms[aTarget], aBondLabel, aBondStereo, aBondNotUsed, aBondTopology, aBondReactionCenter, resetSSSR);
			return( forwardBond );
		}catch( CError e){
			throw( e );
		}

	}
}

Bond* Molecule::linkAtoms(string aSource, string aTarget, int aBondLabel, int aBondStereo, int aBondNotUsed, int aBondTopology, int aBondReactionCenter, bool resetSSSR ) throw( CError ){
	// atoms must exists otherwise cannot link

	int a1 = 0;
	int a2 = 0;
	int i = 0;
	bool foundSource = false;
	bool foundTarget = false;

	// find aSource and aTarget

	vector<Atom*>::iterator atom;
	for( atom = beginAtom(); atom != endAtom(); atom++ ){

		if( (*atom)->getName() == aSource) {
			a1 = i;
			foundSource = true;
		}
		if( (*atom)->getName() == aTarget) {
			a2 = i;
			foundTarget = true;
		}
		i++;
	}

	if( foundSource == false ){
		stringstream out;
		out << "Molecule::linkAtoms: atom " << aSource << " has not been found in molecule " << getName();
		CError e = CError(NOTFOUND, out.str() );
		e.describe();
		throw(e);
	}
	if( foundTarget == false ){
		stringstream out;
		out << "Molecule::linkAtoms: atom " << aTarget << " has not been found in molecule " << getName();
		CError e = CError(NOTFOUND, out.str() );
		e.describe();
		throw(e);
	}

	#ifdef DEBUG
		//cout << "linking " << atoms[a1]->getName() << " to " << atoms[a2]->getName() << endl;
	#endif

	Bond* forwardBond = linkAtoms( atoms[a1], atoms[a2], aBondLabel, aBondStereo, aBondNotUsed, aBondTopology, aBondReactionCenter, resetSSSR);
	return( forwardBond );
}

/** returns a string description of the molecule */
string Molecule::toString(){
	stringstream out;
	//out << "Molecule " << getName() << " (id = " << getId() << ")";
	//if( isSelected() ){
	//	out << " SELECTED";
	//}
	out << toStringShort();
	try{
	  out << " at " << this << ", has " << numAtoms() << " atoms ( + " << numHiddenAtoms() << " hidden atoms ) and " << numBonds() << " bonds ( + " << numHiddenBonds() << " hidden bonds ), activity: " << getActivity();
	} catch( CError ){
	  out << " at " << this << ", has " << numAtoms() << " atoms ( + " << numHiddenAtoms() << " hidden atoms ) and " << numBonds() << " bonds ( + " << numHiddenBonds() << " hidden bonds ), no activity";
	}

	out << " / " << sssr.size() << " rings, size: ";

	vector<Ring*>::iterator ri;
	for( ri = sssr.begin(); ri != sssr.end(); ri++ ){
		out << (*ri)->size() << " ";
	}


	return( out.str() );

}

string Molecule::toStringShort(){
	stringstream out;
	out << "Molecule " << getName() << " (id = " << getId() << ") ";
	if( isSelected() ){
		out << " SELECTED";
	}
	return( out.str() );
}

string Molecule::toStringLong(){
	stringstream out;

	vector<Atom*>::iterator atom;
	map<Atom*, Bond*>::iterator bond;

	out << toString() << endl;
	try{
		if( numRings() == 0 ) {
			out << " no rings ";
		}
	}catch( CError ){
		out << " sssr not detected ";
	}
	// for each atom
	out << "ATOMS: " << endl;
	for( atom = beginAtom(); atom != endAtom(); atom++ ){
		out << (*atom)->toString();

		try{
			out << " " << (*atom)->getMorganLabel(true);
		}catch( CError e ){
			out << " morganLabel NOT SET";
		}

		try{
			out << " (ps=" << (*atom)->getKashimaPS(true) << ", pq=" << (*atom)->getKashimaPQ(true) << ")";
		}catch( CError e ){
			out << " (ps=NOT SET, pq=NOT SET)";
		}

		/*try{

			if( hasRing() ){
				if( (*atom)->numRings() > 0 ){
					out << " " << (*atom)->numRings();

				}else{
					out << " no rings " << endl;
				}
			}
		}catch( CError e ){
			out << " rings not detected" << endl;
		}*/

		out << endl;

		// and for each bond
		//		for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
		//			cout << "    " << (*bond).second->toString() << " (pt=" << (*bond).second->getKashimaPT() << ")" <<endl;
		//		}
	}

	out << "BONDS: " << endl;
	for( atom = beginAtom(); atom != endAtom(); atom++ ){
		for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
			out << (*bond).second->toString() << endl;
		}
	}


	out << "Number of Distinct Morgan Indices: " << getNumberOfDistinctMorganIndices( 1 ) << " " << getNumberOfDistinctMorganIndices( 2 ) << " " << getNumberOfDistinctMorganIndices( 3 ) << " " << endl;
	out << "Maximum number of distinct Morgan Indices reached after " << getMaxMorganIteration() << " iterations with " <<  getNumberOfDistinctMorganIndices( getMaxMorganIteration() ) << " distinct connectivity values" << endl;

	return( out.str() );
}


/** convenience function to retrieve the id as a string */
string Molecule::getIdString(){
	stringstream out;
	out << getId();
	return( out.str() );
}

float Molecule::getActivity( bool silentMode ) throw( CError ){
    if( hasActivity() ){
	return( activity );
    }else{
	stringstream out;
	out << toStringShort() << " has no activity";
	CError e( ERRORNA, out.str() );
	if( silentMode != true ){
	    e.describe();
	}
	throw(e);
    }
}


/** write a description of the molecule to stdout */
void Molecule::describe(){
	cout << toString() << endl;
	//DataContainer::describe();
}
void Molecule::describeShort(){
	cout << toStringShort() << endl;
}
void Molecule::describeLong(){
	cout << toStringLong() << endl;
}

void Molecule::describeEachAtom(){
	vector<Atom*>::iterator atom;
	map<Atom*, Bond*>::iterator bond;

	cout << "Molecule " << getId() << " has " << numAtoms() << " atoms: " << endl;
	// for each atom
	for( atom = beginAtom(); atom != endAtom(); atom++ ){
		(*atom)->describe();
		// describe bonds
		cout << "    " << "this atom has " << (*atom)->numBonds() << " bonds:" << endl;
		for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
			(*bond).second->describe();
		}
	}
}


/** Sets the start, stop and transition probabilities for the calculation of the Kashima kernel as published */
void Molecule::setKashimaKernelProb( double aPq, bool skipSkeleton ){
	// aPq must be a probability between 0 and 1 if not, this is an error

	if(aPq > 1 || aPq < 0){
		stringstream out;
		out << aPq << " is not a valid stop probability (should be between 0 and 1) ";
		CError e(VALUENOTALLOWED, out.str() );
		e.describe();
		throw(e);
	}

	if( skipSkeleton == false ){

		int numatoms = numAtoms();

		if(numatoms>0){
			vector<Atom*>::iterator atom;
			map<Atom*, Bond*>::iterator bond;
			int numbonds = 0;
			for( atom = beginAtom(); atom != endAtom(); atom++ ){

				// for each atom set the start probability to 1/numAtoms
				(*atom)->setKashimaPS( 1/(double) numatoms );

				// for each bond set the transition probability to 1/numbonds
				numbonds = (*atom)->numBonds();
				for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
					// For each bond set the transition probability
					(*bond).second->setKashimaPT( (1.0 - aPq) / numbonds );
				}

  	  	// for each atom set the stop probability to (1-aPq)/numBonds
				(*atom)->setKashimaPQ( aPq );
			}
		}
	}else{
		cout << "setting rdmwk parameters for molecule " << getName() << endl;

		int numatoms = numAtomsNonCSkeleton();

		cout << "found " << numatoms << " non skeleton atoms " << endl;

		if(numatoms>0){
			vector<Atom*>::iterator atom;
			map<Atom*, Bond*>::iterator bond;
			int numbonds = 0;
			for( atom = beginAtom(); atom != endAtom(); atom++ ){

				// for each atom set the start probability to 1/numAtoms
				if( (*atom)->isCSkeleton() ){
					(*atom)->setKashimaPS( 0.0 );
				}else{
					(*atom)->setKashimaPS( 1/(double) numatoms );
				}

				// for each bond set the transition probability to 1/numbonds
				numbonds = (*atom)->numBonds();
				for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
					// For each bond set the transition probability
					(*bond).second->setKashimaPT( (1.0 - aPq) / numbonds );
				}

  	  	// for each atom set the stop probability to (1-aPq)/numBonds
				(*atom)->setKashimaPQ( aPq );
			}
		}


	}
}

double Molecule::computeKernel(
				Molecule* anotherMolecule,
				double (*pt2GraphKernel)
				(
					Molecule* mol1,
					Molecule* mol2,
					double (*pt2AtomKernel)( Atom*, Atom* ),
					double (*pt2BondKernel)( Bond*, Bond* ),
					int parameter1, int parameter2
				),
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)( Bond*, Bond* ),
				int parameter1, int parameter2 ){

				//cout << "START" << endl;
				double k = pt2GraphKernel( this, anotherMolecule, pt2AtomKernel, pt2BondKernel, parameter1, parameter2 );
				//cout << "END" << endl;

				return ( k );
}


double Molecule::calculateSelfKernel(
				double (*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*),
				double (*pt2BondKernel)(Bond*, Bond*), int, int ),
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)( Bond*, Bond* ),
				int parameter1, int parameter2 ){

	//selfKernel = this->moleculeKernel(this, pt2AtomKernel, pt2BondKernel, convergenceCondition);

	//cout << "NOW computing self kernel of molecule " << getName() << endl;

	selfKernel = pt2GraphKernel( this, this, pt2AtomKernel, pt2BondKernel, parameter1, parameter2 );
  //            moleculeKernel( this, anotherMolecule, pt2AtomKernel, pt2BondKernel, convergenceCondition ) );

	selfKernelCalculated = true;
	return( selfKernel );
}



/** move the hydrogen atoms to the hiddenAtoms container, so that further algorithms will run without hydrogens. Returns the number of hydrogens hidden */
int Molecule::hideHydrogens(){

	int i = 0;
	vector<Atom*>::iterator atom;

	for( atom = beginAtom(); atom != endAtom(); atom++ ){

		// hide hidrogen atoms
		if( (*atom)->getAN() == 1 ){
			//cout << "Hiding H atom" << endl;
			hiddenAtoms.push_back( (*atom) );
			//cout << "erasing atom" << endl;
			atoms.erase( atom );
			//cout << "done" << endl;
			if( atom != beginAtom() ){
				atom--;
			}
			i++;

		}
		//cout << "hiding bonds" << endl;
		// hide bonds to hydrogens
		(*atom)->hideHydrogenBonds();
		//cout << "done" << endl;


	}
	return(i);
}

/** restores the hidden atoms (for example hydrogens) */
int Molecule::restoreHiddenAtoms( bool flagRestoreBonds ){
	//cout << "Molecule::restoreHiddenAtoms: restoring hidden atoms" << endl;

	int i = 0;
	vector<Atom*>::iterator atom;

	for( atom = hiddenAtoms.begin(); atom != hiddenAtoms.end(); atom++ ){
		//cout << "    " << i << endl;
		// restore hydrogen atoms
		atoms.push_back( (*atom) );
		//hiddenAtoms.erase(atom);
		i++;
	}
	hiddenAtoms.clear();

	if( flagRestoreBonds == true ){
		for( atom = atoms.begin(); atom != atoms.end(); atom++ ){
			(*atom)->restoreHiddenBonds();
		}
	}

	return(i);
}

int Molecule::restoreHiddenBonds(){
	//cout << "Molecule::restoreHiddenAtoms: restoring hidden atoms" << endl;
	int i = 0;
	vector<Atom*>::iterator atom;
	for( atom = atoms.begin(); atom != atoms.end(); atom++ ){
		i += (*atom)->restoreHiddenBonds();
	}
	return(i);
}

/** This function reads the molecule description from a Mol file */
void Molecule::readMOL( string aFilename, bool genericAtomType ) throw( CError ) {
  // erase all eventually existing atoms in the molecule
	erase();

	// open molfile
	ifstream inFile;
 	inFile.open( aFilename.c_str(), ios::in );
	if( !inFile.good() ){
		CError e = CError(FILENOTFOUND, aFilename + " file not found");
		e.describe();
		throw(e);
	}
	#ifdef DEBUG
		cout << " reading readMDLHeaderBlock " << endl;
	#endif
	MoleculeUtils::readMDLHeaderBlock( *this, inFile, aFilename );
	#ifdef DEBUG
		cout << " reading readMDLCtabBlock " << endl;
	#endif
	MoleculeUtils::readMDLCtabBlock( *this, inFile, genericAtomType );

	inFile.close();

	location = aFilename;

	compute();
}

/** This function reads the molecule description from a Mol file */
void Molecule::readMOLOld( string aFilename ) throw( CError ) {

	stringstream out;
	out << "Molecule::readMOLOld is deprecated use Molecule::readMol instead";
	CError e(NOTIMPLEMENTED, out.str() );
	e.describe();
	throw(e);

	// erase all eventually existing atoms in the molecule
	erase();

	// open molfile
	ifstream inFile;
	inFile.open( aFilename.c_str(), ios::in );
	if( !inFile.good() ){
		CError e = CError(FILENOTFOUND, aFilename + " file not found");
		e.describe();
		throw(e);
	}

	int linesize = 256;
	char *line = new char[linesize];
	string stringLine;

	int nbAtoms = 0;
	int nbBonds = 0;
	int nbAtomlist = 0;
	int nbStext = 0;
	int chiral = 0;
	int nbProperties = 0;
	string version = "";

	int bondType = 0;

	Atom* anAtom = NULL;

	int lineNumber = 0;
	int j = 0;
	// for each line in file
	while( !inFile.eof() ){
		lineNumber++;

		//read line
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;

		if( lineNumber == 1 ){
			// name
			setStringDescriptor( "comment", stringLine, "", "", true, true );

			//vector< string > pathWords;
			//StringUtils::Split(aFilename, "/", pathWords);
			//char fileNameOnly[256] = "";
			//strcpy( pathWords[0].c_str(), fileNameOnly );
			//cout << "MMMM setting name to " << fileNameOnly << endl;
			//setName( fileNameOnly );
			setName( aFilename );

		}else if( lineNumber == 2 ){
			// comment (skipping line 3)
		}else if( lineNumber == 4 ){
			// counts line
			//cout << lineNumber << ": " << stringLine << endl;

			nbAtoms = atoi( stringLine.substr(0,3).c_str() );
			nbBonds = atoi( stringLine.substr(3,3).c_str() );
			nbAtomlist = atoi( stringLine.substr(6,3).c_str() );
			nbStext = atoi( stringLine.substr(15,3).c_str() );
			chiral = atoi( stringLine.substr(9,3).c_str() );
			nbProperties = atoi( stringLine.substr(30,3).c_str() );
			version = stringLine.substr(33,6);

			//cout << "found " << nbAtoms << " atoms" << " and " << nbBonds << " bonds "<< endl;

		}else if( lineNumber > 4 ){
			// data
			if( nbAtoms > 0 && j < nbAtoms ){
				// atoms
				// cout << "atom " << lineNumber << ": " << stringLine << endl;
				anAtom = addAtom( StringUtils::rmSpace( stringLine.substr(30,3) ) );
				anAtom->setCoordinates( atof( stringLine.substr(0,10).c_str() ), atof( stringLine.substr(10,10).c_str() ), atof( stringLine.substr(20,10).c_str() ) );

				j++;
			}else if( j >= nbAtoms && j < nbAtoms + nbBonds ) {
				// bonds
        			//cout << "bond " << lineNumber << ": " << stringLine << endl;
				bondType = atoi( stringLine.substr(6,3).c_str() );
				//bondTopology = atoi( stringLine.substr(15,3) )
				linkAtoms( atoi( stringLine.substr(0,3).c_str() ) - 1, atoi( stringLine.substr(3,3).c_str() ) - 1, bondType );
				j++;
			}else if( j >= nbAtoms + nbBonds && j < nbAtoms + nbBonds + nbAtomlist + nbStext ){
				// atomlist block and stext block (SKIPPING)
				//cout << "skipping line " << lineNumber << " in file " << aFilename << ": " << stringLine << endl;
				j++;
			}else if( j >= nbAtoms + nbBonds + nbAtomlist + nbStext && j < nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties ){
				// properties block (SKIPPING)
				//cout << "skipping line " << lineNumber << " in file " << aFilename << ": " << stringLine << endl;
				j++;
			}else{
  			//cout << "skipping line " << lineNumber << " in file " << aFilename << ": " << stringLine << endl;
				j++;
			}
		}else{
			//cout << "skipping line " << lineNumber << " in file " << aFilename << ": " << stringLine << endl;
		}
	}
	delete[] line;


}


/** Erases all atoms in the current molecule. Returns the molecule in the state of construction by the default constructor.
*/
void Molecule::erase(){
#ifdef DEBUG
		cout << "Molecule::erase()" << endl;
#endif

#ifdef DEBUG
		cout << " deleting fastPQ" << endl;
#endif
	fastPQ.clear();
#ifdef DEBUG
		cout << " deleting fastPQ DONE" << endl;
		cout << " deleting fastPS" << endl;
#endif
	fastPS.clear();
#ifdef DEBUG
		cout << " deleting fastPS DONE" << endl;
		cout << " deleting fastPT" << endl;
#endif

	map<Atom*, map<Atom*, double>* >::iterator j;
	for( j = fastPT->begin(); j != fastPT->end(); j++ ){
		delete (*j).second; // delete map
	}
	fastPT->clear();

#ifdef DEBUG
		cout << " deleted fastPT DONE" << endl;
#endif

	for( j = fastPTSave->begin(); j != fastPTSave->end(); j++ ){
		delete (*j).second; // delete map
	}
	fastPTSave->clear();


	if( fastPTNext != NULL){
		for( j = fastPTNext->begin(); j != fastPTNext->end(); j++ ){
			delete (*j).second; // delete map
		}
		fastPTNext->clear();
	}

#ifdef DEBUG
		cout << " deleted fastPTSave DONE" << endl;
#endif

	// delete all atoms
#ifdef DEBUG
	  cout << "DELETING " << atoms.size() << " ATOMS" << endl;
#endif
	vector<Atom*>::iterator i;
	for( i = atoms.begin(); i != atoms.end(); i++ ){
#ifdef DEBUG
			cout << "deleting " << (*i)->toStringShort() << endl;
#endif
		delete (*i); // delete atom
	}
	atoms.clear();

	eraseHiddenAtoms(); // delete hidden atoms

	eraseRings();

	selfKernelCalculated = false;
#ifdef DEBUG
		cout << "Molecule::erase() done" << endl;
#endif

}



void Molecule::eraseAdjacency(){
 
  for(int i = 0 ; i < numAtoms() ; i++){ // NB : dimension(adjacency) = numAtoms()
    (*adjacency)[i].clear();
  }
  delete adjacency;
  
}

void Molecule::eraseWalks(){

  for(int i = 0 ; i < numAtoms() ; i++){
    (*walks)[i].clear();
  }
  delete walks;
  
}



void Molecule::eraseHiddenAtoms(){
	vector<Atom*>::iterator i;

	// delete hidden atoms
	for( i = hiddenAtoms.begin(); i != hiddenAtoms.end(); i++ ){
		#ifdef DEBUG
			cout << "deleting " << (*i)->toStringShort() << endl;
		#endif
		delete (*i); // delete atom
	}
	hiddenAtoms.clear();
}

void Molecule::eraseRings(){
	// erase rings;
	vector<Ring*>::iterator ri;
	for( ri = sssr.begin(); ri != sssr.end(); ri++ ){
		delete (*ri);
	}
	sssr.clear();
}


/** unset all bonds flags */
void Molecule::unsetBondFlags(){
	vector<Atom*>::iterator ait;
	for( ait = beginAtom(); ait != endAtom(); ait++ ){
		(*ait)->unsetBondFlags();
	}
}

/** unset all bonds flags */
void Molecule::unsetBondFlagsOriginal(){
	vector<Atom*>::iterator ait;
	for( ait = beginAtom(); ait != endAtom(); ait++ ){
		(*ait)->unsetBondFlagsOriginal();
	}
}

/** returns the number of bonds in the molecule */
int Molecule::numBonds(){
	int res = 0;

	if( numAtoms() > 1 ){
		vector<Atom*>::iterator ait;
		for( ait = beginAtom(); ait != endAtom(); ait++ ){
			res = res + (*ait)->numBonds();
		}
		return( res/2 );
	}else{
		return( 0 );
	}
}

long Molecule::bondSum(){
	long res = 0;

	if( numAtoms() > 1 ){
		vector<Atom*>::iterator ait;
		for( ait = beginAtom(); ait != endAtom(); ait++ ){
			res += (*ait)->bondSum();
		}
		return( res );
	}else{
		return( 0 );
	}
}

int Molecule::numHiddenBonds(){
	int res = 0;

	if( numAtoms()>1 ){
		vector<Atom*>::iterator ait;
		for( ait = atoms.begin(); ait != atoms.end(); ait++ ){
			res = res + (*ait)->numHiddenBonds();
		}
		for( ait = hiddenAtoms.begin(); ait != hiddenAtoms.end(); ait++ ){
			res = res + (*ait)->numBonds();
		}

		return( res/2 );
	}else{
		return( 0 );
	}
}


/** write the molecule definition in a mdl mol file */
void Molecule::writeMOL( string aFilename ){
	ofstream outFile;
 	outFile.open( aFilename.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFilename + " could not open file");
		e.describe();
		throw(e);
	}
	MoleculeUtils::writeMDLHeaderBlock( *this, outFile );
	MoleculeUtils::writeMDLCtabBlock( *this, outFile );
	outFile.close();
}

/** writes a dot file representing the molecule. The molecule can be converted to graphioc representationj using Graphviz's dot program */
void Molecule::writeDOT( string aFilename, bool perretLabels ){
	ofstream outFile;
 	outFile.open( aFilename.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFilename + " could not open file");
		e.describe();
		throw(e);
	}
	MoleculeUtils::writeDOTGraph( *this, outFile, perretLabels );
	outFile.close();
}


/** write a MDL structure data (SD) file with a single molecule */
void Molecule::writeSD( string aFilename ){
	ofstream outFile;
 	outFile.open( aFilename.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFilename + " could not open file");
		e.describe();
		throw(e);
	}
	MoleculeUtils::writeMDLHeaderBlock( *this, outFile );
	MoleculeUtils::writeMDLCtabBlock( *this, outFile );
	MoleculeUtils::writeMDLNSDBlock( *this, outFile );

	outFile.close();
}



int Molecule::markFragments(){

	// add componentIndex intDescriptor to all atoms
	//cout << "Molecule::markFragments: adding descriptor COMPONENTINDEX" << endl;
	Descriptor< int >* aDescriptor;
	vector< Atom* >::iterator i;
	for( i = atoms.begin(); i != atoms.end(); i++ ){
		aDescriptor = (*i)->setIntDescriptor( COMPONENTINDEX, -1, "NA", "temporary variable", true, true );
		aDescriptor->setEmpty();
		(*i)->unsetVisited();
	}

	// search for components using Depth First Search algorithm on graphs
	// (marking componentIndex) and save the number of atoms in each component in componentSizes

	//cout << "Molecule::markFragments: search for components with DFS" << endl;


	bool flagAllMarked = false;
	int componentIndex = 0;

	// as long as not all atoms where marked:
	while( flagAllMarked == false ){
		Atom* startAtom = NULL;
		// find start atom
		flagAllMarked = true;
		vector< Atom* >::iterator i;
		for( i = atoms.begin(); i != atoms.end(); i++ ){
			if( (*i)->getIntDescriptor( COMPONENTINDEX )->isEmpty() ){
				startAtom = (*i);
				//cout << "  --> atom " << startAtom->toStringShort() << " was not marked by DFS" << endl;
				flagAllMarked = false;
			}
		}

		if( flagAllMarked == false){
			int componentSize = 0; // contains the number of atoms in the component at the end of DFS
			componentSize = DFS( (Atom*) startAtom, COMPONENTINDEX, componentIndex );
			//                           start    ,desriptor for marking, number to mark
			//cout << "Molecule::markFragments: componentSize :" << componentSize << endl;
			//cout << componentSizes[0] << endl;

			componentSizes[ componentIndex ] = componentSize;

			componentIndex ++;
		}
	}
	//cout << "Molecule::markFragments: done" << endl;
	return( componentIndex );
}

void Molecule::unmarkFragments(){
	// remove componentIndex intDescriptor of all visible and hidden atoms

	vector< Atom* >::iterator i;

	for( i = atoms.begin(); i != atoms.end(); i++ ){
		(*i)->deleteDescriptor( COMPONENTINDEX, false );
	}
	for( i = hiddenAtoms.begin(); i != hiddenAtoms.end(); i++ ){
		(*i)->deleteDescriptor( COMPONENTINDEX, false );
	}

	componentSizes.clear();
}

void Molecule::hideAllFragmentsBut( int aFragmentNumber ){
	//cout << "Molecule::hideAllFragmentsBut fragment " << aFragmentNumber << endl;
	//cout << "   this molecule has " << componentSizes.size() << " fragments" << endl;

	map< int, int >::iterator j;
	for( j = componentSizes.begin(); j != componentSizes.end(); j++ ){
		if( (*j).first != aFragmentNumber ){
			hideAtomsByIntDescriptor( COMPONENTINDEX, (*j).first, false );
			//componentSizes.erase( j );
			//j--;
		}
	}

}

void Molecule::writeFragments( ofstream* outStream ){

	//cout << "Molecule::writeFragments" << endl;
	map< int, int >::iterator j;
	int i = 1;
	for( j = componentSizes.begin(); j != componentSizes.end(); j++ ){
		//describe();

		hideAllFragmentsBut( (*j).first );
		string molName = getName();
		stringstream newName;
		newName << molName << "." << i;
		setName( newName.str() );
		MoleculeUtils::writeMDLHeaderBlock( *this, *outStream );
		MoleculeUtils::writeMDLCtabBlock( *this, *outStream );
		MoleculeUtils::writeMDLNSDBlock( *this, *outStream );
		setName( molName );
		restoreHiddenAtoms( false );
		i++;
	}
}


/** hide all but the biggest (in atom numbers) connected atoms in the molecule.
	Returns the number of components removed
	WARNING: if two compounds have identical maximum size.
		 it keep the compounds with the highest number of Carbons + Nitrogens
	WARNING: deletes all hidden atoms
	*/
int Molecule::hideSalts( stringstream* out ){

	eraseHiddenAtoms();	// delete hidden atoms.
	markFragments();	// mark fragments

	// find the two biggest components
	map< int, int >::iterator j;
	int biggestComponent = 0;
	int secondBiggestComponent = 0;
	int biggestComponentValue = 0;
	int secondBiggestComponentValue = 0;

	for( j = componentSizes.begin(); j != componentSizes.end(); j++ ){
		//cout << (*j).first << ", " << (*j).second << endl;
		if( (*j).second >= biggestComponentValue ){
			secondBiggestComponent = biggestComponent;
			secondBiggestComponentValue = biggestComponentValue;
			biggestComponent = (*j).first;
			biggestComponentValue = (*j).second;
			//cout << " this one";
		}
	}

	// if the two biggest compounds have an identical number of atoms
	// keep the one with the highest number of carbons + nitrogens (hydrogens might be missing)

	if( biggestComponentValue == secondBiggestComponentValue ){
		if( getNumCarbonsOfComponent( COMPONENTINDEX, biggestComponent )
			+ getNumNitrogensOfComponent( COMPONENTINDEX, biggestComponent )
			>
			getNumCarbonsOfComponent( COMPONENTINDEX, secondBiggestComponent )
			+ getNumNitrogensOfComponent( COMPONENTINDEX, secondBiggestComponent )
			){
		}else{
			int toto = biggestComponent;
			int totoValue = biggestComponentValue;

			biggestComponent = secondBiggestComponent;
			biggestComponentValue = secondBiggestComponentValue;

			secondBiggestComponent = toto;
			secondBiggestComponentValue = totoValue;
		}
	}

	// hide all other components and associated bonds
	hideAllFragmentsBut( biggestComponent );

	if( componentSizes.size() - 1 > 0 ){
		if( biggestComponentValue - secondBiggestComponentValue > 7 || secondBiggestComponentValue == 1 ){
			*out << getName() << ";" << componentSizes.size() - 1 << ";"
			<< biggestComponentValue << ";" << secondBiggestComponentValue << ";"
			<< biggestComponentValue - secondBiggestComponentValue << ";"
			<< numHiddenAtoms() << ";" << endl;

			setStringDescriptor( "saltFilterWarning", "FILTERED", "NA", "", true, true );
		}else{
			*out << getName() << ";" << componentSizes.size() - 1 << ";"
			<< biggestComponentValue << ";" << secondBiggestComponentValue << ";"
			<< biggestComponentValue - secondBiggestComponentValue << ";"
			<< numHiddenAtoms() << ";WARNING" << endl;

			setStringDescriptor( "saltFilterWarning", "WARNING", "NA", "", true, true );
		}
	}else{
 			setStringDescriptor( "saltFilterWarning", "UNFILTERED", "NA", "", true, true );
	}

	// document actions
	setIntDescriptor( "removedCompounds", componentSizes.size() - 1, "NA", "", true, true );

	unmarkFragments();

	return( componentSizes.size() - 1 );
}


/** This function performs a DFS search on the graph starting at startAtom,
		marking atoms with markValue in descriptor intDescriptorName,
		and returning the number of nodes in the connected component.
*/
int Molecule::DFS( Atom* startAtom, string intDescriptorName, int markValue ){

	int markedAtoms = 0;

	// mark node as visited
	startAtom->setVisited();

	//cout << "  now at " << startAtom->toStringShort();

	// find next unvisited Atom


	int end = 0;
	Atom* nextAtom;
	while( end == 0 ){
		try{
			nextAtom = startAtom->nextUnvisitedAtom();
			markedAtoms = markedAtoms + DFS( nextAtom, intDescriptorName, markValue );
			//nextAtom = startAtom->nextUnvisitedAtom();
		}catch( CError ){
			end = 1;
			break;
		}
	}
	//cout << " all atoms of this branch visited " << endl;

	startAtom->setIntDescriptor(intDescriptorName, markValue, "", "", true, true);
	return( 1 + markedAtoms );

}

/** hide all atoms (and associated edges) with intDescriptor aDescriptorName == aValue
Returns the number of hidden atoms */
int Molecule::hideAtomsByIntDescriptor( string aDescriptorName, int aValue, bool flagRefreshBonds ){

	//cout << "Molecule::hideAtomsByIntDescriptor: hiding all atoms with " << aDescriptorName << " = " << aValue << endl;

	vector< Atom* >::iterator i;

	int c = 0;

	for( i = atoms.begin(); i != atoms.end(); i++ ){
		if( (*i)->getIntDescriptor( aDescriptorName )->getValue() == aValue ){

			//cout << " found atom " << (*i)->toStringShort() << endl;

			hideAtom( i );

			//hiddenAtoms.push_back( *i );
			//atoms.erase( i );

			c++;
			i--;
		}
	}

	if( flagRefreshBonds == true ){
		//cout << "REFRESHING BONDS" << endl;
		refreshBonds();
	}

	return( c );

}

/** returns the number of carbons in a connected component graph defined with intDescriptor
		aDescriptorName having value aValue
*/
int Molecule::getNumCarbonsOfComponent( string aDescriptorName, int aValue ){

	vector< Atom* >::iterator i;
	int c = 0;

	for( i = atoms.begin(); i != atoms.end(); i++ ){
		if( (*i)->getIntDescriptor( aDescriptorName )->getValue() == aValue ){
			if( (*i)->getElementSymbol() == "C" ){
				c++;
			}
		}
	}
	return( c );
}

/** returns the number of nitrogens in a connected component graph defined with intDescriptor
		aDescriptorName having value aValue
*/
int Molecule::getNumNitrogensOfComponent( string aDescriptorName, int aValue ){

	vector< Atom* >::iterator i;
	int c = 0;

	for( i = atoms.begin(); i != atoms.end(); i++ ){
		if( (*i)->getIntDescriptor( aDescriptorName )->getValue() == aValue ){
			if( (*i)->getElementSymbol() == "N" ){
				c++;
			}
		}
	}
	return( c );
}


/** hides an atom (but not the edges to this node). call refreshBonds after hiding all atoms.
*/
void Molecule::hideAtom( vector< Atom* >::iterator anAtomI ){
//	cout << "Molecule::hideAtom " << (*anAtomI)->toStringShort() << endl;
	hiddenAtoms.push_back( *anAtomI );
	atoms.erase( anAtomI );
}

void Molecule::hideAtomAndToFromBonds( vector< Atom* >::iterator anAtomI ){
	cout << "Molecule::hideAtomAndToFromBonds (iterator), hiding " << (*anAtomI)->toStringShort() << endl;
	cout << "1" << endl;
	(*anAtomI)->hideAllToFromBonds();
	cout << "2" << endl;
	hiddenAtoms.push_back( *anAtomI );
	cout << "3" << endl;
	//atoms.erase( anAtomI );
	cout << "4" << endl;
}

void Molecule::hideAtomAndToFromBonds( Atom* anAtom ){
//	cout << "Molecule::hideAtom " << (*anAtomI)->toStringShort() << endl;
	anAtom->hideAllToFromBonds();
	hiddenAtoms.push_back( anAtom );
	eraseAtom( anAtom );
}

void Molecule::eraseAtom( Atom* anAtom ) throw( CError ){
	vector<Atom*>::iterator ai;
	for( ai = atoms.begin(); ai != atoms.end(); ai++){
		if( (*ai) == anAtom ){
			atoms.erase( ai );
			return;
		}
	}
	stringstream out;
	out << "Atom " << anAtom->toString() << " does not exist in molecule " << toString();
	CError e( ATOMNOTFOUND, out.str() );
	e.describe();
	throw(e);
}


/** hides all bonds to and from hidden atoms.
		Call this function after hiding all atoms in a molecule.
		Returns the number of edges hidden */
int Molecule::refreshBonds(){
	// iterate over all bonds of all atoms. If the destination atom is hidden then remove this bond

	int c = 0;
	vector< Atom* >::iterator n;
	map< Atom*, Bond* >::iterator e;

	for( n = atoms.begin(); n != atoms.end(); n++ ){
		for( e = (*n)->beginBond(); e != (*n)->endBond(); e++){
			if( isHiddenAtom( (*e).second->getTarget() ) && ! isHiddenAtom( (*e).second->getSource() ) ){
				// hide the edges only from non hidden atoms (so in case of restore the edges can be restored

				(*e).first->hideBond( e );
				e--;

				c++;
			}
		}
	}

	return( c );
}

/** returns true if an atom is hidden
*/
bool Molecule::isHiddenAtom( Atom* anAtom ){
	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		if( anAtom == *n ){
			return( false );
		}
	}
	return( true );
}
/** resets the calculations of the morgan indices. CALL THIS FUNCTION EACH TIME THE MOLECULE IS MODIFIED!
*/
void Molecule::resetMorganIndex(){

	maxMorganIteration = -1;

	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		(*n)->resetMorganIndex();
	}
}
/** returns the number of distinct morgan indices (of order anOrder)  present in the molecule */
int Molecule::getNumberOfDistinctMorganIndices( int anOrder ){
	map< int, int > morganIs;

	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		morganIs[ (*n)->getMorganIndex( anOrder ) ] += 1;
	}

	return( morganIs.size() );
}
/** returns the iteration of the Morgan index at which the different connectivity values have reached a maximum */
int Molecule::getMaxMorganIteration(){
	if( maxMorganIteration == -1 ){

		int i = 1;
		bool found = false;
		int lastN = 0;
		int N = 0;
		while( found == false ){
			N = getNumberOfDistinctMorganIndices( i );

			if( N == lastN ){
				found = true;
				return( i-1 );
			}else{
				lastN = N;
			}

			i++;
		}

	}else{
		return( maxMorganIteration );
	}

}

/** sets the morganLabel of atoms to the concatenation of the atomic symbol and the Morgan index of anOrder iteration */
void Molecule::setMorganLabels( int anOrder ){
	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		(*n)->setMorganLabel( anOrder );
	}
}

void Molecule::setPerretLabels(){
	#ifdef DEBUGPERRETLABEL
		cout << "Molecule::setPerretLabels for " << toStringShort() << endl;
		cout << "  molecule has " << sssr.size() << " rings" << endl;
	#endif
	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		(*n)->setPerretLabel();
	}
}

/** sets the uniqueMorganIndex of each atom to the Morgan index having
		the maximum of different connectivity values for the molecule
*/
int Molecule::setUniqueMorganIndices(){

	int anOrder = getMaxMorganIteration();

	vector< Atom* >::iterator n;
	for( n = atoms.begin(); n != atoms.end(); n++ ){
		(*n)->setUniqueMorganIndex( anOrder );
	}
	return( anOrder );
}

double Molecule::getSelfKernel() throw( CError ){
	if( selfKernelCalculated == false ){
		// kernel was not computed, emit error
		stringstream out;
		out << "Molecule::getSelfKernel::self kernel was not calculated in " << toString();
		CError e = CError(NOTCALCULATED, out.str() );
		e.describe();
		throw(e);
	}
	return( selfKernel );
}

/** returns the sum of transition probabilities. Should be used with fused graph. */
double Molecule::sumPT(){
	// add the PT of each bond of each atom
	vector<Atom*>::iterator ai;
	map<Atom*, Bond*>::iterator bi;

	double total = 0;

	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		// source atom is (*ai)
		for( bi = (*ai)->beginBond(); bi != (*ai)->endBond(); bi++){
			// target atom is (*bi).first
			total += (*bi).second->getKashimaPT();
		}
	}
	return( total );
}

/** function used to compute the graph kernel using the fused graph / sum of the powers of the transition matrix approach. See MoleculeUtils::powerKernel */
double Molecule::sumProbabilities(){

	// add the PT of each bond of each atom
	vector<Atom*>::iterator ai;
	map<Atom*, Bond*>::iterator bi;

	double total = 0;

	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		// source atom is (*ai)
		for( bi = (*ai)->beginBond(); bi != (*ai)->endBond(); bi++){
			// target atom is (*bi).first
			total += ( (*ai)->getKashimaPS() * (*bi).first->getKashimaPQ() * (*bi).second->getKashimaPT() );
		}
	}
	return( total );

}

double Molecule::sumProbabilitiesFast(){
	// add the PT of each bond of each atom
	//cout << "Molecule::sumProbabilitiesFast" << endl;

	map<Atom*, map<Atom*, double>* >::iterator ptli;
	map< Atom*, double>::iterator pti;

	double total = 0;

	//int i = 0;
	//int j = 0;
	for( ptli = fastPT->begin(); ptli != fastPT->end(); ptli++ ){
		//cout << i << endl;
		//i++;
		// source atom is (*ai)
		//j = 0;
		for( pti = (*ptli).second->begin(); pti != (*ptli).second->end(); pti++ ){
			//cout << "  " << j << endl;
			//cout << "    ps: " << fastPS[ (*ptli).first ] << endl;
			//cout << "    pq: " << fastPQ[ (*pti).first ] << endl;
			//cout << "    pt: " << (*pti).second << endl;
			//j++;
			// target atom is (*bi).first
			total += ( fastPS[ (*ptli).first ] * fastPQ[ (*pti).first ] * (*pti).second );
			//cout << "    total: " << total << endl;
		}
	}
	return( total );
}

double Molecule::sumPQPS(){

	// add the PT of each bond of each atom
	vector<Atom*>::iterator ai;
	map<Atom*, Bond*>::iterator bi;

	double total = 0;

	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		// source atom is (*ai)
		total += ( (*ai)->getKashimaPS() * (*ai)->getKashimaPQ() );
	}
	return( total );

}

double Molecule::sumPQPSFast(){
	// add the PT of each bond of each atom
	vector<Atom*>::iterator ai;
	double total = 0;
	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		// source atom is (*ai)
		total += ( fastPQ[ *ai ] * fastPS[ *ai ] );
	}
	return( total );

}



void Molecule::raisePowerFast(){
	//cout << "Molecule::raisePowerFast" << endl;

	vector<Atom*>::iterator i;
	vector<Atom*>::iterator j;
	vector<Atom*>::iterator ci;
	double aPt = 0;

	map<Atom*, double>* sources;
	map<Atom*, double>* savedSources;

	double a = 0;
	double b = 0;

	fastPTNext = new map<Atom*, map<Atom*, double>* >;

	for( i = beginAtom(); i != endAtom(); i++ ){
		sources = (*fastPT)[*i];


		if( sources != NULL ){
			//MoleculeUtils::describeMap( sources );

			//cout << "Calculating line of " << (*i)->toStringShort() << " " << sources->size() << " bonds " << endl;

			for( j = beginAtom(); j != endAtom(); j++ ){

				aPt = 0;
				for( ci = beginAtom(); ci != endAtom(); ci++){
					savedSources = (*fastPTSave)[*ci];

					if( savedSources != NULL ){
						a = 0;
						a = (*sources)[*ci];
						if( a != 0 ){
							b = 0;
							b = (*savedSources)[*j];
							aPt += a*b;
							//cout << (*i)->toStringShort() << " - " << (*ci)->toStringShort() << " * " << (*ci)->toStringShort() << " - " << (*j)->toStringShort() << ": " << a << " * " << b << " = " << a*b << " " << aPt <<endl;
						}
					}
					if( aPt != 0 ){
						map<Atom*, double>* aMap = (*fastPTNext)[*i];
						if( aMap == NULL ){
							aMap = new map<Atom*, double>;
							(*aMap)[*j] = aPt;
							(*fastPTNext)[*i] = aMap;
						}else{
							(*aMap)[*j] = aPt;
						}
					}
				}
			}
		}
	}

	// now the sparse matrix fastPT sould contain the data stored in fastPTNext

	// erase the content of fastPT
	map<Atom*, map<Atom*, double>* >::iterator k;
	for( k = fastPT->begin(); k != fastPT->end(); k++ ){
		delete (*k).second; // delete bonds
	}
	fastPT->clear();
	delete fastPT;

	// fastPT is the newly computed version
	fastPT = fastPTNext;
	fastPTNext = NULL;

}



/** delete all bonds in all atoms */
void Molecule::deleteBonds(){
	vector<Atom*>::iterator ai;
	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		(*ai)->deleteBonds();
	}
}


float Molecule::getMW( bool silentError ) throw( CError ){
	float mw = 0;
	vector<Atom*>::iterator ai;
	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		try{
			mw += (*ai)->getFloatDescriptor( AM )->getValue( true );
		}catch( CError e ){
			if( e.getType() == MISSINGDESCRIPTOR || e.getType() == ERRORNA ){
				if( silentError == false) {
					cerr << "Cannot compute molecular weight in " << getName() << endl;
					//cerr << "Molecule::getMW() : WARNING atom " << (*ai)->toStringShort() << " has no descriptor for atomic mass" << endl;
					throw( e );
				}
			}
		}
	}
	return( mw );
}





/*!
    \fn Molecule::deleteHiddenAtoms()
 */
void Molecule::deleteHiddenAtoms()
{
	// delete all hiddenbonds first
	vector<Atom*>::iterator ai;

	for( ai = atoms.begin(); ai != atoms.end(); ai++ ){
		(*ai)->deleteHiddenBonds();
	}

	for( ai = hiddenAtoms.begin(); ai != hiddenAtoms.end(); ai++ ){
		delete (*ai);
	}
	hiddenAtoms.clear();
}




// fill in vector<int> with new kind atoms
void Molecule::atomsLabelsListing(vector<string> *atomLabels ){

  vector<Atom*>::iterator anAtom;
  // for each atom :
  for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ ){
    bool present = false;
    for(int i = 0 ; i < atomLabels->size() ; i++){
      if( (*anAtom)->getMorganLabel( true ) == (*atomLabels)[i] ){
	present = true;
      }
    }
    if( present == false ){
	(*atomLabels).push_back( (*anAtom)->getMorganLabel(true) );
    }
  }
  
}



// fill in vector<int> with new kind atoms
void Molecule::atomsSymbolsListing(vector<string> *atomSymbols ){

  vector<Atom*>::iterator anAtom;
  // for each atom :
  for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ ){
    bool present = false;
    for(int i = 0 ; i < atomSymbols->size() ; i++){
      if( (*anAtom)->getSymbol() == (*atomSymbols)[i] ){
	present = true;
      }
    }
    if( present == false ){
      (*atomSymbols).push_back( (*anAtom)->getSymbol() );
    }
  }
  
}




// fill in vector<int> with new kind of bonds
void Molecule::bondsListing(vector<int> *bondTypes ){
  vector<Atom*>::iterator anAtom;
  // for each atom
  for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ )
    {
      map<Atom* , Bond*>::iterator aBond;
      // for each bond of the atom :
      for( aBond = (*anAtom)->beginBond() ; aBond != (*anAtom)->endBond() ; aBond++ ){
	bool present = false;
	for( int i = 0 ; i < bondTypes->size() ; i++ ){
	  if ( aBond->second->getLabel() == (*bondTypes)[i] ) present = true;
	}
	if( present == false ) (*bondTypes).push_back( aBond->second->getLabel() );
      }
    }
}


// add a value to the self kernel
void Molecule::addToSelfKernel(double value){ 
  selfKernel += value;
  selfKernelCalculated = true;
}


// substract a value to the self kernel
void Molecule::substractToSelfKernel(double value){ 
  selfKernel -= value;
  selfKernelCalculated = true;
}




/*!
    \fn Molecule::numAtomsNonCSkeleton()
 */
int Molecule::numAtomsNonCSkeleton()
{
	int n = 0;
	vector<Atom*>::iterator ai;
	for( ai = atoms.begin(); ai != atoms.end(); ai++ ){
		if( !(*ai)->isCSkeleton() ){
			n++;
		}
	}
	return( n );
}

void Molecule::binClassifyFromDescriptor( string descriptorName, float value, bool smallerOrEqual){
	float dvalue = 0;
	bool biggerthan = false;
	bool smallerthan = false;
	try{

		string svalue = getStringDescriptor( descriptorName, true )->getValue();

		//check if the first character of the descriptor value is > or < which is enough to set the threshold

		if( svalue.substr(0,1) == ">" ){
			biggerthan = true;
			dvalue = atof( svalue.substr(1,svalue.length()-1).c_str() );
		}else if( svalue.substr(0,1) == "<" ){
			smallerthan = true;
			dvalue = atof( svalue.substr(1,svalue.length()-1).c_str() );
		}else{
			dvalue = atof( svalue.c_str() );
		}

	}catch( CError e ){
		try{
			dvalue = getFloatDescriptor( descriptorName, true )->getValue();
		}catch( CError e ){
			try{
				dvalue = getIntDescriptor( descriptorName, true )->getValue();
			}catch( CError e ){
				return;
			}
		}
	}

	cout << "comparing ";
	if( biggerthan == true){
		cout << ">";
	}
	if( smallerthan == true){
		cout << "<";
	}
	cout << dvalue << " with " << value << " setting activity to ";

	if( smallerOrEqual == true ){
		if( dvalue  <= value ){
			if( biggerthan == false ){
				cout << "true" << endl;
				setActivity( 1 );
			}else{
				if( dvalue == value ){
					cout << "false" << endl;
					setActivity( 0 );
				}else{
					cout << "nothing" << endl;
					unsetActivity();
				}
			}
		}else{
			cout << "false" << endl;
			setActivity( 0 );
		}
	}else{
		if(dvalue >= value ){
			if( smallerthan == false ){
				cout << "true" << endl;
				setActivity( 1 );
			}else{
				if( dvalue == value ){
					cout << "false" << endl;
					setActivity( 0 );
				}else{
					cout << "nothing" << endl;
					unsetActivity();
				}
			}
		}else{
			cout << "false" << endl;
			setActivity( 0 );
		}
	}
	return;
}


void Molecule::setSortDescriptor( string aName, int aType ){
	if( aType == INTEGER ){
		if( !hasIntDescriptor( aName ) ){
			stringstream out;
			out << "Molecule has no INT descriptor called " << aName ;
			CError e( MISSINGDESCRIPTOR, out.str() );
			e.describe();
			throw(e);
		}
	}else if( aType == FLOAT ){
		if( !hasFloatDescriptor( aName ) ){
			stringstream out;
			out << "Molecule has no FLOAT descriptor called " << aName ;
			CError e( MISSINGDESCRIPTOR, out.str() );
			e.describe();
			throw(e);
		}
	}else{
		if( !hasFloatDescriptor( aName ) ){
			stringstream out;
			out << "Molecule has no STRING descriptor called " << aName ;
			CError e( MISSINGDESCRIPTOR, out.str() );
			e.describe();
			throw(e);
		}
	}

	sortDescriptorType = aType;
	sortDescriptorName = aName;
}


int Molecule::detectSSSR(){

	if( hasSSSRDetected() ){
		return( sssr.size() );
	}else{

	cout << toStringShort() << "::detectSSSR()" << endl;

	#ifdef DEBUGSSSR
		cout << "a3h detecting SSSR" << endl;
	#endif

	deleteHiddenAtoms();

	// fullset is atoms
	vector<Atom*> trimset;

	// store broken bonds
	Bond* brokenBond = NULL;

	#ifdef DEBUGSSSR
	cout << "a3h 1" << endl;
	#endif

	// add nodes of degree 0 to trimset

	vector<Atom*> substractSet;
	vector<Atom*>::iterator ai;
	for ( ai = beginAtom(); ai != endAtom(); ai++ ){
		substractSet.push_back( (*ai) );
	}

	long finishedAt = atoms.size();

	while( trimset.size() != finishedAt ){
		// while trimset != fullset

		#ifdef DEBUGSSSR
		cout << "a3h trimming nodes with degree 0" << endl;
		#endif
		int i = 0;
		//for ( ai = beginAtom(); ai != endAtom(); ai++ ){
		for ( ai = substractSet.begin(); ai != substractSet.end(); ai++ ){
			if( (*ai)->degree() == 0 ){
				#ifdef DEBUGSSSR
				cout << "a3hb " << (*ai)->toStringShort() << " has degree 0, adding to trimset" << endl;
				#endif

				trimset.push_back( (*ai) );
				//hideAtomAndToFromBonds( (*ai) );
				i++;
			}
		}
		#ifdef DEBUGSSSR
		cout << "a3h trimmed " << i << " nodes " << endl;
		#endif

		// add nodes of degree N2 in fullset but not in trimset to N2
		vector<Atom*> nodesN2;
		vector<Atom*> nodesNLarge;

		MoleculeUtils::substractSet( &atoms, &trimset, &substractSet );
		#ifdef DEBUGSSSR
		cout << "a3hb atoms is " << atoms.size() << " nodes " << endl;
		cout << "a3hb trimset is " << trimset.size() << " nodes " << endl;
		for ( ai = trimset.begin(); ai != trimset.end(); ai++ ){
			cout << "a3hb     " << (*ai)->toStringShort() << endl;
		}
		cout << "a3hb substractSet is " << substractSet.size() << " nodes " << endl;
		for ( ai = substractSet.begin(); ai != substractSet.end(); ai++ ){
			cout << "a3hb     " << (*ai)->toStringShort() << endl;
		}
		#endif

		if( substractSet.size() == 0 ){
			#ifdef DEBUGSSSR
			cout << "a3hb BREAKING" << endl;
			#endif
			break;
		}


		int minDegree = 10000;
		for ( ai = substractSet.begin(); ai != substractSet.end(); ai++ ){
			int degree = (*ai)->degree();
			if( degree == 2 ){
				nodesN2.push_back((*ai));
			}else if( degree > 2 ){
				nodesNLarge.push_back((*ai));
			}
			if( degree < minDegree ){
				minDegree = degree;
			}
		}

		#ifdef DEBUGSSSR
		cout << "a3h found " << nodesN2.size() << " N2 nodes, mindegree is " << minDegree << endl;
		#endif

		// find a node init in fullset but not in trimset having minimum degree
		Atom* init;
		for ( ai = substractSet.begin(); ai != substractSet.end(); ai++ ){
			if( (*ai)->degree() == minDegree ){
				init = (*ai);
				break;
			}
		}

		#ifdef DEBUGSSSR
		cout << "a3h initialising search with " << init->toStringShort() << endl;
		#endif

		if( minDegree == 0 ){
			setIntDescriptor( "sssrSize", sssr.size(), "NA", "", true, true );
			return( sssr.size() );
		}else if( minDegree == 1 ){
			/// if minimum degree == 1 trim(init) and add to trimset
			#ifdef DEBUGSSSR
				cout << "trimming atom " << init->toStringShort() << endl;
			#endif
			init->hideAllToFromBonds();
			//trimset.push_back( init );  // not needed, will be added as node of degree 0.
		}else if( minDegree == 2){
		// elsif minimum degree == 2
			//vector<Atom*>::iterator ci;
			//int gg = 0;
			//for( ci = nodesNLarge.begin(); ci != nodesNLarge.end(); ci++){

				//(*ci)->hideAllToFromBonds();
				//#ifdef DEBUGSSSR
				//	cout << "rrrty    hiding all to from bonds in atom " << (*ci)->toStringShort() << endl;
				//	cout << gg << endl;
				//	gg++;
				//#endif

			vector<Bond*> dBonds;

			bool dring = false;
			do{


				for ( ai = nodesN2.begin(); ai != nodesN2.end(); ai++ ){
					//if( (*ai) != (*ci) ){
				// for each node i in nodesN2
					#ifdef DEBUGSSSR
					cout << "a3hb   examining N2 atom " << (*ai)-> toStringShort() << endl;
					#endif

					if( (*ai)->getBFSVectorSize() == 0 ){
						#ifdef DEBUGSSSR
						cout << "a3h    starting getRingBFS" << endl;
						#endif
						Ring* ring;

						vector<Atom*>* toVisit = new vector<Atom*>;
						vector<Bond*>* toVisitBond = new vector<Bond*>;

						try{

							ring = (*ai)->getRingBFS( toVisit, toVisitBond );  // ring = getRing(i)

							delete toVisit;
							delete toVisitBond;

							// no error, so this atom is part of a ring


							// if ring.size > 0 then
							// check sssr for a duplicate, if no duplicate add ring to sssr
							#ifdef DEBUGSSSR
							cout << "a3h    -> found ring of size " << ring->size() << endl;
							#endif


							if( sssr.size() == 0 ){
								ring->setID( 1 );
								sssr.push_back( ring );										//rememberNodes.push_back( (*ai) );

								#ifdef DEBUGSSSR
									cout << "a3hc describing bonds in new ring: (" << ring->getBonds()->size() << ")" << endl;
									vector<Bond*>::iterator bb;

									for( bb = ring->getBonds()->begin(); bb != ring->getBonds()->end(); bb++ ){
									cout <<  "a3hc  > " << (*bb)->toStringShort() << endl;
									}

									cout << "a3hb found " << ring->toString() << ", added to sssr set " << endl;
								#endif
							}else{
								if( !hasRing( ring, true ) ){
									ring->setID( sssr.size() + 1 );
									sssr.push_back( ring );

									//rememberNodes.push_back( (*ai) );

									#ifdef DEBUGSSSR
										cout << "a3hc describing bonds in new ring: " << endl;
										vector<Bond*>::iterator bb;

										for( bb = ring->getBonds()->begin(); bb != ring->getBonds()->end(); bb++ ){
											cout <<  "  > " << (*bb)->toString() << endl;
										}

										cout << "a3hb found " << ring->toString() << ", added to sssr set " << endl;
									#endif
								}else{
									#ifdef DEBUGSSSR
									cout << "a3hb found " << ring->toString() << ", but this ring already exists in sssr set" << endl;
									#endif
									delete ring;
								}
							}


						}catch( CError e ){
							// This N2 atom is part of no ring

							delete toVisit;
							delete toVisitBond;

							if( e.getType() == ATOMNOTINRING ){
								// this atom is not member of a ring
								#ifdef DEBUGSSSR
								cout << "a3h    getRing returned no ring" << endl;
								#endif
							}else{
								throw( e );
							}
						}


						// resetBFSSearch
						vector<Atom*>::iterator ai;
						for ( ai = beginAtom(); ai != endAtom(); ai++ ){
							(*ai)->resetBFSVector();
						}


					}else{
						#ifdef DEBUGSSSR
						cout << "a3h     this node was already visited, skipping" << endl;
						#endif
					}//}
				}

				if( dring == false ){

				//if( dBonds.size() == 0 && dring == false){

					#ifdef DEBUGSSSR
						cout << "DBOND: first visit " << dring << endl;
					#endif

					// just found the easy rings, search for the difficult now
					// find one bond for each N2 chain and store it in dBonds
					// for simplicity we will hide one bond of all N2 nodes in all easy rings
					// for all rings
					vector<Ring*>::iterator dri;
					for( dri = sssr.begin(); dri != sssr.end(); dri++ ){
						// for all n2 atoms
						vector<Bond*>* aRingBonds = (*dri)->getBonds();
						vector<Bond*>::iterator dbi;
						for( dbi = aRingBonds->begin(); dbi != aRingBonds->end(); dbi++ ){


							if( (*dbi)->getSource()->numBonds() == 2 || (*dbi)->getTarget()->numBonds() == 2 ){
								if( !((*dbi)->getSource()->numBonds() == 2 && (*dbi)->getTarget()->numBonds() == 2 )){



									if( !MoleculeUtils::atomVectorHas( &trimset, (*dbi)->getSource() ) && !MoleculeUtils::atomVectorHas( &trimset, (*dbi)->getTarget() ) ){
										dBonds.push_back(*dbi);
									}
								}
							}

						}
					}

					#ifdef DEBUGSSSR
						cout << "DBOND: dBonds now has size " << dBonds.size() << endl;
					#endif

					if( dBonds.size() != 0 ){

						#ifdef DEBUGSSSR
							cout << "DBOND: hiding bond " <<
							 (*dBonds.begin())->toStringShort()
							 << endl;
						#endif

						// hide the first dBond and repeat the ring search
						(*dBonds.begin())->hideToFrom();
					}

					dring = true;

				}else{
					// restore the hidden dBond
					#ifdef DEBUGSSSR
						cout << "DBOND: restoring " << (*dBonds.begin())->toStringShort() << endl;
					#endif
					(*dBonds.begin())->restoreToFrom();
					dBonds.erase( dBonds.begin() );
					// hide the next dBond and repeat the ring search
					(*dBonds.begin())->hideToFrom();
					#ifdef DEBUGSSSR
						cout << "DBOND: dBonds now has size " << dBonds.size() << endl;
					#endif

				}


			}while( dBonds.size() > 0 );

			#ifdef DEBUGSSSR
				cout << "DBOND:  FINISHED" << endl;
			#endif


				//if( rememberNodes.size() == 0 ){
				//	rememberNodes.push_back( (*nodesN2.begin()) );
				//}

				/*vector<Atom*>::iterator rai;
				for( rai = rememberNodes.begin(); rai != rememberNodes.end(); rai++ ){
					#ifdef DEBUGSSSR
						cout << "hidingToFromBonds in atom " << (*rai)->toStringShort();
					#endif
					(*rai)->hideAllToFromBonds();
				}*/


			/*	// resetBFSSearch
				for ( ai = beginAtom(); ai != endAtom(); ai++ ){
					(*ai)->resetBFSVector();
				}



				#ifdef DEBUGSSSR
					(*ci)->describe();
					cout << "rrrty: restoring broken bond " << (*ci)->toStringShort() << endl;
				#endif
				(*ci)->restoreHiddenBonds();

				#ifdef DEBUGSSSR
					cout << "rrrty: done" << endl;
				#endif*/


			//}

			if( brokenBond != NULL ){
				#ifdef DEBUGSSSR
					cout << "rrrtz: restoring broken bond " << brokenBond->toStringShort() << " ";
				#endif
				brokenBond->getSource()->restoreHiddenBond( brokenBond->getTarget() );
					cout << "2" << endl;
				brokenBond->getTarget()->restoreHiddenBond( brokenBond->getSource() );
				brokenBond = NULL;
			}


			// all N2 atoms were tested for their ring membership, now break
			// one bond for each N2 chain.
			// for each chain of N2 nodes isolate one N2 node and break one bond
			#ifdef DEBUGSSSR
				cout << "a3h     now checking borders " << endl;
			#endif
			for ( ai = nodesN2.begin(); ai != nodesN2.end(); ai++ ){
				#ifdef DEBUGSSSR
					cout << "a3h    removing " << (*ai)->toStringShort() << endl;
				#endif
				hideAtomAndToFromBonds( (*ai) );
				//reduceDegreeToOne( (*ai) );
			}


		}else if( minDegree > 2 ){
		// elsif minimum degree == 3
			vector<Atom*>* toVisit = new vector<Atom*>;
			vector<Bond*>* toVisitBond = new vector<Bond*>;

			try{
				Ring* ring = (*ai)->getRingBFS( toVisit, toVisitBond );  // ring = getRing(i)
				delete toVisit;
				delete toVisitBond;

				if( !hasRing( ring, true ) ){
					ring->setID( sssr.size() + 1 );
					sssr.push_back( ring );
				}

				// resetBFSSearch
				vector<Atom*>::iterator ai;
				for ( ai = beginAtom(); ai != endAtom(); ai++ ){
					(*ai)->resetBFSVector();
				}

				//checkedges( init, ringSet)
				trimset.push_back( init );
				brokenBond = checkEdges( init );

			}catch( CError e ){

				delete toVisit;
				delete toVisitBond;

				// resetBFSSearch
				vector<Atom*>::iterator ai;
				for ( ai = beginAtom(); ai != endAtom(); ai++ ){
					(*ai)->resetBFSVector();
				}

				if( e.getType() == ATOMNOTINRING ){
					// this atom is not member of a ring
					#ifdef DEBUGSSSR
					cout << "a3h    getRing returned no ring" << endl;
					#endif
				}else{
					throw( e );
				}
			}

		}

	}

	#ifdef DEBUGSSSR
	cout << "a3h     restoring hidden atoms and bonds " << endl;
	#endif
	restoreHiddenAtoms( true ); // restore hidden atoms and bonds

	// add the ring information to member atoms and member bonds
	vector<Ring*>::iterator ri;
	for( ri = sssr.begin(); ri != sssr.end(); ri++){
		vector<Atom*>::iterator mi;
		for( mi = (*ri)->begin(); mi != (*ri)->end(); mi++ ){
			(*mi)->addRing( (*ri) );
		}
		vector<Bond*>::iterator bi;
		for( bi = (*ri)->getBonds()->begin(); bi != (*ri)->getBonds()->end(); bi++ ){
			(*bi)->addRing( (*ri) );
			(*bi)->getReverse()->addRing( (*ri) );
		}

	}

	// resetBFSSearch
	for ( ai = beginAtom(); ai != endAtom(); ai++ ){
		(*ai)->resetBFSVector();
	}

	flagHasSSSRDetected = true;

	//#ifdef DEBUGSSSR
	cout << "detectSSSR() found " << sssr.size() << " rings" << endl;
	//#endif

	}

	setIntDescriptor( "sssrSize", sssr.size(), "NA","", true, true );
	return( sssr.size() );
}


Bond* Molecule::checkEdges( Atom* anAtom ){
	// /todo implement Molecule::checkEdges! used at the end of detectSSSR.

	// this varies from the original algorithm in Figueras 1996
	// the first bond of the atom is hidden.

	return( anAtom->hideToFromFirstBond() );

}

bool Molecule::hasRing( Ring* aRing, bool detectingRing ) throw( CError ) {
	if( !flagHasSSSRDetected && detectingRing == false ){
		//detectSSSR();
		CError e( SSSRNOTDETECTED, "Smallest Set of Smallest Rings was not detected before calling Molecule::hasRing( Ring* )" );
		e.describe();
		throw(e);
	}

	vector<Ring*>::iterator i;
	for( i = sssr.begin(); i != sssr.end(); i++ ){
		if( (*i)->equals( aRing ) ){
			return( true );
		}
	}
	return( false );
}

bool Molecule::hasRing() throw( CError ){
	if( !flagHasSSSRDetected ){
		//detectSSSR();
		CError e( SSSRNOTDETECTED, "Smallest Set of Smallest Rings was not detected before calling Molecule::hasRing()" );
		e.describe();
		throw(e);
	}

	if( sssr.size() > 0 ){
		return( true );
	}
	return( false );
}

Ring* Molecule::getRingWithID( int anID, bool createIfMissing) throw( CError ){
	if( !flagHasSSSRDetected ){
		//detectSSSR();
		CError e( SSSRNOTDETECTED, "Smallest Set of Smallest Rings was not detected before calling Molecule::hasRing()" );
		e.describe();
		throw(e);
	}

	if( sssr.size() > 0 ){
		// check if this ring already exist
		vector<Ring*>::iterator ri;
		for( ri = sssr.begin(); ri != sssr.end(); ri++ ){
			if( (*ri)->getID() == anID ){
				return( (*ri) );
			}
		}
		if( createIfMissing ){
			Ring* newRing = new Ring();
			newRing->setID( anID );
			sssr.push_back( newRing );
			return( newRing );
		}else{
			stringstream out;
			out << "There is no ring with id " << anID << " in molecule " << toStringShort();
			CError e( MISSINGRING, out.str() );
			e.describe();
			throw(e);
		}

	}else{
		if( createIfMissing ){
			Ring* newRing = new Ring();
			newRing->setID( anID );
			sssr.push_back( newRing );
			return( newRing );
		}else{
			stringstream out;
			out << "There is no ring with id " << anID << " in molecule " << toStringShort();
			CError e( MISSINGRING, out.str() );
			e.describe();
			throw(e);
		}
	}
}





// no totters transformation
void Molecule::noTottersTransform(){

  vector< vector<int> > coord;
  vector<Atom*>::iterator anAtom;
  map<Atom*, Bond*>::iterator aBond;
  int i = -1;

  Molecule oldMol(*this, false);
  this->erase();

  // DEBUG :
  //cout << "Nbre Atoms initiaux:  " << oldMol.numAtoms() << endl;
  //cout << "Nbre laisons initiales : " << oldMol.numBonds() << endl;


  // 1st step : create new atoms : 1 per original atom + 1 per original bond
   // for each atom in the moldeule --> create an atom
   for( anAtom = oldMol.atoms.begin(); anAtom != oldMol.atoms.end(); anAtom++ )
     {
       Atom* newAtom = new Atom( *(*anAtom) );
       newAtom->setMorganLabel((*anAtom)->getMorganLabel() );

       this->addAtom(newAtom, false, false);
       i++;
       coord.push_back(vector<int>());
       coord[i].push_back(0);  // coord of the 1st atom of the bond (=0 in the case of an atom)
       coord[i].push_back((*anAtom)->getId()); // coord of the second atom of the bond
       coord[i].push_back(0);  // bond label

       // for all the atom's bonds --> create a node
       for( aBond = (*anAtom)->beginBond(); aBond != (*anAtom)->endBond(); aBond++)
	 {
	   Atom* newAtom2 = new Atom(aBond->first->getLabel());
	   newAtom2->setAN(aBond->first->getAN());
	   newAtom2->setType(aBond->first->getType());
	   newAtom2->setMorganLabel(aBond->first->getMorganLabel() );

	   newAtom2->setKashimaPS(0.0);
	   if( aBond->first->numBonds() == 1)
	     {
	       newAtom2->setKashimaPQ(1.0);
	       // cout <<"   --> set a PQ to 1.0" << endl;
	     }
	   else
	     {
	       newAtom2->setKashimaPQ(aBond->first->getKashimaPQ());
	     }


	   this->addAtom(newAtom2,false,false);
	   i++;
	   coord.push_back(vector<int>());
	   coord[i].push_back((*anAtom)->getId());
	   coord[i].push_back(aBond->first->getId());
      	   coord[i].push_back(aBond->second->getLabel());
	 }

     }

   
   // DEBUG
   //cout << "Nbr atoms created : " << numAtoms() << endl;
   //for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ )
   //  {
   //    cout << "\t\t (morgan label after atoms creation = "<<  (*anAtom)->getMorganLabel() << " )"<<  endl;
   //  }



   // 2nd step : create the bonds from the coordinates vector
   vector<int> coord1, coord2;
   int nbLinks = 0;

   for(i = 0 ; i < coord.size() ; i++)
     {
       coord1 = coord[i];
       for(int j = 0 ; j < coord.size() ; j++)
	 {
	   coord2 = coord[j];
	   if(coord1[1] == coord2[0] && coord2[1] != coord1[0])
	     {
	       linkAtomsNoReturn(i,j,coord2[2]);
	       nbLinks++;
	     }
	 }
     }

  // cout << "Nbr bonds created : " << nbLinks << endl;


   // 3rd step : compute the transition probabilities
   for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ )
     {
       //cout << "numbre of bonds per atom = "  << (*anAtom)->numBonds() << endl;
       //cout << "\t\t (morgan label = "<<  (*anAtom)->getMorganLabel() << " )"<<  endl;
       for( aBond = (*anAtom)->beginBond(); aBond != (*anAtom)->endBond(); aBond++)
	 {
	   aBond->second->setKashimaPT( (1 - (*anAtom)->getKashimaPQ()) / (*anAtom)->numBonds() );
	 }

     }
  // cout << "\t in MOLECULE.cpp : nb atoms = " << numAtoms() << endl;
  // cout << "\t in MOLECULE.cpp : nb bonds = " << numBonds() << endl;

  // exit(0);
}


 





void Molecule::linkAtomsNoReturn( int aSource, int aTarget, int aBondLabel ) throw( CError ){
  // atoms must exists otherwise cannot link
  
  if( max(aSource, aTarget) > numAtoms()-1 ){
    stringstream out;
    out << "molecule " << toString() << " has only " << numAtoms() << " atoms, so cannot link atom " << aSource << " to " << aTarget ;
    CError e = CError(NOTENOUGHATOMSINMOLECULE, out.str() );
    e.describe();
    throw(e);
  }else{
    linkAtomsNoReturn( atoms[aSource], atoms[aTarget], aBondLabel);
  }
  
  //	resetMorganIndex();
  // WARNING !!!! 
  // --> modified by PM ??
  //   NB : consistant with linkAtomsNoReturn( Atom* aSource, Atom* aTarget, int aBondLabel) 
}




// 3D KERNEL UTILITIES FUNCTIONS
// --> matrix manipulations

void Molecule::setAdjacency(int i, int j, double value){
  (*adjacency)[i][j] = value;
}


double Molecule::getAdjacency(int i, int j){
  return ( (*adjacency)[i][j] );
}

void Molecule::setWalks(int i, int j, double value){
  (*walks)[i][j] = value;
}


double Molecule::getWalks(int i, int j){
  return ( (*walks)[i][j] );
}


void Molecule::raisePowerAdjacency(){
  //copy old walks matrix in a temp matrix
  vector< vector<double> >* oldWalks;
  oldWalks = new vector< vector<double> >;
  
  for(int i = 0 ; i < numAtoms() ; i++){
    oldWalks->push_back(vector<double>() );
    for(int j = 0 ; j < numAtoms() ; j++){
      (*oldWalks)[i].push_back( getWalks(i,j) );
    }
  }
  
  // matrix product
  for(int i = 0 ; i < numAtoms() ; i++){
    for(int j = 0 ; j < numAtoms() ; j++){
      double temp= 0.0;
      
      for(int k = 0 ; k < numAtoms() ; k++){
	temp += (*oldWalks)[i][k] * getAdjacency(k,j);
      }
      
      setWalks(i,j,temp);
    }
  }

  // free temp matrix
  for(int i = 0 ; i < numAtoms() ; i++){
    (*oldWalks)[i].clear();
  }
  //oldWalks->clear();
  delete oldWalks;

}


double Molecule::traceWalks(){
  double res = 0.0;

  for(int i = 0 ; i < numAtoms() ; i++){
    res += getWalks(i,i);
  }

  return(res);

}



void Molecule::setSelfKernel(double value){
  selfKernel = value;
  selfKernelCalculated = true;
       
}



double Molecule::traceDiagWalks(){
  double res = 0.0;
  
  for(int i = 0 ; i < numAtoms() ; i++){
    double tmp = 0.0;
    for(int j = 0 ; j < numAtoms() ; j++){
      tmp += getWalks(i,j) * getAdjacency(j,i);
    }
    res += tmp;
  }
  
  return res;
}



// transform the 2D molecule into a 3D molecule
// i.e., completely connected graphs, where edges labels are taken in {1,2,...,nBins}
// where labels correspond to a discretization of the distance range into nBins bins
void Molecule::threeDtransform(int nBins, double distMin, double distMax){

  vector<Atom*>::iterator anAtom1, anAtom2;
  double atomDist;
  int bondLabel;
  double binSize;

  deleteBonds();

  binSize = (1.0001*distMax - distMin) / nBins; 
  // NB : trick to ensure that atomDist = distMax falls into the last bin.
  // this is equivalent to use bins of the form ]A,B], where A and B are the bounds of the bins
  // (except for the 1st bin : if atomDist = distMin, atomDist falls in the 1st bin)

  for(anAtom1 = beginAtom() ; anAtom1 != endAtom() - 1 ; anAtom1 ++){
    for(anAtom2 = anAtom1 + 1 ; anAtom2 != endAtom() ; anAtom2 ++){
      atomDist = atomicDistance(*anAtom1, *anAtom2);
      bondLabel = (int)( (atomDist - distMin)/binSize ) + 1; // nb : x = (int)y --> x is the truncation the float y
      //cout << "bondLabel = "<< bondLabel << endl;

      // OLD VERSION: based on distMin and distMax computed from input datasets
      		//if(bondLabel <= 0 || bondLabel > nBins){
		//	cout << "ERROR : Molecule::threeDtransform : invalid bond label : "  << endl;
		//	cout << "aBondLabel = " << bondLabel << " (must be > 0 and <= " << nBins << ")" << endl;
      		//}
      		//linkAtomsNoReturn(*anAtom1, *anAtom2, bondLabel);
      		//linkAtomsNoReturn(*anAtom2, *anAtom1, bondLabel);
      		// // NB : linkAtoms =  link forward + backward / linkAtomsNoReturn = just link forward
      		// // HOWEVER, linAtomsNoReturn does not reset the Morgan Indices 
      		// // WARNING : need to check if this fact is due to a (previous) modification of PM (--> where?)

      // NEW VERSION: distMin and distMax defined a priori (i.e., parameters)
      //        (NB: here, not all atoms are linked but only those having a distance within the distance range)
      if(bondLabel >= 0 && bondLabel <= nBins){
	linkAtomsNoReturn(*anAtom1, *anAtom2, bondLabel);
	linkAtomsNoReturn(*anAtom2, *anAtom1, bondLabel);
	// NB : linkAtoms =  link forward + backward / linkAtomsNoReturn = just link forward
	// HOWEVER, linAtomsNoReturn does not reset the Morgan Indices 
	// WARNING : need to check if this fact is due to a (previous) modification of PM (--> where?)
      }

    }
  }

}



void Molecule::readPartialCharges(string charges){

 vector<string> atomCharges;
 StringUtils::Split(charges, ";", atomCharges);

  if(atomCharges.size() != numAtoms()){
    cout << "ERROR : Molecule.setPartialCharges" << endl;
    cout << "   --> number of charges != number of atoms" << endl;
    cout << "     - numAtoms = " << numAtoms() << " ; numCharges = "  << atomCharges.size()  <<  endl;
  }

  int i = 0;
  vector< Atom* >::iterator anAtom;
  for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ ){
    (*anAtom)->setPartialCharge(atof(atomCharges[i].c_str()));
    i++;
  }

}


void Molecule::setMorganChargesLabels(double threshold){

  vector< Atom* >::iterator anAtom;
  for( anAtom = atoms.begin(); anAtom != atoms.end(); anAtom++ ){
    (*anAtom)->setMorganChargeLabel(threshold);
  }

}




 // ******************************* //
 // **** DEPRECATED FUNCTIONS **** //
 // ****************************** //


/** function used to compute the graph kernel using the fused graph / sum of the powers of the transition matrix approach. See MoleculeUtils::powerKernel

		this function multiplies the adjacency matrix (containing the transition probabilities
		instead of 0 and 1) whith itself
*/
/*void Molecule::raisePower(){
  vector<Atom*>::iterator i;
  vector<Atom*>::iterator j;
  vector<Atom*>::iterator ci;
	Bond* newBond;
	float aPt = 0;

	hideAllBonds();

	// multiplier la matrice d'adjacence par elle-meme. L'ennui c'est qu'on a
	// pas de matrice d'adjacence. Faisons le avec les listes de pointeurs

	// remplissons les lignes de la pseudo matrice resultats
	//cout << endl;
	//int aci = 0;
	float toto = 0;
	for( i = beginAtom(); i != endAtom(); i++ ){
		//cout << aci << endl;
		//aci++;
		// remplissons les colonnes de la pseudo matrice resultats
		if( (*i)->numHiddenBonds() > 0 ){  // les noeuds sans liaisons n'ont pas de pt...
			//int acj = 0;
			for( j = beginAtom(); j != endAtom(); j++ ){
				//cout << "  " << acj << endl;
				//acj++;
				// chaque entree de la nouvelle pseudo matrice est:
				if( (*j)->numHiddenBonds() > 0 ){  // les noeuds sans liaisons n'ont pas de pt...
					aPt = 0;
					for( ci = beginAtom(); ci != endAtom(); ci++){
						toto = (*i)->getHiddenKashimaPT( (*ci) );
						if( toto != 0 ){
							aPt += toto * (*ci)->getSavedKashimaPT( (*j) );
						}
					}

					if( aPt > 0 ){
						newBond = linkAtomsNoReturn( (*i), (*j), 1 );
						newBond->setKashimaPT( aPt );
						//cout << " " << aPt;
					}
				}
			}
		}
		//cout << endl;
	}

}*/




/** hide all bonds in the molecule */
/*void Molecule::hideAllBonds(){
	vector<Atom*>::iterator ai;
	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		(*ai)->clearHiddenBonds();
		(*ai)->hiddeAllBonds();
	}
}

void Molecule::saveAllBonds(){
	vector<Atom*>::iterator ai;
	for( ai = beginAtom(); ai != endAtom(); ai++ ){
		(*ai)->saveAllBonds();
	}
} */






/*void Molecule::exportFragments( string sdFileName, int minAtoms )
{
	MoleculeSet* fragSet = new MoleculeSet();
	fragSet->addFragmentsToSet( this, minAtoms );
	fragSet->writeSD( sdFileName, false );
	delete fragSet;
}*/








