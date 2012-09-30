/****************************************************************************************
					  moleculeset.cpp 
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

//#define DEBUG 1

#include "moleculeset.h"
#include "constant.h"

#include <locale>
#include <string>
#include <algorithm>
//#include <cctype>      // old <ctype.h>

#include <moleculeutils.h>


namespace std {
     template <class charT>
       charT
       toupper (charT c, const locale& loc) ;
     template <class charT>
       charT
       tolower (charT c, const locale& loc) ;
}



struct AscendingOrder
{
     bool operator()(Molecule* const & m1, Molecule* const & m2)
     {

				string aName = m1->getSortDescriptorName();
				int aType = m1->getSortDescriptorType();


				if( aType == DEFAULT ){
					return( m1->getId() < m2->getId() );
				}else if( aType == INTEGER ){
					int v1, v2;

					try{
  					v1 = m1->getIntDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getIntDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getIntDescriptor( aName )->getValue() < m2->getIntDescriptor( aName )->getValue() );

				}else if( aType == FLOAT ){
					float v1, v2;

					try{
  					v1 = m1->getFloatDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getFloatDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getFloatDescriptor( aName )->getValue() < m2->getFloatDescriptor( aName )->getValue() );

				}else if( aType == STRING ){
					string v1, v2;

					try{
  					v1 = m1->getStringDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getStringDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getStringDescriptor( aName )->getValue() < m2->getStringDescriptor( aName )->getValue() );
				}
				return( false );
     }
};

struct DescendingOrder
{
     bool operator()(Molecule* const & m2, Molecule* const & m1)
     {

				string aName = m1->getSortDescriptorName();
				int aType = m1->getSortDescriptorType();


				if( aType == DEFAULT ){
					return( m1->getId() < m2->getId() );
				}else if( aType == INTEGER ){
					int v1, v2;

					try{
  					v1 = m1->getIntDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getIntDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getIntDescriptor( aName )->getValue() < m2->getIntDescriptor( aName )->getValue() );

				}else if( aType == FLOAT ){
					float v1, v2;

					try{
  					v1 = m1->getFloatDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getFloatDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getFloatDescriptor( aName )->getValue() < m2->getFloatDescriptor( aName )->getValue() );

				}else if( aType == STRING ){
					string v1, v2;

					try{
  					v1 = m1->getStringDescriptor( aName )->getValue();
					}catch( CError e ){
						return( false );
					}

					try{
  					v2 = m2->getStringDescriptor( aName )->getValue();
					}catch( CError e ){
						return( true );
					}

					return( v1 < v2 );
					//return( m1->getStringDescriptor( aName )->getValue() < m2->getStringDescriptor( aName )->getValue() );
				}
				return( false );
     }
};



struct AscendingMW
{
	bool operator()(Molecule* const & m1, Molecule* const & m2){
		return( m1->getMW() < m2->getMW() );
	}
};

struct AscendingNumAtoms
{
	bool operator()(Molecule* const & m1, Molecule* const & m2){
		int n1 = m1->getNumAtoms();
		int n2 = m2->getNumAtoms();
		if( n1 == n2){
			return( m1->getMW() < m2->getMW() );
		}
		return( n1 < n2 );
	}
};




MoleculeSet::MoleculeSet(){
	activitySet = false;
	gram = new vector< vector<double> >;
	gramNormal = new vector< vector<double> >;
	gramCalculated = false;
	setKashimaKernelParam( 0.1, 1000, 1 );
}


/** Copy constructor (TODO: check carrefully!)
*/
MoleculeSet::MoleculeSet(const MoleculeSet& aSet){
	activitySet = false;
	gramCalculated = false;
	pq = aSet.pq;
	gram = new vector< vector<double> >;
	gramNormal = new vector< vector<double> >;

/*	vector<Molecule*>::iterator m;
	for( m = aSet.begin(); m != aSet.end(); m++ ){	
		addMoleculeCopy( (*m) );
	}
*/	
}




MoleculeSet::~MoleculeSet(){
	delete gram;
	delete gramNormal;
}

//void MoleculeSet::addMolecule( Molecule* aMolecule, bool updateGram ){
void MoleculeSet::addMolecule( Molecule* aMolecule ){
	push_back( aMolecule );
	aMolecule->setKashimaKernelProb( pq );
	//if(gramCalculated == true && updateGram == true){
	//	completeGramMatrix();
	//}//else{
	//	gramCalculated = false;
	//}
}

/** creates a new molecule in the dataset and reads its definition from a Mol file
*/
Molecule* MoleculeSet::addSingleMOL( string aMolFile, bool genericAtomType ) {
	Molecule* molecule = new Molecule();
	molecule->readMOL( aMolFile, genericAtomType );
	this->addMolecule( molecule );
	return( molecule );
}

/** creates a new molecule in the dataset and reads its definition from a Kcf file
14 June 2004: RENAMED addKcfFile to addSingleKCF
*/
Molecule* MoleculeSet::addSingleKCF( string aKcfFile ) {
	//cout << "adding kcf file: " << aKcfFile << endl;

	KCFMolecule* molecule = new KCFMolecule();
	molecule->readKCF( aKcfFile );
	this->addMolecule( molecule );
	return( molecule );
}

//void MoleculeSet::addMutag(string aFileName, string rFileName, int numMolToRead, bool updateGram ){
void MoleculeSet::addMutag(string aFileName, string rFileName, uint numMolToRead ){

	#ifdef DEBUG
		//cout << "reading atoms and bonds from mutag in " << aFileName << endl;
	#endif

	ifstream inFile;
 	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	int linesize=512;
	int lineNumber = 0;

	char *line = new char[linesize];
	vector<string> typeWords;
	vector<string> valueWords;
	string stringLine;

	string capSymbol;

	Atom* atom;
	Molecule* molecule;

	int atomC = 0;

	uint i=0;
	while( !inFile.eof() && i < numMolToRead + 1 ){
	//while( lineNumber < 55 ){
		lineNumber++;
		inFile.getline(line, linesize-1,'\n');
		stringLine = line;

		// if line is empty or starts with // or # then skip this line
 		if(stringLine.size() == 0){
			continue;
		}
		if(stringLine.substr(0,1) == "#" || stringLine.substr(0,2) == "//"){
			continue;
		}

		StringUtils::Split(line, "(", typeWords);
		StringUtils::Split(typeWords[1], ",", valueWords);

		if(typeWords[0] == "atm"){

				try{
					// try to find the molecule in the moleculeSet
					molecule = (*this)[ valueWords[0] ];  // IMPLEMENTER L'OPERATEUR NA
				}
				catch(CError e){
					// if molecule not found in the dataset then add it
					//#ifdef DEBUG
						//cout << "adding molecule " << valueWords[0] << " to the molecule set"<< endl;
					//#endif
					molecule = new Molecule();
					molecule->setStringDescriptor( "name", valueWords[0], "", "", true, true );
					molecule->setKashimaKernelProb( pq, 1 );

					atomC = 0;

					addMolecule( molecule );
					i++;
				}


				//ToUpper      up(std::locale::classic());
				//ToLower      down(std::locale::classic());

				capSymbol = valueWords[2];
				//std::transform (capSymbol.begin(), capSymbol.begin()+1, capSymbol.begin(), up);

				#ifdef DEBUG
				//cout << "adding atom " << valueWords[1] << " of type " << capSymbol << " to molecule " << valueWords[0] << endl;
				#endif

				atomC++;
				atom = molecule->addAtom( capSymbol );
				atom->setStringDescriptor( "name", valueWords[1], "", "", true, true );
				atom->setIdInMolecule( atomC );

				//molecule->describeShort();

			}else if(typeWords[0] == "bond"){

				// check if molecule exists, eventually create it
				#ifdef DEBUG
				//cout << "adding bond of type " << valueWords[3][0] << " between atoms " << valueWords[1] << " and " << valueWords[2] << " in molecule " << valueWords[0] << endl;
				#endif

				try{
					// try to find the molecule in the moleculeSet
					molecule = (*this)[ valueWords[0] ];
				}
				catch(CError e){
					cout << "ERROR in MoleculeSet::addMutag reading file " << aFileName << endl;
					cout << "CANNOT ADD A BOND TO A NON EXISTING MOLECULE " << endl;
					exit(1);
				}

				try{
					if(atoi(&valueWords[3][0]) == 7 ){
						molecule->linkAtoms( valueWords[1], valueWords[2], 4 );
					}else{
						molecule->linkAtoms( valueWords[1], valueWords[2], atoi(&valueWords[3][0]) );
					}
				}
				catch( CError e ) {
					cout << "Error in line " << lineNumber << " of file " << aFileName << ": ";
				exit(1);
			}

		}else{
			delete[] line;
			stringstream out;
			out << "MoleculeSet::addMutag: Expecting atm() or bond() at line " << lineNumber << " in file " << aFileName;
			CError e(BADFILE, out.str() );
			e.describe();
			throw(e);
		}
		#ifdef DEBUG
			//cout << "found " << typeWords.size() << " words in line " << lineNumber << ": " << typeWords[0]<< " " << typeWords[1] << endl;
		#endif

		typeWords.clear();
		valueWords.clear();

	}

	inFile.close();

	if( numMolecules() > numMolToRead ){
		cout << "poping back " << endl; //(*end()-1)->getName() << endl;
		pop_back();
	}

	subsetStart = 0;
	subsetSize = numMolecules();

//	#ifdef DEBUG
		//cout << "reading mutag biological activity in " << rFileName << endl;
//	#endif
	if( rFileName != "" ){
		ifstream inFile2;
		inFile2.open( rFileName.c_str(), ios::in );
		if(!inFile2.good()){
			delete[] line;
			CError e = CError(FILENOTFOUND, rFileName + " file not found");
			e.describe();
			throw(e);
		}

		i=0;
		lineNumber = 0;
		while( !inFile2.eof() ){
			lineNumber++;
			inFile2.getline(line, linesize-1,'\n');
			stringLine = line;

			//cout << line << endl;

			// if line is empty or starts with // or # then skip this line
 			if(stringLine.size() == 0){
				continue;
			}
			if(stringLine.substr(0,1) == "#" || stringLine.substr(0,2) == "//"){
				continue;
			}

			StringUtils::Split(line, "(", typeWords);
			StringUtils::Split(typeWords[1], ")", valueWords);

			//cout << lineNumber << " " << typeWords[0] << ": " << valueWords[0] << endl;

			if( typeWords[0] == "active" ){
				//cout << "active" << endl;
				try{
			  	//setIntDescriptor(ACTIVITY, valueWords[0], 1);
				setActivity( valueWords[0], 1.0 );
				}
				catch(CError e){
					//e.describe();
				}
			}else{
				//cout << "inactive" << endl;
				try{
					//setIntDescriptor(ACTIVITY, valueWords[0], -1);
					setActivity( valueWords[0], -1.0 );
				}
				catch(CError e){
					//e.describe();
				}
			}

			typeWords.clear();
			valueWords.clear();

			i++;
		}

		inFile2.close();
	}

	delete[] line;

//	if(gramCalculated == true && updateGram == true){
//		completeGramMatrix();
//	}//else{
	//	gramCalculated = false;
	//}
}


bool MoleculeSet::nameExists(string aName){
	#ifdef DEBUG
		//cout << "looking for " << aName << endl;
	#endif

	vector<Molecule*>::iterator molecule;
	// for each molecule
	for( molecule = begin(); molecule != end(); molecule++ ){

		if( (*molecule)->getStringDescriptor("name")->getValue() == aName){
			return(true);
		}
	}
	return(false);
}


Molecule* MoleculeSet::operator[](int anInd) throw( CError ){
	return( getMolByIndex( anInd ) );
}

Molecule* MoleculeSet::getMolByIndex(int anInd) throw( CError ){

  Molecule *res;
  //res = this[anInd];
  res = this->at(anInd);
  return( res );
}



Molecule* MoleculeSet::operator[](string aName) throw( CError ){
	return( getMolByName( aName ) );
}

Molecule* MoleculeSet::getMolByName(string aName) throw( CError ){

	vector<Molecule*>::iterator molecule;
	int numRes = 0;
	Molecule* res;
	// for each molecule
	for( molecule = begin(); molecule != end(); molecule++ ){
		if( (*molecule)->getStringDescriptor("name")->getValue() == aName){
			res = *molecule;
			numRes++;
		}
	}

	if( numRes == 0 ){
		stringstream out;
		out << "MoleculeSet::[]: Error in [] operator for class MoleculeSet: No molecule with name " << aName << " found in the moleculeSet";
		CError e(NOTFOUND, out.str() );
		//e.describe();
		throw(e);
	}else if( numRes > 1 ){
		stringstream out;
		out << "MoleculeSet::[]: Error in [] operator for class MoleculeSet: More than one molecule with name " << aName << " were found in the moleculeSet (" << numRes << " molecules found)";
		CError e(DUPLICATEENTRIES, out.str() );
		//e.describe();
		throw(e);
	}

	return( res );
}


void MoleculeSet::setKashimaKernelParam( double aPq, int aConvergenceCondition, bool skipSkeleton ){

	//cout << "########### " << aPq << " " << pq << " " << aConvergenceCondition << " " << getConvergenceCondition() << endl;

	//if( aPq != pq || aConvergenceCondition != getConvergenceCondition() ){
		resetGramMatrix();
		resetSelfKernels();

		vector<Molecule*>::iterator molecule;
		// for each molecule
		for( molecule = begin(); molecule != end(); molecule++ ){

			//cout << "SETTING K PARAM OF " << (*molecule)->getName() << endl;
			//aPq = 1/(*molecule)->numAtoms();

			(*molecule)->setKashimaKernelProb( aPq, skipSkeleton );
		}
		pq = aPq;
		convergenceCondition = aConvergenceCondition;
	//}

}



void MoleculeSet::gramCompute( double aPq, double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ), double(*pt2AtomKernel)(Atom*, Atom*) = MoleculeUtils::atomKernelExternalMatrix, double(*pt2BondKernel)(Bond*, Bond*) = MoleculeUtils::bondKernelType, int aParameter1, string aReportFileName, int nbThreadsWanted, bool silentMode , bool filterTotters ){

	gramCompute( this, aPq, pt2GraphKernel, pt2AtomKernel, pt2BondKernel , aParameter1, 1, aReportFileName, nbThreadsWanted,  silentMode, filterTotters );
}

void MoleculeSet::gramCompute( double aPq, double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ), double(*pt2AtomKernel)(Atom*, Atom*) = MoleculeUtils::atomKernelExternalMatrix, double(*pt2BondKernel)(Bond*, Bond*) = MoleculeUtils::bondKernelType, int aParameter1, int aParameter2, string aReportFileName, int nbThreadsWanted, bool silentMode , bool filterTotters ){

	gramCompute( this, aPq, pt2GraphKernel, pt2AtomKernel, pt2BondKernel , aParameter1, aParameter2, aReportFileName, nbThreadsWanted, silentMode, filterTotters );
}


void MoleculeSet::gramCompute(	MoleculeSet* anotherSet,
				double aPq, double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ), double(*pt2AtomKernel)(Atom*, Atom*) = MoleculeUtils::atomKernelExternalMatrix, double(*pt2BondKernel)(Bond*, Bond*) = MoleculeUtils::bondKernelType, int aParameter1, string aReportFileName, int nbThreadsWanted, bool silentMode , bool filterTotters ){
	gramCompute( anotherSet, aPq, pt2GraphKernel, pt2AtomKernel, pt2BondKernel , aParameter1, 1, aReportFileName, nbThreadsWanted ,silentMode, filterTotters);

}


void MoleculeSet::gramCompute(	MoleculeSet* anotherSet,
				double aPq, double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ), double(*pt2AtomKernel)(Atom*, Atom*) = MoleculeUtils::atomKernelExternalMatrix, double(*pt2BondKernel)(Bond*, Bond*) = MoleculeUtils::bondKernelType, int aParameter1, int aParameter2, string aReportFileName, int nbThreadsWanted, bool silentMode, bool filterTotters){


	this->comparisonSet = anotherSet;

	subsetStart = 0;
	subsetSize = numMolecules();

	bool noActivity = false; // set to true uf one of the molecule has no activity

	if( !silentMode ){
		cout << "calculating the Kashima gram matrix for " << numMolecules() << " x "
		<< anotherSet->numMolecules() << " molecules" << endl;
	}

	// clear the gram matrix
	gram->clear();
	gramNormal->clear();

	// set the probabilities according to the article by Kashima et al.
	setKashimaKernelParam( aPq, aParameter1 );
	if( anotherSet != this ){
		anotherSet->setKashimaKernelParam( aPq, aParameter1 );
	}


	if(filterTotters)
	  {
// DEBUG (PM)
// 	    vector<Molecule*>::iterator mtest,mtest2;
// 	    int i = 0;
// 	    for( mtest = begin(); mtest != end(); mtest++ )
// 	      {
// 		i++;
// 		cout << "molecule no before transfo " << i << " ; nbr atoms =  "  << (*mtest)->numAtoms() << "  ; nbr edges = "  <<(*mtest)->numBonds() <<  endl;
// 	      }

	    this->noTottersTransform();
	    if( anotherSet != this){
	      anotherSet->noTottersTransform();
	    }
// DEBUG (PM)
// 	    i = 0;
// 	    for( mtest2 = begin(); mtest2 != end(); mtest2++ )
// 	      {
// 		i++;
// 		cout << "molecule no after transfo " << i << " ; nbr atoms =  "  << (*mtest2)->numAtoms()  << "  ; nbr edges = "  <<(*mtest2)->numBonds() << endl;
// 	      }
// 	    exit(1);

	  }



	// calculate Gram Matrix
	vector<Molecule*>::iterator m1;
	vector<Molecule*>::iterator m2;

	int i = 0;
	int j = 0;		// counter

	// initialize the gram matrix and calculate all Selfkernels

	i = 0;
	for( m1 = begin(); m1 != end(); m1++ ){
		gram->push_back( vector<double>() );
		gramNormal->push_back( vector<double>() );

		j = 0;
		for( m2 = anotherSet->begin(); m2 != anotherSet->end(); m2++ ){
			(*gram)[i].push_back( -1 );
			(*gramNormal)[i].push_back( -1 );
			j++;
		}
		i++;
	}

	// for each pair of molecules calculate the Kashima Kernel
	if(nbThreadsWanted > 1) {
		// MULTITHREAD CODE COMES HERE
		stringstream out;
		out << "Atom::atomKernel function is not implemented yet";
		CError e( NOTIMPLEMENTED, out.str() );
		e.describe();
		throw(e);
	}else{
		// no threads

		if( !silentMode ){
			cout << "Using a single thread" << endl;
		}

		double sk1 = -1000.0;		// self kernel mol1
		double sk2 = -1000.0;		// self kernel mol2
		double k = 0;	// kernel
		double nk = 0;						// normalised kernel

		int nbUndistinct = 0;  // number of indistinct molecules
		int nbUndistinctDiffAct = 0; // number of indistinct molecules with different activity status
		int nbUndistinctUnknownAct = 0; // number of indistinct molecules with unknown activity
		stringstream confoundedMolecules; // to report the confounded molecules

		int nbOrthogonal = 0;  // number of orthogonal molecules
		int nbOrthogonalDiffAct = 0; // number of orthogonal molecules with different activity status
		int nbOrthogonalUnknownAct = 0; // number of orthogonal molecules with unknown activity
		stringstream orthogonalMolecules; // to report the orthogonal molecules


		if( !silentMode ){
			cout << "calculating " << this->numMolecules() << " x "
			<< anotherSet->numMolecules() << endl;
		}

		if(this == anotherSet){
			// this is a symetric gram Matrix therefore only one triangular matrix needs to be computed

			i = 0;
			for( m1 = this->begin(); m1 != this->end(); m1++ ){
				if( !silentMode ){
					cout << i << " / " << subsetSize << endl;
				}

				#ifdef DEBUG
					cout << "  sk1" << endl;
				#endif
				sk1 = (*m1)->getSelfKernel(	pt2GraphKernel,
								pt2AtomKernel,
								pt2BondKernel,
								aParameter1, aParameter2 );

				j = 0;
				for( m2 = anotherSet->begin(); m2 != anotherSet->end() - ( anotherSet->numMolecules() - i - 1) ; m2++ ){
					#ifdef DEBUG
						cout << "  " << j << endl;
					#endif
          //cout << i << "/" << j << ": " << (*m1)->getName() << " " << (*m2)->getName() << endl;

					if( (*m1) == (*m2) ){
						#ifdef DEBUG
							cout << "   sk = k" << endl;
						#endif
						k = sk1; //(*m1)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, aParameter1, aParameter2 );
						sk2 = k;
					}else{
						#ifdef DEBUG
							cout << "   k" << endl;
						#endif
						k = (*m1)->computeKernel(*m2, pt2GraphKernel, pt2AtomKernel, pt2BondKernel, aParameter1, aParameter2 );
						#ifdef DEBUG
							cout << "   sk2" << endl;
						#endif
						sk2 = (*m2)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, aParameter1, aParameter2 );

					}


					double denom = sqrt( sk1 * sk2 );
					nk = k / denom;   //sqrt( sk1 * sk2 );

					if( !( denom > 0 )  ){
						if( (*m1) == (*m2) ){
							nk = 1;
						}else{
							nk = 0;
						}
					}

					//cout << sk1 << " " << sk2 << " " << nk << " " << endl;

					(*gram)[i][j] = k;
					(*gram)[j][i] = k;
					(*gramNormal)[i][j] = nk;
					(*gramNormal)[j][i] = nk;

					j++;

					if( nk > MAXSIMILARITYWARNING && m1 != m2 ){
						nbUndistinct++;
						try{
							//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == (*m2)->getIntDescriptor(ACTIVITY)->getValue() ){

							float a1 = 0;
							float a2 = 0;
							//cout << "c1h b" << endl;
							a1 = (*m1)->getActivity( true );
							a2 = (*m2)->getActivity( true );
							//cout << "c2h c" << endl;

							if( a1 == a2 ){
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " cannot be distinguished. Normalized kernel = " << nk << endl;
								}
							}else{
								nbUndistinctDiffAct++;
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " cannot be distinguished and they have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
								}
								confoundedMolecules << (*m1)->getName() << " - " << (*m2)->getName() << endl;

							}
						} catch( CError e ){
							noActivity = true;
							nbUndistinctUnknownAct++;
							continue;
						}
					}
					if( nk < MINSIMILARITYWARNING ){
						nbOrthogonal++;
						try{

							float a1 = 0;
							float a2 = 0;
							//cout << "c2h b" << endl;
							a1 = (*m1)->getActivity( true );
							a2 = (*m2)->getActivity( true );
							//cout << "c2h c" << endl;



							//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == (*m2)->getIntDescriptor(ACTIVITY)->getValue() ){
							if( a1 == a2 ){
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " are orthogonal. Normalized kernel = " << nk << endl;
								}
							}else{
								nbOrthogonalDiffAct++;
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " are orthogonal but have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
								}
								orthogonalMolecules << (*m1)->getName() << " - " << (*m2)->getName() << endl;
							}
						} catch( CError e ){
							nbOrthogonalUnknownAct++;
							noActivity = true;
							continue;
						}
					}

				}
				i++;
			}

		}else{
  		// this is not a symetric matrix, compute all

			i = 0;
			for( m1 = this->begin(); m1 != this->end(); m1++ ){
			  if( !silentMode ){
			    cout << i << " / " << subsetSize << endl;
			  }
			  
			  sk1 = (*m1)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );



				j = 0;
	  			for( m2 = anotherSet->begin(); m2 != anotherSet->end() ; m2++ ){
					if( (*m1) == (*m2) ){
						k = (*m1)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
					}else{
						k = (*m1)->computeKernel(*m2, pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
					}
          				sk2 = (*m2)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
					nk = k / sqrt(sk1 * sk2);

					(*gram)[i][j] = k;
					(*gramNormal)[i][j] = nk;

					j++;

					if( nk > MAXSIMILARITYWARNING && m1 != m2 ){
						nbUndistinct++;
						try{


							float a1 = 0;
							float a2 = 0;
							a1 = (*m1)->getActivity( true );
							a2 = (*m2)->getActivity( true );

							if( a1 == a2 ){
							//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == (*m2)->getIntDescriptor(ACTIVITY)->getValue() ){
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " cannot be distinguished. Normalized kernel = " << nk << endl;
								}
							}else{
								nbUndistinctDiffAct++;
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " cannot be distinguished and they have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
								}
								confoundedMolecules << (*m1)->getName() << " - " << (*m2)->getName() << endl;
							}
						} catch( CError e ){
							nbUndistinctUnknownAct++;
							noActivity = true;
							continue;
						}
					}

					if( nk < MINSIMILARITYWARNING ){
						nbOrthogonal++;
						try{
							//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == (*m2)->getIntDescriptor(ACTIVITY)->getValue() ){
							if( (*m1)->getActivity( true ) == (*m2)->getActivity( true ) ){
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " are orthogonal. Normalized kernel = " << nk << endl;
								}
							}else{
								nbOrthogonalDiffAct++;
								if( !silentMode ){
									cout << "WARNING: molecule " << (*m1)->getName() << " and " << (*m2)->getName() << " are orthogonal but have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
								}
								orthogonalMolecules << (*m1)->getName() << " - " << (*m2)->getName() << endl;
							}
						} catch( CError e ){
							nbOrthogonalUnknownAct++;
							noActivity = true;
							continue;
						}
					}

				}
				i++;
			}

		}

		if( !silentMode ){
			cout << endl<< "Report:" << endl;
			cout << nbUndistinct << " molecules pairs could not be distinguished using the graph kernel" << endl;
			cout << nbUndistinctDiffAct << " of them had a different biological activity" << endl;
			cout << nbUndistinctUnknownAct << " of them had unknown biological activity" << endl<<endl;

			cout << nbOrthogonal << " molecules pairs were orthogonal" << endl;
			cout << nbOrthogonalDiffAct << " of them had a different biological activity" << endl;
			cout << nbOrthogonalUnknownAct << " of them had unknown biological activity" << endl;
		}

		if( nbUndistinctDiffAct>0 && !silentMode ){
			cout << "here is a list of the undistinguished molecules with different activities:" << endl;
			cout << confoundedMolecules.str();
		}

		if( nbOrthogonalDiffAct>0 && !silentMode ){
			cout << "here is a list of the orthogonal molecules with different activities:" << endl;
			cout << orthogonalMolecules.str();
		}


		if( aReportFileName != "" ){
		// write a report in aReportFile

			fstream oFile;
 			oFile.open( aReportFileName.c_str(), ios::out );
			if(!oFile.good()){
				CError e = CError(COULDNOTOPENFILE, "MoleculeSet::kashimaGramMatrix: could not write file " + aReportFileName);
				e.describe();
				throw(e);
			}

			fstream oFileShort;
			string aReportFileNameShort = aReportFileName + ".short";
 			oFileShort.open( aReportFileNameShort.c_str(), ios::out );
			if(!oFileShort.good()){
				CError e2 = CError(COULDNOTOPENFILE, "MoleculeSet::kashimaGramMatrix: could not write file " + aReportFileName + ".short");
				e2.describe();
				throw(e2);
			}

			oFile << "Kashima Gram Kernel Matrix with pq = " << aPq << endl;
			oFile << "convergence until rlk does not change by more than 1/" << aParameter1 << endl;
			oFile << "atom kernel matrix loaded from " << elements.getAtomKernelName() << endl;
			oFile << nbUndistinct << " molecule pairs could not be distinguished using the graph kernel " << endl;
			oFile << nbUndistinctUnknownAct << " of them had unknown biological activity" << endl;

			if( hasActivity() ){
				oFileShort << nbUndistinct << ";" << nbUndistinctDiffAct << ";" << nbUndistinctUnknownAct << ";" << nbOrthogonal << ";" << nbOrthogonalDiffAct << ";" << nbOrthogonalUnknownAct << endl;

		  	oFile << nbUndistinctDiffAct << " of them had a different biological activity" << endl;
				if( nbUndistinctDiffAct>0 ){
					oFile << "here is a list of the undistinguished molecules with different activities:" << endl;
					oFile << confoundedMolecules.str();
				}
			}else{
				oFileShort << nbUndistinct << ";NA;" << nbUndistinctUnknownAct << ";" << nbOrthogonal << ";NA;" << nbOrthogonalUnknownAct << endl;
			}

			oFile << nbOrthogonal << " molecule pairs were orthogonal using the graph kernel" << endl;
			oFile << nbOrthogonalUnknownAct << " of them had unknown biological activity" << endl;
			if( hasActivity() ){
		  	oFile << nbOrthogonalDiffAct << " of them had a different biological activity" << endl;
			}

			oFile.close();
			oFileShort.close();
		}

	}

	gramCalculated = true;

}





void MoleculeSet::gramCompute3D( 
		   double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, 
					    double(*pt2AtomKernel)(Atom*, Atom*), 
					    double(*pt2BondKernel)(float,float,float), float ), 
		   double(*pt2AtomKernel)(Atom*, Atom*), 
		   double(*pt2BondKernel)(float, float, float),
		   float edgeKernelparameter,
		   bool silentMode ){

  gramCompute3D(this, pt2GraphKernel, pt2AtomKernel, pt2BondKernel, edgeKernelparameter, silentMode);


}





void MoleculeSet::gramCompute3D( 
		   MoleculeSet* anotherSet,
		   double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, 
					    double(*pt2AtomKernel)(Atom*, Atom*), 
					    double(*pt2BondKernel)(float,float,float), float ), 
		   double(*pt2AtomKernel)(Atom*, Atom*), 
		   double(*pt2BondKernel)(float, float, float),
		   float edgeKernelParameter,
		   bool silentMode ){

	this->comparisonSet = anotherSet;
	// clear the gram matrix
	gram->clear();
	gramNormal->clear();
	
	// WARNING : you MUST set these parameters
	//           otherwise, problems reading the gram matrix . WHY ??
	subsetStart = 0;
	subsetSize = numMolecules();


	//// set the probabilities according to the article by Kashima et al.
	//setKashimaKernelParam( aPq, aParameter1 );
	//if( anotherSet != this ){
	//	anotherSet->setKashimaKernelParam( aPq, aParameter1 );
	//}


	// calculate Gram Matrix //
	// ----------------------//
	vector<Molecule*>::iterator m1;
	vector<Molecule*>::iterator m2;

	int i = 0;
	int j = 0;		// counter


	// initialize the gram matrix and calculate all Selfkernels
	i = 0;
	for( m1 = begin(); m1 != end(); m1++ ){
		gram->push_back( vector<double>() );
		gramNormal->push_back( vector<double>() );

		j = 0;
		for( m2 = anotherSet->begin(); m2 != anotherSet->end(); m2++ ){
			(*gram)[i].push_back( -1 );
			(*gramNormal)[i].push_back( -1 );
			j++;
		}
		i++;
	}


	// for each pair of molecules calculate the3D Kernel
	double sk1 = -1000.0;		// self kernel mol1
	double sk2 = -1000.0;		// self kernel mol2
	double k = 0;	                // kernel
	double nk = 0;			// normalised kernel
	

	if( !silentMode ){
	  cout << "calculating " << this->numMolecules() << " x "
	       << anotherSet->numMolecules() << endl;
	}
	
	if(this == anotherSet){
	  // this is a symetric gram Matrix therefore only one triangular matrix needs to be computed
	  
	  i = 0;
	  for( m1 = this->begin(); m1 != this->end(); m1++ ){
	    if( !silentMode ){
	      cout << " i = " << i+1 << " / " << this->numMolecules() << endl;
	    }
#ifdef DEBUG
	    cout << "molecule " << i+1 << " , self kernel" << endl;
	    cout << "--------------------------"<< endl;
#endif	    
	    sk1 =  pt2GraphKernel( *m1, *m1, pt2AtomKernel, pt2BondKernel, edgeKernelParameter );
	    (*m1)->setSelfKernel(sk1);

	    j = 0;
	    for( m2 = this->begin(); m2 != this->end() - ( this->numMolecules() - i - 1) ; m2++ ){
	      
	      //cout <<"\t j = " << j+1 << endl;

	      if( (*m1) == (*m2) ){
		k = sk1; 
		sk2 = k;
	      }else{
#ifdef DEBUG
		cout << "molecule " << i+1 << " X molecule " << j+1 << endl;
		cout << "--------------------------"<< endl;
#endif	    
		k = pt2GraphKernel( *m1, *m2, pt2AtomKernel, pt2BondKernel, edgeKernelParameter );

		// NB : core modification to avoid REPETED CALCULATIONS OF SELF KERNELS
		//#ifdef DEBUG
		//cout << "molecule " << j+1 << " , self kernel" << endl;
		//cout << "--------------------------"<< endl;
		//sk2 = pt2GraphKernel( *m2, *m2, pt2AtomKernel, pt2BondKernel, atomKernelParameter, edgeKernelParameter );
		//#endif	    
		sk2 = (*m2)->getSelfKernel();
	      }
	     
	      double denom = sqrt( sk1 * sk2 );
	      nk = k / denom;
	      
	      if( !( denom > 0 )  ){
		if( (*m1) == (*m2) ){
		  nk = 1;
		}else{
		  nk = 0;
		}
	      }
	      
	      (*gram)[i][j] = k;
	      (*gram)[j][i] = k;
	      (*gramNormal)[i][j] = nk;
	      (*gramNormal)[j][i] = nk;
	      
	      j++;
	    }
	    i++;
	  }
	  

	}else{
	  // this is not a symetric matrix, compute all
	  i = 0;
	  for( m1 = this->begin(); m1 != this->end(); m1++ ){
	    //sk1 = (*m1)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
	    sk1 =  pt2GraphKernel( *m1, *m1, pt2AtomKernel, pt2BondKernel, edgeKernelParameter );
	    (*m1)->setSelfKernel(sk1);

	    if( !silentMode ){
	      cout << " i = " << i+1 << " / " << this->numMolecules() << endl;
	    }
	    j = 0;
	    for( m2 = anotherSet->begin(); m2 != anotherSet->end() ; m2++ ){

	      //cout <<"\t j = " << j+1 << endl;

	      //k = (*m1)->computeKernel(*m2, pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
	      k = pt2GraphKernel( *m1, *m2, pt2AtomKernel, pt2BondKernel, edgeKernelParameter );
	      //sk2 = (*m2)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, getConvergenceCondition(), aParameter2 );
	      sk2 = pt2GraphKernel( *m2, *m2, pt2AtomKernel, pt2BondKernel, edgeKernelParameter );
	      (*m2)->setSelfKernel(sk2);


	      nk = k / sqrt(sk1 * sk2);
	      (*gram)[i][j] = k;
	      (*gramNormal)[i][j] = nk;
	      
	      j++;
	      
	    }
	    i++;
	  }
	  
	}
	
	gramCalculated = true;

}











void MoleculeSet::kernelCompute( Molecule* m2, double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ), double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), vector<double>* resultsRaw, vector<double>* resultsNormal, int convergenceCondition, int aParameter2, bool silentMode ){
	vector<Molecule*>::iterator m1;

	double sk1 = -1000.0;		// self kernel mol1
	double sk2 = -1000.0;		// self kernel mol2
	double k = 0;	// kernel
	double nk = 0;						// normalised kernel

	int nbUndistinct = 0;  // number of indistinct molecules
	int nbUndistinctDiffAct = 0; // number of indistinct molecules with different activity status
	int nbUndistinctUnknownAct = 0; // number of indistinct molecules with unknown activity
	stringstream confoundedMolecules; // to report the confounded molecules

	int nbOrthogonal = 0;  // number of orthogonal molecules
	int nbOrthogonalDiffAct = 0; // number of orthogonal molecules with different activity status
	int nbOrthogonalUnknownAct = 0; // number of orthogonal molecules with unknown activity
	stringstream orthogonalMolecules; // to report the orthogonal molecules


	for( m1 = this->begin(); m1 != this->end(); m1++ ){

		sk1 = (*m1)->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, convergenceCondition, aParameter2 );
		sk2 = m2->getSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, convergenceCondition, aParameter2 );
		k = (*m1)->computeKernel(m2, pt2GraphKernel, pt2AtomKernel, pt2BondKernel, convergenceCondition );
		nk = k / sqrt(sk1 * sk2);





		resultsRaw->push_back( k );
		resultsNormal->push_back( nk );


		if( nk > MAXSIMILARITYWARNING ){
			nbUndistinct++;
			try{
				if( (*m1)->getActivity( true ) == m2->getActivity( true ) ){
				//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == m2->getIntDescriptor(ACTIVITY)->getValue() ){
					if( !silentMode ){
						cout << "WARNING: molecule " << (*m1)->getName() << " and " << m2->getName() << " cannot be distinguished. Normalized kernel = " << nk << endl;
					}
				}else{
					nbUndistinctDiffAct++;
					if( !silentMode ){
						cout << "WARNING: molecule " << (*m1)->getName() << " and " << m2->getName() << " cannot be distinguished and they have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
					}
					confoundedMolecules << (*m1)->getName() << " - " << m2->getName() << endl;
				}
			} catch( CError e ){
				nbUndistinctUnknownAct++;
				continue;
			}
		}

		if( nk < MINSIMILARITYWARNING ){
			nbOrthogonal++;
			try{
				if( (*m1)->getActivity( true ) == m2->getActivity( true ) ){
				//if( (*m1)->getIntDescriptor(ACTIVITY)->getValue() == m2->getIntDescriptor(ACTIVITY)->getValue() ){
					if( !silentMode ){
						cout << "WARNING: molecule " << (*m1)->getName() << " and " << m2->getName() << " are orthogonal. Normalized kernel = " << nk << " but have identical biological activity" << endl;
					}
				}else{
					nbOrthogonalDiffAct++;
					//cout << "WARNING: molecule " << (*m1)->getName() << " and " << m2->getName() << " are orthogonal but have DIFFERENT BIOLOGICAL ACTIVITY. Normalized kernel = " << nk << endl;
					orthogonalMolecules << (*m1)->getName() << " - " << m2->getName() << endl;
				}
			} catch( CError e ){
				nbOrthogonalUnknownAct++;
				continue;
			}
		}
	}

}


void MoleculeSet::writeGramMatrix( string aFileName, bool normal, bool self, bool silentMode ){

	// open output file ofile

	if(gramCalculated == true) {
		//cout << "completing matrix" << endl;
//		completeGramMatrix();
	}else{
		//cout << "calculating matrix" << endl;
		//cout << "with param " << getPq() << " " << getConvergenceCondition() << endl;
		// kashimaGramMatrix( getPq(), getConvergenceCondition() );
	}

	string baseName = aFileName;

	if(normal == false){
		aFileName = aFileName + ".gramRaw";
	}else{
		aFileName = aFileName + ".gramNormal";
	}

	fstream oFile;
	oFile.open( aFileName.c_str(), ios::out );
	if(!oFile.good()){
		CError e = CError(COULDNOTOPENFILE, "MoleculeSet::writeGramMatrix: could not write file " + aFileName);
		e.describe();
		throw(e);
	}

	oFile.precision(GRAMDECIMALPRECISION);

	int i = 0;
	int j = 0;

	vector<Molecule*>::iterator m;
	vector< vector<double> >::iterator c1;
	vector<double>::iterator c2;

	//cout << "writing " << subsetStart << " + " << subsetSize << endl;

	oFile << "name\t";
	//for( m = comparisonSet->begin() + subsetStart; m != comparisonSet->begin() + subsetStart + subsetSize - 1; m++ ){
	bool first=true;


	//cout << "  MoleculeSet::writeGramMatrix ; comparisonSet.size() = " << comparisonSet->numMolecules() << endl;
	//cout << "  MoleculeSet::writeGramMatrix ; size() = " << numMolecules() << endl;

	//cout << "  MoleculeSet::writeGramMatrix ; subsetstart = " << subsetStart << endl;

	for( m = comparisonSet->begin(); m != comparisonSet->end(); m++ ){
		//cout  << (*m)->getName() << "\t";
		if( first ){
			oFile << (*m)->getName();
			//cout << "\t - getName = " << (*m)->getName() << endl;
			first = false;
		}else{
			oFile << "\t" << (*m)->getName();
			//cout << "\t - getName = " << (*m)->getName() << endl;
		}
	}
	oFile << endl;
	//cout << endl;
//	oFile << (*m)->getName() << endl;

	if(normal == false){
		// write the raw gram matrix

		for( c1 = gram->begin(), m = begin() + subsetStart; c1 != gram->end(); c1++, m++ ){
			oFile << (*m)->getName() << "\t";
			j = 0;
			for( c2 = (*c1).begin(); c2 != (*c1).end()-1; c2++ ){
	 			oFile << *c2 << "\t";
				//oFile << (*gram)[i][j];
				j++;
			}
			oFile << *c2 << endl;
			i++;
		}
	}else{
		// write the normalised gram matrix

		for( c1 = gramNormal->begin(), m = begin() + subsetStart; c1 != gramNormal->end(); c1++, m++ ){
			oFile << (*m)->getName() << "\t";
			j = 0;
			for( c2 = (*c1).begin(); c2 != (*c1).end()-1; c2++ ){
				oFile << *c2 << "\t";
				j++;
			}
			oFile << *c2 << endl;
			i++;
		}
	}

	if( self == true ){
		writeSelfKernelList( baseName, silentMode );
	}

}

void MoleculeSet::resetGramMatrix(){
	gram->clear();
	gramCalculated = false;
}

/** hides all hydrogens in all molecules of the set
*/
void MoleculeSet::hideHydrogens(){
	cout << "HIDING HYDROGENS in the " << numMolecules() << " molecules of the molecule set" << endl;

	vector<Molecule*>::iterator m;
	int i = 1;
	for( m = begin(); m != end(); m++ ){
		(*m)->hideHydrogens();
		i++;
	}
}


/** hides all but the largest connected graph in all molecules of the set
*/
void MoleculeSet::hideSalts( string aReportFileName ){

	stringstream out;

	cout << "HIDING SALTS COUNTERIONS in the " << numMolecules() <<
	" molecules of the molecule set" << endl;

	out<<"name;graphRemoved;biggestComponent;secondBiggestComponent;difference;numHiddenAtoms;warning"
	<< endl;


	int numRemoved = 0;
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		//cout << "MoleculeSet: calling hideSalts() " << endl;
		numRemoved = (*m)->hideSalts( &out );
	}

	if( aReportFileName != "" ){
		ofstream outFile;
 		outFile.open( aReportFileName.c_str(), ios::out );
		if( !outFile.good() ){
			CError e = CError(COULDNOTOPENFILE, aReportFileName + " could not open file");
			e.describe();
			throw(e);
		}

		outFile << out.str();

		outFile.close();
	}

}

/** restores the hidden atoms for all molecules */
void MoleculeSet::restoreHiddenAtoms(){
	vector<Molecule*>::iterator m;

  for( m = begin(); m != end(); m++ ){
		(*m)->restoreHiddenAtoms();
	}
}
/** sets the value of intDescriptor aLabel of molecule aLabel to aValue
*/
void MoleculeSet::setIntDescriptor( string aLabel, string aMolecule, int aValue ){
	(*this)[aMolecule]->setIntDescriptor( aLabel, aValue, "", "", true, true );
}

/** sets the avtivity of aMolecule to aValue
*/
void MoleculeSet::setActivity( string aMolecule, float aValue ){
	(*this)[aMolecule]->setActivity( aValue );
}
/** writes a file with the biological activity of molecules in a format compatible with GIST */
void MoleculeSet::writeActivityFile(string aFileName, bool addActivityExtension, string activityDescriptor ){

	if( addActivityExtension ){
		aFileName = aFileName + ".activity";
	}

	//cout << "writing activity file " << aFileName << " with descriptor " << activityDescriptor << endl;

	fstream oFile;
 	oFile.open( aFileName.c_str(), ios::out );
	if(!oFile.good()){
		CError e = CError(COULDNOTOPENFILE, "MoleculeSet::writeActivityFile: could not write file " + aFileName);
		e.describe();
		throw(e);
	}

	//if(subsetSize < 1){
	//	subsetSize = numMolecules() - subsetStart;
	//}

	//cout << "for molecules " << subsetStart  << " to " << subsetSize << endl;

	oFile << "label" << "\t" << "class" << endl;

	vector<Molecule*>::iterator m;
	//for( m = begin() + subsetStart; m != begin() + subsetStart + subsetSize; m++ ){
	//cout << "NOW" << endl;
	//cout << toString() << endl;
	for( m = begin(); m != end(); m++ ){
		//cout << (*m)->toStringLong() << endl;
		//Descriptor<int> *D = (*m)->getIntDescriptor( activityDescriptor, false );
		//cout << "\t" << D->getValue() << endl;

		try{
			if( activityDescriptor == "activity"){
				oFile << (*m)->getName() << "\t" << (*m)->getActivity() << endl;
			}else{
				oFile << (*m)->getName() << "\t" << (*m)->getIntDescriptor( activityDescriptor, false )->getValue() << endl;
			}
		}catch( CError e ){
			continue;
		}

	}

	oFile.close();

}

/** reads descriptors for the molecules from a ; separated file
	the first line in the file indicates the Descriptor comment text
	the second line in the file indicates the Descriptor units
	the third line in the file indicates the Descriptor types
	the forth line indicates the Descriptor names

	data start at fifth line.

	the first column should contain the same molecule names as those obtained
	by getName() for the molecules in the dataset

	example:

	Full name;Log octanol partition coefficient;boiling point;special comment
	NA;NA;K;NA
	string;float;float;string
	name;logP;Bp;acomment
	ethane;10;200;junk data


	WARNING: does not check for duplicate name entries in the descriptor file
	WARNING: does not add anything to the molecules not contained in the descriptor file

*/
void MoleculeSet::readDescriptorFile( string aFileName, string separator ){
	ifstream inFile;

	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	int linesize = 512;
	int lineNumber = 0;
	char *line = new char[linesize];
	string stringLine;
	vector<string> words;

	vector<string> comments;
	vector<string> units;
	vector<string> types;
	vector<string> names;

	vector<Molecule*> modMol;
	int numD = 0;

	while( !inFile.eof() ){

		lineNumber++;
		//cout << lineNumber << endl;

		// read line
		inFile.getline(line, linesize-1,'\n');
		stringLine = line;

		//cout << stringLine << endl;

		if( lineNumber == 1 ){
			// comment line
			numD = StringUtils::Split(line, separator, comments);
			cout << "FOUND " << numD << " columns in first line" << endl;
		}else if( lineNumber == 2 ) {
			// unit line
			int a = StringUtils::Split(line, separator, units);
			if( a != numD ){
				delete[] line;
				stringstream out;
				out  << "MoleculeSet::readDescriptorFile: error in line " << lineNumber << " of file " << aFileName << " found " << a << " columns while " << numD << " are required";
				CError e(BADFILE, out.str() );
				e.describe();
				throw(e);
			}

		}else if( lineNumber == 3 ) {
			// type line
			int a = StringUtils::Split(line, separator, types);
			if( a != numD ){
				delete[] line;
				stringstream out;
				out  << "MoleculeSet::readDescriptorFile: error in line " << lineNumber << " of file " << aFileName << " found " << a << " columns while " << numD << " are required";
				CError e(BADFILE, out.str() );
				e.describe();
				throw(e);
			}

		}else if( lineNumber == 4 ) {
			// names line
			int a = StringUtils::Split(line, separator, names);
			if( a != numD ){
				delete[] line;
				stringstream out;
				out  << "MoleculeSet::readDescriptorFile: error in line " << lineNumber << " of file " << aFileName << " found " << a << " columns while " << numD << " are required";
				CError e(BADFILE, out.str() );
				e.describe();
				throw(e);
			}

		}else if( lineNumber > 4) {
			// data line
			if(stringLine!=""){

				//cout << "spliting" << endl;
				int a = StringUtils::Split(line, separator, words);
				if( a != numD ){
					delete[] line;
					stringstream out;
					out  << "MoleculeSet::readDescriptorFile: error in line " << lineNumber << " of file " << aFileName << " found " << a << " columns while " << numD << " are required";
					CError e(BADFILE, out.str() );
					e.describe();
					throw(e);
				}


				string molName = (*words.begin());

				vector<string>::iterator itw = words.begin()+1;
				vector<string>::iterator itc = comments.begin()+1;
				vector<string>::iterator itu = units.begin()+1;
				vector<string>::iterator itt = types.begin()+1;
				vector<string>::iterator itn = names.begin()+1;

				//cout << "looking for molecule " << molName << endl;

				try{
					Molecule* m = findFirstMoleculeWithName( molName );
					Descriptor<int>* di = NULL;
					Descriptor<float>* df = NULL;
					Descriptor<string>* ds = NULL;




					//cout << "found" << endl;

					modMol.push_back( m );

					for(; itw!= words.end(); itw++, itc++, itu++, itt++, itn++ ){
						//cout << "creating descriptor " << *itn << " for molecule " << molName << endl;
						// look for the molecule
						// set descriptor for the molecule
						if( *itt == "integer" ){
							if( *itw == "NA" || *itw == "" ){
								//cout << "setting empty INT" << endl;
								di = m->setIntDescriptor(*itn, 0,
								 *itu, *itw, true, true );

								di->setEmpty();
							}else{
								//cout << "setting INT" << endl;
								m->setIntDescriptor(*itn,
								 StringUtils::toInt( *itw ),
								 *itu, *itw, true, true );
							}
						}else if( *itt == "float" ){
							if( *itw == "NA" || *itw == "" ){
								//cout << "setting empty FLOAT" << endl;
								df = m->setFloatDescriptor(*itn, 0,
								 *itu, *itw, true, true );
								df->setEmpty();
							}else{
								//cout << "setting FLOAT" << endl;
								m->setFloatDescriptor(*itn,
								 StringUtils::toFloat( *itw ), *itu,
								  *itw, true, true );
							}
						}else if( *itt == "string" ){
							if( *itw == "NA" || *itw == "" ){
								//cout << "setting empty STRING" << endl;
								ds = m->setStringDescriptor(*itn, "NA",
								 *itu, *itw, true, true );
								ds->setEmpty();
							}else{
								//cout << "setting STRING" << endl;
								m->setStringDescriptor(*itn, *itw,
								 *itu, *itw, true, true );
							}
						}else{
							delete[] line;
							//unknown data type emit error
							stringstream out;
						  out
						  << "MoleculeSet::readDescriptorFile: Invalid type "
						  << *itt << " found in descriptor in line 3 of file "
						  << aFileName << " (expecting string, integer, or float) ";
							CError e(BADFILE, out.str() );
							e.describe();
							throw(e);
						}

					}
					//m->describe();
				} catch( CError ){
					// molecule not found
					cerr << "WARNING MoleculeSet::readDescriptorFile: molecule "
					<< molName << " present in " << aFileName
					<< ", but is not present in the dataset " << endl;
				}

				//cout << "clear" << endl;

				words.clear();
			}
		}
	}

	//cout << "finish" << endl;
	//cout <<  modMol.size() << ", " <<  numMolecules() << endl;

	if( modMol.size() < numMolecules() ){
		cerr << "WARNING MoleculeSet::readDescriptorFile: descriptor file " << aFileName << " only contains descriptors for " << modMol.size() << " while dataset contains " << numMolecules() << " molecules " << endl;
	}

	modMol.clear();
	comments.clear();
	units.clear();
	names.clear();
	words.clear();

	delete[] line;

	inFile.close();
}



/** reads the activity of molecules from an activity file */
void MoleculeSet::readActivityFile( string aFileName ){

	ifstream inFile2;

 	inFile2.open( aFileName.c_str(), ios::in );
	if(!inFile2.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	int linesize = 512;
	int lineNumber = 0;
	char *line = new char[linesize];
	vector<string> words;
	string stringLine;
	int i = 0;

	while( !inFile2.eof() ){
		lineNumber++;
		inFile2.getline(line, linesize-1,'\n');
		stringLine = line;

		//cout << line << endl;

		// if line is empty or starts with // or # then skip this line
 		if(stringLine.size() == 0){
			continue;
		}
		if(stringLine.substr(0,1) == "#" || stringLine.substr(0,2) == "//"){
			continue;
		}

		StringUtils::Split(line, ";", words);

		if( words[1] == "1" ){
			//cout << "active" << endl;
			try{
			  //setIntDescriptor(ACTIVITY, words[0], 1);
			  setActivity( words[0], 1 );
			}
			catch(CError e){
				//e.describe();
			}
		}else{
			//cout << "inactive" << endl;
			try{
				//setIntDescriptor(ACTIVITY, words[0], -1);
				setActivity( words[0], -1 );
			}
			catch(CError e){
				//e.describe();
			}
		}

		words.clear();

		i++;
	}
	inFile2.close();

	delete[] line;

	activitySet = true;
}

/** loads all .mol files in a directory and adds the molecules to the set */
void MoleculeSet::readMolDirectory(string dataDir, bool genericAtomType, long beginMolecule, long endMolecule){

		// add all files in the directory
	vector<string> fileVector;
	vector<string> pathVector;
	DIR* pDir = opendir ( dataDir.c_str() );
	if ( !pDir ){
		cerr << "could not read directory " << dataDir << endl;
		exit (2);
	}
	dirent* pEntry;
	long j = 0;
	while ( ( pEntry = readdir ( pDir ) ) ) {
		string toto = pEntry->d_name;
		//cout << "found " << toto;
		if( toto != "." && toto != ".." && toto.length()>2){
			if( toto.substr(toto.length()-3,3) == "mol" ){
				//cout << " added" << endl;

			    if( ( j >= beginMolecule || beginMolecule < 0 ) && ( j <= endMolecule || endMolecule < 0 )){
				string sFound = dataDir + string ( "/" ) + string ( pEntry->d_name );
				pathVector.push_back ( sFound );
				fileVector.push_back ( toto.substr(0,toto.length()-4) );
			    }
				
				j++;

			}else{
				//cout << " not added" << endl;
			}
		}else{
			//cout << " not added" << endl;
		}
	}

	cout << "adding mol files" << endl;

	vector<string>::iterator pit;
	vector<string>::iterator fit;
	Molecule* aMolecule;
	for( pit = pathVector.begin(), fit = fileVector.begin(); pit != pathVector.end(); pit++, fit++ ){
		cout << (*pit) << endl;
		aMolecule = this->addSingleMOL( (*pit), genericAtomType );
		aMolecule->setName( (*fit) );
	}
}

/** loads all .kcf files in a directory and adds the molecules to the set */
void MoleculeSet::readKcfDirectory(string dataDir, long beginMolecule, long endMolecule){

		// add all files in the directory
	vector<string> fileVector;
	vector<string> pathVector;
  DIR* pDir = opendir ( dataDir.c_str() );
  if ( !pDir ){
		cerr << "could not read directory " << dataDir << endl;
		exit (2);
	}
  dirent* pEntry;
  int j = 0;
	while ( ( pEntry = readdir ( pDir ) ) ) {
		string toto = pEntry->d_name;
		//cout << "found " << toto;
		if( toto != "." && toto != ".." && toto.length()>2){
			if( toto.substr(toto.length()-3,3) == "kcf" ){
			    if( ( j >= beginMolecule || beginMolecule < 0 ) && ( j <= endMolecule || endMolecule < 0 )){
				//cout << " added" << endl;
	 	  	        string sFound = dataDir + string ( "/" ) + string ( pEntry->d_name );
			        pathVector.push_back ( sFound );
				fileVector.push_back ( toto.substr(0,toto.length()-4) );
			    }
			    j++;
			}else{
				//cout << " not added" << endl;
			}
		}else{
			//cout << " not added" << endl;
		}
  }

	cout << "adding kcf files" << endl;

	vector<string>::iterator pit;
	vector<string>::iterator fit;
	Molecule* aMolecule;
	for( pit = pathVector.begin(), fit = fileVector.begin(); pit != pathVector.end(); pit++, fit++ ){
		cout << (*pit) << endl;
		aMolecule = this->addSingleKCF( (*pit) );
		aMolecule->setName( (*fit) );
	}
}

/** Writes all self kernel values in a file
*/
void MoleculeSet::writeSelfKernelList( string aFileName, bool silentMode ){
	aFileName = aFileName + ".self";

	//cout << silentMode << endl;

	if( !silentMode ){
		cout << "writing self kernel file " << aFileName << endl;
	}

	fstream oFile;
 	oFile.open( aFileName.c_str(), ios::out );
	if(!oFile.good()){
		CError e = CError(COULDNOTOPENFILE, "MoleculeSet::writeSelfKernelList: could not write file " + aFileName);
		e.describe();
		throw(e);
	}

	oFile << "name" << "\t" << "selfKernel" << endl;

	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		oFile << (*m)->getName() << "\t" << (*m)->getSelfKernel() << endl;
	}

	oFile.close();

}
/** reset all calculated selfkernels for the molecules in the dataset */
void MoleculeSet::resetSelfKernels(){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->resetSelfKernel();
	}

}

void MoleculeSet::deleteAll(){
	// write all molecules
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		delete (*m);
	}
	clear();
}



/** write a MDL structure data (SD) file with the whole moleculeSet */
void MoleculeSet::writeSD( string aFileName, bool selectedOnly ){

	ofstream outFile;
	outFile.open( aFileName.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFileName + " could not open file");
		e.describe();
		throw(e);
	}

	if( selectedOnly == true ){
		// write only selected molecules
		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			if( (*m)->isSelected() ){
				MoleculeUtils::writeMDLHeaderBlock( *(*m), outFile );
				MoleculeUtils::writeMDLCtabBlock( *(*m), outFile );
				MoleculeUtils::writeMDLNSDBlock( *(*m), outFile );
			}
		}
	}else{
		// write all molecules
		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			//cout << "writing header" << endl;
			MoleculeUtils::writeMDLHeaderBlock( *(*m), outFile );
			//cout << "writing ctab" << endl;
			MoleculeUtils::writeMDLCtabBlock( *(*m), outFile );
			//cout << "writing NSD" << endl;
			MoleculeUtils::writeMDLNSDBlock( *(*m), outFile );
		}
	}
		outFile.close();
}


void MoleculeSet::writeSubsetSD( string aFileName, vector<string>* anOrder ){

	ofstream outFile;
	outFile.open( aFileName.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFileName + " could not open file");
		e.describe();
		throw(e);
	}

	// write all molecules
	vector<string>::iterator m;
	Molecule* aMol = NULL;
	for( m = anOrder->begin(); m != anOrder->end(); m++ ){
		try{
			aMol = getMolByName( (*m) );

			MoleculeUtils::writeMDLHeaderBlock( *aMol, outFile );
			MoleculeUtils::writeMDLCtabBlock( *aMol, outFile );
			MoleculeUtils::writeMDLNSDBlock( *aMol, outFile );

		}catch( CError e ){
			continue;
		}

	}

	outFile.close();
}

void MoleculeSet::writeSubsetKCF( string aFileName, vector<string>* anOrder ){

	ofstream outFile;
	outFile.open( aFileName.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(COULDNOTOPENFILE, aFileName + " could not open file");
		e.describe();
		throw(e);
	}

	// write all molecules
	vector<string>::iterator m;
	Molecule* aMol = NULL;
	for( m = anOrder->begin(); m != anOrder->end(); m++ ){
		try{
			aMol = getMolByName( (*m) );

			MoleculeUtils::writeKCF( *aMol, outFile );

		}catch( CError e ){
			continue;
		}

	}

	outFile.close();
}


/** write a KCF file with the whole moleculeSet */
void MoleculeSet::writeKCF( string aFileName, bool selectedOnly ){

  ofstream outFile;
 		outFile.open( aFileName.c_str(), ios::out );
		if( !outFile.good() ){
			CError e = CError(COULDNOTOPENFILE, aFileName + " could not open file");
			e.describe();
			throw(e);
		}

	if( selectedOnly == true ){
		// write only selected molecules
		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			if( (*m)->isSelected() ){
				MoleculeUtils::writeKCF( *(*m), outFile );
			}
		}
	}else{
		// write all molecules
		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			MoleculeUtils::writeKCF( *(*m), outFile );
		}
	}
		outFile.close();
}


/** write a description of the moleculeSet to cout */
void MoleculeSet::describeShort(){
	cout << toStringShort() << endl;
	//DataContainer::describe();
}
void MoleculeSet::describe( bool selectedOnly ){
	cout << toString() << endl;
	//DataContainer::describe();
}
void MoleculeSet::describeLong(){
	cout << toStringLong() << endl;
	//DataContainer::describe();
}


/** returns a short string description of the MoleculeSet
		(number of molecules in the set and short description of each molecule)
*/
string MoleculeSet::toStringShort(){
	stringstream out;
	//cout << "MoleculeSet has " << numMolecules() << " molecules"<< endl;
	out << "MoleculeSet has " << numMolecules() << " molecules"<< endl;

	vector<Molecule*>::iterator m;
	//int i = 1;
	for( m = begin(); m != end(); m++ ){
		//cout << i << endl;
		out << (*m)->toStringShort() << endl;
		//cout << (*m)->toStringShort() << endl;
		//i++;
	}

	return( out.str() );
}

string MoleculeSet::toString( bool selectedOnly ){
	stringstream out;

	if( selectedOnly = false ){

		out << "MoleculeSet has " << numMolecules() << " molecules"<< endl;

		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			out << (*m)->toString() << endl;
		}

	}else{

		out << "MoleculeSet has " << numMolecules() << " molecules out of which: "<< endl;

		vector<Molecule*>::iterator m;
		int c = 0;
		for( m = begin(); m != end(); m++ ){
			if( (*m)->isSelected() ){
				out << (*m)->toString() << endl;
				c++;
			}
		}
		out << c << " are selected" << endl;
	}


	return( out.str() );
}

string MoleculeSet::toStringLong(){
	stringstream out;
	out << "MoleculeSet has " << numMolecules() << " molecules"<< endl;

	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		out << (*m)->toStringLong() << endl;
	}

	return( out.str() );
}


/** unselect all molecules in the dataset */
void MoleculeSet::unSelectAll(){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->unSelect();
	}
}
/** select all molecules in the dataset */
void MoleculeSet::selectAll(){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->select();
	}
}


void MoleculeSet::sortByMW(){
	sort( begin(), end(), AscendingMW() );
}

void MoleculeSet::sortByNumAtoms(){
	sort( begin(), end(), AscendingNumAtoms() );
}



/** sorts the molecule collection according to the molecule Descriptor descriptorName of type descriptorType */
void MoleculeSet::sortByDescriptor( string aDescriptorName, int aDescriptorType, bool reverse ){
	setSortDescriptor( aDescriptorName, aDescriptorType );
	//cout << "MoleculeSet: sorting by " << aDescriptorName << endl;

	if( reverse ){
	cout << "MoleculeSet: sorting by " << getSortDescriptorName() << " descending " << endl;
		sort( begin(), end(), DescendingOrder() );
	}else{
		cout << "MoleculeSet: sorting by " << getSortDescriptorName() << " ascending" << endl;
		sort( begin(), end(), AscendingOrder() );
	}
}

/** sorts the molecule collection according to the molecule Descriptor descriptorName.
			descriptor type will be read from descriptor name if name is of type
			******.integer, or *******.float, set to string otherwise*/
void MoleculeSet::sortByDescriptor( string aDescriptorName, bool reverse ){

	string ext = StringUtils::getExtension( aDescriptorName );

	//cout << "TYPE: " << ext << endl;

	if( ext == "integer" ){
		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 8 );
		//cout << "sorting by " << aDescriptorName << endl;
		sortByDescriptor( aDescriptorName, INTEGER, reverse );
	}else if( ext == "float" ){
		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 6 );
		//cout << "sorting by " << aDescriptorName << endl;
		sortByDescriptor( aDescriptorName, FLOAT, reverse );
	}else if( ext == "string" ){
		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 7 );
		//cout << "sorting by " << aDescriptorName << endl;
		sortByDescriptor( aDescriptorName, STRING, reverse );
	}else{
		//cout << "sorting by " << aDescriptorName << endl;
		sortByDescriptor( aDescriptorName, STRING, reverse );
	}
}

/** sets the type and name of descriptor to be used when sorting molecules */
void MoleculeSet::setSortDescriptor( string aName, int aType, bool reverse ) throw( CError ){
	// if aType == default, then molecules will be sorted by id whatever aName is
	if( aType == INTEGER || aType == FLOAT || aType == STRING || aType == DEFAULT ){
		vector<Molecule*>::iterator m;
		for( m = begin(); m != end(); m++ ){
			(*m)->setSortDescriptor( aName, aType );
		}
	}else{
		stringstream out;
		out << "MoleculeSet::setSortDescriptor: aType is none of integer, float, string or default";
		CError e( UNKNOWNDATATYPE, out.str() );
		e.describe();
		throw(e);
	}
}

string MoleculeSet::getSortDescriptorName(){
	vector<Molecule*>::iterator m = begin();

	return( (*m)->getSortDescriptorName() );
}

/** returns a pointer to the molecule with name aName in the MoleculeSet */
Molecule* MoleculeSet::findFirstMoleculeWithName( string aName ) throw( CError ){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		if( (*m)->getName() == aName ){
			return( (*m) );
		}
	}
	CError e( MOLECULENOTFOUND, "No molecule with name " + aName + " in dataset " + toStringShort() );
	throw( e );
}

/** writes a ; separated file containing all descriptors for all molecules */
void MoleculeSet::writeDescriptors( string aFileName, bool selectedOnly ) throw( CError ) {
	
	vector<Molecule*>::iterator m;
	map< string, Descriptor< string >* >::iterator ds;
	map< string, Descriptor< int >* >::iterator di;
	map< string, Descriptor< float >* >::iterator df;

	// get list of all existing descriptors with their type
	map< string, int > vs; // string descriptors
	map< string, int >::iterator v;
	bool found = false;

	for( m = begin(); m != end(); m++ ){
		for( ds = (*m)->beginStringDescriptor(); ds != (*m)->endStringDescriptor(); ds++ ){
			found = false;
			for( v = vs.begin(); v != vs.begin(); v++ ){
				if( (*v).first == (*ds).first ){
					found = true;
					break;
				}
			}
			if( found == false ){
				vs[ (*ds).first ] = STRING;
			}
		}

		for( di = (*m)->beginIntDescriptor(); di != (*m)->endIntDescriptor(); di++ ){
			found = false;
			for( v = vs.begin(); v != vs.begin(); v++ ){
				if( (*v).first == (*di).first ){
					found = true;
					break;
				}
			}
			if( found == false ){
				vs[ (*di).first ] = INTEGER;
			}
		}

		for( df = (*m)->beginFloatDescriptor(); df != (*m)->endFloatDescriptor(); df++ ){
			found = false;
			for( v = vs.begin(); v != vs.begin(); v++ ){
				if( (*v).first == (*df).first ){
					found = true;
					break;
				}
			}
			if( found == false ){
				vs[ (*df).first ] = FLOAT;
			}
		}
	}

	std::ostream* out;
	if( aFileName == "cout" ){
    out = &std::cout;
  }else{
	  out = new std::ofstream( aFileName.c_str() );

		if( !out->good() ){
			delete out;
			CError e = CError(COULDNOTOPENFILE, aFileName + " could not open file");
			e.describe();
			throw(e);
		}
  }

  //cout << "found " << vs.size() << " descriptors" << endl;

	// write list of descriptors:
	*out << "name;";
	for( v = vs.begin(); v != vs.end(); v++ ){
		//cout << "@" << endl;
		*out << (*v).first << "." << (*v).second << ";";
  }
	*out << endl;
	

	// write descriptors
	for( m = begin(); m != end(); m++ ){
		if( ( selectedOnly && (*m)->isSelected() ) || !selectedOnly ){
		  *out << (*m)->getName() << ";";

			for( v = vs.begin(); v != vs.end(); v++ ){
				if( (*v).second == INTEGER ){
					try{
					  if( (*m)->getIntDescriptor( (*v).first )->isEmpty() ){
						  *out << "NA;";
						}else{
							*out << (*m)->getIntDescriptor( (*v).first )->getValue() << ";";
						}
					}catch( CError e ){
						if( e.getType() == MISSINGDESCRIPTOR ){
							*out << "NA;";
						}else{
							throw( e );
						}
					}
				}else if( (*v).second == FLOAT ){
					try{
						if( (*m)->getFloatDescriptor( (*v).first )->isEmpty() ){
							*out << "NA;";
						}else{
		  	   		*out << (*m)->getFloatDescriptor( (*v).first )->getValue() << ";";
						}
					}catch( CError e ){
						if( e.getType() == MISSINGDESCRIPTOR ){
							*out << "NA;";
						}else{
							throw( e );
						}
					}
				}else{
					try{
						if( (*m)->getStringDescriptor( (*v).first )->isEmpty() ){
							*out << "NA;";
						}else{
							*out << (*m)->getStringDescriptor( (*v).first )->getValue() << ";";
						}
					}catch( CError e ){
						if( e.getType() == MISSINGDESCRIPTOR ){
							*out << "NA;";
						}else{
							throw( e );
						}
					}
				}
			}

			*out << endl;
		}
	}

  if (out != &std::cout){
		//(*out).close();
	  delete out;
	}

}
/** selects the molecules with the names provided as arguments in a vector< string >. Unselect all others. Returns the number of selected molecules. */
int MoleculeSet::select( vector< string >* aSubset ){
	int nbselected = 0;

	unSelectAll();

	vector< string >::iterator it;
	vector< Molecule* >::iterator m;

	for( m = begin(); m != end(); m++ ){
		for( it = aSubset->begin(); it != aSubset->end(); it++ ){
			//cout << (*m)->getName() << "?=" << (*it) << endl;
			if( (*m)->getName() == (*it) ){
				//cout << "  MATCH" << endl;
				(*m)->select();
				nbselected++;
			}
		}
	}

	return( nbselected );
}

/** unselects the molecules with the names provided as arguments in a vector< string >. Returns the number of selected molecules. */
int MoleculeSet::unSelect( vector< string >* aSubset ){
	int nbselected = 0;

	vector< string >::iterator it;
	vector< Molecule* >::iterator m;

	for( m = begin(); m != end(); m++ ){
		for( it = aSubset->begin(); it != aSubset->end(); it++ ){
			//cout << (*m)->getName() << "?=" << (*it) << endl;
			if( (*m)->getName() == (*it) ){
				//cout << "  MATCH" << endl;
				(*m)->unSelect();
				nbselected++;
			}
		}
	}

	return( nbselected );
}



/** write a mol file for each molecule in the set to aDirName. if selectedOnly == true then only selected molecules are written */
long MoleculeSet::writeMolToDir( string aDirName, bool selectedOnly ){
	vector<Molecule*>::iterator m;

	int i = 0;
	string fileName = "";
	if( selectedOnly == false ){
		for( m = begin(); m != end(); m++ ){

			//cout << "next" << endl;
			//cout << (*m)->getName() << endl;

			if(  StringUtils::right( (*m)->getName(), 4 ) == ".mol" ){
				fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() );
			}else{
				fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() ) + ".mol";
			}

			//cout << "writing " << fileName << endl;
			(*m)->writeMOL( fileName );
			//cout << "done" << endl;
			i++;
		}
	}else{
		for( m = begin(); m != end(); m++ ){
			if( (*m)->isSelected() ){
				if(  StringUtils::right( (*m)->getName(), 4 ) == ".mol" ){
					fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() );
				}else{
					fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() ) + ".mol";
				}
				//cout << "writing " << fileName << endl;
				(*m)->writeMOL( fileName );
				//cout << "done" << endl;
				i++;
			}
		}
	}
	return i;
}

/** reads an sd file and returns the number of created molecules
	beginMolecule and endMolecule specify the index of the first and last
	molecule to include in the training set (starting count from 0). values of -1 (default)
	means no limit.
	*/
int MoleculeSet::addSD( string aFileName, bool genericAtomTypeFlag, long beginMolecule, long endMolecule ){

	ifstream inFile;
	inFile.open( aFileName.c_str(), ios::in );
	if( !inFile.good() ){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	Molecule* m;
	int j=0;
	while ( !inFile.eof() ){

		if( j >= beginMolecule && ( j <= endMolecule || endMolecule == -1 ) ){
		    //#ifdef DEBUG
			  cout << "reading molecule " << j+1 << endl;
			  //#endif

			m = new Molecule;
			#ifdef DEBUG
		  	  cout << "reading Header" << endl;
			#endif
			try{
				MoleculeUtils::readMDLHeaderBlock( *m, inFile );
				//cout << m->getStringDescriptor("comment")->isEmpty() << endl;
			}catch( CError e ) {
				delete m;
				if( e.getType() == EOFERROR ){
					//e.describe();
					inFile.close();
					return(j);
				}else{
					cout << "MoleculeSet::addSD: CAUGHT UNKNOWN ERROR:" << endl;
					e.describe();
					throw(e);
				}
			}
			


			#ifdef DEBUG
			  cout << "reading CTAB" << endl;
			#endif
			bool hasStructure = true;
			try{
				MoleculeUtils::readMDLCtabBlock( *m, inFile, genericAtomTypeFlag );
  			}catch( CError e ) {
				if( e.getType() == EOFERROR ){
					//e.describe();
					delete m;
					inFile.close();
					return(j);
				}else if( e.getType() == NOSTRUCTURE ){
					hasStructure = false;
				}else{
					delete m;
					cout << "MoleculeSet::addSD: CAUGHT UNKNOWN ERROR:" << endl;
					e.describe();
					throw(e);
				}
			}

			#ifdef DEBUG
				cout << "reading SD block" << endl;
			#endif
			try{
				MoleculeUtils::readMDLNSDBlock( *m, inFile );
			}catch( CError e ) {
				if( e.getType() == EOFERROR ){
					delete m;
					//e.describe();
					inFile.close();
					return(j);
				}else{
					delete m;
					cout << "MoleculeSet::addSD: CAUGHT UNKNOWN ERROR:" << endl;
					e.describe();
					throw(e);
				}
			}

			//m->describe();
			#ifdef DEBUG
			  cout << m->getName() << endl;
			#endif

			if( hasStructure == true ){
				addMolecule( m );
				//cout << " MW: " << m->getMW() << endl;
			}else{
				cerr << "WARNING: molecule " << m->getName() << " was not loaded from " << aFileName << " because it has no structure information " << endl;
				m->describe();
				delete m;
			}
		}else{
			// skip entry
			try{
			  #ifdef DEBUG
			  cout << "skipping molecule " << j+1 << endl;
			  #endif
				MoleculeUtils::skipMDLEntry( *m, inFile );
  			}catch( CError e ) {
				if( e.getType() == EOFERROR ){
					inFile.close();
					return(j);
				}else{
					cout << "MoleculeSet::addSD: CAUGHT UNKNOWN ERROR:" << endl;
					e.describe();
					throw(e);
				}
			}

		}

		if( j > endMolecule && endMolecule > 0 ){
			inFile.close();
			return(j);
		}

		j++;
	}

	inFile.close();

	return( j );
}


/** reads an KCF file and returns the number of created molecules */
int MoleculeSet::addKCF( string aFileName, long beginMolecule, long endMolecule ){

	ifstream inFile;
 	inFile.open( aFileName.c_str(), ios::in );
	if( !inFile.good() ){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}


	KCFMolecule* m;
	bool valid = false;
	int j=0;
	while ( !inFile.eof() ){

	    if( j >= beginMolecule && ( j <= endMolecule || endMolecule == -1 ) ){

		m = new KCFMolecule;

		//cout << " READING MOLECULE " << j << endl;

		try{
			valid = MoleculeUtils::readKCFMolecule( *m, inFile );
		}catch( CError e ) {
			if( e.getType() == EOFERROR ){
				delete m;
				//e.describe();
				return(j);
			}else{
				delete m;
				cout << "MoleculeSet::addKCF: CAUGHT UNKNOWN ERROR:" << endl;
				e.describe();
				throw(e);
			}
		}

		if( valid) {
			//m->describe();
			addMolecule( m );
		}
	    }

	    if( j > endMolecule && endMolecule > 0 ){
			inFile.close();
			return(j);
		}

		j++;
	}

	inFile.close();

	return( j );
}

/** add all molecules in aSet to the current set */
int MoleculeSet::add( MoleculeSet* aSet ){
	vector<Molecule*>::iterator m;
	int i = 0;
	for( m = aSet->begin(); m != aSet->end(); m++ ){
		addMoleculeCopy( (*m) );
		i++;
	}
	return( i );
}
/** add a copy of the molecule in argument. Used by add to merge to datasets */
Molecule* MoleculeSet::addMoleculeCopy( Molecule* aMolecule ){
	Molecule* aNewMolecule = new Molecule( *aMolecule );
	//aNewMolecule = aMolecule;

	push_back( aNewMolecule );
	aNewMolecule->setKashimaKernelProb( pq );
	return( aNewMolecule );
}
/** set integer aName to aValue for all compounds in the dataset */
void MoleculeSet::setIntDescriptor( string aName, int aValue ){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->setIntDescriptor( aName, aValue, "", "", true, true );
	}
}


 
/** returns the mean distance of all molecules to the moleculeSet barycenter */
double MoleculeSet::diversityBaryMean(){

	//cout << "MoleculeSet::diversityBaryMean: " << endl;

	if(!gramCalculated){
		// gram matrix was not calculated, emit error
		CError e = CError(MISSINGDATA, "MoleculeSet::diversityBaryMean: no gram matrix");
		e.describe();
		throw(e);
	}

	double sumSelf = 0; // we use normal kernel matrix so average of self kernels is 1
	double sumNonSelf = 0;


	vector< vector<double> >::iterator l;
	vector<double>::iterator c;

	long lc = 0;
	long cc = 0;
	long numMol = gramNormal->size();
	//numMol = ((numMol * numMol) - numMol)/2;

	double kernelNormalValue = 0;

	for( l = gramNormal->begin(); l != gramNormal->end(); l++ ){
		//cout << lc << ": ";
		cc = lc;
		for( c = (*l).begin()+lc; c != (*l).end() ; c++ ){
			kernelNormalValue = (*c);
			//cout << lc << ";" << cc << ";" << kernelNormalValue << endl;
			if( cc == lc ){
				sumSelf += kernelNormalValue;
			}
			else{
				sumNonSelf += kernelNormalValue;
			}
			cc++;
		}
		lc++;
	}

	sumNonSelf = ( sumNonSelf * 2 ) + sumSelf;

	sumSelf = sumSelf / numMol;
	sumNonSelf = sumNonSelf / (numMol*numMol);

	//sumNonSelf = sumNonSelf/ ((( numMol*numMol) - numMol) / 2);
	//sumNonSelf = sumNonSelf / ( ( ( ( numMol * numMol) - numMol) / 2 ) + numMol );
	//cout << numMol << " " << sumNonSelf << " " << sumSelf << " " << endl;

	return( sumSelf - sumNonSelf );

}



/** this function reads a gram matric matching the dataset. Emits an error if the dimension of the read gram matrix read does not match the number of compounds in the dataset. */
void MoleculeSet::readGramNormal( string aFileName ){
	readGram( aFileName, gramNormal );
}
void MoleculeSet::readGramRaw( string aFileName ){
	readGram( aFileName, gram );
}


void MoleculeSet::readGram( string aFileName, vector< vector<double> >* aGram ){

	// clear potentially existing matrix
	bool readingGramNormal = true;
	if( aGram == gram){
		readingGramNormal = false;
	}

	delete gram;
	gram = new vector< vector<double> >;

	delete gramNormal;
	gramNormal = new vector< vector<double> >;

	if( readingGramNormal == true ){
		aGram = gramNormal;
	}else{
		aGram = gram;
	}

	// clear MoleculeSet
	clear();

	ifstream inFile;
 	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}
	int linesize = 30000;
	char *line = new char[linesize];
	vector<string> words;
	string stringLine;
	int lineNumber = 1;

	inFile.getline(line, linesize-1,'\n');
	lineNumber++;

	//cout << "MoleculeSet::readGram 2 " << endl;


	//cout << "MoleculeSet::readGram  reading names" << endl;
	// add molecules (no structure) with their names
	StringUtils::Split(line, "	", words);
	vector< string >::iterator wi2;

	for( wi2 = words.begin(); wi2 != words.end(); wi2++ ){
		Molecule* aMol = new Molecule();
		aMol->setName( (*wi2) );
		push_back( aMol );
	}

	//cout << "MoleculeSet::readGram 3 " << endl;


	//cout << "MoleculeSet::readGram  reading values" << endl;
	// read matrix values
	while( !inFile.eof() ){

		//cout << lineNumber << endl;

		inFile.getline(line, linesize-1,'\n');
		stringLine = line;

		if(stringLine.size() == 0){
			continue;
		}
		if(stringLine.substr(0,1) == "#" || stringLine.substr(0,2) == "//"){
			continue;
		}

		words.clear();
		//cout << "splitting line " << lineNumber << ": " << line << " into words" << endl;
		StringUtils::Split(line, "	", words);
		//cout << "splitted in " << words.size() << endl;
		if( words.size() < 2 ){
			// bad line emit error
			stringstream out;
			out << aFileName << " error at line " << lineNumber;
			CError e = CError( BADFILE, out.str() );
			e.describe();
			throw(e);

		}else{

			vector< string >::iterator wi;
			aGram->push_back( vector<double>() );
			vector< vector<double> >::iterator li = gramNormal->begin();
			vector<double>::iterator ci = (*li).begin();

			for( wi = words.begin()+1, ci = (*li).begin() ; wi != words.end(); wi++, ci++){
				//cout << atof( (*wi).c_str() ) << " " << lineNumber-2 << endl;
				(*aGram)[lineNumber-2].push_back( atof( (*wi).c_str() ) );
			}
			li++;
		}
		words.clear();
		lineNumber++;
	}

	delete[] line;

	inFile.close();
	gramCalculated = true;

	//cout << "gramNormal has now size " << gramNormal->size() << ";" << aGram->size() << endl;
}




/** sets the morganLabels of each molecule to the anOrder iteration of
		the Morgan index calculation process
*/
void MoleculeSet::setMorganLabels( int anOrder ){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->setMorganLabels( anOrder );
	}
}

/** sets the uniqueMorganIndex of each atom to the Morgan index having
		the maximum of different connectivity values for the molecule
*/
void MoleculeSet::setUniqueMorganIndices(){
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->setUniqueMorganIndices();
	}
}


/*!
    \fn MoleculeSet::selectByMW( float minmw, float maxmw, bool addDescriptor )
    selects all molecules in the dataset with mw >= minmw and < maxmw
    if maxmw = -1 then there is no maximum limit
 */
int MoleculeSet::selectByMW( float minmw, float maxmw, bool addD){
	int i = 0;
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		float mw = (*m)->getMW();
		if( mw >= minmw ){
			if( mw <= maxmw || maxmw == -1 ){
				(*m)->select();
				i++;
				if( addD == true ){
					(*m)->setFloatDescriptor(
						"mw", mw, "", "Molecular weight", true, true
					);
				}
			}
		}
	}
	return( i );
}

int MoleculeSet::selectByNumAtoms( float minAtoms, float maxAtoms, bool addD){
	int i = 0;
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		long numAtoms = (*m)->getNumAtoms();
		if( numAtoms >= minAtoms ){
			if( numAtoms <= maxAtoms || maxAtoms == -1 ){
				(*m)->select();
				i++;
				if( addD == true ){
					(*m)->setIntDescriptor(
						"numAtoms", numAtoms, "", "Number of atoms", true, true
					);
				}
			}
		}
	}
	return( i );
}


void MoleculeSet::addFragmentsToSet( Molecule* aMol, int minAtoms ){

	string molName = aMol->getName();
	string fragName = "";

	aMol->unsetBondFlags();
	vector<Atom*>::iterator ai;

	vector<Atom*>* atomList = new vector<Atom*>; 	// a vector containing all atoms
							// to be successively examined
	for( ai = aMol->beginAtom(); ai != aMol->endAtom(); ai++ ){
		//cout << " adding atom " << (*ai)->toStringShort() << endl;
		atomList->push_back( (*ai) );
	}

	for( ai = atomList->begin(); ai!= atomList->end(); ai++){
		//cout << "####################################" << endl;
		//(*ai)->describeShort();

		map<Atom*, Bond*>::iterator bi;
		vector<Bond*> bondList;  // a vector containing all bonds to be successively removed
		for( bi = (*ai)->beginBond(); bi != (*ai)->endBond(); bi++ ){
			//cout << "  # examining bond " << (*bi).second->toStringShort();

			if( (*bi).second->getLabel() != AROMATICBOND && !(*bi).second->hasFlag() ){
			//if( (*bi).second->getLabel() != AROMATICBOND ){
				//cout << " stored " << endl;
				bondList.push_back( (*bi).second );
				try{
					Bond* oppositeBond = (*bi).first->getBondWithTarget( (*ai) );
					oppositeBond->setFlag();
				}catch( CError ){
				}

			}//else{
				//cout << " ignored " << endl;
			//}
		}

		vector<Bond*>::iterator bli;
		for( bli = bondList.begin(); bli != bondList.end(); bli++ ){
			//cout << "##" << " treating bond " << (*bli)->toStringShort() << endl;
			aMol->setName( molName + "." + (*bli)->toStringShort() );

			//describe();
			//(*bli)->describeShort();
			(*bli)->getSource()->hideBond( (*bli) );
			(*bli)->getTarget()->hideBond( (*bli)->getSource() );

			aMol->markFragments();
			pushFragments( aMol, minAtoms );
			aMol->unmarkFragments();

			(*bli)->getSource()->restoreHiddenBonds();
			(*bli)->getTarget()->restoreHiddenBonds();
			//describe();
		}
	}

	aMol->setName( molName );
	aMol->unsetBondFlags();

	atomList->clear();
	delete atomList;

	//removeDuplicates();
}

void MoleculeSet::pushFragments( Molecule* aMol, int minAtoms ){

	map< int, int >::iterator j;
	int i = 1;
	for( j = aMol->beginComponentSizes(); j != aMol->endComponentSizes(); j++ ){
		//describe();
		aMol->hideAllFragmentsBut( (*j).first );

		if( aMol->getNumAtoms() >= minAtoms ){

			string molName = aMol->getName();
			stringstream newName;
			newName << molName << "." << i;
			aMol->setName( newName.str() );

			addMoleculeCopy( aMol );
			aMol->setName( molName );

			i++;
		}

		aMol->restoreHiddenAtoms( false );

	}
}

void MoleculeSet::removeDuplicates()
{
	// for each molecule filter based on
	// 1) number of atoms
	// 2) molecular weight
	// 3) graph kernel

	vector<Molecule*>::iterator m;
	vector<Molecule*>::iterator n;

	int na1 = 0;
	int na2 = 0;
	float mw1 = 0;
	float mw2 = 0;

	for( m = begin(); m != end(); m++ ){
		na1 = (*m)->getNumAtoms();
		for( n = begin(); n != end(); n++ ){
			// do not compare a molecule with itself
			if( (*m) != (*n) ){
				// same number of atoms?
				na2 =(*n)->getNumAtoms();

				//cout << "Num atoms: " << na1-na2 << endl;

				if( na1 == na2 ){
					// was n not already marked as a duplicate?
					if(
					!(*n)->hasIntDescriptor( "markForDelete" ) &&
					!(*n)->hasIntDescriptor( "markForKeep" )
					){


						// if it has not the same sum of bond orders

						long bo1 = (*m)->bondSum();
						long bo2 = (*n)->bondSum();

						//cout << "Bond Orders: " << bo1 << " " << bo2 << endl;

						if( bo1 == bo2 ){

							//cout << "MW" << endl;

							mw1 = (*m)->getMW();
							mw2 = (*n)->getMW();

							//cout << "MW: " << mw1-mw2 << endl;

							// same molecular weight?
							if( fabs( mw1 - mw2 ) < DUPLICATEMWIDENTITY ){
								// same graphKernel?

								//cout << "calculating kernel bewteen " << (*n)->getName() << " and " << (*m)->getName() << endl;

								(*m)->setKashimaKernelProb( 0.01 );
								(*n)->setKashimaKernelProb( 0.01 );

						 		double sk1 = (*m)->getSelfKernel(
								MoleculeUtils::moleculeKernel,
								MoleculeUtils::atomKernelSymbol,
								MoleculeUtils::bondKernelType, 1000
								);

								double sk2 = (*n)->getSelfKernel(
								MoleculeUtils::moleculeKernel,
								MoleculeUtils::atomKernelSymbol,
								MoleculeUtils::bondKernelType, 1000
								);

								double k = (*m)->computeKernel(
								(*n),
								MoleculeUtils::moleculeKernel,
								MoleculeUtils::atomKernelSymbol,
								MoleculeUtils::bondKernelType, 1000
								);

								double nk = k / sqrt(sk1 * sk2);

								//cout << "kernel: " << nk << endl;

								if( nk >= 1){
									//(*m) and (*n) are duplicates, mark (*n) for delete
									(*n)->setIntDescriptor(
									"markForDelete", 1, "", "", true, true
									);

									(*m)->setIntDescriptor(
									"markForKeep", 1, "", "", true, true
									);
									//cout << "DUPLICATES" << endl;
								}
							}
						}
					}
				}
			}
		}
	}

	// delete marked compounds
	for( m = begin(); m != end(); m++ ){
		if( (*m)->hasIntDescriptor( "markForDelete" ) ){
			delete (*m);
			erase( m );
			m--;
		}else{
			(*m)->deleteDescriptor( "markForKeep" );
		}
	}
}

void MoleculeSet::deleteHiddenAtoms()
{
	vector<Molecule*>::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->deleteHiddenAtoms();
	}

}


void MoleculeSet::readGistClassifyFile( string aFileName )
{
	ifstream inFile;

	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	int linesize = 512;
	int lineNumber = 0;
	char *line = new char[linesize];
	string stringLine;
	vector<string> words;

	vector<Molecule*> modMol;

	bool readingComments = true;

	while( !inFile.eof() ){

		lineNumber++;
		// read line
		inFile.getline(line, linesize-1,'\n');
		stringLine = line;

		// data line
		if( stringLine != "" && stringLine.substr( 0, 1 ) != "#" ){

			if( readingComments == true ){
				if(  stringLine.substr(0, 4) == "name" ){
					readingComments = false;
				}
			}else{

				StringUtils::Split(line, "\t", words);
				string molName = (*words.begin());

				//cout << "looking for molecule " << molName << endl;

				try{
					Molecule* m = findFirstMoleculeWithName( molName );
					modMol.push_back( m );

					m->setStringDescriptor( "gistname", words[0],
					 "", "name in gist classify file", true, true );

					m->setIntDescriptor( "gistclass",
					 StringUtils::toInt( words[1] ),"", "activity class",
					 true, true );

					m->setFloatDescriptor( "gistdiscr",
					 StringUtils::toFloat( words[2] ),"", "activity class",
					 true, true );

					//m->describe();

				} catch( CError e ){

					// molecule not found
					cerr << "WARNING MoleculeSet::readDescriptorFile: molecule "
					<< molName << " present in " << aFileName
					<< ", but is not present in the dataset " << endl;

				}

				//cout << "clear" << endl;
				words.clear();
			}
		}
	}


	//cout << "finish" << endl;
	//cout <<  modMol.size() << ", " <<  numMolecules() << endl;

	if( modMol.size() < numMolecules() ){
		cerr << "WARNING MoleculeSet::readDescriptorFile: descriptor file " << aFileName << " only contains descriptors for " << modMol.size() << " while dataset contains " << numMolecules() << " molecules " << endl;
	}

	delete[] line;

	modMol.clear();
	words.clear();

	inFile.close();
}

void MoleculeSet::readGistActivityFile( string aFileName, string aDescriptor ){
	ifstream inFile;

	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	int linesize = 512;
	int lineNumber = 0;
	char *line = new char[linesize];
	string stringLine;
	vector<string> words;

	vector<Molecule*> modMol;

	bool readingFirstLine = true;

	while( !inFile.eof() ){

		lineNumber++;
		// read line
		inFile.getline(line, linesize-1,'\n');
		stringLine = line;

		// data line
		if( stringLine != "" && stringLine.substr( 0, 1 ) != "#" ){

			if( readingFirstLine == true ){
				//if(  stringLine.substr(0, 4) == "name" ){
					readingFirstLine = false;
				//}
			}else{

				StringUtils::Split(line, "\t", words);

				if( words.size() < 2 ){
					// bad file emit error
					delete[] line;
					stringstream s;
					s << aFileName << " error reading line " << lineNumber << endl;
					CError e = CError(BADFILE, s.str() );
					e.describe();
					throw(e);
				}

				string molName = ( *words.begin() );
				string activity = ( words[1] );

				//cout << "looking for molecule " << molName << endl;

				try{
					Molecule* m = findFirstMoleculeWithName( molName );

					modMol.push_back( m );

					m->setIntDescriptor( aDescriptor,
					 StringUtils::toInt( activity ),"", "activity class",
					 true, true );

					//m->describe();

				}catch( CError e ){
					// molecule not found
					cerr << "WARNING MoleculeSet::readDescriptorFile: molecule "
					<< molName << " present in " << aFileName
					<< ", but is not present in the dataset " << endl;
				}

				//cout << "clear" << endl;
				words.clear();
			}
		}
	}


	//cout << "finish" << endl;
	//cout <<  modMol.size() << ", " <<  numMolecules() << endl;

	if( modMol.size() < numMolecules() ){
		cerr << "WARNING MoleculeSet::readDescriptorFile: descriptor file " << aFileName << " only contains descriptors for " << modMol.size() << " while dataset contains " << numMolecules() << " molecules " << endl;
	}

	delete[] line;

	modMol.clear();
	words.clear();

	inFile.close();
}


long MoleculeSet::selectByFloatDescriptor( string aName, float aValue )
{
	vector< Molecule* >::iterator m;
	long c = 0;
	for( m = begin(); m != end(); m++ ){
		Descriptor< float >* d = NULL;

		try{
			d = (*m)->getFloatDescriptor( aName );
		}catch( CError e ){
			d = NULL;
		}

		if( d != NULL ){
			if( d->getValue( false ) == aValue ){
				(*m)->select();
				c++;
			}else{
				(*m)->unSelect();
			}
		}else{
			(*m)->unSelect();
		}
	}
	return( c );
}

long MoleculeSet::selectByIntDescriptor( string aName, int aValue )
{

	if( aName.substr( aName.length() - min( (int) aName.length(),8), 8 ) == ".integer"){
		aName = aName.substr( 0, aName.size() - 8 );
	}

	vector< Molecule* >::iterator m;
	long c = 0;
	for( m = begin(); m != end(); m++ ){
		Descriptor< int >* d = NULL;

		try{
			d = (*m)->getIntDescriptor( aName );
		}catch( CError e ){
			d = NULL;
		}

		if( d != NULL ){
			if( d->getValue( false ) == aValue ){
				(*m)->select();
				c++;
			}else{
				(*m)->unSelect();
			}
		}else{
			(*m)->unSelect();
		}
	}
	return( c );
}

long MoleculeSet::selectByActivity( float aValue )
{
	vector< Molecule* >::iterator m;
	long c = 0;
	for( m = begin(); m != end(); m++ ){
		float activityValue = -1;

		try{
			if( (*m)->getActivity() == aValue ){
				(*m)->select();
				c++;
			}else{
				(*m)->unSelect();
			}

		}catch( CError e ){
			(*m)->unSelect();
		}
	}
	return( c );
}

long MoleculeSet::selectHasActivity()
{
	vector< Molecule* >::iterator m;
	long c = 0;
	for( m = begin(); m != end(); m++ ){
		if( (*m)->hasActivity() ){
			(*m)->select();
		}//else{
			//(*m)->unSelect();
		//}
	}
	return( c );
}


long MoleculeSet::getPossibleValuesInIntDescriptor( string aDescriptorName, vector< int >* p ){


	long res = 0;

	if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),8), 8 ) == ".integer"){
		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 8 );
	}



	vector< Molecule* >::iterator m;
	for( m = begin(); m != end(); m++ ){
		int newValue = -1;
		try{
			if( aDescriptorName == "activity" ){
				newValue = (int) (*m)->getActivity();
			}else{
				newValue = (*m)->getIntDescriptor( aDescriptorName )->getValue();
			}

			vector<int>::iterator vi;
			long nfound = 0;
			for( vi = p->begin(); vi != p->end(); vi++){
				if( (*vi) == newValue ){
					nfound += 1;
					break;
				}
			}
			if( nfound == 0 ){
				p->push_back( newValue );
				res += 1;
			}
		}catch( CError e ){
			continue;
		}


	}
	return( res );


}


void MoleculeSet::binClassifyFromDescriptor( string descriptorName, float value, bool smallerOrEqual ){
	vector< Molecule* >::iterator m;
	for( m = begin(); m != end(); m++ ){
		(*m)->binClassifyFromDescriptor( descriptorName, value, smallerOrEqual );
	}
}


// initialize every element gram matrix to the given value
void MoleculeSet::initializeGram( double value ){
	gram->clear();
	gramNormal->clear();
	vector<Molecule*>::iterator m1;
	vector<Molecule*>::iterator m2;
	int i = 0;

	for( m1 = begin(); m1 != end(); m1++ ){
		gram->push_back( vector<double>() );
		gramNormal->push_back( vector<double>() );
		for( m2 = comparisonSet->begin(); m2 != comparisonSet->end(); m2++ ){
			(*gram)[i].push_back( value );
			(*gramNormal)[i].push_back( value );
		}
		i++;
	}
}


// initialize every self-kernelto the given value
void MoleculeSet::initializeSelfKernel( double value ){
  vector<Molecule*>::iterator m;
  for( m = begin(); m != end(); m++ ){
    (*m)->resetSelfKernel();
    (*m)->setSelfKernel(0.0);
    // (*m)->addToSelfKernel(0.0);
    // --> bug fixed by PM (november 2007)
  }
}



/*
// normalize the gram matrix
void MoleculeSet::normalizeGram(){
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;

  for(mol1 = begin(); mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin() ; mol2 != comparisonSet->end() ; mol2++){
      //(*gramNormal)[i][j]  = (*gram)[i][j]/sqrt( (*mol1)->getSelfKernel() * (*mol2)->getSelfKernel() );
      (*gramNormal)[i][j]  = (*gram)[i][j]/sqrt( (*gram)[i][i] * (*gram)[j][j]  );
      j++;
    }
    i++;
  }
}

*/


// normalize the gram matrix
void MoleculeSet::normalizeGram(){
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;

  for(mol1 = begin(); mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin() ; mol2 != comparisonSet->end() ; mol2++){

      if( ( (*mol1)->getSelfKernel() == 0 ) ||  ( (*mol2)->getSelfKernel() == 0 )  ){
	(*gramNormal)[i][j]  = 0;	
      }
      else{
	(*gramNormal)[i][j]  = (*gram)[i][j]/sqrt( (*mol1)->getSelfKernel() * (*mol2)->getSelfKernel() );
      }
      j++;
    }
    i++;
  }

}


// normalize the gram matrix BASED ON THE RAW GRAM MATRIX (instead of self-kernels)
// WARNING: can only work in "train set mode", i.e., with a square gram matrix
void MoleculeSet::normalizeGram_raw(){
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;

  for(mol1 = begin(); mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin() ; mol2 != comparisonSet->end() ; mol2++){

      //if( ( (*mol1)->getSelfKernel() == 0 ) ||  ( (*mol2)->getSelfKernel() == 0 )  ){
      if( ((*gram)[i][i] == 0) ||  ((*gram)[j][j] == 0)  ){
	(*gramNormal)[i][j]  = 0;	
      }
      else{
	//(*gramNormal)[i][j]  = (*gram)[i][j]/sqrt( (*mol1)->getSelfKernel() * (*mol2)->getSelfKernel() );
	(*gramNormal)[i][j]  = (*gram)[i][j]/sqrt( (*gram)[i][i] * (*gram)[j][j] );
      }
      j++;
    }
    i++;
  }

}


// normalize the tanimoto coefficient
void MoleculeSet::normalizeTanimoto(){
  
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;
  
  for(mol1 = begin() ; mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin(); mol2 != comparisonSet->end() ; mol2++){
	if( ( (*mol1)->getSelfKernel() == 0 ) ||  ( (*mol2)->getSelfKernel() == 0 )  ){
		(*gramNormal)[i][j]  = 0;	
      	}
      	else{
 	     (*gramNormal)[i][j] = (*gram)[i][j]/ ( (*mol1)->getSelfKernel()  + (*mol2)->getSelfKernel() - (*gram)[i][j] );
	}
      j++;
    }
    i++;
  }

}



// normalize the tanimoto coefficient BASED ON THE RAW GRAM MATRIX (instead of self-kernels)
// WARNING: can only work in "train set mode", i.e., with a square gram matrix
void MoleculeSet::normalizeTanimoto_raw(){
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;

  for(mol1 = begin(); mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin() ; mol2 != comparisonSet->end() ; mol2++){

      //if( ( (*mol1)->getSelfKernel() == 0 ) ||  ( (*mol2)->getSelfKernel() == 0 )  ){
      if( ( (*gram)[i][i] == 0 ) ||  ( (*gram)[j][j] == 0 )  ){
	(*gramNormal)[i][j]  = 0;	
      }
      else{
	//(*gramNormal)[i][j]  = (*gram)[i][j]/( (*mol1)->getSelfKernel() + (*mol2)->getSelfKernel() - (*gram)[i][j]);
	(*gramNormal)[i][j]  = (*gram)[i][j]/( (*gram)[i][i] + (*gram)[j][j] - (*gram)[i][j]);
      }
      j++;
    }
    i++;
  }

}




// normalize the Min/Max tanimoto coefficient
void MoleculeSet::normalizeTanimotoMinMax(){
  vector<Molecule*>::iterator mol1;
  vector<Molecule*>::iterator mol2;
  int i = 0;
  int j = 0;
  
  for(mol1 = begin() ; mol1 != end() ; mol1++){
    j = 0;
    for(mol2 = comparisonSet->begin(); mol2 != comparisonSet->end() ; mol2++){
	if( (*gramNormal)[i][j] != 0){
	      (*gramNormal)[i][j] = (*gram)[i][j]/ (*gramNormal)[i][j];
	}
      j++;
    }
    i++;
  }
  
}


// add a value to the gram matrix at the given coordinates
void MoleculeSet::addToGram( int row, int col, double value ){
	(*gram)[row][col] += value;
}

// add a value to the normal gram matrix  at the given coordinates
void MoleculeSet::addToGramNormal( int row, int col, double value ){
	(*gramNormal)[row][col] += value;
}

// substract a value to the gram at the given coordinates
void MoleculeSet::substractToGram( int row, int col, double value ){
	(*gram)[row][col] -= value;
}


// return a gram matrix value
double MoleculeSet::getGramValue( int row, int col ){
	return (*gram)[row][col];
}


// set parameters to use only training data....
void  MoleculeSet::setComparisonSet(MoleculeSet *aSet){
	comparisonSet = aSet;
	subsetStart = 0;
}

// list the different kinds of atoms present in the set
vector<string> MoleculeSet::atomsLabelsListing(){

	vector<string> atomLabels;

	vector<Molecule*>::iterator molecule;
	// for each molecule : update atomLabels
	for( molecule = begin(); molecule != end(); molecule++ ){
		(*molecule)->atomsLabelsListing( &atomLabels );
	}

	for(molecule = comparisonSet->begin(); molecule != comparisonSet->end() ; molecule++){
		(*molecule)->atomsLabelsListing( &atomLabels );
	}

	return atomLabels;
}


// list the different kinds of atoms present in the set
vector<string> MoleculeSet::atomsSymbolsListing(){

	vector<string> atomSymbols;

	vector<Molecule*>::iterator molecule;
	// for each molecule : update atomSymbols
	for( molecule = begin(); molecule != end(); molecule++ ){
		(*molecule)->atomsSymbolsListing( &atomSymbols );
	}

	for(molecule = comparisonSet->begin(); molecule != comparisonSet->end() ; molecule++){
		(*molecule)->atomsSymbolsListing( &atomSymbols );
	}

	return atomSymbols;
}





// list the different kinds of bonds present in the set
vector<int> MoleculeSet::bondsListing(){

	vector<int> bondTypes;

	vector<Molecule*>::iterator molecule;
	// for each molecule : update bondTypes
	for( molecule = begin(); molecule != end(); molecule++ ){
		(*molecule)->bondsListing( &bondTypes );
	}

	for(molecule = comparisonSet->begin(); molecule != comparisonSet->end() ; molecule++){
		(*molecule)->bondsListing( &bondTypes );
	}


	return bondTypes;
}




long MoleculeSet::writeDotsToDir( string aDirName, bool selectedOnly, bool perretLabels )
{
	vector<Molecule*>::iterator m;

	int i = 0;
	string fileName = "";
	if( selectedOnly == false ){
		for( m = begin(); m != end(); m++ ){

			//cout << "next" << endl;
			//cout << (*m)->getName() << endl;

			if(  StringUtils::right( (*m)->getName(), 4 ) == ".dot" ){
				fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() );
			}else{
				fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() ) + ".dot";
			}

			//cout << "writing " << fileName << endl;
			(*m)->writeDOT( fileName, perretLabels );
			//cout << "done" << endl;
			i++;
		}
	}else{
		for( m = begin(); m != end(); m++ ){
			if( (*m)->isSelected() ){
				if(  StringUtils::right( (*m)->getName(), 4 ) == ".dot" ){
					fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() );
				}else{
					fileName = aDirName + "/" + StringUtils::slashToUnderscore( (*m)->getName() ) + ".dot";
				}
				//cout << "writing " << fileName << endl;
				(*m)->writeDOT( fileName );
				//cout << "done" << endl;
				i++;
			}
		}
	}
	return i;
}





void MoleculeSet::noTottersTransform(){  
  vector<Molecule*>::iterator molecule;
  //for each molecule
  for( molecule = begin(); molecule != end(); molecule++){
    (*molecule)->noTottersTransform();
  }


/*
// DEBUG (PM): 
int i = 0;
int j = 0;
  for( molecule = begin(); molecule != end(); molecule++){
    i++;
    cout << "Molecule no : " << i << " ; " << (*molecule)->numAtoms() << " atom.s" << endl;
    vector<Atom*>::iterator anAtom;
    j = 0;
    for( anAtom = (*molecule)->beginAtom(); anAtom != (*molecule)->endAtom(); anAtom++ )
      {
        j++;
        cout << "\t atom no " << j << " ; morgan label = " << (*anAtom)->getMorganLabel() << " ; ps =  " << (*anAtom)->getKashimaPS() << "  ; pq = "  << (*anAtom)->getKashimaPQ() << "  ;  AN = " <<  (*anAtom)->getAN()  <<endl;

         map<Atom*, Bond*>::iterator aBond;
        // for all the atom's bonds --> create a node
        for( aBond = (*anAtom)->beginBond(); aBond != (*anAtom)->endBond(); aBond++)
 	 {
 	   cout << "\t\t bond with prob :" <<  aBond->second->getKashimaPT() << endl;
 	 }

      }
  }
exit(142);
*/

}



// find smallest and largest distances between atoms
void MoleculeSet::minMaxDistances(double *distMin, double* distMax){
 
  // step1 : in moleculeSet
  vector<Molecule*>::iterator aMol;
  for( aMol = begin(); aMol != end(); aMol++){
    vector<Atom*>::iterator anAtom1, anAtom2;
    for(anAtom1 = (*aMol)->beginAtom() ; anAtom1 != (*aMol)->endAtom() - 1 ; anAtom1++ ){
      for(anAtom2 = anAtom1 + 1 ; anAtom2 != (*aMol)->endAtom() ; anAtom2++ ){
	float atomDist;
	atomDist = (*aMol)->atomicDistance(*anAtom1,*anAtom2);
	if(atomDist < *distMin)
	  *distMin = atomDist;
	if(atomDist > *distMax)
	  *distMax = atomDist;
      }
    }
  }

  // step2 : in comparisonSet
  if(this != comparisonSet){
    for( aMol = comparisonSet->begin(); aMol != comparisonSet->end(); aMol++){
      vector<Atom*>::iterator anAtom1, anAtom2;
      for(anAtom1 = (*aMol)->beginAtom() ; anAtom1 != (*aMol)->endAtom() - 1 ; anAtom1++ ){
	for(anAtom2 = anAtom1 + 1 ; anAtom2 != (*aMol)->endAtom() ; anAtom2++ ){
	  float atomDist;
	  atomDist = (*aMol)->atomicDistance(*anAtom1,*anAtom2);
	  if(atomDist < *distMin)
	    *distMin = atomDist;
	  if(atomDist > *distMax)
	  *distMax = atomDist;
	}
      }
    } 
  }

}


// transform the 2D molecule set into a 3D molecule set
// i.e., completely connected graphs, where edges labels are taken in {1,2,...,nBins}
// where labels correspond to a discretization of the distance range into nBins bins
void MoleculeSet::threeDtransform(int nBins, double distMin, double distMax){
  
  vector<Molecule*>::iterator aMol;
  //  transform all molecules
  for( aMol = begin(); aMol != end(); aMol++){
    (*aMol)->threeDtransform(nBins, distMin, distMax);
  }
  
}




void MoleculeSet::readPartialCharges(string fileName){
  
  vector<string> molCharges;
  ifstream fileInCharges;

  fileInCharges.open(fileName.c_str(), ios::in);
  if(!fileInCharges.good()){
    CError e = CError(FILENOTFOUND, fileName + " file not found");
    e.describe();
    throw(e);
  }
   
  int linesize = 1024;
  char *line = new char[linesize];
  
  while( !fileInCharges.eof() ){
    fileInCharges.getline(line, linesize-1,'\n');
    molCharges.push_back(line);
  }


if( molCharges[molCharges.size()-1] == ""){
	//cout << "EMPTY LINE!!!"<< endl;
	molCharges.pop_back();
}

  if(molCharges.size() != numMolecules()){
    cout << "ERROR : MoleculeSet::setPartialCharges" << endl;
    cout << "  --> number of lines read in charges file != numMolecules"  << endl;
    cout << "   - numMolecules = " << numMolecules() << " ; number of lines = " << molCharges.size()  << endl;
    cout << " last entry :  " << molCharges[molCharges.size()-1] << endl; 
   exit(12);
  }

  int i = 0;
  vector<Molecule*>::iterator m;
  for( m = begin(); m != end(); m++ ){
    (*m)->readPartialCharges(molCharges[i]);
    i++;
  }
 

  fileInCharges.close();
 
}


void MoleculeSet::setMorganChargesLabels(double threshold){
  vector<Molecule*>::iterator m;
  for( m = begin(); m != end(); m++ ){
    (*m)->setMorganChargesLabels(threshold);
  }



}


// ********************************* //
// ***** DEPRECATED FUNCTIONS ***** //
// ******************************** //


//normalize the gram matrix, i.e. fill in gramNormal
//void MoleculeSet::normalizeGram_self(){
//	for( int i = 0 ; i < numMolecules() ; i++ ){
//		for( int j = i ; j < numMolecules() ; j++ ){
//			(*gramNormal)[i][j] = (*gram)[i][j] / sqrt( (*gram)[i][i] * (*gram)[j][j] );
//			(*gramNormal)[j][i] = (*gram)[i][j] / sqrt( (*gram)[i][i] * (*gram)[j][j] );
//		}
//	}
//}


// // normalize the gram matrix, i.e. fill in gramNormal WITH RBF NORMALIZATION
// void MoleculeSet::normalizeGram_self(){
// 	double dist2, sigma2, ker;
// 	sigma2 = 1.0;
// 	for( int i = 0 ; i < numMolecules() ; i++ ){
// 		for( int j = i ; j < numMolecules() ; j++ ){
// 			dist2 = (*gram)[i][i] + (*gram)[j][j] - 2 * (*gram)[i][j];
// 			ker = exp( -dist2/sigma2 );
// 			(*gramNormal)[i][j] = ker;
// 			(*gramNormal)[j][i] = ker;
// 		}
// 	}
// }




// normalize the gram matrix, i.e. fill in gramNormal
//void MoleculeSet::normalizeGram_test(){
//vector<Molecule*>::iterator mol1;
//  vector<Molecule*>::iterator mol2;
//  int i = 0;
//  int j = 0;
  
//  for(mol1 = begin() ; mol1 != end() ; mol1++){
//    j = 0;
//    for(mol2 = comparisonSet->begin(); mol2 != comparisonSet->end() ; mol2++){
//      (*gramNormal)[i][j] = (*gram)[i][j]/sqrt( (*mol1)->getSelfKernel() * (*mol2)->getSelfKernel() );
//      j++;
//    }
//    i++;
//  }
//}


// // normalize the gram matrix, i.e. fill in gramNormal WITH RBF NORMALIZATION
// void MoleculeSet::normalizeGram_test(){
// 	vector<Molecule*>::iterator mol1;
// 	vector<Molecule*>::iterator mol2;
// 	int i = 0;
// 	int j = 0;
// 	double dist2, sigma2, ker;
// 	sigma2 = 1.0;
	
// 	for(mol1 = begin() ; mol1 != end() ; mol1++){
// 		j = 0;
// 		for(mol2 = comparisonSet->begin(); mol2 != comparisonSet->end() ; mol2++){
// 			dist2 = (*mol1)->getSelfKernel() + (*mol2)->getSelfKernel() - 2 * (*gram)[i][j];
// 			ker = exp( -dist2/sigma2 );
// 			(*gramNormal)[i][j] = ker;
// 			j++;
// 		}
// 		i++;
// 	}
// }



/*struct ToLower
{
  char operator() (char c) const  { return std::tolower(c); }
}
struct ToUpper
{
  char operator() (char c) const  { return std::toupper(c); }
}*/




/*void* kashimaKernelThread(void* arg){

	kashimaThreadStruct *str = (kashimaThreadStruct *) arg;
	*str->kernel = (str->molecule1)->newKernel( str->molecule2, str->convergence);

	float sk1 = (str->molecule1)->getSelfNewKernel();
	float sk2 = (str->molecule2)->getSelfNewKernel();

	*str->normalKernel = *str->kernel / sqrt(sk1 * sk2);

//	cout << "   comparing " << (str->molecule1)->getName() << " with " << (str->molecule2)->getName() << " " << "K = " << *str->kernel;
//	cout << " storing result in " << str << endl;

	str->noErrors = true;

	pthread_exit(0);

}  */


