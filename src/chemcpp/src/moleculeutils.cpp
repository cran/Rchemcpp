/****************************************************************************************
					  moleculeutils.cpp 
					---------------------
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
//#define VERBOSE 1

#include "moleculeutils.h"
#include <math.h>

//MoleculeUtils::MoleculeUtils(){
//}
//MoleculeUtils::~MoleculeUtils(){
//}


/**
	read the ctab block of a stream (connection table),
	and add atoms and bonds to aMolecule
*/
void MoleculeUtils::readMDLCtabBlock( Molecule& aMolecule, ifstream& inFile, bool genericAtomTypeFlag ) throw( CError ) {

	//cout << "reading connection table " << endl;

	int linesize = 256;
	char *line = new char[linesize];
	string stringLine;

	int nbAtoms = 0;
	int nbBonds = 0;
	int nbAtomlist = 0;
	int nbStext = 0;
	int chiral = 0;
	int nbProperties = 0;
	string versionString = "";
	int version = NAVALUE;

	int bondType = 0;

	Atom* anAtom = NULL;

	int lineNumber = 0;
	int j = 0;

	//read block line
	if( inFile.eof() ){
		stringstream out;
		out << "MoleculeUtils::readMDLCtabBlock: eof" << endl;
		CError e(EOFERROR, out.str() );
		//e.describe();
		throw(e);
	}
	inFile.getline( line, linesize-1,'\n' );
	stringLine = line;
	//cout << "# " << stringLine << endl;

	nbAtoms = atoi( stringLine.substr(0,3).c_str() );
	nbBonds = atoi( stringLine.substr(3,3).c_str() );
	nbAtomlist = atoi( stringLine.substr(6,3).c_str() );
	nbStext = atoi( stringLine.substr(15,3).c_str() );
	chiral = atoi( stringLine.substr(9,3).c_str() );
	nbProperties = atoi( stringLine.substr(30,3).c_str() );
	versionString = StringUtils::rmSpace( stringLine.substr( 33, stringLine.length() - 33 ) );
	#ifdef DEBUG
	  cout << "FOUND version " << versionString << endl;
	#endif
	if( versionString == "V2000" || versionString == "v2000" ){
		aMolecule.setOriginalFormat( MDL2000 );
		version = MDL2000;
	}else if( versionString == "V3000" || versionString == "v3000" ){
		aMolecule.setOriginalFormat( MDL3000 );
		version = MDL3000;
	}else if( versionString == "V2000JLP" || versionString == "v2000JLP" ){
		aMolecule.setOriginalFormat( MDL2000JLP );
		version = MDL2000JLP;
	}else{
	    /*stringstream out;
		out << "MoleculeUtils::readMDLCtabBlock: Bad file, cannot read " << versionString << " files" << endl;
		CError e(BADFILE, out.str() );
		e.describe();*/
		cout << "could not find MDL version number, infering this is a MDL 2000 file and continuing" << endl;
		//throw(e);
		aMolecule.setOriginalFormat( MDL2000 );
		version = MDL2000;
	}


	int bondStereo = 0;
	int bondNotUsed = 0;
	int bondTopology = 0;
	int bondReactionCenter = 0;
	int numRings = 0;

	if(nbProperties == 999){
		nbProperties = 0;
	}else{
		nbProperties--; // nbProperties usually includes the M  END line
	}

	//version = stringLine.substr(33,6);
	aMolecule.addStringDescriptor( "formatVersion", versionString, "", "mol file version" );

	//cout << nbAtoms << " atoms (56567)" << endl;

	aMolecule.setHasSSSR();

	//while( !inFile.eof() ){
	int atomC = 0;
	bool finished = false;
	while( lineNumber < nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties && finished == false ){
		lineNumber++;

		//read line
		if( inFile.eof() ){
			stringstream out;
			out << "MoleculeUtils::readMDLCtabBlock: eof" << endl;
			CError e(EOFERROR, out.str() );
			//e.describe();
			throw(e);
		}
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << "  " << stringLine << endl;

		// data
		if( nbAtoms > 0 ){

			if( j < nbAtoms ){
				// atoms
			 	//cout << "atom " << lineNumber << ": " << stringLine << endl;
				atomC++;
				if( genericAtomTypeFlag == false ){
					try{
						//cout << "MoleculeUtils::readMDLCtabBlock: " << stringLine.substr(30,3) << endl;
						anAtom = aMolecule.addAtom( StringUtils::rmSpace( stringLine.substr(30,3) ), false );
					}catch( CError E ){
						cerr << "->REPLACED WITH H" << endl;
						anAtom = aMolecule.addAtom( "H", false );
					}
				}else{
					//cout << "MoleculeUtils::readMDLCtabBlock: " << stringLine.substr(30,3) << endl;
					anAtom = new Atom( stringLine.substr(30,3) );
					aMolecule.addAtom( anAtom, false , true);
				}
				anAtom->setIdInMolecule( atomC );
				anAtom->setCoordinates( atof( stringLine.substr(0,10).c_str() ), atof( stringLine.substr(10,10).c_str() ), atof( stringLine.substr(20,10).c_str() ) );

				//j++;
			}
			if( j >= nbAtoms && j < nbAtoms + nbBonds ) {
				// bonds
        //cout << "bond " << lineNumber << ": " << stringLine << endl;
				bondType = atoi( stringLine.substr(6,3).c_str() );
				bondStereo = atoi( stringLine.substr(9,3).c_str() );
				bondNotUsed = atoi( stringLine.substr(12,3).c_str() );
				bondTopology = atoi( stringLine.substr(15,3).c_str() );
				bondReactionCenter = atoi( stringLine.substr(18,3).c_str() );

				Bond* newForward = aMolecule.linkAtoms( atoi( stringLine.substr(0,3).c_str() ) - 1, atoi( stringLine.substr(3,3).c_str() ) - 1, bondType, bondStereo, bondNotUsed, bondTopology, bondReactionCenter, false );
				Bond* newBackward = newForward->getReverse();

				if( bondTopology == 1 ){
					// check if the number of rings and the ring ids are appended at the end of line
					if( version == MDL2000JLP && stringLine.size() > 24 ){
						// read the number of rings it is member of
						numRings = atoi( stringLine.substr( 21, 3 ).c_str() );
						// read each ring id, adding the bond membership
						for( int j = 0; j != numRings; j++ ){
							if( stringLine.size() >= 24+j*3 ){
								int newID = atoi( stringLine.substr( 24+(j*3),3 ).c_str() );
								Ring* aRing = aMolecule.getRingWithID( newID, true );
								aRing->addBond( newForward, true );
								aRing->addBond( newBackward, true );
								aRing->addAtom( newForward->getSource(), true );
								aRing->addAtom( newForward->getTarget(), true );

								newForward->addRing( aRing );
								newBackward->addRing( aRing );
							}else{
								#ifdef DEBUG
								cout << "error reading ring membership, will recompute sssr" << endl;
								#endif
								aMolecule.resetSSSR();
							}

						}
					}else{
						#ifdef DEBUG
						cout << "no ring membership info, will recompute sssr" << endl;
						#endif

						aMolecule.resetSSSR();
					}
				}else if( bondTopology == 0){
					// sssr needs to be recomputed
					#ifdef DEBUG
						cout << "one bond's topology is 0, will recompute sssr" << endl;
					#endif

					aMolecule.resetSSSR();
				}

				//j++;
			}
			if( j >= nbAtoms + nbBonds && j < nbAtoms + nbBonds + nbAtomlist + nbStext ){
				// atomlist block and stext block (SKIPPING)
			}
			if( j >= nbAtoms + nbBonds + nbAtomlist + nbStext && j < nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties ){
				// properties block (SKIPPING)
			}
			//if( j >= nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties && j < nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties + 1){
			if( j >= nbAtoms + nbBonds + nbAtomlist + nbStext + nbProperties ){
				if(stringLine.substr(0,6) != "M  END" ){
					//cout << "MoleculeUtils::readConnectionTable: WARNING END LINE OF CONNECTION TABLE IS \"" << stringLine.substr(0,6) << "\" WHILE IT SHOULD BE: \"M  END\" " << endl;
					cout << "skiping " << stringLine.substr(0,6) << endl;
				}else{
					//cout << "found M  END" << endl;
					finished = true;
				}
			}
			j++;
		}
	}

	// there can be properties even when nbProperties is == 0 read until meeting M  END
	if( finished == false ){
		do{
			inFile.getline( line, linesize-1,'\n' );
			stringLine = line;
			if(stringLine.substr(0,6) == "M  END" ){
				//cout << "Prop finished " << endl;
				finished = true;
			}else{
				#ifdef VERBOSE
					cout << "  skipping: " << stringLine << endl;
				#endif
			}
		}while( finished == false );
	}

	delete[] line;

	if( nbAtoms < 1 ){
		// this molecule has no structure information, emit error
		stringstream out;
		out << "MoleculeUtils::readMDLCtabBlock: no structure information" << endl;
		CError e(NOSTRUCTURE, out.str() );
		//e.describe();
		throw(e);
	}

	aMolecule.compute();

}



void MoleculeUtils::skipMDLEntry( Molecule& aMolecule, ifstream& inFile ) throw( CError ) {

	//cout << "MoleculeUtils::skipMDLEntry: skiping entry" << endl;

	int linesize = 256;
	char *line = new char[linesize];
	string stringLine;

	int lineNumber = 0;
	int j = 0;

	bool finished = false;

	//read block line
	if( inFile.eof() ){
		stringstream out;
		out << "MoleculeUtils::readMDLCtabBlock: eof" << endl;
		CError e(EOFERROR, out.str() );
		//e.describe();
		throw(e);
	}

	while( finished == false ){
		lineNumber++;

		//read line
		if( inFile.eof() ){
			stringstream out;
			out << "MoleculeUtils::readMDLCtabBlock: eof" << endl;
			CError e(EOFERROR, out.str() );
			//e.describe();
			throw(e);
		}
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << "  " << stringLine << endl;

		if( stringLine.length() > 0){
			if( StringUtils::rmSpace( stringLine.substr(0,4) ) == "$$$$"){
				finished = true;
			}
		}
	}

	delete[] line;
}





/**
	read the header block of a stream,
	and set molecule name to the comment content or to aName if specified
*/

void MoleculeUtils::readMDLHeaderBlock( Molecule& aMolecule, ifstream& inFile, string aName ) throw( CError ){
	int linesize = 256;
	char *line = new char[linesize];
	string stringLine;

	Descriptor<string>* ds;
	//Descriptor<int>* di;
	//Descriptor<float>* df;



	// read comment block
	if( inFile.eof() ){
		stringstream out;
		out << "MoleculeUtils::readMDLHeaderBlock: eof" << endl;
		CError e(EOFERROR, out.str() );
		//e.describe();
		throw(e);
	}
	inFile.getline( line, linesize-1,'\n' );
	stringLine = line;

	#ifdef DEBUG
		cout << " adding comment descriptor " << endl;
	#endif

	ds = aMolecule.setStringDescriptor( "comment", stringLine, "", "", true, true );


	//cout << "comment: " << ds->getValue() << endl;

	//if( stringLine == "" ){
	//	cout << "comment is empty " << endl;
	//	ds->setEmpty();
	//}

	#ifdef DEBUG
		cout << "done, setting mol name" << endl;
	#endif

	if( aName != "COMMENT" ){
		#ifdef DEBUG
			cout << "filename: " << StringUtils::getFileName( aName ) << endl;
			//cout << "setting name to " << StringUtils::getNoExtension( StringUtils::getFileName( aName ) ) << endl;
		#endif
		aMolecule.setName( StringUtils::getNoExtension( StringUtils::getFileName( aName ) ) );
		//aMolecule.setName( StringUtils::getNoExtension( aName ) );
	}else{
		//cout << "setting name to " << stringLine << endl;
		aMolecule.setName( stringLine );
	}

	#ifdef DEBUG
		cout << " adding comment2 descriptor " << endl;
	#endif


	if( inFile.eof() ){
		stringstream out;
		out << "MoleculeUtils::readMDLHeaderBlock: eof" << endl;
		CError e(EOFERROR, out.str() );
		//e.describe();
		throw(e);
	}
	inFile.getline( line, linesize-1,'\n' );
	stringLine = line;

	//cout << "comment2: " << stringLine << endl;
	ds = aMolecule.setStringDescriptor( "comment2", stringLine, "", "", true, true );
	//if( stringLine == "" ){
		//cout << "comment2 is empty " << endl;
		//ds->setEmpty();
	//}

	#ifdef DEBUG
		cout << "done, adding comment 3 descriptor" << endl;
	#endif

	if( inFile.eof() ){
		stringstream out;
		out << "MoleculeUtils::readMDLHeaderBlock: eof" << endl;
		CError e(EOFERROR, out.str() );
		//e.describe();
		throw(e);
	}
	inFile.getline( line, linesize-1,'\n' );
	stringLine = line;

	if( stringLine.substr( 0, 21 ) == "Copyright by the U.S." ){
		aMolecule.setStringDescriptor( "comment3", "", "", "", true, true );
		//cout << "comment3 is empty " << endl;
		//ds->setEmpty();
	}else{
		aMolecule.setStringDescriptor( "comment3", stringLine, "", "", true, true );
		//if( stringLine == "" ){
		//	ds->setEmpty();
		//	cout << "comment3 is empty " << endl;
		//}
	}


	// aMolecule.describe();
	#ifdef DEBUG
		cout << "finished reading header " << endl;
	#endif

	// comment (skipping line 3)
	delete[] line;
}

void MoleculeUtils::writeMDLCtabBlock( Molecule& aMolecule, ofstream& outFile ){
	outFile << StringUtils::preFill( StringUtils::toString( aMolecule.numAtoms() ), 3, " " );
	outFile << StringUtils::preFill( StringUtils::toString( aMolecule.numBonds() ), 3, " " );
	outFile << StringUtils::preFill( "0", 3, " " ); // atom lists
	outFile << StringUtils::preFill( " ", 3, " " ); // obsolete
	if( aMolecule.isChiral() ){
		outFile << StringUtils::preFill( "1", 3, " " ); // chiral
	}else{
		outFile << StringUtils::preFill( "0", 3, " " ); // not chiral
	}
	outFile << StringUtils::preFill( "0", 3, " " ); // stext
	outFile << StringUtils::preFill( "1", 15, " " ); // stext
	outFile << " V2000JLP" << endl;

	// write atom block
	vector<Atom*>::iterator ait;
	int id = 1;
	for( ait = aMolecule.beginAtom(); ait != aMolecule.endAtom(); ait++ ){
		(*ait)->setIdInMolecule( id );
		outFile <<  setw(10) << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << (*ait)->getX();
		outFile <<  setw(10) << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << (*ait)->getY();
		outFile <<  setw(10) << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << (*ait)->getZ();

		outFile << " ";
		outFile << StringUtils::fill( (*ait)->getElementSymbol(), 3, " "); // symbol
		//outFile << StringUtils::fill( (*ait)->getSymbol(), 3, " "); // symbol

		outFile << StringUtils::preFill( "0", 2, " "); // mass diff
		outFile << StringUtils::preFill( "0", 3, " "); // charge diff
		outFile << StringUtils::preFill( "0", 3, " "); // stereo parity
		outFile << StringUtils::preFill( "0", 3, " "); // hydrogen counts
		outFile << StringUtils::preFill( "0", 3, " "); // stereocare box
		outFile << StringUtils::preFill( "0", 3, " "); // valence
		outFile << StringUtils::preFill( " ", 3, " "); // H0 designator
		outFile << StringUtils::preFill( " ", 3, " "); // reaction component type
		outFile << StringUtils::preFill( " ", 3, " "); // reaction component number
		outFile << StringUtils::preFill( "0", 3, " "); // atom-atom mapping number
		outFile << StringUtils::preFill( "0", 3, " "); // inversion/retention flag
		outFile << StringUtils::preFill( "0", 3, " "); // exact change flag
		outFile << endl;
		id++;
	}

	// write bond block
	map<Atom*, Bond*>::iterator bit;
	for( ait = aMolecule.beginAtom(); ait != aMolecule.endAtom(); ait++ ){
		for( bit = (*ait)->beginBond(); bit != (*ait)->endBond(); bit++ ){
			if( (*bit).second->hasFlagOriginal() ){
				outFile << StringUtils::preFill( (*bit).second->getSource()->getIdInMolecule(), 3, " " );
				outFile << StringUtils::preFill( (*bit).first->getIdInMolecule(), 3, " " );
				outFile << StringUtils::preFill( (*bit).second->getLabel(), 3, " " );
				outFile << StringUtils::preFill( (*bit).second->getStereo(), 3, " "); // bond stereo

				outFile << StringUtils::preFill( (*bit).second->getNotUsed(), 3, " "); // notUsed field

				if( aMolecule.hasSSSRDetected() ){
					// the molecule has ring information, save it
					if( (*bit).second->hasRing() ){
						outFile << "  1"; // bond topology: 1: ring
					}else{
						outFile << "  2"; // bond topology: 2: chain
					}
				}else{
					outFile << "  0"; // bond topology: 0 either bond or chain
				}
				outFile << StringUtils::preFill( (*bit).second->getReactionCenter(), 3, " "); // reaction center

				// add the ring information at the end of the line

				outFile << MoleculeUtils::getRingString( (*bit).second );

				outFile << endl;
			}
		}
	}

	outFile << "M  END" << endl;

}

string MoleculeUtils::getRingString( Bond* aBond ) throw( CError ){
	// NNNRRR... where NNN is the number of rings this bond is part of, RRR is the first ring id. other rings membership are appended totalising NNN*RRR columns. Ring id should start with 1.
	stringstream out;

	//cout << "OUTPUT " << aBond->numRings() << " RINGS " << endl;

	try{
		out << StringUtils::preFill( aBond->numRings(), 3, " " );

		if( aBond->numRings() > 0 ){
			vector<Ring*>::iterator ri;
			for( ri = aBond->beginRing(); ri != aBond->endRing(); ri++ ){
				out << StringUtils::preFill( (*ri)->getID(), 3, " " );
			}
		}
		//cout << out.str() << endl;
		return( out.str() );
	}catch( CError e ){
		e.describe();
		cout << "MoleculeUtils::getRingString: Error, sssr was not detected, please call Molecule::detectSSSR() before calling me" << endl;
		throw( e );
	}

}

/** read the non structural data block of a compound in a MDL SDfile
	WARNING only data entry with a structure > <aLabel> will be considered. entries
		missing the <aLabel> entry will be ignored
	*/
void MoleculeUtils::readMDLNSDBlock( Molecule& aMolecule, ifstream& inFile ) throw( CError ){

	int linesize = 256;
	char *line = new char[linesize];
	string stringLine = "";

	string aLabel = "";
	string aValue = "";
	int aType = STRING;

	bool nameSet = false;

	int lc = 0;

	while( stringLine != "$$$$" ){

		lc++;

		if( inFile.eof() ){
			delete[] line;
			stringstream out;
			out << "MoleculeUtils::readMDLNSDBlock: eof" << endl;
			CError e(EOFERROR, out.str() );
			//e.describe();
			throw(e);
		}
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << "1: " << stringLine << endl;

		if( stringLine.substr( 0, 1 ) != ">" ){
			delete[] line;
			if( stringLine.substr( 0, 4 ) == "$$$$" ){
				return;
			}
			// bad formated file or 2 values for one descriptor
			stringstream out;
			out << "MoleculeUtils::readMDLNSDBlock: > not found in descriptor " << aLabel;
			out << " in molecule " << aMolecule.getName() << "; line: ";
			out << " (the sd import function does not support descriptors with multiple values)" << endl;
			CError e(BADFILE, out.str() );
			e.describe();
			throw(e);

		}else{
			// look for <> content
			uint dls = 0;
			uint dle = 0;
			dls = stringLine.find("<", 1);
			dle = stringLine.find(">", 1);

			if( dls == string::npos || dls == string::npos ){
				stringstream out;
  				out << "MoleculeUtils::readMDLNSDBlock: <aLabel> not found in line " << stringLine << endl;
				delete[] line;
				CError e(BADFILE, out.str() );
				e.describe();
				throw(e);
			}else{
				aLabel = stringLine.substr( dls + 1, (dle - dls) - 1 );
				//cout << "found descriptor " << aLabel << endl;

				if( aLabel.substr( aLabel.length() - min( (int) aLabel.length(), 7), 7 ) == "integer" || aLabel.substr( aLabel.length() - min( (int) aLabel.length(), 4), 4 ) == ".int"){
					aLabel = aLabel.substr( 0, aLabel.length() - 8 );
					aType = INTEGER;
				}else if( aLabel.substr( aLabel.length() - min( (int) aLabel.length(), 5), 5 ) == "float" || aLabel.substr( aLabel.length() - min( (int) aLabel.length(), 4), 4 ) == ".flo"){
					aLabel = aLabel.substr( 0, aLabel.length() - 6 );
					aType = FLOAT;
				}else if( aLabel.substr( aLabel.length() - min( (int) aLabel.length(),6), 6 ) == "string" || aLabel.substr( aLabel.length() - min( (int) aLabel.length(),4), 4 ) == ".str" ){
					aLabel = aLabel.substr( 0, aLabel.length() - 7 );
					aType = STRING;
				}else{
					aType = STRING;
				}
			}

			// look for () content
			uint pls = 0;
			uint ple = 0;
			pls = stringLine.find("(", 1);
			ple = stringLine.find(")", 1);

			// there is a registry number in parenthesis. Use this name for the compound name
			if( aType == STRING && pls != string::npos && ple != string::npos && pls != ple ){
				aMolecule.setName( stringLine.substr( pls + 1, (ple - pls) - 1 ) );
				nameSet = true;
			}

		}

		// read value lines

		if( inFile.eof() ){
			delete[] line;
			stringstream out;
			out << "MoleculeUtils::readMDLNSDBlock: eof" << endl;
			CError e(EOFERROR, out.str() );
			//e.describe();
			throw(e);
		}
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << "2: " << stringLine << endl;

		if( aType == STRING ){
			Descriptor<string>* d = NULL;
			if( aLabel == "name" ){
			    if( nameSet == false ){
				aMolecule.setName( stringLine );
				nameSet = true;
			    }
			}else if( aLabel == ACTIVITY ){
			    aMolecule.setActivity( stringLine );
			}else{
			  d = aMolecule.setStringDescriptor( aLabel, stringLine, "", "", true, true );
			  if( stringLine == "" ){
			    d->setEmpty();
			  }
			}
		}else if( aType == INTEGER ) {
			Descriptor<int>* d = NULL;
			if( aLabel == ACTIVITY ){
			    aMolecule.setActivity( stringLine );
			}else{
			    //cout << "adding " << aLabel << endl;
			    d = aMolecule.setIntDescriptor( aLabel, atoi( stringLine.c_str() ), "", "", true, true );

			    if( stringLine == "" ){
				d->setEmpty();
			    }
			}
		}else if( aType == FLOAT ) {
			Descriptor<float>* d = NULL;
			if( aLabel == ACTIVITY ){
			    aMolecule.setActivity( stringLine );
			}else{
			    d = aMolecule.setFloatDescriptor( aLabel, atof( stringLine.c_str() ), "", "", true, true );
			    if( stringLine == "" ){
				d->setEmpty();
			    }
			}
		}

		// read blank line
		if( inFile.eof() ){
			delete[] line;
			stringstream out;
			out << "MoleculeUtils::readMDLNSDBlock: eof" << endl;
			CError e(EOFERROR, out.str() );
			//e.describe();
			throw(e);
		}
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << "3: " << stringLine << endl;


	}
	delete[] line;

}


void MoleculeUtils::writeMDLHeaderBlock( Molecule& aMolecule, ofstream& outFile ){
	//int flagNameWritten = 0;

	string s = "";
	//aMolecule.describe();

	try{
		s = aMolecule.getStringDescriptor( "comment" )->getValue( true );
	}catch( CError e ){

		//e.describe();
		//cerr << "error accessing comment" << endl;
		// comment is empty
		//try{
		//	// try get the name
		//	s = aMolecule.getStringDescriptor( "name" )->getValue();
		//	cout << " writing name instead " << endl;
		//	flagNameWritten++;
		//}catch( CError e ){
		//	// name is empty or does not exist
		//	cout << " writing nothing instead " << endl;
		s = "";
		//}

	}
	outFile << s << endl;

	try{
		s = aMolecule.getStringDescriptor( "comment2" )->getValue( true );
	}catch( CError e ){
		//e.describe();
		//cerr << "error accessing comment2" << endl;
		// comment is empty
		//if( flagNameWritten == 0){

		//	try{
		//		// try get the name
		//		s = aMolecule.getStringDescriptor( "name" )->getValue();
		//		cout << " writing name instead " << endl;
		//		flagNameWritten++;
		//	}catch( CError e ){
		//		// name is empty or does not exist
		//		s = "";
		//		cout << " writing nothing instead " << endl;
		//	}
		//}else{
		//	cout << " writing nothing instead " << endl;
		s = "";
		//}
	}
	outFile << s << endl;

	try{
		s = aMolecule.getStringDescriptor( "comment3" )->getValue( true );
	}catch( CError e ){
		//e.describe();
		//cerr << "error accessing comment3" << endl;
		// comment is empty
		//if( flagNameWritten == 0){

		//	try{
		//		// try get the name
		//		s = aMolecule.getStringDescriptor( "name" )->getValue();
		//		cout << " writing name instead " << endl;
		//		flagNameWritten++;
		//	}catch( CError e ){
		//		// name is empty or does not exist
		//		cout << " writing nothing instead " << endl;
		//		s = "";
		//	}
		//}else{
		//	cout << " writing nothing instead " << endl;
		 s = "";
		//}
	}
	outFile << s << endl;




	//if( aMolecule.hasStringDescriptor( "comment" ) ){
	//	outFile << aMolecule.getStringDescriptor( "comment" )->getValue() << endl;
	//}else{
	//	outFile << aMolecule.getStringDescriptor( "name" )->getValue() << endl;
	//	flagNameWritten++;
	//}

	//if( aMolecule.hasStringDescriptor( "comment2" ) ){
	//	outFile << aMolecule.getStringDescriptor( "comment2" )->getValue() << endl;
	//}else{
	//	if( flagNameWritten == 0){
	//		outFile << aMolecule.getStringDescriptor( "name" )->getValue() << endl;
	//		flagNameWritten++;
	//	}else{
	//		outFile << endl;
	//	}
	//}
	//if( aMolecule.hasStringDescriptor( "comment3" ) ){
	//	outFile << aMolecule.getStringDescriptor( "comment3" )->getValue() << endl;
	//}else{
	//	if( flagNameWritten == 0){
	//		outFile << aMolecule.getStringDescriptor( "name" )->getValue() << endl;
	//		flagNameWritten++;
	//	}else{
	//		outFile << endl;
	//	}
	//}

	//outFile << endl << endl;
}

/** writes the non structural data block to a stream, for aMolecule
  	this function writes all stringDescriptors, floatDescriptors and IntDescriptors
		of a molecule in appropriate format
	*/
void MoleculeUtils::writeMDLNSDBlock( Molecule& aMolecule, ofstream& outFile ){
	map< string, Descriptor< int >* >::iterator iit;
	map< string, Descriptor< float >* >::iterator fit;
	map< string, Descriptor< string >* >::iterator sit;

	Descriptor< int >* di;
	Descriptor< float >* df;
	Descriptor< string >* ds;


	// if the molecule has activity write the activity descriptor
	if( aMolecule.hasActivity() ){
		outFile << "> <" << "activity>" << endl << aMolecule.getActivity() << endl << endl;
	}

	int i;
	float f;
	string s;

	for( iit = aMolecule.beginIntDescriptor(); iit != aMolecule.endIntDescriptor(); iit++ ){
		di = (*iit).second;
		try{
			i = di->getValue( true ); // silentError
			outFile << "> <" << di->getLabel();
			outFile << ".integer";
			outFile << ">" << endl;
			outFile << i << endl;
			outFile << endl;
		}catch( CError e ){
			// descriptor is empty, do nothing
		}
	}

	for( fit = aMolecule.beginFloatDescriptor(); fit != aMolecule.endFloatDescriptor(); fit++ ){
		df = (*fit).second;

		try{
			f = df->getValue( true ); // silentError
			outFile << "> <" << df->getLabel();
			outFile << ".float";
			outFile << ">" << endl;
			outFile << f << endl;
			outFile << endl;
		}catch( CError e ){
			// descriptor is empty, do nothing
		}
	}

	for( sit = aMolecule.beginStringDescriptor(); sit != aMolecule.endStringDescriptor(); sit++ ){
		ds = (*sit).second;
		if( ds->getLabel() != "comment" &&
			ds->getLabel() != "comment2" &&
			ds->getLabel() != "comment3" &&
			ds->getLabel() != "formatVersion" ){

			try{
				s = ds->getValue( true ); // silentError
				outFile << "> <" << ds->getLabel();
				outFile << ".string";
				outFile << ">" << endl;
				outFile << s << endl;
				outFile << endl;
			}catch( CError e ){
				// descriptor is empty, do nothing
			}
		}
	}
	outFile << "$$$$" << endl;
}


void MoleculeUtils::writeKCF( Molecule& aMolecule, ofstream& outFile ){
	outFile << StringUtils::fill( "ENTRY", 12 ) <<
	StringUtils::fill( aMolecule.getName(), 28 ) << "Compound  #chiral" << endl;

	// write ATOMS
	outFile << StringUtils::fill( "ATOM", 12) << aMolecule.numAtoms() << endl;
	vector<Atom*>::iterator ait;
	int i = 1;
	stringstream istream;
	string istring = "";
	string xstring = "";
	string ystring = "";

	int minId = (*aMolecule.beginAtom())->getId()-1;

	for( ait = aMolecule.beginAtom(); ait != aMolecule.endAtom(); ait++ ){
		istream << (*ait)->getId() - minId;
		istream >> istring;
		istream.clear();

		istream << (*ait)->getX();
		istream >> xstring;
		istream.clear();

		istream << (*ait)->getY();
		istream >> ystring;
		istream.clear();

		outFile << StringUtils::fill(" ", 12) << StringUtils::fill( istring , 4) << StringUtils::fill( (*ait)->getSymbol() , 4) <<
			StringUtils::fill( (*ait)->getStringDescriptor("ElementSymbol")->getValue() , 5) << StringUtils::preFill( xstring , 7) << StringUtils::preFill( ystring , 10) << endl;

		i++;
	}

	// write BONDS

	outFile << StringUtils::fill( "BOND", 12) << aMolecule.numBonds() << endl;

	map<Atom*, Bond*>::iterator aib;
	i = 1;
	string source = "";
	string target = "";
	string label = "";
	for( ait = aMolecule.beginAtom(); ait != aMolecule.endAtom(); ait++ ){
		for( aib = (*ait)->beginBond(); aib != (*ait)->endBond(); aib++ ){

			istream << i;
			istream >> istring;
			istream.clear();

			istream << (*aib).second->getSource()->getId() - minId;
			istream >> source;
			istream.clear();

			istream << (*aib).first->getId() - minId;
			istream >> target;
			istream.clear();

			istream << (*aib).second->getLabel();
			istream >> label;
			istream.clear();

			// set flag on back bond


			if( !(*aib).second->hasFlagOriginal() ){
				outFile << StringUtils::fill(" ", 12) << StringUtils::fill( istring , 3) <<
					StringUtils::preFill( source , 4) <<
					StringUtils::preFill( target , 4) <<
					StringUtils::preFill( label , 2) << endl;
			i++;
			}
		}
	}

	// write other descriptors
	MoleculeUtils::writeKCFNSDBlock( aMolecule, outFile );
	outFile << "///" << endl;

}


/** writes the non structural data block to a stream, for aMolecule in KCF format
	this function writes all stringDescriptors, floatDescriptors and IntDescriptors
	of a molecule in appropriate format
*/
void MoleculeUtils::writeKCFNSDBlock( Molecule& aMolecule, ofstream& outFile ){
	map< string, Descriptor< int >* >::iterator iit;
	map< string, Descriptor< float >* >::iterator fit;
	map< string, Descriptor< string >* >::iterator sit;

	Descriptor< int >* di = NULL;
	Descriptor< float >* df = NULL;
	Descriptor< string >* ds = NULL;

	int i;
	float f;
	string s;
	string outLabel;

	for( iit = aMolecule.beginIntDescriptor(); iit != aMolecule.endIntDescriptor(); iit++ ){
		di = (*iit).second;
		try{
			i = di->getValue();
			outLabel = di->getLabel().substr( 0, 8 ) + ".int";
			outFile << StringUtils::fill( outLabel, 12, " " );
			outFile << i << endl;
		}catch( CError e ){
			// descriptor is empty, do nothing
		}
	}

	for( fit = aMolecule.beginFloatDescriptor(); fit != aMolecule.endFloatDescriptor(); fit++ ){
		df = (*fit).second;
		try{
			f = df->getValue();
			outLabel = df->getLabel().substr( 0, 8 ) + ".flo";
			outFile << StringUtils::fill( outLabel, 12, " " );
			outFile << f << endl;
		}catch( CError e ){
			// descriptor is empty, do nothing
		}
	}

	for( sit = aMolecule.beginStringDescriptor(); sit != aMolecule.endStringDescriptor(); sit++ ){
		ds = (*sit).second;
		if( ds->getLabel() != "comment" &&
			ds->getLabel() != "comment2" &&
			ds->getLabel() != "comment3" &&
			ds->getLabel() != "formatVersion" ){

			try{
				s = ds->getValue();
				outLabel = ds->getLabel().substr( 0, 8 ) + ".str";
				outFile << StringUtils::fill( outLabel, 12, " " );
				outFile << s << endl;
			}catch( CError e ){
				// descriptor is empty, do nothing
			}
		}
	}
	//outFile << "$$$$" << endl;
}

/** read one molecule from a kcf file. returns true if a valid molecule was found (containing
at least one atom.
*/
bool MoleculeUtils::readKCFMolecule( KCFMolecule& m, ifstream& inFile ) throw( CError ){
	// COMPLETE ***

	// erase all eventually existing atoms in the molecule
	m.erase();

	bool valid = false;

	int linesize = 256;
	char *line = new char[linesize];
	string stringLine;
	string lineLabel;
	string lineContent;
	int readerStatus = FINDENTRY;

	int nbAtoms = 0;
	int nbBonds = 0;
	string KCFType = "";

	m.unsetBondFlagsOriginal();

	Atom* anAtom = NULL;

	int lineNumber = 0;
	// for each line in file

	int atomC = 0;

	while( !inFile.eof()){

		lineNumber++;
		//cout << lineNumber << endl;

		//read line
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;
		//cout << stringLine.length() << " " << stringLine << endl;
		if(stringLine.length()>11){

			// parse line
			lineLabel = StringUtils::rmSpace( stringLine.substr( 0, 11 ) );
			lineContent = stringLine.substr( 12, stringLine.length() - 12 );

			//cout << "line: " << stringLine << endl;
			//cout << "label  :" << lineLabel << "." << readerStatus << endl;
			//cout << "content:" << lineContent << "." << readerStatus << endl;

			if(lineLabel == "ENTRY"){
			    //cout << "FOUND ENTRY:" << StringUtils::rmTailSpace( lineContent ) << endl;
			    //cout << " setting name to " << StringUtils::rmSpace( StringUtils::field( lineContent, 0, 18 ) ) << endl;
				m.setName( StringUtils::rmSpace( StringUtils::field( lineContent, 0, 18 ) ) );
				KCFType = StringUtils::field( stringLine, 30, 8 );
			}

			//cout << "checking node" << endl;
			if(lineLabel == "NODE" || lineLabel == "ATOM"){
				//cout << "FOUND NODE"<< endl;
				valid = true;
				readerStatus = READINGNODES;
				nbAtoms = atoi( stringLine.substr(12,stringLine.length()-12).c_str() );
				//cout << "this molecule has " << nbAtoms << " atoms" << endl;
			}else if(lineLabel == "EDGE" || lineLabel == "BOND"){
				//cout << "FOUND EDGE" << endl;
				readerStatus = READINGEDGES;
				nbBonds = atoi( stringLine.substr(12,stringLine.length()-12).c_str() );
				//cout << "this molecule has " << nbBonds << " bonds" << endl;

			}else if( lineLabel == "" ){
				//cout << "label is empty" << endl;

				if( readerStatus == READINGNODES ){
						//cout << "p1" << endl;

						if(lineLabel!=""){
							readerStatus = FINDENTRY;
						}else{
							// read node
							//cout << "MoleculeUtils::readKCFMolecule:  adding atom " << StringUtils::rmSpace( stringLine.substr(16,3) ) << endl;
							atomC++;
							anAtom = m.addAtom( StringUtils::rmSpace( stringLine.substr(16,3) ) );
							anAtom->setIdInMolecule( atomC );
							//cout << " coordinates: " << atof( stringLine.substr(23,9).c_str() ) << "; " << atof( stringLine.substr(33,9).c_str() ) << endl;
							anAtom->setCoordinates( atof( stringLine.substr(23,9).c_str() ), atof( stringLine.substr(33,9).c_str() ), 0 );
							anAtom->setElementSymbol( stringLine.substr(20,2) );

							//cout << anAtom->getX() << " " << anAtom->getY() << endl;
						}
				}else if( readerStatus == READINGEDGES ){
					//cout << "p2" << endl;

						if(lineLabel!=""){
							readerStatus = FINDENTRY;
						}else{
							// read edge
							//cout << "linking " << atoi( stringLine.substr(15,3).c_str() ) - 1 << " -" << atoi( stringLine.substr(23,1).c_str() ) << "- " << atoi( stringLine.substr(19,3).c_str() ) - 1 << endl;
							m.linkAtoms( atoi( stringLine.substr(16,3).c_str() ) - 1, atoi( stringLine.substr(20,3).c_str() ) - 1, atoi( stringLine.substr(24,1).c_str() ) );
							//cout << "done" << endl;
						}
				}else{
					//cout << "p3" << endl;
					//cout << "doing nothing" << endl;
				}

				//cout << "p4" << endl;
			}else{
				// unknown tag met
				if( lineLabel != "" && lineLabel != "ENTRY" ){
					//cout << "unknown tag: " << lineLabel << endl;
					//cout << "adding descriptor " << lineLabel << " with value " << lineContent << endl;
					m.addUnknownTypeDescriptor( lineLabel, lineContent );
					readerStatus = FINDENTRY;
				}else if( lineLabel != "ENTRY" ){
					cout << "skipping unlabeled descriptor with value " << lineContent << endl;
				}
			}
		}else{
			if(stringLine.substr(0,3) == "///"){
				//cout << "END OF MOLECULE " << endl;
				if( valid ){
					int i = m.hideHydrogens();
					if( i > 0 ){
						cout << i << " H atoms hidden" << endl;
					}
				}
				return( valid );
			}else if( stringLine.substr(0,2) == "//" ){
			}else if( stringLine.substr(0,4) == "$$$$" ){
			}else{
				if( stringLine!="" ){
					cout << "skipping short line: " << stringLine << endl;
				}
			}
		}

	}
	if( valid ){
		cout << m.hideHydrogens() << " H atoms hidden" << endl;
	}

	delete[] line;

	m.compute();

  return( valid );
}

double MoleculeUtils::moleculeKernel(
					Molecule* mol1, Molecule* mol2,
					double (*pt2AtomKernel)(Atom*, Atom*),
					double (*pt2BondKernel)(Bond*, Bond*),
					int convergenceCondition, int parameter2
					){

	// we need an array containing values calculated for each atom pairs among the two molecules

	#ifdef DEBUG
		cout << "MoleculeUtils::moleculeKernel " << convergenceCondition << endl;
	#endif
	if( convergenceCondition <= 0 ){
		stringstream out;
		out << "MoleculeUtils::moleculeKernel: Invalid convergence condition: should be > 0 and it is = " << convergenceCondition << endl;
		CError e( VALUENOTALLOWED, out.str() );
		e.describe();
		throw(e);
	}

	vector< vector<double> >* r = new vector< vector<double> >;
	vector< vector<double> >* rwork = new vector< vector<double> >;
	vector< vector<double> >* rstart = new vector< vector<double> >;  // contains the values at length parameter2

	//vector< vector<float> >* rFinal;

	vector<Atom*>::iterator atom;
	vector<Atom*>::iterator otherAtom;

	// initial Array filled with the probability to die at the first step
	int i = 0;
	int j = 0;

	//cout << "initializing matrix" << endl;

	for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
		// add each line
		//cout << "MoleculeUtils::moleculeKernel: i = " << i << endl;
		r->push_back( vector<double>() );
		rwork->push_back( vector<double>() );
		rstart->push_back( vector<double>() );
		j = 0;
		for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){
			//cout << "MoleculeUtils::moleculeKernel: j = " << j << endl;
		  // fill line with the values for each column
			try{
				//if( parameter2 < 2 ){
					(*r)[i].push_back( (*atom)->getKashimaPQ() * (*otherAtom)->getKashimaPQ() );
				//}else{
					//cout << "gghgh 0 " << endl;
					//(*r)[i].push_back( 0 );
				//}
			}catch(CError e){
				cout << "MoleculeUtils::graphKernel(): random walk probability was not set" <<endl;
				cout << "use Molecule::setKashimaKernelProb( int ) before calculating the kernel" << endl << endl;
				throw( e );
			}

      //cout << "MoleculeUtils::moleculeKernel: F1" << endl;
			(*rwork)[i].push_back( 1000.0 );
			(*rstart)[i].push_back( 0.0 );

			//cout << "MoleculeUtils::moleculeKernel: F2" << endl;
			(*otherAtom)->setRPosition( j );
			//cout << "MoleculeUtils::moleculeKernel: F3" << endl;

			j++;
		}

		(*atom)->setRPosition( i );
		i++;
	}

	//cout << "MoleculeUtils::moleculeKernel: start rlk" << endl;

  // calcul de la somme des probabilites de tous les chemins possibles communs aux deux molecules partant de tous les noeuds
	MoleculeUtils::rlk( r, rwork, rstart, mol1, mol2, convergenceCondition, parameter2, pt2AtomKernel, pt2BondKernel, 1 );
	//cout << "MoleculeUtils::moleculeKernel: rlk done" << endl;

	double s = 0.0;
	double k = 0.0;
	i = 0;

	//int localAN = 0;
	double localPS = 0;

	for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
		j = 0;
		//localAN = (*atom)->getAN()-1;
		//localAN = (*atom)->getType();
		localPS = (*atom)->getKashimaPS();
		for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){
			s = pt2AtomKernel( (*atom), (*otherAtom) );
			if( s != 0 ){
				if( parameter2 < 2 ){

					s = localPS * (*otherAtom)->getKashimaPS() * s;
					k = k + ( s * (*rwork)[i][j] );
				}else{
					k = k + ( (*rwork)[i][j] );
				}
			}
			j++;
		}
		i++;
	}

	delete r;
	delete rstart;
	delete rwork;

	return(k);
}

vector< vector<double> >* MoleculeUtils::rlk(
						vector< vector<double> >* r,
						vector< vector<double> >* rwork,
						vector< vector<double> >* rstart,
						Molecule* mol1, Molecule* mol2,
						int convergenceCondition,
						int parameter2,
						double (*pt2AtomKernel)(Atom*, Atom*),
						double (*pt2BondKernel)(Bond*, Bond*), int depth ){

	//cout << "MOLECULE RLK" << endl;

	vector<Atom*>::iterator atom;
	vector<Atom*>::iterator otherAtom;

	map<Atom*, Bond*>::iterator bond;
	map<Atom*, Bond*>::iterator otherBond;

	double q = 0;
	double t = 0;

	int i = 0;
	int j = 0;

	Bond* bondSecond = NULL;
	Atom* bondFirst = NULL;
	double bondSecondPT = 0;
	int bondFirstRP = 0;

	#ifdef DEBUG
		cout << "MoleculeUtils::rlk( depth = " << depth << " )" << endl;
	#endif


	//if( depth >= parameter2 ){

		//cout << "MoleculeUtils::rlk 1" << endl;

		//cout << "molecules have " << mol1->numAtoms() << " and " << mol2->numAtoms() << " atoms " << endl;
		// pour tous les atomes deux a deux
		for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
			j = 0;
			for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){
				//if( parameter2 < 2 ){
					q = (*atom)->getKashimaPQ() * (*otherAtom)->getKashimaPQ();
				//}else{
					//if( depth >= parameter2 ){
					//	q = (*atom)->getKashimaPQ() * (*otherAtom)->getKashimaPQ();
					//}else{
				//		q = 0;
					//}
				//}
				//cout << "A" << endl;
				t = 0;

				// pour toutes les liaisons des atomes deux a deux
				for( bond = (*atom)->beginBond(); bond != (*atom)->endBond(); bond++ ){
					//cout << "B" << endl;
					bondSecond = (*bond).second;
					bondFirst = (*bond).first;
					//cout << bondSecond << endl;
					//cout << "C" << endl;
					bondSecondPT = bondSecond->getKashimaPT();
					//cout << "D" << endl;
					bondFirstRP = bondFirst->getRPosition();
					//cout << "E" << endl;


					for( otherBond = (*otherAtom)->beginBond(); otherBond != (*otherAtom)->endBond(); otherBond++ ){


						t = t + bondSecondPT *
						(*otherBond).second->getKashimaPT() *
						pt2BondKernel( bondSecond, (*otherBond).second ) *
						pt2AtomKernel( bondFirst, (*otherBond).first ) *
						(*r)[ bondFirstRP ][ (*otherBond).first->getRPosition() ];
					}
				}
				(*rwork)[i][j] = q + t;
				j++;
			}
			i++;
		}

		if( parameter2 > 1 && depth == parameter2-1 ){
			//cout << "SUBSTR RLK at depth " << depth << endl;
			i = 0;
			for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
				j = 0;
				for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){

					(*rstart)[i][j] = (*rwork)[i][j];
					j++;
				}
				i++;
			}
		}


		//cout << "MoleculeUtils::rlk 2" << endl;


		// if r and rwork did not converge yet, make another recursion
		//if( depth >= parameter2 ){
			if( MoleculeUtils::converge( r, rwork, mol1, mol2, convergenceCondition )  || depth >= 50 ) {
			}else{
			//cout << "MoleculeUtils::rlk 3" << endl;
			  //if( t > 0 ){
				MoleculeUtils::rlk( rwork, r, rstart, mol1, mol2, convergenceCondition, parameter2, pt2AtomKernel, pt2BondKernel, depth + 1 );
			  //}else{
			//	cout << "i stopped at depth " << depth << endl;
			  //}
			}//else{
			//	cout << "stopped at depth " << depth << endl;
			//}
		//}else{
		//		MoleculeUtils::rlk( rwork, r, mol1, mol2, convergenceCondition, parameter2, pt2AtomKernel, pt2BondKernel, depth + 1 );
		//}
	//}else{
	//	MoleculeUtils::rlk( rwork, r, mol1, mol2, convergenceCondition, parameter2, pt2AtomKernel, pt2BondKernel, depth + 1 );
	//}

	if( parameter2 > 1 && depth == 1){
		i = 0;
		for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
			j = 0;
			for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){
				//cout <<	(*rwork)[i][j] << " - " << (*rstart)[i][j] << endl;

				(*rwork)[i][j] = (*rwork)[i][j] - (*rstart)[i][j];
				j++;
			}
			i++;
		}
	}

	return( rwork );
}

inline
bool MoleculeUtils::converge(
				vector< vector<double> >* r1,
				vector< vector<double> >* r2,
				Molecule* mol1,
				Molecule* mol2,
				int convergenceCondition ) {

	vector<Atom*>::iterator atom;
	vector<Atom*>::iterator otherAtom;
	int i = 0;
	int j = 0;
	for( atom = mol1->beginAtom(); atom != mol1->endAtom(); atom++ ){
		j = 0;
		for( otherAtom = mol2->beginAtom(); otherAtom != mol2->endAtom(); otherAtom++ ){
			if( fabs ( (*r2)[i][j] - (*r1)[i][j] ) > (*r1)[i][j] / convergenceCondition ){
				#ifdef DEBUG
					cout << "  MoleculeUtils::converge: " << fabs ( (*r2)[i][j] - (*r1)[i][j] ) << " >? " << (*r1)[i][j] / convergenceCondition << endl;
				#endif
				return( false );
			}
			j++;
		}
		i++;
	}
	return( true );
}

double MoleculeUtils::atomKernelSymbol( Atom* a1, Atom* a2 ){
	if( a1->getAN() == a2->getAN() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}

double MoleculeUtils::atomKernelMorganLabel( Atom* a1, Atom* a2 ){
	if( a1->getMorganLabel() == a2->getMorganLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}

double MoleculeUtils::atomKernelPerretLabel( Atom* a1, Atom* a2 ){
	//#ifdef DEBUG
	//	cout << "MoleculeUtils::atomKernelPerretLabel" << endl;
	//#endif
	if( a1->getPerretLabel() == a2->getPerretLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}

double MoleculeUtils::atomKernelPerretLabelExternalMatrix( Atom* a1, Atom* a2 ){
	string l1 = a1->getPerretLabel();
	string l2 = a2->getPerretLabel();
	if( l1 == l2 ) {
		return 1.0;
	}else{
		if( l1 == "CJ" || l1 == "CK" || l2 == "CJ" || l2 == "CK" ){
			return 0.0;
		}else{
			return( elements.gramAtom[ a1->getType() ][ a2->getType() ] );
		}
	}
}


double MoleculeUtils::atomKernelExternalMatrix( Atom* a1, Atom* a2 ){
	return( elements.gramAtom[ a1->getType() ][ a2->getType() ] );
}

/*!
    \fn MoleculeUtils::atomKernelLabel( Atom* a1, Atom*, a2 )
    returns 1 if two atoms have the same label, 0 othwerwise
 */
double MoleculeUtils::atomKernelLabel( Atom* a1, Atom* a2 ){
	if( a1->getLabel() == a2->getLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}

double MoleculeUtils::bondKernelType( Bond* b1, Bond* b2 ){
	if( b1->getLabel() == b2->getLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}

double MoleculeUtils::bondKernelRotable( Bond* b1, Bond* b2 ){
	if( b1->getLabel() == SINGLEBOND || b2->getLabel() == SINGLEBOND ){
		if( b1->getLabel() == b2->getLabel() ) {
			return 1.0;
		}else{
			return 0.0;
		}
	}else{
		return 1.0;
	}
}

double MoleculeUtils::bondKernelPerretLabelStrict( Bond* b1, Bond* b2 ){
	if( b1->getPerretLabel() == b2->getPerretLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}
}


double MoleculeUtils::bondKernelPerretLabel( Bond* b1, Bond* b2 ){

	//#ifdef DEBUG
	//	cout << "MoleculeUtils::bondKernelPerretLabel(...)" << endl;
	//#endif

	bool b1IsCycle = b1->hasRing();
	bool b2IsCycle = b2->hasRing();

	int b1Label = b1->getPerretLabel();
	int b2Label = b2->getPerretLabel();

	//cout << "rrh comparing bond " << b1->toStringShort() << " with " << b2->toStringShort() << endl;

	double result = -1;

	if( b1IsCycle != b2IsCycle ){
		// one bond is part of a cycle, the other not, return 0
		//return( 0 );
		result = 0;
	}else{
		if( b1IsCycle == true ){
			// both bonds are member of a cycle
			switch( b1Label ){
				case AROMATICBOND:
					switch( b2Label ){
						case AROMATICBOND:
							result = 1.0;
							break;
						case SINGLECYCLEBOND:
							result = 0.99;
							break;
						case DOUBLECYCLEBOND:
							result = 0.25;
							break;
						case TRIPLECYCLEBOND:
							result = 0.125;
							break;
					}
					break;
				case SINGLECYCLEBOND:
					switch( b2Label ){
						case AROMATICBOND:
							result = 0.99;
							break;
						case SINGLECYCLEBOND:
							result = 1.0;
							break;
						case DOUBLECYCLEBOND:
							result = 0.5;
							break;
						case TRIPLECYCLEBOND:
							result = 0.25;
							break;
					}
					break;
				case DOUBLECYCLEBOND:
					switch( b2Label ){
						case AROMATICBOND:
							result = 0.25;
							break;
						case SINGLECYCLEBOND:
							result = 0.5;
							break;
						case DOUBLECYCLEBOND:
							result = 1.0;
							break;
						case TRIPLECYCLEBOND:
							result = 0.5;
							break;
					}
					break;
				case TRIPLECYCLEBOND:
					switch( b2Label ){
						case AROMATICBOND:
							result = 0.125;
							break;
						case SINGLECYCLEBOND:
							result = 0.25;
							break;
						case DOUBLECYCLEBOND:
							result = 0.5;
							break;
						case TRIPLECYCLEBOND:
							result = 1.0;
							break;
					}
					break;
			}

		}else{
			// both bonds are aliphatic
			switch( b1Label ){
				case SINGLEBOND:
					switch( b2Label ){
						case SINGLEBOND:
							result = 1.0;
							break;
						case DOUBLEBOND:
							result = 0.5;
							break;
						case TRIPLEBOND:
							result = 0.25;
							break;
					}
					break;
				case DOUBLEBOND:
					switch( b2Label ){
						case SINGLEBOND:
							result = 0.5;
							break;
						case DOUBLEBOND:
							result = 1.0;
							break;
						case TRIPLEBOND:
							result = 0.5;
							break;
					}
					break;
				case TRIPLEBOND:
					switch( b2Label ){
						case SINGLEBOND:
							result = 0.25;
							break;
						case DOUBLEBOND:
							result = 0.5;
							break;
						case TRIPLEBOND:
							result = 1.0;
							break;
					}
					break;
			}
		}
	}

	//cout << "   = " << result << endl;

	/*if( result == -1 ){
		cout << "rrh comparing bond " << b1->toStringShort() << " with " << b2->toStringShort() << endl;
		cout << b1IsCycle << " " << b2IsCycle << endl;
		exit(1);
	}*/

	return( result );

	/*if( b1->getPerretLabel() == b2->getPerretLabel() ) {
		return 1.0;
	}else{
		return 0.0;
	}*/
}




void MoleculeUtils::writeDOTGraph( Molecule& aMolecule, ofstream& outFile, bool perretLabels ){

	//cout << "MoleculeUtils::writeDOTGraph" << endl;

	//cout << "digraph " << aMolecule.getName() << " {" << endl;
	//cout << "  size = \"4,4\";" << endl;


	outFile << "digraph \"" << aMolecule.getName() << "\" {" << endl;
	outFile << "  size = \"4,4\";" << endl;
	Atom* source;
	Bond* bond;
	Atom* target;

	// write atoms
	vector<Atom*>::iterator ai;
	map<Atom*, Bond*>::iterator bi;
	for( ai = aMolecule.beginAtom(); ai != aMolecule.endAtom(); ai++ ){
		source = (*ai);
		outFile << "  " << source->getElementSymbol() << "_" << source->getId();
		//outFile << " [ label = \"" << source->getElementSymbol();
		outFile << " [ label = \"" << source->getMorganLabel();

		if( perretLabels == true ){
			outFile << "/" << source->getPerretLabel();
		}

		//try{
			//outFile << "/" << source->getKashimaPS();
			//outFile << "/" << source->getMorganLabel();
		//} catch( CError e){
		//	continue;
		//}
		//try{
		//	outFile << "/" << source->getKashimaPQ();
		//} catch( CError e){
		//	continue;
		//}

		outFile << "\"]" << endl;

	}
	// write bonds
	for( ai = aMolecule.beginAtom(); ai != aMolecule.endAtom(); ai++ ){
		source = (*ai);
		for( bi = source->beginBond(); bi != source->endBond(); bi++ ){
			target = (*bi).first;
			bond = (*bi).second;
			outFile << "  " << source->getElementSymbol() << "_" << source->getId();
			outFile << "->" << target->getElementSymbol() << "_" << target->getId();
			outFile << " [label = \"" << bond->getLabel();
			//try{
			//	outFile << "/" << bond->getKashimaPT();
			//} catch( CError e){
			//	continue;
			//}
			if( perretLabels == true ){
				outFile << "/" << bond->getPerretLabel();
			}
			outFile << "\"]" << endl;
		}
	}

	outFile << "}" << endl;
}



/** returns the graph kernel value between two molecules using the Ecole des Mine
			approach, based on the fused graph sum of matrix power.
			The kernel computes the sum of te probabilities of finite length
			paths, as opposed to infinite length path as in Kashima's original paper
*/
double MoleculeUtils::powerKernelUntilN(
				Molecule* mol1, Molecule* mol2,
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)(Bond*, Bond*),
				int maxPower, int minLength ) throw( CError ){

	//cout << "MoleculeUtils::powerKernelUntilN " << minLength << " -> " << maxPower << endl;

	if( maxPower < 0 ){
		stringstream out;
		out << "MoleculeUtils::powerKernelUntilN: bad number of itterations: " << maxPower << " should be integer > 0 " << endl;
		CError e(BADVALUE, out.str() );
		e.describe();
		throw(e);
	}

	// make the fused graph
	Molecule* fused = new Molecule( *mol1, *mol2, pt2AtomKernel, pt2BondKernel );

	// compute the sum of the probabilities for the first maxPower powers of the transition matrix
        double kernel = 0;
	if( minLength < 2 ){
		kernel = fused->sumPQPSFast();
		//cout << "START PROB" << endl;
	}

	if( maxPower > 0 ){
		int i = 0;
		//cout << "rrr " << endl;
		//cout << "---" << endl;
		i++;
		while( i < maxPower ){
			//cout << "power " << i << endl;
			if( i >= minLength ){
				//cout << "USING LENGTH " << i << endl;
				kernel += fused->sumProbabilitiesFast();
			}
			fused->raisePowerFast();
			i++;
		}
		kernel += fused->sumProbabilitiesFast();
		//cout << "USING FINAL LENGTH " << i << endl;

	}

	delete fused;

	return( kernel );
}


double MoleculeUtils::powerKernelOrderN(
				Molecule* mol1, Molecule* mol2,
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)(Bond*, Bond*),
				int anOrder, int parameter2 ) throw( CError ) {

	//cout << "MoleculeUtils::powerKernelOrderN " << anOrder << endl;

	if( anOrder < 0 ){
		stringstream out;
		out << "MoleculeUtils::powerKernelOrderN: bad order: " << anOrder << " should be integer > 0 " << endl;
		CError e(BADVALUE, out.str() );
		e.describe();
		throw(e);
	}

	double kernel = 0;

	// make the fused graph
	Molecule* fused = new Molecule( *mol1, *mol2, pt2AtomKernel, pt2BondKernel );

	if( anOrder == 0 ){
		kernel = fused->sumPQPSFast();
	}else{

		// compute the transition probabilities for paths of order anOrder
		int i = 0;
		//fused->saveAllBonds();
		//kernel = fused->sumProbabilities();
		//kernel = fused->sumPQPS();

		i++;
		while( i < anOrder ){
			//cout << "rrr" << endl;
			fused->raisePowerFast();
			//cout << "---" << endl;
			i++;
		}
		kernel += fused->sumProbabilitiesFast();
	}
	delete fused;
	return( kernel );
}






double MoleculeUtils::powerKernelConverge(
				Molecule* mol1, Molecule* mol2,
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)(Bond*, Bond*),
				int converge, int minLength ) throw( CError ){

	//cout << "MoleculeUtils::powerKernelConverge " << converge << endl; //<< mol1->toStringShort() << " / " << mol2->toStringShort() << endl;

	if( converge <= 0 ){
		stringstream out;
		out << "MoleculeUtils::powerKernelUntilN: bad number of itterations: " << converge << " should be > 0 " << endl;
		CError e(VALUENOTALLOWED, out.str() );
		e.describe();
		throw(e);
	}

	double okernel = 0;
	double condition = 1.0/converge;

	// make the fused graph
	Molecule* fused = new Molecule( *mol1, *mol2, pt2AtomKernel, pt2BondKernel );
	//Molecule* fused = new Molecule( *mol1 );

	//fused->describeLong();

	// compute the sum of the probabilities for the first maxPower powers of the transition matrix
	//fused->saveAllBonds();

        double kernel = 0;
	if( minLength < 2 ){
		kernel = fused->sumPQPSFast();
	}


	int i = 0;
	//cout << "rrr " << endl;
	//cout << "---" << endl;
	i++;
	while( kernel - okernel > condition || i < minLength+1 ){
		//cout << "power " << i << endl;
		okernel = kernel;
		if( i >= minLength ){
			//cout << "USING LENGTH " << i << endl;
			kernel += fused->sumProbabilitiesFast();
		}
		fused->raisePowerFast();
		i++;
	}
	//cout << "stopped at " << i << " " << minLength+1 << endl;
	kernel += fused->sumProbabilitiesFast();
	//cout << "USING FINAL LENGTH " << i << endl;

	delete fused;
	return( kernel );
}



/** No descriptions */
void MoleculeUtils::describeMap( map<Atom*, float>* aMap ){
	map<Atom*, float>::iterator i;
	for( i = aMap->begin(); i != aMap->end(); i++ ){
		cout << (*i).first->toStringShort() << " " << (*i).second << endl;
	}
}


/** helper function for detectSSSR. returns in result a vector containing atom pointer
	for atoms present in full but not in exclude
*/
void MoleculeUtils::substractSet( vector<Atom*>* full, vector<Atom*>* exclude, vector<Atom*>* result){

	result->clear();

	vector<Atom*>::iterator ai;
	for ( ai = full->begin(); ai != full->end(); ai++ ){
		vector<Atom*>::iterator aj;
		bool found = false;
		for ( aj = exclude->begin(); aj != exclude->end(); aj++ ){
			if( (*aj) == (*ai) ){
				#ifdef DEBUG
				cout << "a10h   excluding " << (*ai)->toStringShort() << endl;
				#endif
				found = true;
				break;
			}
		}
		if( found == false ){
			result->push_back( (*ai) );
		}
	}
}

void MoleculeUtils::mergeSet( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result){
	vector<Atom*>::iterator ai;
	for ( ai = v1->begin(); ai != v1->end(); ai++ ){
		result->push_back( (*ai) );
	}

	for ( ai = v2->begin(); ai != v2->end(); ai++ ){
		vector<Atom*>::iterator aj;
		bool found = false;
		for ( aj = result->begin(); aj != result->end(); aj++ ){
			if( (*aj) == (*ai) ){
				found = true;
				break;
			}
		}
		if( found == false ){
			result->push_back( (*ai) );
		}

	}
}

void MoleculeUtils::selectRingMemberBonds( vector<Bond*>* v1, vector<Atom*>* a, vector<Bond*>* v2 ){
	vector<Bond*>::iterator bi;
	for ( bi = v1->begin(); bi != v1->end(); bi++ ){

		#ifdef DEBUG
			cout << "ffh verifying bond " << (*bi)->toStringShort() << endl;
		#endif

		Atom* source = (*bi)->getSource();
		Atom* target = (*bi)->getTarget();

		bool foundSource = false;
		bool foundTarget = false;

		vector<Atom*>::iterator ai;
		for( ai = a->begin(); ai != a->end(); ai++ ){
			#ifdef DEBUG
			cout << "ffh   with " << (*ai)->toStringShort();
			#endif
			if( (*ai) == source ){
				//cout << " found source ";
				foundSource = true;
			}
			if( (*ai) == target ){
				//cout << " found target ";
				foundTarget = true;
			}
			if( foundSource == true && foundTarget == true ){
				v2->push_back( (*bi) );
				#ifdef DEBUG
					cout << " member " << endl;
				#endif
				break;
			}
			//cout << endl;
		}
		#ifdef DEBUG
			if( foundSource == false || foundTarget == false ){
				cout << " non member" << endl;
			}
		#endif

	}
}

void MoleculeUtils::mergeBondSet( vector<Bond*>* v1, vector<Bond*>* v2, vector<Bond*>* result){
	#ifdef DEBUG
	cout << "MoleculeUtils::mergeBondSet" << endl;
	#endif
	vector<Bond*>::iterator ai;
	for ( ai = v1->begin(); ai != v1->end(); ai++ ){

		#ifdef DEBUG
		cout << (*ai)->toStringShort() << endl;
		#endif


		vector<Bond*>::iterator aj;
		bool found = false;
		for ( aj = result->begin(); aj != result->end(); aj++ ){
			if( (*aj)->getSource() == (*ai)->getSource() && (*aj)->getTarget() == (*ai)->getTarget() ){
				found = true;
				break;
			}
		}
		if( found == false ){
			result->push_back( (*ai) );
			#ifdef DEBUG
			cout << " added 1" << endl;
			#endif
		}
		#ifdef DEBUG
			else{

			cout << " skip" << endl;
			}
		#endif
	}

	for ( ai = v2->begin(); ai != v2->end(); ai++ ){
		vector<Bond*>::iterator aj;
		bool found = false;

		#ifdef DEBUG
		cout << (*ai)->toStringShort() << endl;
		#endif

		for ( aj = result->begin(); aj != result->end(); aj++ ){
			if( (*aj)->getSource() == (*ai)->getSource() && (*aj)->getTarget() == (*ai)->getTarget() ){
				found = true;
				break;
			}
		}
		if( found == false ){
			result->push_back( (*ai) );
			#ifdef DEBUG
			cout << " added 2" << endl;
			#endif
		}
		#ifdef DEBUG
			else{

			cout << " skip" << endl;
			}
		#endif
	}
}

bool MoleculeUtils::atomVectorHas( vector<Atom*>* atomVector, Atom* anAtom ){
	vector<Atom*>::iterator ai;
	for( ai = atomVector->begin(); ai != atomVector->end(); ai++ ){
		if( (*ai) == anAtom ){
			return( true );
		}
	}
	return( false );
}










// THREE-D KERNELS //
//-----------------//

double MoleculeUtils::threeDkernel(
				Molecule* mol1, Molecule* mol2,
				double (*pt2AtomKernel)( Atom*, Atom*),
				double (*pt2BondKernel)(float, float, float),
				float edgeKernelParameter ) {

  double kernel = 0;

  // make the fused graph
  Molecule* fused = new Molecule( *mol1, *mol2, pt2AtomKernel, pt2BondKernel, edgeKernelParameter);


#ifdef DEBUG
  cout << "fused->numAtoms() : " << fused->numAtoms() << endl;
  cout << "ADJACENCY MATRIX";
  for(int i = 0 ; i < fused->numAtoms() ; i++){
    cout << endl;
    for(int j = 0 ; j < fused->numAtoms() ; j++){
      cout << fused->getAdjacency(i,j) << "\t";
    }
  }
  cout << endl << endl;
#endif


  // 1st possibility : compute Ax^3 and get its trace
  //      --> NOT OPTIMUM : we don't need all the cubic matrix
  //fused->raisePowerAdjacency();
         // --> fused->raisePowerAdjacency_sparse()
  //fused->raisePowerAdjacency();
         // --> fused->raisePowerAdjacency_sparse()
  //kernel = fused->traceWalks();
         // --> fused->traceWalks_sparse()


  // 2nd possibility : compute Ax^2 and only the diagonal of Ax^3
  fused->raisePowerAdjacency();
  kernel = fused->traceDiagWalks();



#ifdef DEBUG
  cout << "CUBE MATRIX";
  for(int i = 0 ; i < fused->numAtoms() ; i++){
    cout << endl;
    for(int j = 0 ; j < fused->numAtoms() ; j++){
      cout << fused->getWalks(i,j) << "\t";
    }
  }
  cout << endl << endl << endl;
#endif
  


  fused->eraseAdjacency();
  fused->eraseWalks();
  delete fused;


  return(kernel);
}



double MoleculeUtils::threeDedgeKernelRBF( float dist1, float dist2 , float param ){
  // RBF kernel based on the difference between the edges distances : exp( -(d1-d2)^2 / 2*sigma^2)
  // param = sigma (RBF bandwidth)
  double res;
  res = exp( -0.5 * (dist1-dist2) * (dist1-dist2) * (1/(param*param)) ); 

  return res;
}


double MoleculeUtils::threeDedgeKernelTriangle( float dist1, float dist2 , float param ){
  // triangular kernel : k(d1,d2) = max( 0 , (C - |d1-d2|)/C )
  // Rq : 0 <= k(d1,d2) <= 1
  // param = C = maximum difference tolerated

  double tmp;

  tmp = (param - fabs(dist1 - dist2)) / param;

  if (tmp > 0.0)
    return tmp;
  else 
    return 0.0;

}





