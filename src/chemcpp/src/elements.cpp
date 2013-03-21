/****************************************************************************************
					  elements.cpp 
					---------------
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


#include "elements.h"
#include "constant.h"



//#define ELEMENTS_DEBUG 1


// instanciate Elements (Atoms)
Elements elements( ELEMENTSFILENAME, GRAMATOMSIDENTITYFILENAME );
Elements KEGGelements( KEGGATOMS, GRAMKEGGATOMSIDENTITYFILENAME );



// load gram matrices for comparing atoms
void Elements::loadGramAtoms( string aFileName ){

	//periodicTable.clear();

	#ifdef ELEMENTS_DEBUG
		cout << "loading Atom Kernel matrix from file: " << aFileName << endl;
	#endif

	ifstream inFile;
	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError( FILENOTFOUND, aFileName + " file not found" );
		e.describe();
		throw(e);
	}

	int linesize=2048;
	int lineNumber = 0;
	string delimiter = ";";


	char *line = new char[linesize];

	vector<string> value;

	// read and instanciate each new element

	string stringLine;
	while( !inFile.eof() ){
		lineNumber++;
		inFile.getline( line, linesize-1,'\n' );
		stringLine = line;

		// if line is empty or starts with // or # then skip this line
 		if( stringLine.size() == 0 ){
			continue;
		}
		if( stringLine.substr( 0, 1 ) == "#" || stringLine.substr(0,2) == "//" ){
			continue;
		}

		StringUtils::Split( line, delimiter, value );

		if( value.size() < numElements() ){
			stringstream out;
			out << aFileName << " line " << lineNumber << ": has " << value.size() << " values while " << numElements() << " are required ";
			CError e = CError( BADFILE, out.str() );
			e.describe();

			cout << "found:" << endl;

			cout << line << endl;

			for( uint c = 0; c < value.size(); c++ ){
				cout << c+1 << " " << value[c] << endl;
			}

			throw(e);
		}

		// PROCESS LINE
    		for( uint i = 0; i < numElements(); i++ ){
			gramAtom[lineNumber - 1][i] = atof( value[i].c_str() );
		}

		value.clear();
	}

	delete[] line;

	gramAtomName = aFileName;

	#ifdef ELEMENTS_DEBUG
		cout << "loading Atom Kernel matrix from file: " << aFileName << " DONE " << endl;
	#endif

}


/** loads the elements from a file */
Elements::Elements( string aFileName, string aGramAtomFileName ){

	//cout << "loading element from file " << aFileName << "... " << endl;

	loadDefinition( aFileName );
	loadGramAtoms( aGramAtomFileName );

}

Elements::~Elements(){
	map<string, Atom*>::iterator i;
	for( i = periodicTable.begin(); i != periodicTable.end(); i++ ){
  		delete (*i).second; // delete bonds
	}
	periodicTable.clear();
}

/** returns an atom pointer to an element designed using its chemical symbol */
Atom* Elements::operator[]( string aSymbol ){
	// return( *periodicTable[aSymbol] );
	return( getElement( aSymbol ) );
}

/** returns a pointer to the element designed by aSymbol */
Atom* Elements::getElement( string aSymbol ){
	map<string, Atom*>::iterator anElement = periodicTable.find( aSymbol );

	if( anElement == periodicTable.end() ) {
		// element does not exist, emit exception
		stringstream out;
		out << "Error in Elements::getElement: element " << aSymbol << " not found in the Elements Set";
		CError e( NOTFOUND, out.str() );
		e.describe();
		throw( e );
	}

	return( (*anElement).second );
}




/** returns the filename of the atom kernel matrix loaded. returns "default" if no kernel matrix was loaded. */
string Elements::getAtomKernelName(){
	return gramAtomName;
}



/** loads the definition of elements
 */
void Elements::loadDefinition( string aFileName ){

	#ifdef ELEMENTS_DEBUG
		cout << "Elements::loadDefinition loading element from file " << aFileName << "... " << endl;
	#endif

	// empty periodic table
	periodicTable.clear();
	Atom::resetCounter();

	Atom* atomP;
	gramAtomName = "default";

	ifstream inFile;
 	inFile.open( aFileName.c_str(), ios::in );
	if(!inFile.good()){
		CError e = CError( FILENOTFOUND, aFileName + " file not found" );
		e.describe();
		throw(e);
	}

	int linesize=512;
	int lineNumber = 0;
	string delimiter( ";" );

  	char *line = new char[linesize];

	vector<string> header;
	vector<string> fullName;
	vector<string> unit;
	vector<string> type;
	vector<string> value;


	Descriptor<string>* newDs;
	Descriptor<int>* newDi;
	Descriptor<float>* newDf;



	//cout << "found " << SYMBOLNAME << " at position " << symbolPosition << endl;

	// read fullName line and split fields
	lineNumber++;
	inFile.getline( line, linesize-1, '\n' );
	//cout << "read fullname line, found: " << line << endl;
	StringUtils::Split( line, delimiter, fullName );   // LINE 1: full name

	// read units line and split fields
	lineNumber++;
	inFile.getline( line, linesize-1, '\n' );
	//cout << "read units line, found: " << line << endl;
	StringUtils::Split( line, delimiter, unit );   		 // LINE 2: unit

	// read descriptor types line and split fields
	lineNumber++;
	inFile.getline( line, linesize-1, '\n' );
	//cout << "read types line, found: " << line << endl;
	StringUtils::Split( line, delimiter, type );			 // LINE 3: type


		// read header line and split fields
	lineNumber++;
	inFile.getline( line, linesize-1, '\n' );
	StringUtils::Split( line, delimiter, header );		// LINE 4: header

	// find symbol position
	int symbolPosition = 0;
	bool found = false;
	uint i;
	for( i = 0; i < header.size(); i++ ){
		// cout << "C:"<<header[i] << endl;
		if( header[i] == SYMBOLNAME ){
			found = true;
			symbolPosition = i;
		}
	}
	if ( found == false ){
		// ERROR no symbol defined
		CError e = CError(BADFILE, "No Symbol field found in element file " + aFileName );
		e.describe();
		throw(e);
	}


	// read and instanciate each new element
	int typeC = 0;
	#ifdef DEBUG
		cout << "CREATING ELEMENTS: ";
	#endif
	string stringLine;
	while( !inFile.eof() ){
		lineNumber++;
		inFile.getline( line, linesize-1, '\n' );
		stringLine = line;

		// if line is empty or starts with // or # then skip this line
 		if(stringLine.size() == 0){
			continue;
		}
		if(stringLine.substr( 0, 1 ) == "#" || stringLine.substr( 0, 2 ) == "//"){
			continue;
		}

		StringUtils::Split( line, delimiter, value );
		
		if (header.size() != value.size())
		{
				stringstream out2;
				out2 << aFileName << ": header/values mismatch";
				CError e( UNKNOWNDATATYPE, out2.str() );
				e.describe();
				throw( e );
		}

    		atomP = new Atom;
		periodicTable[ value[symbolPosition] ] = atomP;
		atomP->setType( typeC );
		typeC++;
		#ifdef DEBUG
			cout << value[symbolPosition] << ", ";
		#endif

		// add all descriptors to the new Element
		for( i = 0; i < header.size(); i++ ){
			if( type[i] == "String" ) {
				if( value[i] == "Nil" || value[i] == "NA" ){
					newDs = atomP->addKindStringDescriptor( header[i], "", unit[i], fullName[i] );
					newDs->setEmpty();
				}else{
					newDs = atomP->addKindStringDescriptor( header[i], value[i], unit[i], fullName[i] );
				}
			}else if( type[i] == "Integer" ){
				if(header[i] == "AN" || header[i] == "An" ){
					atomP->setAN( atoi(value[i].c_str()) );
				}else{
        				if( value[i] == "Nil" || value[i] == "NA" ){
						newDi = atomP->addKindIntDescriptor( header[i], 0, unit[i], fullName[i] );
						newDi->setEmpty();
					}else{
						newDi = atomP->addKindIntDescriptor( header[i], atoi(value[i].c_str()), unit[i], fullName[i] );
					}
				}
			}else if( type[i] == "Float" ){
				if( value[i] == "Nil" || value[i] == "NA" ){
					newDf = atomP->addKindFloatDescriptor( header[i], 0.0, unit[i], fullName[i] );
					newDf->setEmpty();
				}else{
					atomP->addKindFloatDescriptor( header[i], atof( value[i].c_str() ), unit[i], fullName[i] );
				}
			}else{
				// unknown data type
				stringstream out;
				out << "type " << type[i] << " in line " << lineNumber-1 << " of file " << aFileName << " is not allowed. String, Float and Integer only are allowed ";
				CError e( UNKNOWNDATATYPE, out.str() );
				e.describe();
				throw( e );
			}
		}
		atomP->setLabel( atomP->getSymbol() );
		//cout << "atom " << atomP->getSymbol() << " now has label " << atomP->getLabel() << endl;

		//atomP->describe();
		value.clear();

	}



	#ifdef DEBUG
	  cout << endl;
	#endif

	delete[] line;
	inFile.close();

	// by default define an identity matrix to compare atoms

	uint j = 0;
	for( i = 0 ; i < numElements() ; i++ ){
		for( j = 0 ; j < numElements() ; j++ ){
			if(i == j){
				gramAtom[i][j] = 1;
			}else{
				gramAtom[i][j] = 0;
			}
		}
	}
	#ifdef ELEMENTS_DEBUG
		cout << "Elements::loadDefinition done" << endl;
	#endif
}


/** write a description of the elements to cout */
//void Elements::describe(){
//	map<string, Atom*>::iterator i;

	//for(i = periodicTable.begin(); i!=periodicTable.end(); i++){
		//cout << i->second << endl;
//	}
//}    //



