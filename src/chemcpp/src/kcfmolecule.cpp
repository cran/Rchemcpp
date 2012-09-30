/****************************************************************************************
					  kcfmolecule.cpp 
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



#include "kcfmolecule.h"
#include "moleculeutils.h"



extern Elements KEGGelements; 


KCFMolecule::KCFMolecule(){
	//cout << "KCFMolecule constructor" << endl;
}
KCFMolecule::~KCFMolecule(){
}

/** reads the molecule description from a KCF file. The internal representation of the molecule is the Atom KEGG atoms, so the function Elements:loadDefinition should have previously loaded the KEGG atom definition types).
It is not possible to mix KCF and MOL files.
A function readKCFconverted should be written to read KCF file with an internal representation of elements. */
void KCFMolecule::readKCF( string aFileName ) throw( CError ){

  // erase all eventually existing atoms in the molecule
	erase();

	// open molfile
	ifstream inFile;
 	inFile.open( aFileName.c_str(), ios::in );
	if( !inFile.good() ){
		CError e = CError(FILENOTFOUND, aFileName + " file not found");
		e.describe();
		throw(e);
	}

	MoleculeUtils::readKCFMolecule( *this, inFile );

	inFile.close();

}
/** write a kcf file */
void KCFMolecule::writeKCF( string aFileName ){

	// open kcf file
	ofstream outFile;
	outFile.open( aFileName.c_str(), ios::out );
	if( !outFile.good() ){
		CError e = CError(FILENOTFOUND, aFileName + " could not be created");
		e.describe();
		throw(e);
	}

	MoleculeUtils::writeKCF( *this, outFile );

	outFile.close();

}


Atom* KCFMolecule::addAtom(string aSymbol) throw( CError ){

	/*ToUpper      up(std::locale::classic());
 	ToLower      down(std::locale::classic());

  std::transform (aSymbol.begin(), aSymbol.begin()+1, aSymbol.begin(), up);
	if( aSymbol.length() > 1 ){
		std::transform (aSymbol.begin()+1, aSymbol.end(), aSymbol.begin()+1, down);
	}   */

	string h = aSymbol.substr(0, 1);
	string r = "";
	if( aSymbol.length() > 1 ){
		r = aSymbol.substr( 1, aSymbol.length()-1 );
	}

	h = StringUtils::toUpper( h );
	r = StringUtils::toLower( r );
	aSymbol = h + r;



	Atom* anAtom = new Atom( *KEGGelements[aSymbol] );

	atoms.push_back(anAtom);

	#ifdef DEBUG
		//cout << "created atom " << anAtom->toString() << endl;
	#endif
	return( anAtom );

}
