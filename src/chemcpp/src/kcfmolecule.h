/****************************************************************************************
					  kcfmolecule.h 
					-----------------
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



#ifndef KCFMOLECULE_H
#define KCFMOLECULE_H

#include <molecule.h>

/**this is a subclass of molecules providing special functions for molecules represented in KCF format
  *@author Jean-Luc Perret
  */

class KCFMolecule : public Molecule  {
public:

	/** class constructor.
	*/
	KCFMolecule();
	
	/** class desctructor.
	*/
	virtual ~KCFMolecule();

	/** reads the molecule description from a KCF file. The internal representation of the molecule is the Atom KEGG atoms, so the function Elements:loadDefinition should have previously loaded the KEGG atom definition types).
It is not possible to mix KCF and MOL files.
A function readKCFconverted should be written to read KCF file with an internal representation of elements.
WARNING this function only reads the first entry in the kcf file.
	*/
	void readKCF( string aFileName ) throw( CError );
	
	/** writes a kcf file.
	*/
	void writeKCF( string aFileName );

	/** adds an atom to the kcf molecule.
	*/
	virtual Atom* addAtom(string aSymbol) throw( CError );

	

//virtual float newKernel( Molecule* anotherMolecule, int convergenceCondition);
  	//virtual vector< vector<float> >* rlk( vector< vector<float> >* r, vector< vector<float> >* rwork, Molecule* anotherMolecule, int convergenceCondition );

};

#endif
