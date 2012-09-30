
/****************************************************************************************
					  moleculeutils.h 
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


#ifndef MOLECULEUTILS_H
#define MOLECULEUTILS_H

//#include <vector>
#include <istream>
#include <locale>
#include <string>
#include <iomanip>

//#include <stringutils.h>
//#include <molecule.h>
#include <moleculeset.h>
//#include <atom.h>

#define FINDENTRY 0
#define READINGNODES 1
#define READINGEDGES 2
#define READINGUNKNOWN 3


using namespace std;


/**static functions to be used to process molecules
  *@author Jean-Luc Perret
  */

class MoleculeUtils {

public:


///@name Input/Ouput utility functions
//@{


	/** reads the ctab block of a stream (connection table),
		and adds atoms and bonds to aMolecule.
		if genericAtomTypeFlag is true then the atoms are not taken from the periodic
		table of elements, but are created from the label provided.
	*/
	static void readMDLCtabBlock(
		Molecule& aMolecule,
		ifstream& inFile,
		bool genericAtomTypeFlag = false
		) throw( CError );

	/** reads the header block of a stream, and sets molecule name to the comment content
		or to aName if specified.
	*/
	static void readMDLHeaderBlock(
		Molecule& aMolecule,
		ifstream& inFile,
		string aName = "COMMENT"
		) throw( CError );

	/** reads the non structural data block of a compound in a MDL SDfile.
  	WARNING only data entry with a structure aLabel will be considered. entries
		missing the aLabel entry will be ignored.
	*/
	static void readMDLNSDBlock( Molecule& aMolecule, ifstream& inFile ) throw( CError );

	/** skips all lines until next entry.
	*/
	static void skipMDLEntry( Molecule& aMolecule, ifstream& inFile ) throw( CError );

	/** writes the header block to a stream, for aMolecule.
	*/
	static void writeMDLHeaderBlock( Molecule& aMolecule, ofstream& outFile );

	/** writes the ctab block to a stream (connection table),
			for aMolecule.
	 */
	static void writeMDLCtabBlock( Molecule& aMolecule, ofstream& outFile );

	/** writes the non structural data block to a stream, for aMolecule.
  	this function writes all stringDescriptors, floatDescriptors and IntDescriptors
		of a molecule in appropriate MDL SD format.
	*/
	static void writeMDLNSDBlock( Molecule& aMolecule, ofstream& outFile );

	/** writes a molecule under the KCF format in the outFile stream
	*/
	static void writeKCF( Molecule& aMolecule, ofstream& outFile );

	/** writes the non structural data block to a stream, for aMolecule.
		this function writes all stringDescriptors, floatDescriptors and IntDescriptors
		of a molecule in appropriate KEGG KCF format.
	*/
	static void writeKCFNSDBlock( Molecule& aMolecule, ofstream& outFile );

	/** reads one molecule from a kcf file.
	*/
	static bool readKCFMolecule( KCFMolecule& aMolecule, ifstream& inFile ) throw( CError );


	/** writes a graph in Graphviz's dot format.
	*/
	static void writeDOTGraph( Molecule& aMolecule, ofstream& outFile, bool perretabels = false );

//@}




///@name Graph kernels related functions
//@{

	/** computes the marginalized graph kernel between two molecules.
			takes two molecules, a function pointer to the atom and the bond
			kernels and the convergence condition as parameter.
			WARNING: random walk parameters must first be set.
	*/
	static double moleculeKernel(
		Molecule* mol1,
		Molecule* mol2,
		double (*pt2AtomKernel)(Atom*, Atom*) = &atomKernelSymbol,
		double (*pt2BondKernel)(Bond*, Bond*) = &bondKernelType,
		int convergenceCondition = 1000, int parameter2 = 1
		);

	/** returns the marginalized graph kernel value between two molecules using the product graph
			approach and the sum of matrix powers.
			The kernel computes the sum of the probabilities of finite length
			paths, until order maxPower.
	*/
	static double powerKernelUntilN( Molecule* mol1, Molecule* mol2, double (*pt2AtomKernel)( Atom*, Atom* ), double (*pt2BondKernel)(Bond*, Bond*), int maxPower, int minLength = 1 ) throw( CError );

	/** returns the marginalized graph kernel value between two molecules using the product graph
			approach and the sum of matrix powers.
			The kernel computes the probabilities of finite length of length anOrder.
	*/
	static double powerKernelOrderN( Molecule* mol1, Molecule* mol2, double (*pt2AtomKernel)( Atom*, Atom* ), double (*pt2BondKernel)(Bond*, Bond*), int anOrder, int parameter2 = 1 ) throw( CError );

	/** returns the marginalized graph kernel value between two molecules using the product graph
			approach and the sum of matrix powers.
			The kernel computes the sum of the probabilities of infinite length
			paths, until the kernel value does not change by more than 1/converge.
	*/
	static double powerKernelConverge( Molecule* mol1, Molecule* mol2, double (*pt2AtomKernel)( Atom*, Atom* ), double (*pt2BondKernel)(Bond*, Bond*), int converge, int parameter2 = 1 ) throw( CError );

	/** computes the "rlk" term involved in the marginalized kernel computation.
	*/
	static vector< vector<double> >* rlk(
		vector< vector<double> >* r,
		vector< vector<double> >* rwork,
		vector< vector<double> >* rstart,
		Molecule* mol1,
		Molecule* mol2,
		int convergenceCondition,
		int parameter2,
		double (*pt2AtomKernel)(Atom*, Atom*),
		double (*pt2BondKernel)(Bond*, Bond*),
		int depth
		);
	
	/** returns true if the marginalized kernel computation has converged according to a convergence condition.
	*/
	static bool converge(
		vector< vector<double> >* r1,
		vector< vector<double> >* r2,
		Molecule* mol1,
		Molecule* mol2,
		int convergenceCondition
		);

	/** atom kernel based on atomic number.
	*/
	static double atomKernelSymbol( Atom* a1, Atom* a2 );

	/** atom kernel based on Morgan label.
			WARNING: labels must first be set.
	*/
	static double atomKernelMorganLabel( Atom* a1, Atom* a2 );

	/** atom kernel based on Perret label.
	*/
	static double atomKernelPerretLabel( Atom* a1, Atom* a2 );

	/** atom kernel based on Perret label and external gram matrix.
	*/
	static double atomKernelPerretLabelExternalMatrix( Atom* a1, Atom* a2 );


	/** atom kernel based on an external gram matrix.
			WARNING: the gram matrix must first be loaded using Elements::loadGramAtoms().
	*/
	static double atomKernelExternalMatrix( Atom* a1, Atom* a2 );

	/** generic atom kernel label.
	*/
	static double atomKernelLabel( Atom* a1, Atom* a2 );


	/** bond kernel.
		returns 1 if two bonds have the same type, 0 otherwise.
	*/
	static double bondKernelType( Bond* b1, Bond* b2 );

	/** bond kernel.
		returns 1 if two bonds have the same perretLabel, 0 otherwise.
	*/
	static double bondKernelPerretLabelStrict( Bond* b1, Bond* b2 );

	/** bond kernel.
		returns 1 if two bonds have the same perretLabel.
		returns 0 if one bond is a cycle while the other is not.
		returns a weight is both bonds are cycles or non cycles.
	*/
	static double bondKernelPerretLabel( Bond* b1, Bond* b2 );


	/** bond kernel.
		Single bonds are rotable, others are not. Returns 1 if two bonds.
		are rotable or two bonds are non rotable, 0 if one of the bonds.
		is rotable and the other is not.
	*/
	static double bondKernelRotable( Bond* b1, Bond* b2 );

//@}




///@name  3D kernels related functions
//@{

	/** (3D) generic pharmacophore kernel between two molecules. defined for a choice of atom an edge kernels, and corresponding parameters.
	*/
	static double threeDkernel(
				   Molecule* mol1, Molecule* mol2,
				   double (*pt2AtomKernel)( Atom*, Atom*),
				   double (*pt2BondKernel)(float, float, float),
				   float edgeKernelParameter );

	/** (3D) RBF edge kernel. param = bandwith.
	*/
	static double threeDedgeKernelRBF( float dist1, float dist2 , float param );

	/** (3D) triangular edge kernel. param = cut-off.
	*/
	static double threeDedgeKernelTriangle( float dist1, float dist2 , float param );

//@}





///@name Miscallenous utility functions
//@{

	/** returns a string with the ring membership of aBond.
	*/
	static string getRingString( Bond* aBond ) throw( CError );

	/** NO DOCUMENTATION
	*/
	static void describeMap( map<Atom*, float>* aMap );

	/** helper function for Molecule::detectSSSR. returns in result a vector containing atom pointer
	for atoms present in full but not in exclude.
	*/
        static void substractSet( vector<Atom*>* full, vector<Atom*>* exclude, vector<Atom*>* result);

	/** returns true if anAtom is present into the atomVector vector.
	*/
	static bool atomVectorHas( vector<Atom*>* atomVector, Atom* anAtom );

	/** helper function for Molecule::detextSSSR. the result vector is the
		concatenation of v1 and v2 without duplicates.
	*/
	static void mergeSet( vector<Atom*>* v1, vector<Atom*>* v2, vector<Atom*>* result);
	
	/** merges two bonds sets.
	*/
	static void mergeBondSet( vector<Bond*>* v1, vector<Bond*>* v2, vector<Bond*>* result);
	
	/** NO DOCUMENTATION
	*/
	static void selectRingMemberBonds( vector<Bond*>* v1, vector<Atom*>* a, vector<Bond*>* v2 );


//@}


	//MoleculeUtils();
	//~MoleculeUtils();

};

#endif
