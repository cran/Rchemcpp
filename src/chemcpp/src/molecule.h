/****************************************************************************************
					  molecule.h 
					--------------
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



#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <map>
#include <sstream>
#include <locale>
#include <string>
#include <algorithm>


#include <datacontainer.h>
#include <atom.h>
#include <elements.h>
#include <constant.h>

//#include <kashimathreadstruct.h>

extern Elements elements;


/** Molecule class which can be compared using Graph Kernels

    @author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
		@version 0.3
		@date 17 Jan 2004

		CLASS NAME: 	Molecule

		FOR:					SNSF SPONSORED PROJECT

		PURPOSE:  		This class implements the notion of molecule. Molecules have

								- a vector of atoms
								- properties

		Properties are implemented by deriving the Molecule from the DataContainer class
		which	takes care of	memory allocation for these descriptors. Therefore descriptors
		can be added and removed at runtime using the DataContainer functions.

		The Molecule class can be seen as a graph containing a list of vertex (Atom), each
		vertex containing a list of edges (Bond) to the neighbour vertex (Atom). Edges
		are directed and have a label. Chemical bonds are not directed so bonds are always
		double, one in each direction.
		So a Molecule can be seen as a directed labeled graph with bidirectional edges.
		... many other changes

*/
class Molecule : public DataContainer  {

/**
		\example molecule_example.cpp
*/

public:

///@name Molecule construction functions
//@{

	/** class constructor.
	*/
	Molecule();

	/** class destructor. deletes all atoms in the molecule.
	*/
	virtual ~Molecule();

	/** copy constructor.
	*/
	Molecule( Molecule& aMolecule , bool bool_resetMorganIndex = true );

	/** product graph constructor. this constructor constructs a fusion molecule from two molecules using a product graph operation.
			two atoms which are connected in both molecules have an edge
			connecting them in the fused graph.
	*/
	Molecule(
				Molecule& molecule1,
				Molecule& molecule2,
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)( Bond*, Bond* )
				);

	/** "3D" constructor.
		Note : BondKernel prototype is different of those used on the "product graph constructor".
	*/    	
	Molecule(
			Molecule& m1,
			Molecule& m2,
			double (*pt2AtomKernel)( Atom*, Atom*),
			double (*pt2BondKernel)( float, float, float ),
			float edgeKernelParameter
			);


	/** assignment operator.
	*/
	Molecule& operator=( const Molecule& aMolecule );

	/** adds an Atom with aSymbol chemical symbol to the molecule.
		aSymbol can be in any case. The first letter of aSymbol will be automatically changed
		to capital and the second letter to normal caps.
		Returns a pointer to the created Atom.
		WARNING this atom will be deleted when the molecule is deleted so don't use this pointer!
		if resetSSSR = true, then this action invalidates previous sssr computations.
	*/
	virtual Atom* addAtom( string aSymbol, bool resetSSSR = true ) throw( CError );

	/** adds an Atom to the atoms vector.
		Returns a pointer to the created Atom.
		WARNING this atom will be deleted when the molecule is deleted so don't use this pointer!
		if resetSSSR = true, then this action invalidates previous sssr computations.
	*/
	Atom* addAtom( Atom*, bool resetSSSR = true, bool resetMorganIndex = true ) throw( CError);


	/** links two atoms poited by aSource and aTarget using aBondLabel bond type.
		If the two atoms are not in the atoms vector they are first added to the molecule.
		returns a pointer to the forward bond just created (the return bond can be retrived with forwardBond->getReverse()).
	*/
	Bond* linkAtoms( Atom* aSource, Atom* aTarget, int aBondLabel, int aBondStereo = 0,  int aBondNotUsed = 0, int aBondTopology = 0, int aBondReactionCenter = 0, bool resetSSSR = true );

	/** links two atoms using their name.
		If no atom with name aSource or aTarget are found a CError error is thrown.
		returns a pointer to the forward bond just created (the return bond can be retrived with forwardBond->getReverse()).
	*/
	Bond* linkAtoms( string aSource, string aTarget, int aBondLabel, int aBondStereo = 0,  int aBondNotUsed = 0, int aBondTopology = 0, int aBondReactionCenter = 0, bool resetSSSR = true ) throw( CError );

	/** links two atoms using their position in the atom container.
		if firstAtom or secondAtom are not between 0 and numAtoms()-1 a CError error is thrown.
		returns a pointer to the forward bond just created (the return bond can be retrived with forwardBond->getReverse()).
	*/
	Bond* linkAtoms( int firstAtom, int secondAtom, int aBondLabel, int aBondStereo = 0,  int aBondNotUsed = 0, int aBondTopology = 0, int aBondReactionCenter = 0, bool resetSSSR = true ) throw( CError );


	/** reads the molecule description from a Mol file.
		MUST STILL BE COMPLETED

		\todo atomlist block is not taken into account
		\todo stext block is not taken into account
		\todo properties block is not taken into account
		\todo chirality is not taken into account (chiral flag is set by reading mol file)
		\todo bondTopology is not taken into account
	*/
	void readMOL( string aFilename, bool genericAtomType = false ) throw(CError);

	/** OLD version of readMOL().
	*/
	void readMOLOld( string aFilename ) throw(CError);

	/** erases all atoms in the current molecule.
		returns the molecule in the state of construction by the default constructor. 
	*/
	void erase();
	
	/** erases hidden atoms of the molecule.
	*/
	void eraseHiddenAtoms();

	/** erases the rings of the molecule.
	*/
	void eraseRings();

	/** erases the adjacency matrix of the molecule.
	*/
	void eraseAdjacency();
	
	/** erases the 'walks' matrix of the molecule (used to compute the pharmacophore kernel, see (Mahe, 2006)).
	*/
	void eraseWalks();

//@}



///@name Accessor functions
//@{
	/** returns the unique id of the molecule.
	*/
	int getId() const { return( id ); }

	/** returns the unique id of the molecule as a string.
	*/
	string getIdString();

	/** returns the molecule name, i.e., the name descriptor value (Nb :not necessarily unique).
	*/
	string getName(){ return( getStringDescriptor( "name" )->getValue() ); }

	/** sets the molecule name.
	*/
	void setName( string aName ) { setStringDescriptor( "name", aName, "", "", true, true ); }

	/** sets the activity of the molecule (a descriptor in float native format for faster retrieval).
	*/
	void setActivity( float aNumber ){ activity = aNumber; flagActivity = true; }

	/** sets the activity of the molecule from a string activity. provided for convenience.
	*/
	void setActivity( string aNumber ){ activity = atof( aNumber.c_str() ); flagActivity = true; };

	/** returns the activity of the molecule.
	*/
	float getActivity( bool silentMode = false ) throw( CError );
	
	/** returns true if the activity of the molecule was initialized.
	*/
	float hasActivity(){ return( flagActivity ); }

	/** turns off the activity flag.
	*/
	void unsetActivity(){ flagActivity = false; }


	/** returns a pointer to the atom's vector.
	*/
	vector<Atom*>& getAtoms(){ return( atoms ); }

	/** returns a pointer to an atom designated by its unique Id.
		raises a CError exception if no atom with this Id is present in the molecule.
	*/
	Atom* getAtom(int anId) throw( CError );

	Atom* getAtomByIndex(int ind);

	/** checks if an atom is part of the molecule. 
		returns true if anAtom is included in the atom vector, false otherwise.
	*/
	bool atomExists(Atom* anAtom);

	/** returns the number of atoms of the molecule.
	*/
	int numAtoms(){ return( atoms.size() ); }

	/** returns the number of hidden atoms of the molecule.
	*/
	int numHiddenAtoms(){ return( hiddenAtoms.size() ); }


	/** returns the number of bonds of the molecule.
	*/
	int numBonds();
	/** returns the number of hidden bonds of the molecule.
	*/
	int numHiddenBonds();

	/** returns the sum of the bond types.
	*/
	long bondSum();

	/** unsets all bonds flags.
	*/
	void unsetBondFlags();
	/** unsets all bonds flags original.
	*/
	void unsetBondFlagsOriginal();

	/** returns true if the molecule has chiral information.
	*/
	bool isChiral(){ return( chiral ); }

	/** sets selectedFlag.
	*/
	void select(){ selectedFlag = true; }
	/** unset selectedFlag.
	*/
	void unSelect(){ selectedFlag = false; }
	/** returns selectedFlag.
	*/
	bool isSelected(){ return( selectedFlag ); }

	/** sets the descriptor type and name used by the < operator to compare molecules.
		if reverse == true then the sorting will be in descending order.
	*/
	void setSortDescriptor( string aName, int aType );


	/** returns the type of the descriptor used to sort molecules.
	*/
	int getSortDescriptorType(){ return( sortDescriptorType ); }

	/** returns the name of the descriptor used to sort molecules.
	*/
	string getSortDescriptorName(){ return( sortDescriptorName ); }

	/** returns the ring of the corresponding Id.
	*/
	Ring* getRingWithID( int anID, bool createIfMissing) throw( CError );

	/** returns true if the molecule already has newRing in sssr.
		If sssr was not detected before this call, detectSSSR() is called.
	*/
	bool hasRing( Ring* newRing, bool detectingRing = false ) throw( CError );

	/** returns true if the molecule has a ring.
		If sssr was not detected before this call, detectSSSR() is called.
	*/
	bool hasRing() throw( CError );

	/** returns the number of rings present in the molecule.
	*/
	int numRings() throw( CError ){
		if( !hasSSSRDetected() ){
			//detectSSSR();
			CError e( SSSRNOTDETECTED, "Smallest Set of Smallest Rings was not detected before calling Molecule::numRings()" );
			//e.describe();
			throw(e);
		}
		return( sssr.size() );
	}


	/** returns the value of the 'hasSSRDetected' flag.
	*/
	bool hasSSSRDetected(){ return( flagHasSSSRDetected ); }

	/** returns the molecular weight.
	*/
	float getMW( bool silentError = false ) throw( CError );

	/** returns the number of atoms in the molecule (ignoring hidden atoms).
	*/
	int getNumAtoms(){
		return( atoms.size() );
	}

        /** returns the location where the molecule is stored.
	*/
	string getLocation(){ return( location ); }

	/** returns the format of the original file.
		one of MDL2000, MDL3000, MDL2000JLP.
	*/
	int getOriginalFormat(){ return( originalFormat ); }
	
	/** sets the molecule's original format.
		one of MDL2000, MDL3000, MDL2000JLP.
	*/
	void setOriginalFormat( int a ){ originalFormat = a; }

//@}


	



///@name Molecule manipulation functions
//@{

	/** hides all atoms (and associated edges) with intDescriptor aDescriptorName == aValue.
			returns the number of hidden atoms 
	*/
	int hideAtomsByIntDescriptor( string aDescriptorName, int aValue, bool refreshBonds = true );

	/** hides an atom (but not the bonds to this atom).
		call refreshBonds() after hiding all atoms.
	*/
	void hideAtom( vector< Atom* >::iterator anAtomI );

	/** hides anAtom and all to / from bonds for this atom.
	*/
	void hideAtomAndToFromBonds( vector< Atom* >::iterator anAtomI );

	/** hides anAtom and all to / from bonds for this atom.
	*/
	void hideAtomAndToFromBonds( Atom* anAtom );

	/** moves the hydrogen atoms to the hiddenAtoms container, so that further algorithms will run without hydrogens.
		returns the number of hydrogens hidden.
		hidden hydrogens can be restored using restoreHiddenAtoms().
	*/
	int hideHydrogens();

	/** returns true if an atom is hidden.
	*/
	bool isHiddenAtom( Atom* anAtom );

	/** restores the hidden atoms (for example hydrogens).
	*/
	int restoreHiddenAtoms( bool flagRestoreBonds = true );

	/** restores the hidden bonds among atoms (does not restore bonds in hidden atoms).
	*/
	int restoreHiddenBonds();

	/** hides all bonds to and from hidden atoms.
		call this function after hiding all atoms in a molecule.
		returns the number of edges hidden.
	*/
	int refreshBonds();

	/** erases anAtom from the list. Does not delete the atom though!
	*/
	void eraseAtom( Atom* anAtom ) throw( CError );

	/** deletes all bonds in all atoms.
	*/
	void deleteBonds();

	/** deletes all hidden atoms.
	*/
	void deleteHiddenAtoms();

	/** sets the morganLabel of atoms to the concatenation of the atomic symbol
		and the Morgan index of anOrder iteration.
	*/
	void setMorganLabels( int anOrder );

	/** sets the perretLabel of atoms: the atomic symbol. 
		Only carbons get an additional J if they have more than two aromatic bonds (aromatic rings joining carbons).
	*/
	void setPerretLabels();

	/** sets the uniqueMorganIndex of each atom to the Morgan index having the maximum of
		different connectivity values for the molecule.
		returns the iteration of the smallest Morgan Index with the
		maximum number of connectivity values.
	*/
	int setUniqueMorganIndices();

	/** resets the calculations of the morgan indices.
		CALL THIS FUNCTION EACH TIME THE MOLECULE IS MODIFIED! 	*/
	void resetMorganIndex();

	/** returns the number of distinct morgan indices (of order anOrder)  present in the molecule.
	*/
	int getNumberOfDistinctMorganIndices( int anOrder );

	/** returns the iteration of the Morgan index at which the different connectivity values
		have reached a maximum.
	*/
	int getMaxMorganIteration();


	/** returns the number of carbons in a connected component graph defined with an Int descriptor named 
		aDescriptorName having value aValue.
	*/
	int getNumCarbonsOfComponent( string aDescriptorName, int aValue );

	/** returns the number of nitrogens in a connected component graph defined with an Int descriptor named
		aDescriptorName having value aValue.
	*/
	int getNumNitrogensOfComponent( string aDescriptorName, int aValue );

	/** returns the number of terminal atoms of the molecule.
	*/
	int numAtomsNonCSkeleton();	

	/** computes the Euclidian distance between two atoms.
	*/
	float atomicDistance(Atom* atom1, Atom* atom2);

	/** this function should be called everytime a molecule is created / modified.
		automatically called at the end of Molecule::readMol() and
		MoleculeUtils::readMDLCtabBlock().
	*/
	void compute();

	/** sets the activity according to the value of a descriptor and a value.
		If true is given as third argument (default) then descriptorName <= value are considered positive.
		If false is given as third argument then molecules with descriptorName >= value are considered positive.
		if a descriptor is missing then the molecule is left without activity.
	*/
	void binClassifyFromDescriptor( string descriptorName, float value, bool smallerOrEqual = true );

	
	/** fills in a vector of String with the distinct labels of the atoms of the molecule.
	*/
	void atomsLabelsListing( vector<string>* );

	/** fills in vector of String with the distinct symbols of the atoms of the molecule.
	*/
	void atomsSymbolsListing( vector<string>* );

	/** fills in vector of Int with the distinct types of bonds of the molecule.
	*/
	void bondsListing( vector<int>* );

	/** transforms the molecule into a graph preventing 'totters' (see (Mahe 2004)).
	*/
	void noTottersTransform();

	/** transforms the molecule into a 3D graph (see (Mahe 2006)).
	*/ 
	void threeDtransform(int nBins, double distMin, double distMax);

	/** reads partial charges of the atoms of the molecule from a ';' separated string.
	*/
	void readPartialCharges(string charges);
	
	/** introduces the sign of the partial charge into the label of the atoms of the molecule.
	*/
	void setMorganChargesLabels(double threshold);

	/** hides all but the biggest (in atom numbers) connected atoms in the molecule.
		returns the number of components removed.
		WARNING: does not test if two compounds have identical maximum size.
			in that case, it will hide all but the first one.
	*/
	int hideSalts( stringstream* out );

	/** NOT DOCUMENTED
	*/
	int markFragments();

	/** NOT DOCUMENTED
	*/
	void unmarkFragments();

	/** NOT DOCUMENTED
	*/
	void hideAllFragmentsBut( int aFragmentNumber );

	/** NOT DOCUMENTED
	*/
	void writeFragments( ofstream* outStream );

	/** performs a DFS search on the graph starting at startAtom,
		marking atoms with markValue in descriptor intDescriptorName,
		and returning the number of nodes in the connected component.
	*/
	int DFS( Atom* startAtom, string intDescriptorName, int markValue );

	/** detects the smallest set of smallest rings and returns the number of rings.
		this procedure is automatically called when a molecule is loaded from a mol or sd file through the function compute. 
		call compute() explicitly after modifying atom or bond composition.
		this function also sets the ring membership of atoms and bonds.
		according to: Ring Perception Using Breadth-First Search, John Figueras 399 Baker s Pond Road, Orleans, Massachusetts 02653, J. Chem. Inf. Comput. Sci. 1996, 36, 986-991.
	*/
	int detectSSSR();

	/** sets the flagHasSSSRDetected to false.
	*/
	void resetSSSR(){ flagHasSSSRDetected = false; }

	/** sets the flagHasSSSRDetected to true.
	*/
	void setHasSSSR(){ flagHasSSSRDetected = true; }

	/** helper function for detectSSSR.
		selects an optimum edge for elimination in structures without N2 nodes.
	*/
	Bond* checkEdges( Atom* );


//@}





///@name Iterators
//@{

	/** returns an iterator to the first Atom.
	*/
	vector<Atom*>::iterator beginAtom(){ return( atoms.begin() ); }
	/** returns an iterator to the last Atom.
	*/
	vector<Atom*>::iterator endAtom(){ return( atoms.end() ); }

	/** returns an iterator to the first Bond of an Atom.
	*/
	map<Atom*, Bond*>::iterator beginBond( Atom* anAtom ){ return( anAtom->beginBond() ); }
	/** returns an iterator to the first Bond of an Atom.
	*/
	map<Atom*, Bond*>::iterator endBond( Atom* anAtom ){ return( anAtom->endBond() ); }

	/** returns an iterator to the first Bond of an Atom designated by its id.
	*/
	map<Atom*, Bond*>::iterator beginBond( int anId ){ return( getAtom( anId )->beginBond() ); }
	/** returns an iterator to the first Bond of an Atom designated by its id.
	*/
	map<Atom*, Bond*>::iterator endBond( int anId ){ return( getAtom( anId )->endBond() ); }

	/** NOT DOCUMENTED
	*/
	map<int, int>::iterator beginComponentSizes(){
		return( componentSizes.begin() );
	}

	/** NOT DOCUMENTED
	*/
	map<int, int>::iterator endComponentSizes(){
		return( componentSizes.end() );
	}

//@}



///@name Marginalized and general kernel related functions
//@{

	/** sets the start, stop and transition probabilities for the calculation
		of the Kashima kernel as published in (Kashima et al. 2003).
		the start probability is set to an identical value of 1/numAtoms() for all atoms.
		the stop probability can be choosen between 0 and 1.
		the transition probability of each bond of an Atom is set to
		(1.0 - stop probability) / numbonds.

	*/
	void setKashimaKernelProb( double aPq, bool skipSkeleton = false );

	/** returns the sum of transition probabilities. should be used with fused graph.
	*/
	double sumPT();

	/** function used to compute the graph kernel using the fused graph + 
		sum the powers of the transition matrix approach.
		see MoleculeUtils::powerKernel().
	*/
	double sumProbabilities();
	
	/** fast version of sumProbabilities().
	*/
	double sumProbabilitiesFast();

	/** function used to compute the graph kernel using the fused graph + 
		sum of the powers of the transition matrix approach.
		see MoleculeUtils::powerKernel()
	*/
	double sumPQPS();

	/** fast version of sumPQPS().
	*/
	double sumPQPSFast();


	/** function used to compute the graph kernel using the fused graph + sum of
		the powers of the transition matrix approach. 
		see MoleculeUtils::powerKernel().
		same as raisePower() but instead of using the Bond class the transition
		probabilities are read from fastPt and fastPTSave.
	*/
	void raisePowerFast();

	/** compares this molecule with anotherMolecule and returns the value defined by
		the function pt2GraphKernel, using the function pt2AtomKernel to compare atoms
		and pt2BondKernel to compare bonds.
	*/
	double computeKernel(
				Molecule* anotherMolecule,
				double (*pt2GraphKernel)(
					Molecule* mol1,
					Molecule* mol2,
					double (*pt2AtomKernel) ( Atom*, Atom* ),
					double (*pt2BondKernel)( Bond*, Bond* ),
					int parameter1, int parameter2
				),
				double (*pt2AtomKernel)( Atom*, Atom* ),
				double (*pt2BondKernel)( Bond*, Bond* ),
				int parameter1, int parameter2 = 1 );


	/** returns the kernel value for the molecule with itself.
	*/
	double getSelfKernel( double(*pt2GraphKernel)( Molecule* mol1,
		Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)( Bond*, Bond* ), int, int ),
		double(*pt2AtomKernel)( Atom*, Atom* ),
		double(*pt2BondKernel)( Bond*, Bond* ) ,
		int parameter1, int parameter2 = 1
		){
			#ifdef DEBUG
				cout << "Molecule::getSelfKernel(...) " << endl;
			#endif

			if(selfKernelCalculated == false){
				#ifdef DEBUG
					cout << "  sk not computed, computing" << endl;
				#endif
				calculateSelfKernel( pt2GraphKernel, pt2AtomKernel, pt2BondKernel, parameter1, parameter2 );
			}
			#ifdef DEBUG
				cout << "  returning " << selfKernel << endl;
			#endif
			return( selfKernel );
		}


	/** calculates the molecule's self Kashima kernel and stores the result in
		the private selfKernel variable.
	*/
	double calculateSelfKernel
	(
				double(*pt2GraphKernel)
				(
						Molecule* mol1,
						Molecule* mol2,
						double(*pt2AtomKernel)(Atom*, Atom*),
						double(*pt2BondKernel)(Bond*, Bond*),
						int, int
				),
				double(*pt2AtomKernel)( Atom*, Atom* ),
				double(*pt2BondKernel)( Bond*, Bond* ),
				int paramter1, int parameter2
	);


	/** emits an error if the kernel was not calculated before.
	*/ 
	double getSelfKernel() throw( CError );

	/** sets the set kernel to the given value.
	*/
	void setSelfKernel(double value);

	/** adds a value to the self kernel.
	*/
	void addToSelfKernel( double );
	
	/** substracts a value to the self kernel.
	*/
	void substractToSelfKernel( double );

	/** resets the self kernel.
	*/
	void resetSelfKernel(){ selfKernelCalculated = false; }

//@}



///@name Matrix utilities functions related to the 3D pharmacophore kernel (Mahe 2006)
//@{
	/** sets an entry of the adjacency matrix of the molecule to a given value.
	*/
	void setAdjacency(int i, int j, double value);

	/** returns the value of an entry of the adjacency matrix of the molecule.
	*/
	double getAdjacency(int i, int j);

	/** sets an entry of the 'walks' matrix (i.e., a weighted adjacency matrix) of the molecule to a given value.
	*/
	void setWalks(int i, int j, double value);

	/** returns the value of an entry of the walks matrix of the molecule.
	*/
	double getWalks(int i, int j);

	/** multiplies the walks matrix by the adjacency matrix.
	*/
	void raisePowerAdjacency();
	
	/** computes the trace of the walks matrix.
	*/
	double traceWalks();

	/** computes the diagonal elements of the walk matrix multiplied by the adjacency matrix, and sums to get the trace.
	*/
	double traceDiagWalks();


//@}



///@name Output functions
//@{

	/** returns a string description of the molecule
		(name, unique Id, memory Adress, number of atoms).
	*/
	string toString();

	/** returns a short string description of the molecule.
		(name and unique Id).
	*/
	string toStringShort();

	/** returns a string description of the molecule
		(name, unique Id, memory Adress, number of atoms and atom descriptions).
	*/
	string toStringLong();


	/** prints a description of the molecule to stdout
		using toString() and DataContainer describe().
	*/
	void describe();

	/** prints a description of the molecule to stdout
		using toString().
	*/
	void describeShort();

	/** prints a description of the molecule to stdout
		using the molecule toString(), all atoms toStringShort(),
		and all bonds toStringShort().
	*/
	void describeLong();

	/** prints a description of each atom (using describe())
		and each bond (using describe())to stdout.
	*/
	void describeEachAtom();

	/** write the molecule definition in a MDL MOL file.
	*/
	void writeMOL( string aFileName );

	/** writes a dot file representing the molecule.
		the molecule can be converted to graphic representation using Graphviz's dot program.
	*/
	void writeDOT( string aFilename, bool perretLabels = false );

	/** write a MDL structure data (SD) file with a single molecule.
	*/
	void writeSD( string aFileName );

//@}


// ******************************* //
// **** DEPRECATED FUNCTIONS **** //
// ****************************** //


	/** save all bonds in the molecule
	*/
	//void saveAllBonds();





	/** returns the graph kernel value between two molecules by default the atom kernel is:
			MoleculeUtils::atomKernelExternalMatrix
			the values used as kernel for atoms are read in a table.
			by default this table returns the same value as the Kashima kernel.
			However the kernel used to compare atoms can be changed by loading a gram
			matrix for elements using elements.loadGramAtoms(...),
			or by specifying another atom kernel function pointer
			as the pt2AtomKernel argument
			and the bond kernel is:
			MoleculeUtils::bondKernelType
	*/
	//virtual float moleculeKernel( Molecule* anotherMolecule, float (*pt2AtomKernel)( Atom*, Atom* ), float (*pt2BondKernel)( Bond*, Bond* ), int convergenceCondition = 1000 );


	/** returns the graph kernel value between two molecules using the Ecole des Mine
			approach, based on the fused graph sum of matrix power.
			The kernel computes the sum of te probabilities of finite length
			paths, as opposed to infinite length path as in Kashima's original paper
	*/
	//virtual float powerKernelUntilN( Molecule* anotherMolecule, float (*pt2AtomKernel)( Atom*, Atom* ), int maxPower = 4 );

	/** returns the graph kernel value between two molecules using the Ecole des Mine
			approach, based on the fused graph sum of matrix power.
			The kernel computes an approximation of the sum of the probabilities of
			infinite finite length paths (until the kernel value does not change by more
			than 1/converg
	*/
	//virtual float powerKernelConverge( Molecule* anotherMolecule, float (*pt2AtomKernel)( Atom*, Atom* ), int converg = 1000 );

	/** returns the graph kernel value between two molecules using the Ecole des Mine
			approach, based on the fused graph sum of matrix power.
			The kernel computes the sum of the probabilities of paths of length N
	*/
	//virtual float powerKernelOrderN( Molecule* anotherMolecule, float (*pt2AtomKernel)( Atom*, Atom* ), int length );




	/** thread version of the Kashima kernel.
		Takes a kashimaThreadStruct argument as input
		(using the molecule1, molecule2, and convergence members to calculate the kernel).

		returns the result of the kernel calculation in the member "result" of the
		kashimaThreadStruct. If calculation occurred without errors than the noErrors
		member of the kashimaThreadStruct is set to true.

		to create a new threads with this function use:

		#define MAXTHREADS 100

	pthread_t thread;
	kashimaThreadStruct args;

	if(error=pthread_create( &thread, NULL, (void*)&kashimaKernelThread, (void*) &args) !=0 )
	{
		switch(error){
			case EAGAIN: cout<<"EAGAIN\n";break;
			case EINVAL: cout<<"EINVAL\n";break;
			case ENOMEM: cout<<"ENOMEM\n";break;
		}
		exit(1);
	}

	// ensure thread is finished before comtinuing
	pthread_join( thread, NULL );

	*/
	//  void* kashimaKernelThread(void* arg);




	/** function used to compute the graph kernel using the fused graph / sum of
	the powers of the transition matrix approach. See MoleculeUtils::powerKernel
	*/
	//void raisePower();





	/** this function outputs fragments of the molecule into a new SDFile
		the fragments are generated by removing one bound. Aromatic
		bounds are not broken.
	*/
	//void exportFragments( string sdFileName, int minAtoms = 1 );



protected:

	/** this function is called internally every time a molecule is changed (atom or bonding), reset flags.
	*/
	void moleculeChanged( bool resetSSSR = true, bool resetMorganIndex = true );


	/** container for atoms.
	*/
	vector<Atom*> atoms;

	/** container for hidden atoms (for example the removed hydrogens are moved
	to this container).
	*/
	vector<Atom*> hiddenAtoms;


	/** container for rings (smallest set of smallest rings).
	*/
	vector<Ring*> sssr;

	/** flag indicating if smallest set of smallest rings was detected.
	*/
	bool flagHasSSSRDetected;

	/** counter of the number of molecules which have been instanciated.
	*/
	static int counter;

	/** molecule unique Id.
	*/
	int id;

	/** flag to select subsets of molecules.
	*/
	bool selectedFlag;

	/** molecule's self Kernel.
	*/
	double selfKernel;

	/** true if selfKernel was calculated, false otherwise.
	*/
	bool selfKernelCalculated;

	/** true if the molecule has chiral information.
	*/
	bool chiral;

	/** true if the molecule has activity information.
	 */
	bool flagActivity;

	/** type and name of descriptor to be used to sort molecules.
		use setSortDescriptor( aType, aName )
			getSortDescriptorType(),
			setSortDescriptorName() to access.
	*/
	int sortDescriptorType;

	/** name of the sorting descriptor.
	*/
	string sortDescriptorName;

	/** contains the iteration number at which the maximum
			number of distinct connectivity values is reached.
	*/
	int maxMorganIteration;

	

private: // Private methods

  	/** hide all bonds in the molecule
			warning clears all previously hidden bonds
  	*/
  	//void hideAllBonds();

	/** links two atoms poited by aSource and aTarget using aBondLabel bond type.
			If the two atoms are not in the atoms vector they are first
			added to the molecule
			WARNING this function adds a single directed bond, only to be used with in
			Molecule::raisePower().
	*/
	Bond* linkAtomsNoReturn( Atom* aSource, Atom* aTarget, int aBondLabel, int aBondStereo = 0,  int aBondNotUsed = 0, int aBondTopology = 0, int aBondReactionCenter = 0 );

	/** links a pair of atoms in a single direction only : firstAtom -> secondAtom.
	*/
	void linkAtomsNoReturn( int firstAtom, int secondAtom, int aBondLabel ) throw( CError );


	/** hash used in the fast computations of the powers of the product graph adjacency matrix.
	*/
	map<Atom*, map<Atom*, double>* >* fastPT;
	/** hash used in the fast computations of the powers of the product graph adjacency matrix.
	*/
	map<Atom*, map<Atom*, double>* >* fastPTNext;
	/** hash used in the fast computations of the powers of the product graph adjacency matrix.
	*/
	map<Atom*, map<Atom*, double>* >* fastPTSave;	
	/** hash used in the fast computations of the powers of the product graph adjacency matrix.
	*/
	map<Atom*, double> fastPQ;
	/** hash used in the fast computations of the powers of the product graph adjacency matrix.
	*/
	map<Atom*, double> fastPS;


	/** NO DOCUMENTATION
	*/
	map< int, int > componentSizes;

	/** activity value of the molecule.
	*/
	float activity;

	/** location where the molecule is stored. automatically set by readMOL.
	*/
	string location; 	
	
	/** original format of the molecule.
	*/
	int originalFormat;

	/** adjacency matrix of the molecule. used to compute the 3D pharmacophore kernel.
	*/
	vector< vector<double> >* adjacency;
	
	/** 'walks' matrix of the molecule. used to compute the 3D pharmacophore kernel.
	*/	
	vector< vector<double> >* walks;

};

#endif
