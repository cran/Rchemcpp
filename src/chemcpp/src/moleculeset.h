/****************************************************************************************
					  moleculeset.h 
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


#ifndef MOLECULESET_H
#define MOLECULESET_H

#include <vector>
#include <fstream>
#include <sstream>
#include <pthread.h>
#include <math.h>

#include <datacontainer.h>
#include <molecule.h>
//#include <moleculeutils.h>
#include <kcfmolecule.h>
#include <jlpioutils.h>
#include <constant.h>

#define VERBOSECALC 1


/**Set of molecules on which virtual experiments may be performed

	@author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
	@version 0.3
	@date 17 Jan 2004

	CLASS NAME: 	MoleculeSet
	FOR:					SNSF SPONSORED PROJECT
	PURPOSE:  		This class implements the notion of molecule set

	It is a set in the mathematical way in the sense that no two
	molecule should have the same name. However a MoleculeSet can
	contain two identical graphs.

	WARNING: No checks are made for throwing error in case there
	are two molecules with the same name, except when calling the
	[] operator.



*/
class MoleculeSet : public std::vector<Molecule*> {

/**
		\example moleculeset_example.cpp
*/

public:

///@name MoleculeSet construction functions
//@{

	/** class constructor.
	*/
	MoleculeSet();
	MoleculeSet(const MoleculeSet& aSet);

	/** destructor for the MoleculeSet.
		deletes gram and gramNormal.
	*/
	~MoleculeSet();

	/** adds all molecules in aSet to the current set.
	*/
	int add( MoleculeSet* aSet );

 	/** adds a molecule to the set.
		molecules added to the molecule set are not deleted when the set is deleted.
		use the deleteAllMolecule() function to do so
		(CAUTION if you reference these molecules from elsewhere).
	*/
	void addMolecule( Molecule* aMolecule );

	/** adds a copy of the molecule in argument. used by add() to merge two datasets.
	*/
	Molecule* addMoleculeCopy( Molecule* aMolecule );

	/** deletes all molecules in dataset (memory deallocation) and clears the vector containing the pointers.
	*/
	void deleteAll();


//@}


///@name Input functions
//@{

	/** adds an sd file content and returns the number of created molecules.
		if genericAtomType is true, then the atoms are not read from the
		periodic table, but are created based only on the label provided.

		beginMolecule and endMolecule specify the index of the first and last
		molecule to include in the training set (starting count from 0).
		values of -1 (default) means no limit.
	*/
	int addSD( string aFileName, bool genericAtomType = false, long beginMolecule = -1, long endMolecule = -1 );

	/** reads an KCF file and returns the number of created molecules.
	*/
	int addKCF( string aFileName, long beginMolecule = -1, long endMolecule = -1  );

	/** creates a new molecule in the dataset and reads its definition from a MDL MOL file.
	*/
	Molecule* addSingleMOL( string aMolFile, bool genericAtomType = false );

	/** creates a new molecule in the dataset and reads its definition from a Kcf file.
	*/
	Molecule* addSingleKCF( string aMolFile );

	/** loads all .mol files in a directory and adds the molecules to the set.
	*/
	void readMolDirectory( string aPath, bool genericAtomType = false, long beginMolecule = -1, long endMolecule = -1 );

	/** loads all .kcf files in a directory and adds the molecules to the set.
	*/
	void readKcfDirectory( string dataDir, long beginMolecule = -1, long endMolecule = -1 );

	/** reads the Mutag dataset.
		Atoms and bonds are read from
		aFilename while the biological activity of the molecule is read
		from file rFilename. numMolToRead allows to specify the number of
		first molecules to read.
	*/
	void addMutag( string aFileName, string rFileName="", uint numMolToRead = 500 );


	/** reads the activity of molecules from an activity file.
		expects a file containing : 
		label	[tab] class
		molname [tab] 1     for active molecules
		molname [tab] -1    for inactive molecules
		molname is the fileName of the Mol file.
	*/
	void readActivityFile( string aFileName );

	/** reads a descriptor file.
	the first line in the file indicates the Descriptor comment text.
	the second line in the file indicates the Descriptor units.
	the third line in the file indicates the Descriptor types.
	the forth line indicates the Descriptor names.
	data start at fifth line.

	the first column should contain the same molecule names as those obtained
	by getName() for the molecules in the dataset.

	example:

	Full name;Log octanol partition coefficient;boiling point;special comment
	NA;NA;K;NA
	string;float;float;string
	name;logP;Bp;acomment
	ethane;10;200;junk data


	WARNING: does not check for duplicate name entries in the descriptor file.
	WARNING: does not add anything to the molecules not contained in the descriptor file.
	*/
	void readDescriptorFile( string aFileName, string separator = ";" );

	/** reads the content of a file produced by gist-classify and adds the information
	matching molecule names in the set as descriptors gistname, gistclass, and gistdiscr.
	*/
	void readGistClassifyFile( string aFileName );

	/** reads the content of a gist activity file and adds the information to
	descriptor aDescriptor for each molecule matching molecule names. 
	*/
	void readGistActivityFile( string aFileName, string aDescriptor );

	/** reads a gram matrix matching the dataset.
		emits an error if the dimension of the read gram matrix read does not
		match the number of compounds in the dataset.
	*/
	void readGram( string aFileName, vector< vector<double> >* gram );

	/** reads a normalized gram matrix matching the dataset.
	*/
	void readGramNormal( string aFileName );
	
	/** reads a raw gram matrix matching the dataset.
	*/	
	void readGramRaw( string aFileName );

	/** reads the partial charges associated to the molecule set from an input file.
	NOTE : the input file has one line per molecule of the molecule set, and within each line, 
	the values of the partial charges are separated by ';'.
	*/
	void readPartialCharges(string fileName);


//@}





///@name Accessor functions
//@{

	/** returns the number of molecules in the MoleculeSet.
	*/
	uint numMolecules(){ return( size() ); }

	/** returns a pointer to the first molecule with name aName in the set.
		Throws a CError exception if no molecule with that name exists in the set.
		WARNING: if more than one molecule have the same name than a CError is also
		thrown.
	*/
	Molecule* operator[]( string aName ) throw( CError );

	/** returns a pointer to the first molecule with name aName in the set.
		Throws a CError exception if no molecule with that name exists in the set.
		WARNING: if more than one molecule have the same name than a CError is also
		thrown.
	*/
	Molecule* getMolByName( string aName ) throw( CError );

	/** returns a pointer to the molecule of index anInd in the set.
	*/
	Molecule* operator[]( int anInd ) throw( CError );

	/** returns a pointer to the molecule of index anInd in the set.
	*/
	Molecule* getMolByIndex( int anInd ) throw( CError );


	/** fills a vector of int with the values an int descriptor can take among all molecules and
		returns the number of such values.
	*/
	long getPossibleValuesInIntDescriptor( string aDescriptorName, vector< int >* p );

	/** returns the stop probability currently in use in the MoleculeSet
			for the calculation of random walk graph Kernel.
	*/
	double getPq(){ return( pq ); }


	/** returns the convergence condition currently in use in the MoleculeSet
			for the calculation of random walk graph Kernel.
	*/
	int getConvergenceCondition(){ return( convergenceCondition ); }


	/** returns true if a molecule with name aName exists in the set, false otherwise.
	*/
	bool nameExists( string aName );

	/** returns true if the activity of molecules in the set was set using readActivityFile.
	WARNING: no check is made to verify if all molecules were set.
	*/
	bool hasActivity(){ return activitySet; }

	/** sets integer aName to aValue for all compounds in the dataset.
	*/
	void setIntDescriptor( string aName, int aValue );

	/** sets the uniqueMorganIndex of each atom to the Morgan index having
		the maximum of different connectivity values for the molecule.
	*/
	void setUniqueMorganIndices();

	/** sets the morganLabels of each molecule to the anOrder iteration of
		the Morgan index calculation process.
	*/
	void setMorganLabels( int anOrder );

	/** sets the value of intDescriptor aLabel of molecule aLabel to aValue.
	*/
	void setIntDescriptor( string aLabel, string aMolecule, int aValue );

	/** sets the avtivity of aMolecule to aValue.
	*/
	void setActivity( string aMolecule, float aValue );
  
	/** sets the start, stop (aPq) and transition probabilities.
		sets the start, stop (aPq) and transition probabilities for all molecules
		according to the article by Kashima et al. using the setKashimaKernelProb(aFloat)
		function of the Molecule class.
		WARNING: a call to this function erases the gram Matrix and it will
		be recalculated.
	*/
	void setKashimaKernelParam( double aPq, int aConvergenceCondition, bool skipSkeleton = false );

	/** sets the comparison set of the moleculeSet.
	*/
	void setComparisonSet( MoleculeSet* );
	
	/** sets the 'morgan charges' labels of the atoms, i.e., the concatenation of the Morgan labels of the atoms 
	and the (+/-) sign of their partial charges.
	*/
	void setMorganChargesLabels(double threshold);


//@}



///@name MoleculeSet manipulation functions
//@{

 	/** selects all molecules in the dataset.
 		WARNING the selection status is stored in the Molecule class. Therefore
		molecules in other datasets will be selected too...
	*/
 	void selectAll();

 	/** unselects all molecules in the dataset.
		WARNING the selection status is stored in the Molecule class. Therefore
		molecules in other datasets will be unselected too...
	*/
	void unSelectAll();

	/** selects the molecules with the names provided as arguments in a vector of string.
		unselect all others. returns the number of selected molecules.
	*/
	int select( vector< string >* aSubset );

	/** unselects the molecules with the names provided as arguments in a vector of string.
		returns the number of selected molecules.
	*/
	int unSelect( vector< string >* aSubset );

	/** selects all molecules having float descriptor aName with value aValue, unselects all others.
	*/
	long selectByFloatDescriptor( string aName, float aValue );

	/** selects all molecules having int descriptor aName with value aValue, unselects all others.
	*/
	long selectByIntDescriptor( string aName, int aValue );

	/** selects all molecules having activity equal to aValue.
	*/
	long selectByActivity( float aValue );

	/** selects all molecules which have a defined activity status.
	*/
	long selectHasActivity();

	/** selects all molecules in the dataset with mw >= minmw and <= maxmw.
		if maxmw = -1 then there is no maximum limit.
		if addMolecularDescriptor is true then a floatDescriptor with label mw is added
		to the molecule.
	*/
	int selectByMW( float minmw, float maxmw = -1 , bool addMolecularDescriptor = false );

	/** selects all molecules in the dataset with number of atoms >= numAtoms and <= numAtoms.
		if maxNumAtoms = -1 then there is no maximum limit.
		if addMolecularDescriptor is true then a intDescriptor with label numAtoms is added
		to the molecule.
	*/
	int selectByNumAtoms( float minNumAtoms, float maxNumAtoms = -1 , bool addMolecularDescriptor = false );

	/** sorts the molecule collection according to the molecule
		Descriptor descriptorName of type descriptorType.
		if reverse = true then the sorting is in decreasing order.
	 */
	void sortByDescriptor( string aDescriptorName, int aDescriptorType, bool reverse = false );

	/** sorts the molecule collection according to the molecule Descriptor descriptorName.
		descriptor type will be read from descriptor name if name is of type
		******.integer, or *******.float, set to string otherwise.
		if reverse = true then the sorting is in decreasing order. 
	*/
	void sortByDescriptor( string aDescriptorName, bool reverse = false );

	/** sorts all compounds in the set by Molecular weight.
	*/
	void sortByMW();

	/** sorts all compounds in the set by their number of atoms.
	*/
	void sortByNumAtoms();

	/** sets the activity of molecules in the set according to the value of a descriptor and a value.
		if true is passed as third argument (default) then descriptorName <= value are considered positive.
		if false is passed as third argument then molecules with descriptorName >= value are considered positive.
		if a descriptor is missing then the molecule is left without activity.
	*/
	void binClassifyFromDescriptor( string descriptorName, float value, bool smallerOrEqual = true );


	/** returns a pointer to the molecule with name aName in the MoleculeSet.
	*/
	Molecule* findFirstMoleculeWithName( string aName ) throw( CError );


	/** removes duplicate entries in the set.
	*/
	void removeDuplicates();

	/** deletes all hidden atoms in all molecules of the set.
	*/
	void deleteHiddenAtoms();

	/** hides all hydrogens in all molecules of the set.
	*/
	void hideHydrogens();

	/** hides all but the largest connected graph in all molecules of the set.
	*/
	void hideSalts( string aReportFileName = "" );

	/** restores the hidden atoms for all molecules.
	*/
	void restoreHiddenAtoms();

	/** NOT DOCUMENTED
	*/
	void addFragmentsToSet( Molecule* aMol, int minAtoms = 1 );

	/** NOT DOCUMENTED
	*/
	void pushFragments( Molecule* aMol, int minAtoms = 1 );

	/** returns the mean distance of all molecules to the moleculeSet barycenter.
	*/
	double diversityBaryMean();

	/** lists the different kinds of atoms present in the set based on their Morgan labels.
	*/
	vector<string> atomsLabelsListing();

	/** lists the different kinds of atoms present in the set based on their symbols.
	*/
	vector<string> atomsSymbolsListing();

	/** lists the different kinds of bonds present in the set.
	*/
	vector<int> bondsListing();

	/** transforms the molecular graphs into graph preventing tottering paths (see (Mahe et al., 2004)).
	*/
	void noTottersTransform();
	
	/** transforms the molecular graphs into '3D complete graphs' with edges labeled by inter atomic distances
	in order to compute the pharmacophore kernel (see (Mahe et al., 2006)).
	*/
	void threeDtransform(int nBins, double distMin, double distMax);

	/** computes the Euclidian distance between the XYZ coordinates of two atoms.
	*/
	void minMaxDistances(double *distMin, double* distMax);


//@}



///@name Kernels and Gram matrices related functions
//@{

	/** calculates the gram matrix of similarity using the marginalized graph kernel for all molecules in
		the MoleculeSet using a specified graph, atom and bond kernel.
		ff aReportFilename is specified a report is saved in that file.
	*/
	void gramCompute(
		double aPq,
		double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2,
					 double(*pt2AtomKernel)(Atom*, Atom*),
					 double(*pt2BondKernel)(Bond*, Bond*), int, int ),
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*),
		int aParameter = 1000,
		string aReportFileName = "",
		int nbThreadsWanted = 1,
		bool silentMode = false,
		bool filterTotters = false);

	/** calculates the gram matrix of similarity using the marginalized graph kernel for all molecules in
		the MoleculeSet using a specified graph, atom and bond kernel.
		ff aReportFilename is specified a report is saved in that file.
	*/
	void gramCompute(
		MoleculeSet* anotherSet,
		double aPq,
		double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2,
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*), int, int ),
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*),
		int aParameter = 1000,
		string aReportFileName = "",
		int nbThreadsWanted = 1,
		bool silentMode = false, 
		bool filterTotters = false);

	/** calculates the gram matrix of similarity using the marginalized graph kernel for all molecules in
		the MoleculeSet using a specified graph, atom and bond kernel.
		ff aReportFilename is specified a report is saved in that file.
	*/
	void gramCompute(
		double aPq,
		double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2,
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*), int, int ),
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*),
		int parameter1 = 1000,
		int parameter2 = 1,
		string aReportFileName = "",
		int nbThreadsWanted = 1,
		bool silentMode = false, 
		bool filterTotters = false);

	/**  "TRUE" gramCompute() function, i.e., the one called by EVERY OTHER GRAMCOMPUTE FUNCTION.
	*/
	void gramCompute(
		MoleculeSet* anotherSet,
		double aPq,
		double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2,
					 double(*pt2AtomKernel)(Atom*, Atom*),
					 double(*pt2BondKernel)(Bond*, Bond*), int, int ),
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*),
		int parameter1 = 1000,
		int parameter2 = 1,
		string aReportFileName = "",
		int nbThreadsWanted = 1,
		bool silentMode = false, 
		bool filterTotters = false);

	/**  3D kernel computation. NOTE : AtomKernel/BondKernel prototyes different from gramCompute.
	*/
	void gramCompute3D( 
			 double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, 
						  double(*pt2AtomKernel)(Atom*, Atom*), 
						  double(*pt2BondKernel)(float,float,float), float), 
			 double(*pt2AtomKernel)(Atom*, Atom*), 
			 double(*pt2BondKernel)(float, float, float),
			 float edgeKernelparameter, 
			 bool silentMode );

	/**  "TRUE" gramCompute3D function, i.e., the one called by EVERY OTHER GRAMCOMPUTE3D FUNCTION.
	*/
	void gramCompute3D( 
			   MoleculeSet* anotherSet,
			   double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, 
						    double(*pt2AtomKernel)(Atom*, Atom*), 
						    double(*pt2BondKernel)(float,float,float), float), 
			   double(*pt2AtomKernel)(Atom*, Atom*), 
			   double(*pt2BondKernel)(float, float, float),
			   float edgeKernelparameter, 
			   bool silentMode );

	/** computes all kernel values between a molecule aMol and all compounds in the set.
	*/
	void kernelCompute(
		Molecule* aMol,
		double(*pt2GraphKernel)( Molecule* mol1, Molecule* mol2, double(*pt2AtomKernel)(Atom*, Atom*), double(*pt2BondKernel)(Bond*, Bond*), int, int ),
		double(*pt2AtomKernel)(Atom*, Atom*),
		double(*pt2BondKernel)(Bond*, Bond*),
		vector<double>* resultsRaw,
		vector<double>* resultsNormal,
		int convergenceCondition = 1000,
		int parameter2 = 1,
		bool silentMode = false
	);


	/** erases the Gram Matrix and sets the gramCalculated flag to false.
	*/
	void resetGramMatrix();

	/** resets all calculated selfkernels for the molecules in the dataset.
	*/
	void resetSelfKernels();

	/** initializes every element gram matrix to the given value.
	*/
	void initializeGram( double value );

	/** initializes the self kernels to the given value.
	*/
	void initializeSelfKernel( double value );

	/** normalizes the gram matrix, i.e. compute gramNormal.
	*/
	void normalizeGram();

	/** normalizes the gram matrix, i.e. compute gramNormal. 
 	WARNING: NORMALIZATION BASED ON THE RAW GRAM MATRIX (instead of self-kernels)
	*/
	void normalizeGram_raw();
	
	/** normalizes the Gram matrix according to the Tanimoto kernel definition (Ralaivola et al. 2005).
	*/
	void normalizeTanimoto();

	/** normalizes the Gram matrix according to the Tanimoto kernel definition, BASED ON THE RAW GRAM MATRIX
	*/
	void normalizeTanimoto_raw();	
	/** normalizes the Gram matrix according to the 'min-max' Tanimoto kernel definition (Ralaivola et al. 2005).
	*/
	void normalizeTanimotoMinMax();

	/** adds a value to a Gram matrix entry.
	*/
	void addToGram( int row, int col, double value );

	/** adds a value to a Gram normal matrix entry.
	*/
	void addToGramNormal( int row, int col, double value );

	/** substracts a value to a Grammatrix entry.
	*/
	void substractToGram( int row, int col, double value );

	/** returns a Gram matrix entry.
	*/
	double getGramValue( int row, int col );


//@}







///@name Output functions
//@{
	/** writes a file with the biological activity of molecules in a format compatible with GIST.
	*/
	void writeActivityFile( string aFilename, bool addActivityExtension = true, string activityDescriptor = ACTIVITY );

	/** writes the gram matrix.
	*/
	void writeGramMatrix( string aFileName, bool normal = false, bool self = false, bool silentMode = false );

	/** writes all self kernel values in a file.
	*/
	void writeSelfKernelList( string aFilename, bool silentMode = false );

	/** writes a MDL structure data (SD) file with the molecules in the moleculeSet
		setting selectedOnly to true outputs only the selectwed compounds.
		a pointer to a vector of strings containing the molecule names for ordered
		output can be specified.
	*/
	void writeSD( string aFileName, bool selectedOnly = false );

	/** writes a MDL structure data (SD) file with the molecules in the moleculeSet
		matching the names given as argument.
	*/
	void writeSubsetSD( string aFileName, vector<string>* anOrder );

	/** writes a KCF file with the molecules in the moleculeSet
		matching the names given as argument.
	*/
	void writeSubsetKCF( string aFileName, vector<string>* anOrder );

	/** writes a KCF file with the whole moleculeSet.
	 */
	void writeKCF( string aFileName, bool selectedOnly = false );

	/** returns a string description of the MoleculeSet.
	*/
	string toString( bool selectedOnly = false );
	/** returns a short string description of the MoleculeSet
		(number of molecules in the set and short description of each molecule).
	*/
	string toStringShort();

	/** returns a long string description of the MoleculeSet.
	*/
	string toStringLong();


	/** writes a description of the moleculeSet to cout.
	*/
	void describe( bool selectedOnly = false );

	/** writes a short description of the moleculeSet to cout.
	*/
	void describeShort();

	/** writes a long description of the moleculeSet to cout.
	*/
	void describeLong();


	/** writes a ';' separated file containing all descriptors for all molecules.
	*/
	void writeDescriptors( string aFileName, bool selectedOnly = false ) throw( CError ) ;

	/** writes a mol file for each molecule in the set to aDirName.
		if selectedOnly == true then only selected molecules are written.
	*/
	long writeMolToDir( string aDirName, bool selectedOnly = false );

	/** writes a dot file for all molecule in the set to aDirectory.
	*/
	long writeDotsToDir( string aDirectory, bool selectedOnly = false, bool perretLabels = false  );


//@}



// DEPRECATED FUNCTIONS : 
//void addMolecule( Molecule* aMolecule, bool updateGram = true );

/** removes the last molecule from the set
*/
//void pop_back();
//void normalizeGram_self();
//void normalizeGram_test();
//void setSelfGram();
	


protected: // Protected attributes

	/** comparison set of the molecule set.
	In test set mode, the comparison set is another set of compounds. In self set mode, the comparison set 
	is the molecule set itself.
 	*/
	MoleculeSet* comparisonSet;

	/** Gram matrix.
	*/
	vector< vector<double> >* gram;

	/** normal Gram matrix.
	*/
	vector< vector<double> >* gramNormal;

	/** contains true if the gram matrix was evaluated for the current MoleculeSet, false otherwise. 
			note that setKashimaKernelProb sets this flag to false.
			{Addition of a new molecule to the set when the flag is set to true induces
			the calculation of a new line to the gram matrix DEPRECATED}
	*/
	bool gramCalculated;

	/** kashima Stop probability for the moleculeSet
			used to set kashimaProb in added molecules.
	*/
	double pq;

	/** convergence condition used in the calculation of the Kashima Kernel.
	*/
	int convergenceCondition;

	/** subset start.
	*/
	int subsetStart;

	/** subset size.
	*/
	int subsetSize;

	/** stores if the activity of the molecules in the set was specified.
	*/
	bool activitySet;


	/** sets the type and name of descriptor to be used when sorting molecules.
		if reverse == true the sorting will be in descending order.
	*/
	void setSortDescriptor( string aName, int aType, bool reverse = false ) throw( CError );

	/** returns the name of the sorting descriptor.
	*/
	string getSortDescriptorName();


	/** vector of molecules
	*/
	//vector<Molecule*> molecules;


};

#endif
