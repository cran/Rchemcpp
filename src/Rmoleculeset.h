
#ifndef RMOLECULESET_H
#define RMOLECULESET_H

#include <Rcpp.h>
#include <RcppClassic.h>
#include <Rinternals.h>

#include "Rmolecule.h"

#include <moleculeset.h>
#include <moleculeutils.h>

class Rmoleculeset : public MoleculeSet {
public:
	

	Rmoleculeset()
	{
		//superclass constructor is called...

		comparisonSet = NULL;
	}	

	~Rmoleculeset()
	{
		MoleculeSet::deleteAll(); //It is assumed that all molecules are a part of the set!

		if (comparisonSet != NULL)
			if (comparisonSet != this)
				delete comparisonSet; //Also the comparisonset is a part of the set!

	}



	Rmoleculeset * getComparisonSet2() //no R accessor
	{return (Rmoleculeset*)comparisonSet;}

	vector< vector<double> > * getGram2() //no R accessor
	{return (gram);}

	vector< vector<double> > * getGramNormal2() //no R accessor
	{return (gramNormal);}

	Rmoleculeset(Rmoleculeset * aSet) //copy constructor (more or less)
	{
		gramCalculated = (*aSet).gramCalculated;
		pq = (*aSet).pq;
		convergenceCondition = (*aSet).convergenceCondition;
		subsetStart = (*aSet).subsetStart;
		subsetSize = (*aSet).subsetSize;
		activitySet = (*aSet).activitySet;

		if (aSet->getComparisonSet2() != NULL)
			comparisonSet = new Rmoleculeset(aSet->getComparisonSet2()); //deep copy (recursive)
		else
			comparisonSet = NULL;
	
		if (aSet->getGram2() != NULL)
			(*gram) = (*aSet->getGram2()); //deep copy
		else
			gram = NULL;

		if (aSet->getGramNormal2() != NULL)
			(*gramNormal) = (*aSet->getGramNormal2()); //deep copy
		else
			gramNormal = NULL;

		//deep copy all contained molecules
		this->clear();
		for( vector<Molecule*>::iterator m1 = aSet->begin(); m1 != aSet->end(); m1++ )
		{
			this->push_back(new Molecule(*(*m1),false));
		}

	}



	//RcppModules in Rcpp 0.9.13 does not allow exposing base class method, so they need to be redefined here
	int addSD2( string aFileName, bool genericAtomTypeFlag = false, long beginMolecule = -1, long endMolecule = -1 )
	{return MoleculeSet::addSD(aFileName, genericAtomTypeFlag, beginMolecule, endMolecule);}

	//R or RcppModules also seems not to allow overloaded methods in Rcpp 0.9.13
	//It also does not seem to allow "Formal argument specification" (default-values) with ".method" member-functions. 
	int addSD( string aFileName, bool genericAtomTypeFlag)
	{return MoleculeSet::addSD(aFileName, genericAtomTypeFlag);}

	//--

	int addKCF2( string aFileName, long beginMolecule = -1, long endMolecule = -1  )
	{return MoleculeSet::addKCF(aFileName, beginMolecule, endMolecule);}

	int addKCF( string aFileName)
	{return MoleculeSet::addKCF(aFileName);}

	//--

	uint numMolecules()
	{return MoleculeSet::numMolecules();}

	void hideHydrogens()
	{MoleculeSet::hideHydrogens();}

	void setMorganLabels(int order)
	{MoleculeSet::setMorganLabels(order);}

	void setKashimaKernelParam( double aPq, int aConvergenceCondition, bool skipSkeleton = false )
	{MoleculeSet::setKashimaKernelParam( aPq, aConvergenceCondition, skipSkeleton );}

	//avoid a segfault...
	void initializeGram( double value )
	{
		if (comparisonSet == NULL)
		{
		  CError e = CError(MOLECULENOTFOUND, "Comparisonset is not yet set");
		  e.describe();
		  throw(e);
		}

		MoleculeSet::initializeGram(value);
	}

	void initializeSelfKernel( double value )
	{MoleculeSet::initializeSelfKernel(value);}

	void setComparisonSetCopy(SEXP s)
	{
		//if (comparisonSet != NULL)
		delete comparisonSet;

		std::string rtypename("Rcpp_Rmoleculeset");

		Rcpp::S4 s4obj( s );
		if ( !s4obj.is( rtypename.c_str() ) ) {
			Rf_error( (std::string("object is not of the type ")+rtypename).c_str() );
		}

		Rcpp::Environment env( s4obj );
		Rcpp::XPtr<Rmoleculeset> xptr( env.get(".pointer") );

		Rmoleculeset *o = static_cast<Rmoleculeset*> (R_ExternalPtrAddr(xptr));
			  
		//Note: it could already be assigned here - however it is cloned to avoid scope/garbage-collection issues
  		//MoleculeSet::setComparisonSet(o);  

  		MoleculeSet::setComparisonSet(new Rmoleculeset(o)); //invoke copy-constructor to clone object  
	}

	void setComparisonSetSelf()
	{
		comparisonSet = this;
	}

	SEXP getComparisonSet()
	{
		Rcpp::XPtr<Rmoleculeset> xp((Rmoleculeset*)comparisonSet, false); //do NOT mark as finalizable!
		Rcpp::Function maker=Rcpp::Environment::Rcpp_namespace()[ "cpp_object_maker"];
		return maker ( typeid(Rmoleculeset).name() , xp );
	}

	void normalizeTanimoto()
	{MoleculeSet::normalizeTanimoto();}

	void normalizeTanimotoMinMax()
	{MoleculeSet::normalizeTanimotoMinMax();}

	void normalizeGram()
	{MoleculeSet::normalizeGram();}

	void noTottersTransform()
	{MoleculeSet::noTottersTransform();}

	void readPartialCharges(string fileName)
	{MoleculeSet::readPartialCharges(fileName);}

	void setMorganChargesLabels(double threshold)
	{MoleculeSet::setMorganChargesLabels(threshold);}

	//R cannot pass C++ function pointers - therfore it is wrapped here bluntly
	//the function now is a litte more intelligent, because it recognizes whether there is a comparison-set
	void gramCompute3D (float edgeKernelParameter, bool silentMode, bool external, bool rbf)
	{

		double(*atomkernel)(Atom*, Atom*);
		double(*bondkernel)(float, float, float);

		if (external)
			atomkernel = MoleculeUtils::atomKernelExternalMatrix;
		else
			atomkernel = MoleculeUtils::atomKernelMorganLabel;

		if (rbf)
			bondkernel = MoleculeUtils::threeDedgeKernelRBF;
		else
			bondkernel = MoleculeUtils::threeDedgeKernelTriangle;


		if (comparisonSet == NULL)
		{
			MoleculeSet::gramCompute3D(
				MoleculeUtils::threeDkernel,
				atomkernel,
				bondkernel,
				edgeKernelParameter,
				silentMode);

			//sorts out a very ugly hack in chemcpp...
			comparisonSet = NULL;
		}
		else if(comparisonSet == this)
		{
			MoleculeSet::gramCompute3D(
				MoleculeUtils::threeDkernel,
				atomkernel,
				bondkernel,
				edgeKernelParameter,
				silentMode);
		}
		else
		{
			MoleculeSet::gramCompute3D(comparisonSet,
				MoleculeUtils::threeDkernel,
				atomkernel,
				bondkernel,
				edgeKernelParameter,
				silentMode);
		}
	}



	void gramCompute (double stopP, int converg, int fromN, int nbThreadsWanted, bool silentMode, bool filterTotters, bool external)
	{
		double(*atomkernel)(Atom*, Atom*);

		if (external)
			atomkernel = MoleculeUtils::atomKernelExternalMatrix;
		else
			atomkernel = MoleculeUtils::atomKernelMorganLabel;

		
		if (comparisonSet == NULL)
		{
			MoleculeSet::gramCompute( 	stopP,
							MoleculeUtils::moleculeKernel,
							atomkernel,
							MoleculeUtils::bondKernelType,
							converg, fromN, "",
							nbThreadsWanted, silentMode, filterTotters );

			//sorts out a very ugly hack in chemcpp...
			comparisonSet = NULL;
		}
		else if(comparisonSet == this)
		{
			MoleculeSet::gramCompute( 	stopP,
							MoleculeUtils::moleculeKernel,
							atomkernel,
							MoleculeUtils::bondKernelType,
							converg, fromN, "",
							nbThreadsWanted, silentMode, filterTotters );			
		}
		else
		{
			MoleculeSet::gramCompute(comparisonSet, stopP,
							MoleculeUtils::moleculeKernel,
							atomkernel,
							MoleculeUtils::bondKernelType,
							converg, fromN, "",
							nbThreadsWanted, silentMode, filterTotters );
		}
	}


	void addMoleculeCopy(SEXP s)
	{
		std::string rtypename("Rcpp_Rmolecule");

		Rcpp::S4 s4obj( s );
		if ( !s4obj.is( rtypename.c_str() ) ) {
			Rf_error( (std::string("object is not of the type ")+rtypename).c_str() );
		}

		Rcpp::Environment env( s4obj );
		Rcpp::XPtr<Rmoleculeset> xptr( env.get(".pointer") );

		Rmolecule *o = static_cast<Rmolecule*> (R_ExternalPtrAddr(xptr));
			    		
		MoleculeSet::addMoleculeCopy(o); 		
	}

	vector< vector<double> > getGram()
	{return (*gram);}

	vector< vector<double> > getGramNormal()
	{return (*gramNormal);}

	//currently not expicitly necessary
	std::vector<std::string> atomsLabelsListing()
	{return MoleculeSet::atomsLabelsListing();}

	std::vector<int> bondsListing()
	{return MoleculeSet::bondsListing();}

	void normalizeTanimoto_raw()
	{MoleculeSet::normalizeTanimoto_raw();}

	void normalizeGram_raw()
	{MoleculeSet::normalizeGram_raw();}

	void writeGramMatrix(string aFileName, bool normal, bool self, bool silentMode)
	{MoleculeSet::writeGramMatrix( aFileName, normal, self, silentMode );}
	 
	void writeSelfKernelList( string aFileName, bool silentMode )
	{MoleculeSet::writeSelfKernelList(aFileName,silentMode);}

private:


};



#endif
