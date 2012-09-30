/****************************************************************************************
					spectrumhelper.cpp 
					------------------
    copyright            : (C) 2012 Mahr Michael
    email                : michael.mahr@gmail.com

rearranged from:
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



#include "spectrumhelper.h"


// -----------------------------------------------
// ------------- UTILITY FUNCTIONS ---------------
// -----------------------------------------------

//--------
//External
//--------

void gramSpectrum_self(SEXP s, int depthMax, int kernelType, double kernelParam, bool onlyN, bool silentMode){

    std::string rtypename("Rcpp_Rmoleculeset");
    Rcpp::S4 s4obj( s );
    if ( !s4obj.is( rtypename.c_str() ) ) {
      Rf_error( (std::string("object is not of the type ")+rtypename).c_str() );
    }
    Rcpp::Environment env( s4obj );
    Rcpp::XPtr<Rmoleculeset> xptr( env.get(".pointer") );
    Rmoleculeset *aSet = static_cast<Rmoleculeset*> (R_ExternalPtrAddr(xptr));

    //---------
    //---------
    //---------

    int depth = -1;

    vector< pathsInMol > molPaths;

    // 2 - get the database description
    // --------------------------------
    // i.e. fill in atom and bond types vectors
    vector<string> atomLabels;
    vector<int> bondTypes;
    
    atomLabels = aSet->atomsLabelsListing();
    bondTypes = aSet->bondsListing();
    
    if( !silentMode ){
      for( int i = 0 ; i < atomLabels.size() ; i++ ){
        Rcpp::Rcout << "atom type no " << i+1 << " ; atomic number = " << atomLabels[i] << endl;
      }
      for( int i = 0 ; i < bondTypes.size() ; i++ ){
        Rcpp::Rcout << "bond type no " << i+1 << " ; bond type = " << bondTypes[i] << endl;
      }
    }

    gramComputeSpectrum_self( aSet, depth, depthMax, kernelType, kernelParam, &molPaths, &atomLabels, &bondTypes, onlyN, silentMode);
    if( !silentMode ){
      Rcpp::Rcout << "gramComputeSpectrum (self) OK" << endl;
    }
}



void gramSpectrum_test(SEXP s, int depthMax, int kernelType, double kernelParam, bool onlyN, bool silentMode){

    std::string rtypename("Rcpp_Rmoleculeset");
    Rcpp::S4 s4obj( s );
    if ( !s4obj.is( rtypename.c_str() ) ) {
      Rf_error( (std::string("object is not of the type ")+rtypename).c_str() );
    }
    Rcpp::Environment env( s4obj );
    Rcpp::XPtr<Rmoleculeset> xptr( env.get(".pointer") );
    Rmoleculeset *aSet = static_cast<Rmoleculeset*> (R_ExternalPtrAddr(xptr));

/*
    std::string rtypename2("Rcpp_Rmoleculeset");
    Rcpp::S4 s4obj2( s2 );
    if ( !s4obj2.is( rtypename2.c_str() ) ) {
      Rf_error( (std::string("object is not of the type ")+rtypename2).c_str() );
    }
    Rcpp::Environment env2( s4obj2 );
    Rcpp::XPtr<Rmoleculeset> xptr2( env2.get(".pointer") );
    Rmoleculeset *aSet2 = static_cast<Rmoleculeset*> (R_ExternalPtrAddr(xptr2));
*/
    Rmoleculeset *aSet2 = aSet->getComparisonSet2();

    //---------
    //---------
    //---------

    int depth = -1;

    vector< pathsInMol > molPaths;
    vector< pathsInMol > molPaths_test;

    // 2 - get the database description
    // --------------------------------
    // i.e. fill in atom and bond types vectors
    vector<string> atomLabels;
    vector<int> bondTypes;
    
    atomLabels = aSet->atomsLabelsListing();
    bondTypes = aSet->bondsListing();
    
    if( !silentMode ){
      for( int i = 0 ; i < atomLabels.size() ; i++ ){
        Rcpp::Rcout << "atom type no " << i+1 << " ; atomic number = " << atomLabels[i] << endl;
      }
      for( int i = 0 ; i < bondTypes.size() ; i++ ){
        Rcpp::Rcout << "bond type no " << i+1 << " ; bond type = " << bondTypes[i] << endl;
      }
    }
    
    // 3 - compute the gram matrix
    // ----------------------------
    gramComputeSpectrum_test( aSet, aSet2,  depth, depthMax, kernelType, kernelParam, &molPaths, &molPaths_test,  &atomLabels, &bondTypes, onlyN, silentMode);
    if( !silentMode ){
      Rcpp::Rcout << "gramComputeSpectrum (test) OK" << endl;
    }

}

//--------
//internal
//--------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void gramComputeSpectrum_self(MoleculeSet* aSet, int depth, int depthMax, int kernelType, double kernelParam, vector< pathsInMol >*  molPaths, vector<string>* atomLabels, vector<int>* bondTypes, bool onlyN, bool silentMode){
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  vector< pathsInMol> molPaths_NEW;
  
  depth = depth + 1;
  
  if (depth == 0){ // depth = 0  --> initialization step
    // for each kind of atoms, find paths starting points and launch the recursion
    for(int i = 0 ;  i < atomLabels->size() ; i++){
      molPaths->clear();
      
      if( !silentMode ){
	Rcpp::Rcout <<  " \t finding paths starting from atoms labeled = " << (*atomLabels)[i] << endl;
      }
      
      init_path_pointers(aSet, molPaths, (*atomLabels)[i], kernelType);
      
      if( !onlyN ){
	updateGram_self(aSet,molPaths, kernelType, kernelParam, depth);
	updateSelfKernel(aSet,molPaths, kernelType, kernelParam, depth);
      }
      
      if(depthMax == 0){  //...then, just include the weights in the gram matrix
	if(onlyN){
	  updateGram_self(aSet,molPaths, kernelType, kernelParam, depth);
	  updateSelfKernel(aSet,molPaths, kernelType, kernelParam, depth);
	}
	continue;
      }else{
	gramComputeSpectrum_self(aSet, depth, depthMax, kernelType, kernelParam, molPaths, atomLabels,  bondTypes, onlyN, silentMode);
      }
    }
    
  } // end "if(depth == 0)"
  else{  // if depth >= 1 :  we can try to extend the paths stored in path_pointers
	for(int j = 0 ; j < atomLabels->size() ; j++){
		for(int k = 0 ; k < bondTypes->size() ; k++){
			// for each atom type/bond type : extend the paths
			updatePaths(aSet, (*atomLabels)[j], (*bondTypes)[k], molPaths, &molPaths_NEW, kernelType, depth);

			int nonEmpty;
			nonEmpty = molPaths_NEW.size();
    			if( nonEmpty > 0 ){  // i.e : if at least 1 vector of paths is non-empty, we go further in the recursion
				if( !onlyN ){
					updateGram_self(aSet, &molPaths_NEW, kernelType, kernelParam, depth);
					updateSelfKernel(aSet,&molPaths_NEW, kernelType, kernelParam, depth);
				}
				if(depth == depthMax){
					if( onlyN ){
						updateGram_self(aSet, &molPaths_NEW, kernelType, kernelParam, depth);
						updateSelfKernel(aSet,&molPaths_NEW, kernelType, kernelParam, depth);
					}
					continue;
				}else{
					gramComputeSpectrum_self(aSet, depth, depthMax, kernelType, kernelParam, &molPaths_NEW, atomLabels, bondTypes, onlyN, silentMode);
				}
			}// end "if(nonEmpty > 0)"

		}
	}
    
  } // end "else if(depth == 0)"
}



//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void gramComputeSpectrum_test(MoleculeSet* aSet, MoleculeSet* aSet2, int depth, int depthMax, int kernelType, double kernelParam, vector< pathsInMol >*  molPaths, vector< pathsInMol >*  molPaths_test, vector<string>* atomLabels, vector<int>* bondTypes, bool onlyN, bool silentMode){
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
  vector< pathsInMol > molPaths_NEW;
  vector< pathsInMol > molPaths_test_NEW;
  
  depth = depth + 1;

  if (depth == 0){ // depth = 0  --> initialization step
    // for each kind of atoms, find paths starting points and launch the recursion
    for(int i = 0 ;  i < atomLabels->size() ; i++){
      molPaths->clear();
      molPaths_test->clear();
      
      if( !silentMode ){
	Rcpp::Rcout << "\t - finding paths starting from atoms labeled = " << (*atomLabels)[i] << endl;
      }
      init_path_pointers(aSet, molPaths, (*atomLabels)[i], kernelType);
      init_path_pointers(aSet2, molPaths_test, (*atomLabels)[i], kernelType);
      
      if( !onlyN ){
	updateGram_test(aSet, aSet2, molPaths, molPaths_test, kernelType, kernelParam, depth);
	updateSelfKernel(aSet,molPaths, kernelType, kernelParam, depth);
	updateSelfKernel(aSet2,molPaths_test, kernelType, kernelParam, depth);
      }
      
      if(depthMax == 0){
	if( onlyN ){
	  updateGram_test(aSet, aSet2, molPaths, molPaths_test, kernelType, kernelParam, depth);
	  updateSelfKernel(aSet,molPaths, kernelType, kernelParam, depth);
	  updateSelfKernel(aSet2,molPaths_test, kernelType, kernelParam, depth);
	}
	continue;
      }else{
	gramComputeSpectrum_test(aSet, aSet2, depth, depthMax, kernelType, kernelParam, molPaths, molPaths_test, atomLabels, bondTypes, onlyN, silentMode);
      }
    }
  } // end "if (depth == 0)"
  else{  // if depth >= 1, we can try to extend the paths previously stored in path_pointers
	for(int j = 0 ; j < atomLabels->size() ; j++){
		for(int k = 0 ; k < bondTypes->size() ; k++){
			// for each atom type/bond type : extend the paths
			updatePaths(aSet, (*atomLabels)[j], (*bondTypes)[k], molPaths, &molPaths_NEW, kernelType, depth);
			updatePaths(aSet2, (*atomLabels)[j], (*bondTypes)[k], molPaths_test, &molPaths_test_NEW, kernelType, depth);
	      
			int nonEmpty1, nonEmpty2;
			nonEmpty1 = molPaths_NEW.size();
			nonEmpty2 = molPaths_test_NEW.size();
			if(nonEmpty1 > 0 || nonEmpty2 > 0){  // i.e : if at least 1 vector of paths is non-empty IN THE TWO SETS, we go further in the recursion
				if( !onlyN ){
					updateGram_test(aSet, aSet2, &molPaths_NEW, &molPaths_test_NEW, kernelType, kernelParam, depth);
					updateSelfKernel(aSet,&molPaths_NEW, kernelType, kernelParam, depth);
					updateSelfKernel(aSet2,&molPaths_test_NEW, kernelType, kernelParam, depth);
				}
	  			if(depth == depthMax){
					if( onlyN ){
						updateGram_test(aSet, aSet2, &molPaths_NEW, &molPaths_test_NEW, kernelType, kernelParam, depth);
						updateSelfKernel(aSet,&molPaths_NEW, kernelType, kernelParam, depth);
						updateSelfKernel(aSet2,&molPaths_test_NEW, kernelType, kernelParam, depth);
					}
					continue;
				}else{
					gramComputeSpectrum_test(aSet, aSet2, depth, depthMax, kernelType, kernelParam, &molPaths_NEW, &molPaths_test_NEW, atomLabels, bondTypes, onlyN, silentMode);
				}
			} // end if(nonEmpty1 > 0 || nonEmpty2 > 0)

		}
	}     
  } // end "else if(depth == 0)"

}



//************************************************************************************************
// -----------------------------------------------------------------------------------------------
void init_path_pointers(MoleculeSet* aSet, vector< pathsInMol >* molPaths, string atomLabel,int kernelType ){
// -----------------------------------------------------------------------------------------------
      // --> fill in "path_pointers" for starting atoms labeled as "atomLabel"

  path aPath;
  pathsInMol aPathInMol;
  int ind = 0;

  vector<Molecule*>::iterator aMol;
  // for each molecule of the set :
  for(aMol = (*aSet).begin() ; aMol != (*aSet).end() ; aMol++ ){
    bool first = true;
    aPathInMol.path_pointers.clear();
    vector<Atom*>::iterator anAtom;
    // for each atom of the molecule
    for(anAtom = (*aMol)->beginAtom() ; anAtom != (*aMol)->endAtom() ; anAtom++){
      // if the atomic number corresponds to the one needed (i.e. "atomNumb")
      if((*anAtom)->getMorganLabel(true) == atomLabel){
	if(first == true){ 
	  if(kernelType == 3){ // i.e., marginalized kernel approximation
	    aPath.weight = (*anAtom)->getKashimaPS();
	  }
	  else{
	    aPath.weight = 1;
	  }
	  
	  if( aPath.weight > 0){
	    aPath.endAtom = *anAtom;
	    aPathInMol.path_pointers.push_back(aPath);
	    aPathInMol.molName = (*aMol)->getName();
	    aPathInMol.molInd = ind;
	    first = false;
	  }
	} // end "if(first == true)"
	else{
	  if(kernelType == 3){
	    aPath.weight = (*anAtom)->getKashimaPS();
	  }
	  else{
	    aPath.weight = 1;
	  }
	  
	  if( aPath.weight > 0){
	    aPath.endAtom = *anAtom;
	    aPathInMol.path_pointers.push_back(aPath);
	  }
	} // end "else"
      } // end "if((*anAtom)->getMorganLabel(true) == atomLabel)"
    }
    if(aPathInMol.path_pointers.size() > 0 ){
      (*molPaths).push_back(aPathInMol);
    } 
    
    ind++;
  }
  
}



// ************************************************************************************************************************************************
// ------------------------------------------------------------------------------------------------------------------------------------------------
void updatePaths(MoleculeSet *aSet, string atomLabel, int bondType, vector< pathsInMol >* molPaths, vector< pathsInMol >* molPaths_NEW, int kernelType, int depth){
// ------------------------------------------------------------------------------------------------------------------------------------------------
  // --> extends path stored in "path_pointers" to atom of type "atomNumb". Stores the results in "path_pointers_NEW"

  path aPath;
  pathsInMol aPathInMol; // object initialized for each molecule present in the input molPaths, possibly inserted in the output molPaths_NEW
  

  molPaths_NEW->clear();

  // for each molecule in molPaths
  for(int i = 0 ; i < molPaths->size() ; i++){
    bool first = true;
    aPathInMol.path_pointers.clear();
    //for each previous path pointer
    for(int j = 0 ; j < (*molPaths)[i].path_pointers.size() ; j++ ){
      map<Atom* , Bond*>::iterator aBond;
      // along each bond of the path's ending atom :
      for(aBond = (*molPaths)[i].path_pointers[j].endAtom->beginBond() ; aBond != (*molPaths)[i].path_pointers[j].endAtom->endBond() ; aBond++ ){
	
	if( (
	     ( depth > 1)
	     &&  ( aBond->second->getLabel() == bondType ) 
	     &&  ( aBond->first->getMorganLabel(true) == atomLabel ) 
	     &&  ( aBond->first->getId() !=  (*molPaths)[i].path_pointers[j].previousAtom->getId() )  
	     ) // if depth > 1 : we have to prevent tottering paths 
	    ||
	    (
	     (depth == 1) 
	     &&  ( aBond->second->getLabel() == bondType ) 
	     &&  ( aBond->first->getMorganLabel(true) == atomLabel )
	     ) // if depth = 1 : path can't totter	   
	    )
	  {
	    
	    if(first == true){
	      if(kernelType == 3){ // i.e., marginalized kernel approximation
		if(depth == 1) {
		  aPath.weight = (*molPaths)[i].path_pointers[j].weight * aBond->second->getKashimaPT();
		} // if depth = 1 : paths can't totter
		else{
		  aPath.weight = (*molPaths)[i].path_pointers[j].weight * ( 1.0 - (*molPaths)[i].path_pointers[j].endAtom->getKashimaPQ() ) / ( (*molPaths)[i].path_pointers[j].endAtom->numBonds() - 1 );
		} // if depth > 1 : modify weight for the non-tottering paths
	      }
	      else{
		aPath.weight = 1;
	      }
	      
	      if(aPath.weight > 0 ){
		aPath.endAtom = aBond->first;
		aPath.previousAtom = (*molPaths)[i].path_pointers[j].endAtom;
		aPathInMol.path_pointers.push_back(aPath);
		aPathInMol.molName = (*molPaths)[i].molName;
		aPathInMol.molInd = (*molPaths)[i].molInd;
		first = false;
	      }
	    } // end "if(first == true)"
	    else{
	      if(kernelType == 3){
		if(depth == 1){
		  aPath.weight = (*molPaths)[i].path_pointers[j].weight * aBond->second->getKashimaPT();
		}
		else{
		  aPath.weight = (*molPaths)[i].path_pointers[j].weight * ( 1.0 - (*molPaths)[i].path_pointers[j].endAtom->getKashimaPQ() ) / ( (*molPaths)[i].path_pointers[j].endAtom->numBonds() - 1 );
		}
	      }
	      else{
		aPath.weight = 1;
	      }
	      
	      if(aPath.weight > 0 ){    
		aPath.endAtom = aBond->first;
		aPath.previousAtom = (*molPaths)[i].path_pointers[j].endAtom;
		aPathInMol.path_pointers.push_back(aPath);
	      }
	    }
	  }
      }
    }
    if(aPathInMol.path_pointers.size() > 0){
      (*molPaths_NEW).push_back(aPathInMol);
    }
  }
  
}



// ***********************************************************************************************************************
// -----------------------------------------------------------------------------------------------------------------------
void updateGram_self(MoleculeSet *aSet, vector< pathsInMol >*  molPaths, int kernelType, double kernelParam, int depth){
// -----------------------------------------------------------------------------------------------------------------------

  vector<double> sumWeight;
  
  if(kernelType == 3){ // i.e., marginalized approximation --> fill in a vector summing the weight of the paths for each molecule
    for(int i = 0 ; i < molPaths->size() ; i++){
      sumWeight.push_back(0.0);
      for(int j = 0 ; j < (*molPaths)[i].path_pointers.size() ; j++ ){
	if( depth >= 1 && (*molPaths)[i].path_pointers[j].endAtom->numBonds() == 1 ) 
	  sumWeight[i] += (*molPaths)[i].path_pointers[j].weight;
	else 
	  sumWeight[i] += (*molPaths)[i].path_pointers[j].weight * (*molPaths)[i].path_pointers[j].endAtom->getKashimaPQ();
      }
    }
  }


  //  update the gram matrix
  double update, update2;
  for(int i = 0 ; i < molPaths->size() ; i++){
    for(int j = i ; j < molPaths->size() ; j++){
      // compute the update factor depending on the kernel type
      switch(kernelType){
      case 0 : // 'spectrum'
	update = (*molPaths)[i].path_pointers.size() * (*molPaths)[j].path_pointers.size();
	break;
      case 1 : // tanimoto
	update = 1.0;
	break;
      case 2 : // min/max tanimoto
	update = min( (*molPaths)[i].path_pointers.size(), (*molPaths)[j].path_pointers.size() );
	// NB : for the min/max kernel we have to keep the maximum values too, for further normalization
	update2 = max( (*molPaths)[i].path_pointers.size(), (*molPaths)[j].path_pointers.size() );
	(*aSet).addToGramNormal( (*molPaths)[i].molInd , (*molPaths)[j].molInd, update2);
	if( i!= j){
	  (*aSet).addToGramNormal( (*molPaths)[j].molInd , (*molPaths)[i].molInd, update2);   
	}
	break;
      case 3 : // marginalized approximation
	update =  sumWeight[i] * sumWeight[j]; 
	break;
      case 4 : // lambda^k
	update =  (*molPaths)[i].path_pointers.size() * (*molPaths)[j].path_pointers.size() *  pow(kernelParam, depth+1);
	// NB : we put depth+1 because 0 edge <-> 1 atom => x lambda
	// Rq : equivalent to case 0 with lambda = 1
	break;	
      }
      
      // add the update factor to the gram matrix
      (*aSet).addToGram( (*molPaths)[i].molInd , (*molPaths)[j].molInd, update);
      if( i!= j){
	(*aSet).addToGram( (*molPaths)[j].molInd , (*molPaths)[i].molInd, update);   
      }
      
   }
  }

}
      


//**************************************************************************
//--------------------------------------------------------------------------
void updateGram_test(MoleculeSet *aSet, MoleculeSet *aSet2, vector< pathsInMol >*  molPaths, vector< pathsInMol >*  molPaths_test, int kernelType, double kernelParam, int depth){
//--------------------------------------------------------------------------
  

  vector<double> sumWeight, sumWeight_test;
  
  if(kernelType == 3){ // i.e., marginalized approximation --> fill in a vector summing the weight of the paths for each molecule (train + test)
    for(int i = 0 ; i < molPaths->size() ; i++){
      sumWeight.push_back(0.0);
      for(int k = 0 ; k < (*molPaths)[i].path_pointers.size() ; k++){
	if( (*molPaths)[i].path_pointers[k].endAtom->numBonds() ==1  && depth >= 1) // NB : no-tottering implementation
	  sumWeight[i] += (*molPaths)[i].path_pointers[k].weight;
	else
	  sumWeight[i] += (*molPaths)[i].path_pointers[k].weight *  (*molPaths)[i].path_pointers[k].endAtom->getKashimaPQ();
      }
    }

    for(int i = 0 ; i < molPaths_test->size() ; i++){
      sumWeight_test.push_back(0.0);
      for(int k = 0 ; k < (*molPaths_test)[i].path_pointers.size() ; k++){
	if( (*molPaths_test)[i].path_pointers[k].endAtom->numBonds() == 1  && depth >= 1) // NB : no-tottering implementation
	  sumWeight_test[i] += (*molPaths_test)[i].path_pointers[k].weight;     
	else
	  sumWeight_test[i] += (*molPaths_test)[i].path_pointers[k].weight *  (*molPaths_test)[i].path_pointers[k].endAtom->getKashimaPQ();
      }
    }
    
  } // end  "if(kernelType == 3)"


  // update the gram matrix
  double update, update2;
  for(int i = 0 ; i < molPaths->size() ; i++){
      for(int j = 0 ; j < molPaths_test->size() ; j++){
	// compute the update factor depending on the kernel type
	switch(kernelType){
	case 0 : // spectrum
	  update = (*molPaths)[i].path_pointers.size() *  (*molPaths_test)[j].path_pointers.size();
	  break;
	case 1 : // tanimoto
	  update = 1.0;
	  break;
	case 2 : // min/max tanimoto
	  update = min( (*molPaths)[i].path_pointers.size(), (*molPaths_test)[j].path_pointers.size() ); 
	  // NB : for the min/max kernel we have to keep the maximum values too, for further normalization
	  update2 = max( (*molPaths)[i].path_pointers.size(), (*molPaths_test)[j].path_pointers.size() ); 
	  (*aSet).addToGramNormal( (*molPaths)[i].molInd , (*molPaths_test)[j].molInd, update2);
	  break;
	case 3 : // marginalized approximation
	  update =  sumWeight[i] * sumWeight_test[j];
	  break;
	case 4 : // lambda^k
	  update = (*molPaths)[i].path_pointers.size() *  (*molPaths_test)[j].path_pointers.size() * pow(kernelParam, depth+1);
	  // NB : we put depth+1 because 0 edge <-> 1 atom => x lambda
	  // Rq : equivalent to case 0 with lambda = 1
	  break;
	}

	// add the update factor to the gram matrix
	(*aSet).addToGram( (*molPaths)[i].molInd , (*molPaths_test)[j].molInd, update);
      }
  }


}





// **************************************************************************
// --------------------------------------------------------------------------
void updateSelfKernel(MoleculeSet* aSet, vector< pathsInMol >* molPaths, int kernelType, double kernelParam, int depth){
// --------------------------------------------------------------------------


  vector<double> sumWeight;
  Molecule *aMol;


  if(kernelType == 3){ // i.e., marginalized approximation --> fill in a vector summing the weight of the paths for each molecule 
    for(int i = 0 ; i < molPaths->size() ; i++){
      sumWeight.push_back(0.0);
      for(int j = 0 ; j < (*molPaths)[i].path_pointers.size() ; j++ ){
	if( (*molPaths)[i].path_pointers[j].endAtom->numBonds() == 1 && depth >= 1) // NB : no-tottering implementation
	  sumWeight[i] += (*molPaths)[i].path_pointers[j].weight;     
	else 
	  sumWeight[i] += (*molPaths)[i].path_pointers[j].weight * (*molPaths)[i].path_pointers[j].endAtom->getKashimaPQ();
      }
    }
  }

 
  // update the self kernel values
  double update; 
 for(int i = 0 ; i < molPaths->size() ; i++){
    //aMol = (*aSet)[(*molPaths)[i].molName];
    //cout << "### MolInd = " << (*molPaths)[i].molInd << endl;
    aMol = (*aSet)[(*molPaths)[i].molInd];
    // compute the update facor depending on the kernel type
    switch(kernelType){
    case 0 : // spectrum
      update =  (*molPaths)[i].path_pointers.size() *  (*molPaths)[i].path_pointers.size();
      break;
    case 1 : // tanimoto
      update = 1.0;
      break;
    case 2 : // min/max tanimoto
      update =  (*molPaths)[i].path_pointers.size();
      break;
    case 3 : // marginalized approximation
      update = sumWeight[i]*sumWeight[i];
      break;
    case 4 : // lambda^k
      update =  (*molPaths)[i].path_pointers.size() *  (*molPaths)[i].path_pointers.size() *  pow(kernelParam, depth+1);
      break;
    }

    // add the update factor to the self kernel values
    aMol->addToSelfKernel(update);
 }


}

RCPP_MODULE(spectrumhelper) {
	Rcpp::function( "gramSpectrum_self", &gramSpectrum_self );
	Rcpp::function( "gramSpectrum_test", &gramSpectrum_test );
}




