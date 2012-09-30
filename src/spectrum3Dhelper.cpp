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



#include "spectrum3Dhelper.h"

// -----------------------------------------------
// ------------- UTILITY FUNCTIONS ---------------
// -----------------------------------------------

//--------
//External
//--------

void gramSpectrum3D_self(SEXP s, int depthMax, int kernelType, int nBins, double distMin, double distMax, bool silentMode){

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

    vector< pathsInMol3D > molPaths;

    // 2 - get the database description
    // --------------------------------
    // i.e. fill in atom and bond types vectors
    vector<string> atomLabels;
    vector<int> bondTypes;

    atomLabels = aSet->atomsLabelsListing();
        
    for(int i = 1 ; i <= nBins; i++){
      bondTypes.push_back(i);
    }

    if( !silentMode ){
      for( int i = 0 ; i < atomLabels.size() ; i++ ){
        Rcpp::Rcout << "atom type no " << i+1 << " ; atomic number = " << atomLabels[i] << endl;
      }
      //for( int i = 0 ; i < bondTypes.size() ; i++ ){
	//cout << "bond type no " << i+1 << " ; bond type = " << bondTypes[i] << endl;
      //}
    }
    
        
    // 3 - transform the dataset from 2D to 3D
    // ---------------------------------------
    	// OLD VERSION: distMin and distMin were computed from the data
	//            --> now: specified as parameters
    	//double distMin = 100000.0;
    	//double distMax = 0.0;
    	//aSet.minMaxDistances(&distMin, &distMax);

    if( !silentMode ){
      Rcpp::Rcout << " - distMin = " << distMin << endl;
      Rcpp::Rcout << " - distMax = " << distMax << endl;
      Rcpp::Rcout << " - nBins = " << nBins << endl;
      Rcpp::Rcout << "   --> binSize = " << (1.0001*distMax - distMin)/nBins << endl;
    }
    
    aSet->threeDtransform(nBins, distMin, distMax);


//     vector<Molecule*>::iterator aMol;
//     for(aMol = aSet.begin() ; aMol != aSet.end() ; aMol++ ){
//       (*aMol)->describe();
//       vector<Atom*>::iterator anAtom;
//       //vector<Bond*>::iterator aBond;
//       map<Atom*, Bond*>::iterator aBond;

//       for(anAtom = (*aMol)->beginAtom(); anAtom != (*aMol)->endAtom() ; anAtom++ ){
// 	//for( aBond = (*anAtom)->beginBond(); aBond != (*anAtom)->endBond(); aBond++)
// 	// cout << "bondLabel = " << aBond->second->getLabel() << endl;
// 	cout << "Atom Label after 3D-transfo ; " << (*anAtom)->getMorganLabel() << endl;
//       }
//     }
//     exit(12);



    // 4 - compute the gram matrix
    // ---------------------------
    gramComputeSpectrum3D_self( aSet, depth, depthMax, kernelType, &molPaths, &atomLabels, &bondTypes, silentMode);

    if( !silentMode ){
      Rcpp::Rcout << "gramComputeSpectrum (self) OK" << endl;
    }
}



void gramSpectrum3D_test(SEXP s, int depthMax, int kernelType, int nBins, double distMin, double distMax, bool silentMode){

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

    vector< pathsInMol3D > molPaths;
    vector< pathsInMol3D > molPaths_test;

    // 2 - get the database description
    // --------------------------------
    // i.e. fill in atomLabels  vectors
    vector<string> atomLabels;
    vector<int> bondTypes;

    atomLabels = aSet->atomsLabelsListing();

    for(int i = 1 ; i <= nBins; i++){
      bondTypes.push_back(i);
    }

    if( !silentMode ){
      for( int i = 0 ; i < atomLabels.size() ; i++ ){
        Rcpp::Rcout << "atom type no " << i+1 << " ; atomic number = " << atomLabels[i] << endl;
      }
      //for( int i = 0 ; i < bondTypes.size() ; i++ ){
	//cout << "bond type no " << i+1 << " ; bond type = " << bondTypes[i] << endl;
      //}
    }
    
    // 3 - transform the dataset from 2D to 3D
    // ---------------------------------------
    	// OLD VERSION: distMin and distMin were computed from the data
	//            --> now: specified as parameters
	//double distMin = 100000.0;
    	//double distMax = 0.0;
    	//aSet.minMaxDistances(&distMin, &distMax);

    if( !silentMode ){
      Rcpp::Rcout << " - distMin = " << distMin << endl;
      Rcpp::Rcout << " - distMax = " << distMax << endl;
      Rcpp::Rcout << " - nBins = " << nBins << endl;
      Rcpp::Rcout << "   --> binSize = " << (1.0001*distMax - distMin)/nBins << endl;
    }
    
    aSet->threeDtransform(nBins, distMin, distMax);
    aSet2->threeDtransform(nBins, distMin, distMax);


    // 4 - compute the gram matrix
    // ----------------------------
    gramComputeSpectrum3D_test( aSet, aSet2, depth, depthMax, kernelType, &molPaths, &molPaths_test, &atomLabels, &bondTypes, silentMode);

    if( !silentMode ){
      Rcpp::Rcout << "gramComputeSpectrum (test) OK" << endl;
    }

}

//--------
//internal
//--------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void gramComputeSpectrum3D_self(MoleculeSet* aSet, int depth, int depthMax, int kernelType, vector< pathsInMol3D >*  molPaths, vector<string>* atomLabels, vector<int>* bondTypes, bool silentMode){
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  vector< pathsInMol3D> molPaths_NEW;
  
  depth = depth + 1;

  if (depth == 0){ // depth = 0  --> initialization step
    // for each kind of atoms, find paths starting points and launch the recursion
    for(int i = 0 ;  i < atomLabels->size() ; i++){
      molPaths->clear();
      
      if( !silentMode ){
	Rcpp::Rcout <<  " \t finding paths starting from atoms labeled = " << (*atomLabels)[i] << endl;
      }
      
      init_path_pointers3D(aSet, molPaths, (*atomLabels)[i]);

      gramComputeSpectrum3D_self(aSet, depth, depthMax, kernelType, molPaths, atomLabels,  bondTypes, silentMode);

    } // end "for(int i = 0 ;  i < atomLabels->size() ; i++)"
    
  } // end "if(depth == 0)"

  else{  // if depth >= 1 :  we can try to extend the paths stored in path_pointers
	for(int j = 0 ; j < atomLabels->size() ; j++){
		for(int k = 0 ; k < bondTypes->size() ; k++){
	  		// for each atom type/bond type : extend the paths
	  		updatePaths3D(aSet, (*atomLabels)[j], (*bondTypes)[k], molPaths, &molPaths_NEW, kernelType, depth);

			int nonEmpty;
			nonEmpty = molPaths_NEW.size();
			if(nonEmpty > 0){  // i.e : if at least 1 vector of paths is non-empty we update the gram matrix or launch the recursion once again
				if(depth == depthMax){
					updateGram3D_self(aSet, &molPaths_NEW, kernelType);
					updateSelfKernel3D(aSet,&molPaths_NEW, kernelType);
					continue;
				}else{
					gramComputeSpectrum3D_self(aSet, depth, depthMax, kernelType, &molPaths_NEW, atomLabels, bondTypes, silentMode);
				}
			} // end "if(nonEmpty > 0)"

		}
	}
    
  } // end "else if(depth == 0)"

}



//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void gramComputeSpectrum3D_test(MoleculeSet* aSet, MoleculeSet* aSet2, int depth, int depthMax, int kernelType, vector< pathsInMol3D >*  molPaths, vector< pathsInMol3D >*  molPaths_test, vector<string>* atomLabels, vector<int>* bondTypes, bool silentMode){
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

  vector< pathsInMol3D > molPaths_NEW;
  vector< pathsInMol3D > molPaths_test_NEW;
  
  depth = depth + 1;


  if (depth == 0){ // depth = 0  --> initialization step
    // for each kind of atoms, find paths starting points and launch the recursion
    for(int i = 0 ;  i < atomLabels->size() ; i++){
      molPaths->clear();
      molPaths_test->clear();
      
      if( !silentMode ){
	Rcpp::Rcout << "\t - finding paths starting from atoms labeled = " << (*atomLabels)[i] << endl;
      }
      init_path_pointers3D(aSet, molPaths, (*atomLabels)[i]);
      init_path_pointers3D(aSet2, molPaths_test, (*atomLabels)[i]);
      
      gramComputeSpectrum3D_test(aSet, aSet2, depth, depthMax, kernelType, molPaths, molPaths_test, atomLabels, bondTypes, silentMode);

    } // end "for(int i = 0 ;  i < atomLabels->size() ; i++)"
  } // end "if (depth == 0)"

  else{  // if depth >= 1, we can try to extend the paths previously stored in path_pointers
	for(int j = 0 ; j < atomLabels->size() ; j++){
		for(int k = 0 ; k < bondTypes->size() ; k++){
			// for each atom type/bond type : extend the paths
			updatePaths3D(aSet, (*atomLabels)[j], (*bondTypes)[k], molPaths, &molPaths_NEW, kernelType, depth);
			updatePaths3D(aSet2, (*atomLabels)[j], (*bondTypes)[k], molPaths_test, &molPaths_test_NEW, kernelType, depth);
			
			int nonEmpty1, nonEmpty2;
			nonEmpty1 = molPaths_NEW.size();
    			nonEmpty2 = molPaths_test_NEW.size();
			if(nonEmpty1 > 0 || nonEmpty2 > 0){  // i.e : if at least 1 vector of paths is non-empty IN THE TWO SETS, we can try to extend the paths (otherwise, we do nothing)
				if(depth == depthMax){
					updateGram3D_test(aSet, aSet2, &molPaths_NEW, &molPaths_test_NEW, kernelType);
					updateSelfKernel3D(aSet,&molPaths_NEW, kernelType);
					updateSelfKernel3D(aSet2,&molPaths_test_NEW, kernelType);
					continue;
				}else{
					gramComputeSpectrum3D_test(aSet, aSet2, depth, depthMax, kernelType, &molPaths_NEW, &molPaths_test_NEW, atomLabels, bondTypes, silentMode);
				}

			} // end if(nonEmpty1 > 0 || nonEmpty2 > 0)

		}
	}

  } // end "else if(depth == 0)"

}




void init_path_pointers3D(MoleculeSet* aSet, vector< pathsInMol3D >* molPaths, string atomLabel){
// -----------------------------------------------------------------------------------------------
      // --> fill in "path_pointers" for starting atoms labeled as "atomLabel"

  path3D aPath;
  pathsInMol3D aPathInMol;
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
	    aPath.weight = 1;
	    aPath.endAtom = *anAtom;
	    aPath.startAtom = *anAtom;
	    aPathInMol.path_pointers.push_back(aPath);
	    aPathInMol.molName = (*aMol)->getName();
	    aPathInMol.molInd = ind;
	    first = false;
	} // end "if(first == true)"
	else{
	    aPath.weight = 1;
	    aPath.endAtom = *anAtom;
	    aPath.startAtom = *anAtom;
	    aPathInMol.path_pointers.push_back(aPath);
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
void updatePaths3D(MoleculeSet *aSet, string atomLabel, int bondType, vector< pathsInMol3D >* molPaths, vector< pathsInMol3D >* molPaths_NEW, int kernelType, int depth){
// ------------------------------------------------------------------------------------------------------------------------------------------------
  // --> extends path stored in "path_pointers" to atom of type "atomNumb". Stores the results in "path_pointers_NEW"

  path3D aPath;
  pathsInMol3D aPathInMol; // object initialized for each molecule present in the input molPaths, possibly inserted in the output molPaths_NEW
  
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
	     depth < 3
	     && ( aBond->second->getLabel() == bondType ) 
	     && ( aBond->first->getMorganLabel(true) == atomLabel )
	     )
	    ||
	    (
	     depth == 3  // (NB : implicit here : kernelType = 0, 1 or 2. kernelType parameter useless for the moment, but may be usefull later...)
	     && ( aBond->second->getLabel() == bondType ) 
	     && ( aBond->first->getMorganLabel(true) == atomLabel )
	     && ( aBond->first->getId() ==  (*molPaths)[i].path_pointers[j].startAtom->getId() ) 
	     )
	    )
	  {

	    if(first == true){
	      aPath.weight = 1;
	      aPath.endAtom = aBond->first;
	      aPath.previousAtom = (*molPaths)[i].path_pointers[j].endAtom;
	      aPath.startAtom = (*molPaths)[i].path_pointers[j].startAtom;
	      aPathInMol.path_pointers.push_back(aPath);
	      aPathInMol.molName = (*molPaths)[i].molName;
	      aPathInMol.molInd = (*molPaths)[i].molInd;
	      first = false;
	    } // end "if(first == true)"
	    else{
	      aPath.weight = 1;
	      aPath.endAtom = aBond->first;
	      aPath.previousAtom = (*molPaths)[i].path_pointers[j].endAtom;
	      aPath.startAtom = (*molPaths)[i].path_pointers[j].startAtom;
	      aPathInMol.path_pointers.push_back(aPath);
	    }

	  } // end "if (depth < 3 .....)"
      }
    }

    if(aPathInMol.path_pointers.size() > 0){
      (*molPaths_NEW).push_back(aPathInMol);
    }
  }
  
}


// ***********************************************************************************************************************
// -----------------------------------------------------------------------------------------------------------------------
void updateGram3D_self(MoleculeSet *aSet, vector< pathsInMol3D >*  molPaths, int kernelType){
// -----------------------------------------------------------------------------------------------------------------------

  double update;

  for(int i = 0 ; i < molPaths->size() ; i++){
    for(int j = i ; j < molPaths->size() ; j++){
      // 1 - compute the update factor depending on the kernel type
      switch(kernelType){
      case 0 : // 3-points spectrum
	update = (*molPaths)[i].path_pointers.size() * (*molPaths)[j].path_pointers.size();
	break;
      case 1 : // 3-points binary
	update = 1.0;
	break;
      case 2 : // 3-points Tanimoto
	update = 1.0;
	break;
      case 3 : // 2-points spectrum
	update = (*molPaths)[i].path_pointers.size() * (*molPaths)[j].path_pointers.size();
	break;
      case 4 : // 2-points binary
	update = 1.0;
	break;
      case 5 : // 2-points Tanimoto
	update = 1.0;
	break;
      }

      // 2 - add the update factor to the gram matrix
      (*aSet).addToGram( (*molPaths)[i].molInd , (*molPaths)[j].molInd, update);
      if( i!= j){
	(*aSet).addToGram( (*molPaths)[j].molInd , (*molPaths)[i].molInd, update);   
      }
      
   }
  }

}
      

//**************************************************************************
//--------------------------------------------------------------------------
void updateGram3D_test(MoleculeSet *aSet, MoleculeSet *aSet2, vector< pathsInMol3D >*  molPaths, vector< pathsInMol3D >*  molPaths_test, int kernelType){
//--------------------------------------------------------------------------
 
  double update;
 
  for(int i = 0 ; i < molPaths->size() ; i++){
      for(int j = 0 ; j < molPaths_test->size() ; j++){
	// 1 - compute the update factor depending on the kernel type
	switch(kernelType){
	case 0 : // 3-points spectrum
	  update = (*molPaths)[i].path_pointers.size() *  (*molPaths_test)[j].path_pointers.size();
	  break;
	case 1 : // 3-points binary
	  update = 1.0;
	  break;
	case 2 : // 3-points Tanimoto
	  update = 1.0;
	  break;
	case 3 : // 2-points spectrum
	  update = (*molPaths)[i].path_pointers.size() *  (*molPaths_test)[j].path_pointers.size();
	  break;
	case 4 : // 2-points binary
	  update = 1.0;
	  break;
	case 5 : // 2-points Tanimoto
	  update = 1.0;
	  break;
	}

	// 2 - add the update factor to the gram matrix
	(*aSet).addToGram( (*molPaths)[i].molInd , (*molPaths_test)[j].molInd, update);
      }
  }


}



// **************************************************************************
// --------------------------------------------------------------------------
void updateSelfKernel3D(MoleculeSet* aSet, vector< pathsInMol3D >* molPaths, int kernelType){
// --------------------------------------------------------------------------

  Molecule *aMol;  
  double update; 

 for(int i = 0 ; i < molPaths->size() ; i++){
    //aMol = (*aSet)[(*molPaths)[i].molName];
    //cout << "### MolInd = " << (*molPaths)[i].molInd << endl;
    aMol = (*aSet)[(*molPaths)[i].molInd];
    // 1 - compute the update facor depending on the kernel type
    switch(kernelType){
    case 0 : // 3-points spectrum
      update =  (*molPaths)[i].path_pointers.size() *  (*molPaths)[i].path_pointers.size();
      break;
    case 1 : // 3-points binary
      update = 1.0;
      break;
    case 2 : // 3-points Tanimoto
      update = 1.0;
      break;
    case 3 : // 2-points spectrum
      update =  (*molPaths)[i].path_pointers.size() *  (*molPaths)[i].path_pointers.size();
      break;
    case 4 : // 2-points binary
      update =  1.0;
      break;
    case 5 : // 2-points Tanimoto
      update = 1.0;
      break;
    }

    // 2 - add the update factor to the self kernel values
    aMol->addToSelfKernel(update);
 }


}


RCPP_MODULE(spectrum3Dhelper) {
	Rcpp::function( "gramSpectrum3D_self", &gramSpectrum3D_self );
	Rcpp::function( "gramSpectrum3D_test", &gramSpectrum3D_test );
}





