/****************************************************************************************
					subtreehelper.cpp 
					-----------------
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



#include "subtreehelper.h"

//--------
//external
//--------

void gramSubtree_self(SEXP s, double lambda, int depthMax, bool filterTotters, bool branchFlag, bool branchKernelUntilN, bool silentMode){

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

	vector< vector< vector< Nextatom> > > edges; // the efficient indexation of the edges
	initialize_extended( aSet, &edges); 

	// 2 - compute the gram matrix
	// ---------------------------

	int depth = -1;

	if (!silentMode)
	{
		Rcpp::Rcout << "Subtree-kernel computation:" << endl;
		Rcpp::Rcout << "\t- depthMax = " << depthMax << endl;
		Rcpp::Rcout << "\t- lambda = " << lambda << endl;
		if(filterTotters)
			Rcpp::Rcout << "\t- no-totters" << endl;
		else
			Rcpp::Rcout << "\t- with-totters" << endl;	
		if(branchKernelUntilN)
			Rcpp::Rcout << "\t- until-N patterns" << endl;
		Rcpp::Rcout << endl;	
	}

	vector<Molecule*>::iterator mol1, mol2;
	int i = 0;
	int j = 0;
	double kernel;

	for( mol1 = aSet->begin(); mol1 != aSet->end(); mol1++ ){
		if( !silentMode)
			Rcpp::Rcout << "\t\t-molecule no " << i+1 << "/" << aSet->numMolecules() << endl;
		j = 0;
		for( mol2 = aSet->begin() + i ; mol2!= aSet->end() ; mol2++){
			kernel = subTreeKernel(*mol1, *mol2, &(edges[i]) , &(edges[i+j]) , depthMax, lambda, filterTotters, branchKernelUntilN, branchFlag);
			aSet->addToGram(i, i+j, kernel);
			if(j == 0){
				(*mol1)->addToSelfKernel(kernel);	
			}
			else{
				aSet->addToGram(i+j, i, kernel);
			}
			j++;
		}
		i++;
	}

}
  

void gramSubtree_test(SEXP s, double lambda, int depthMax, bool filterTotters, bool branchFlag, bool branchKernelUntilN, bool silentMode){

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

	vector< vector< vector< Nextatom> > > edges; // the efficient indexation of the edges
	initialize_extended( aSet, &edges); 

	vector< vector< vector< Nextatom> > > edges2; // the efficient indexation of the edges
	initialize_extended( aSet2, &edges2); 


	// 2 - compute the gram matrix
	// ---------------------------

	int depth = -1;

	if (!silentMode)
	{
		Rcpp::Rcout << "Subtree-kernel computation:" << endl;
		Rcpp::Rcout << "\t- depthMax = " << depthMax << endl;
		Rcpp::Rcout << "\t- lambda = " << lambda << endl;
		if(filterTotters)
			Rcpp::Rcout << "\t- no-totters" << endl;
		else
			Rcpp::Rcout << "\t- with-totters" << endl;	
		if(branchKernelUntilN)
			Rcpp::Rcout << "\t- until-N patterns" << endl;
		Rcpp::Rcout << endl;	
	}


	vector<Molecule*>::iterator mol1, mol2;
	int i = 0;
	int j = 0;
	double kernel;

	for( mol1 = aSet->begin(); mol1 != aSet->end(); mol1++ ){
		if( !silentMode)
			Rcpp::Rcout << "\t\t-molecule no " << i+1 << "/" << aSet->numMolecules() << endl;
		j = 0;
		for( mol2 = aSet2->begin(); mol2!= aSet2->end() ; mol2++){
			kernel = subTreeKernel(*mol1, *mol2, &(edges[i]), &(edges2[j]), depthMax, lambda, filterTotters, branchKernelUntilN, branchFlag);
			aSet->addToGram(i, j, kernel);
			// self kernel: only once!!
			if(i == 0){
				kernel = subTreeKernel(*mol2, *mol2, &(edges2[j]) , &(edges2[j]) , depthMax, lambda, filterTotters, branchKernelUntilN, branchFlag);		
				(*mol2)->addToSelfKernel(kernel);
			}
			j++;
		}
		// self kernel
		kernel = subTreeKernel(*mol1, *mol1, &(edges[i]) , &(edges[i]) , depthMax, lambda, filterTotters, branchKernelUntilN, branchFlag);		
		(*mol1)->addToSelfKernel(kernel);
		i++;
	}

}

//--------
//internal
//--------

void initialize_extended(MoleculeSet* aSet, vector< vector< vector< Nextatom> > > * edges)
{
  		vector<Molecule*>::iterator mol;
		map<Atom*, Bond*>::iterator aBond;
		
		// For each molecule in the set:
		int n=0;
		for( mol = aSet->begin(); mol != aSet->end(); mol++ ){
			edges->push_back(vector< vector< Nextatom> >());

			// Make a map between the atoms and their index
			map< Atom*, int> indAtoms;
			for(int i = 0 ; i < (*mol)->numAtoms() ; i++){
				indAtoms[ (*mol)->getAtomByIndex(i) ] = i;
			}
						
			// For each atom in the molecule:
			for (int i=0 ; i<(*mol)->numAtoms() ; i++){
				Atom *atom = (*mol)->getAtomByIndex(i);
				edges->at(n).push_back( vector< Nextatom>() );
				
				// For each bond leaving this atom
				for(aBond = atom->beginBond() ; aBond != atom->endBond() ; aBond++ ){
					Atom *atom2 = aBond->first;
					int isIn = 0;
					int bondtype = atom->getBondWithTarget(atom2)->getLabel();
					string atomtype = atom2->getMorganLabel();
					int atomindex = indAtoms[atom2];
					// Check whether a neighbor with the same bound/atom type already exists
					for (int k=0; k < edges->at(n)[i].size() ; k++){
						if ( (atomtype == edges->at(n)[i][k].a) && (bondtype == edges->at(n)[i][k].b )){
							// If a similar link is found, add the current atom to the corresponding list
							edges->at(n)[i][k].i.push_back(atomindex);
							isIn = 1;
						}
					}
					if (isIn == 0){
						// If no similar link is found, create a nex list and add the current atom to it
						edges->at(n)[i].push_back(Nextatom());
						int k=edges->at(n)[i].size()-1;
						edges->at(n)[i][k].b = bondtype;
						edges->at(n)[i][k].a = atomtype;
						edges->at(n)[i][k].i.push_back(atomindex);
					}
					
				} // All neighbours of *atom have been processed 
				
				// Sort the list of neighbor links for faster matching later
				sort(edges->at(n)[i].begin(),edges->at(n)[i].end());
			} // The current molecule is done.
			n++;
		}

}




double subTreeKernel(	Molecule* mol1,
						Molecule* mol2,
						vector< vector< Nextatom> >* edges1,
						vector< vector< Nextatom> >* edges2,
						int depthMax,
						double lambda,
						bool filterTotters,
						bool untilFlag,
						bool branchFlag){
	// ------------------------------------------------------------------
	
	
	// 1 - create matrices of k_1(r,s)...k_h(r,s)
	// --> subMatrices[0],...,subMatrices[h-1]
	vector<  vector< vector<double> > > subMatrices;
	
	for(int i = 0; i < depthMax ; i++){
		subMatrices.push_back( vector< vector<double> >() );
		for(int j = 0 ; j < mol1->numAtoms() ; j++){
			subMatrices[i].push_back( vector<double>() );
			for (int k = 0 ; k < mol2->numAtoms() ; k++){
				subMatrices[i][j].push_back(0.0);
			}
		}
	}
	
	
	// 2 - initialize subMatrices[0] with binary kernels based on Morgan Labels
	// NB: MorganLabelKernel can be changed to a different kernel (--> SEE SD2GRAM, OPTION -K)
	vector<int> toupdate1,toupdate2;
	for(int i = 0 ; i < mol1->numAtoms() ; i++){
		for(int j = 0 ; j < mol2->numAtoms() ; j++){
			if(branchFlag){
				subMatrices[0][i][j] = MoleculeUtils::atomKernelMorganLabel( mol1->getAtomByIndex(i) , mol2->getAtomByIndex(j) );
			}
			else{
				subMatrices[0][i][j] = lambda * MoleculeUtils::atomKernelMorganLabel( mol1->getAtomByIndex(i) , mol2->getAtomByIndex(j) );
			}
			if (subMatrices[0][i][j] != 0.0 ) {
				toupdate1.push_back(i);
				toupdate2.push_back(j);
			}
		}
	}
	

	
	// 3 - compute subMatrices[1]...subMatrices[depthMax]
	double update, temp, ker;
	
	for(int depth = 1 ; depth < depthMax ; depth++){
				
		// FOR EACH PAIR OF ATOMS THAT MATCH:
		for (int istep=0 ; istep < toupdate1.size() ; istep++) {
			int atom1 = toupdate1[istep];
			int atom2 = toupdate2[istep];
			vector<double> tempKernel;
			
			// Detect common neighbor labels
			int n1=0,n2=0;
			while ( (n1 < (*edges1)[atom1].size()) && (n2 < (*edges2)[atom2].size() ) ) {
				
				if ((*edges1)[atom1][n1] < (*edges2)[atom2][n2]) {
					n1++;
					continue;
				}
				if ((*edges2)[atom2][n2] < (*edges1)[atom1][n1]) {
					n2++;
					continue;
				}
				
				// A matching pair has been found
				int natoms1 = (*edges1)[atom1][n1].i.size();
				int natoms2 = (*edges2)[atom2][n2].i.size();
				int maxCard = min( natoms1, natoms2); 
				temp = 0.0;
				for(int d = 1; d <= maxCard; d++){
					for(int l = 0 ; l < tuples[natoms1-1][d-1].size() ; l++){
						for(int m = 0 ; m < tuples[natoms2-1][d-1].size() ; m++){
							update = 1.0;
							for (int toto = 0 ; toto<d ; toto++) {

								update = update * subMatrices[depth-1][ (*edges1)[atom1][n1].i[tuples[natoms1-1][d-1][l][toto] ] ][ (*edges2)[atom2][n2].i[tuples[natoms2-1][d-1][m][toto] ] ];
							}
							
							if(branchFlag){
								temp += update * pow(lambda,d-1);
							}
							else{
								temp += update;
							}
						}//loop m
					} //loop l
				}//loop d
				tempKernel.push_back(temp);
				
				n1++;
				n2++;
			} // end of while( (n1<... loop
			
			// Compute kernel update from tempkernel
			if (tempKernel.size() == 0)
				// No common neighbour has been found
				ker = 0.0;
			else {
				ker = tempKernel[0] + 1;
				for (int n = 1 ; n < tempKernel.size() ; n++){
					if(branchFlag){
						ker = ker * (lambda*tempKernel[n] + 1) + tempKernel[n]*(1 - lambda);
					}
					else{
						ker = ker * (tempKernel[n] + 1);
					}
				}
				ker = ker - 1;
			}
			
			// update DP matrix
			subMatrices[depth][atom1][atom2] = ker;
			if(untilFlag){
				subMatrices[depth][atom1][atom2] = 1.0 + subMatrices[depth][atom1][atom2];
			}
			if(!branchFlag){
				subMatrices[depth][atom1][atom2] = lambda * subMatrices[depth][atom1][atom2];
			}
			
		} //loop i
		
	} //loop depth
				
		
	
	// 4 - finally : sum subMatrices[depthMax-1] to get the kernel
	double kernel = 0.0;
	for(int i = 0 ; i < mol1->numAtoms() ; i++){
		for(int j = 0 ; j < mol2->numAtoms() ; j++){
			if( !(filterTotters && (mol1->getAtomByIndex(i)->getKashimaPS() * mol2->getAtomByIndex(j)->getKashimaPS() == 0.0) ) ){ 
				kernel += subMatrices[depthMax-1][i][j];
			}
		}
	}
	

	if(!branchFlag){
		kernel = kernel / pow(lambda,depthMax);
	}
	
	
	/*
	 // DEBUG:
	 // CHECK CORRECTNESS --> print DP matrices on screen
	 for(int k = 0 ; k < depthMax ; k++){
		 Rcpp::Rcout << "####--- subMatrices[" << k << "] ---####"<< endl;
		 Rcpp::Rcout << "corner";
		 for(int i = 0 ; i < mol2->numAtoms() ; i++){
			 if( mol2->getAtomByIndex(i)->getKashimaPS() != 0.0 ){
				 Rcpp::Rcout << "\tno" << i+1;
			 }
		 }
		 Rcpp::Rcout << endl;
		 
		 for(int i = 0 ; i < mol1->numAtoms() ; i++){
			 if( mol1->getAtomByIndex(i)->getKashimaPS() != 0.0 ){
				 Rcpp::Rcout << "no" << i+1;
				 for(int j = 0 ; j < mol2->numAtoms() ; j++){
					 if(mol2->getAtomByIndex(j)->getKashimaPS() != 0.0 ){
						 Rcpp::Rcout << "\t" << subMatrices[k][i][j];
					 }
				 }
				 Rcpp::Rcout << endl;
			 }
		 }
		 Rcpp::Rcout << "#####################" << endl ;
	 }
	 Rcpp::Rcout << endl << endl;
	 */
	
	
	return kernel;
	
}



// Initialize the static variable the contains the list of all p-tuples among n indices for various p and n<nmax
// -> Rchemcpp: Deleted MoleculeSetExtended; converted Member-function to this
void initialize_tuples( int nmax){

	tuples.clear(); //MAHRMIC: Remove all of them first - otherwise the algorithms segfault

	for (int n=0; n<nmax; n++) {
		// (n+1) is the number of points among which we choose
		tuples.push_back( vector< vector< vector<int> > >() );
		
		for (int p=0; p<=n ; p++) {
			// (p+1) is the number of points that we choose
			tuples[n].push_back( vector< vector<int> >() );
			
			if (p==0) {
				// If we choose only one point, then we directly list all 1-tuples (the n singletons)
				for (int i=0; i<=n ; i++) {
					tuples[n][p].push_back( vector<int>() );
					tuples[n][p][i].push_back(i);
				}
			}
			else {
				// If we choose (p+1)>1 points, we deduce the (p+1)-tuples from the p-uples
				for (int i=0 ; i<tuples[n][p-1].size() ; i++) {
					
					for (int newi=0 ; newi<=n ; newi++) {
						// newi is a candidate point to add to the i-th p-uple.
						// Check whether it is already in the p-uple
						int isIn = 0;
						for (int j=0 ; j<tuples[n][p-1][i].size() ; j++) {
							if (tuples[n][p-1][i][j] == newi)
								isIn = 1;
						}
						// If newi is not already in the i-th p-uple, then create a (p+1)-tuple by adding it
						if (isIn == 0) {
							tuples[n][p].push_back(tuples[n][p-1][i]);
							tuples[n][p][tuples[n][p].size()-1].push_back(newi);
						}
					}
				}
			}
		}
	}
}





RCPP_MODULE(subtreehelper) {
	Rcpp::function( "gramSubtree_self", &gramSubtree_self,"kalsfj" );
	Rcpp::function( "gramSubtree_test", &gramSubtree_test,"kkleuh" );
	Rcpp::function( "initialize_tuples", &initialize_tuples, "oashf" );
}





