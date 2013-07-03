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


#ifndef SUBTREEHELPER_H
#define SUBTREEHELPER_H

#include <Rcpp.h>
#include <Rinternals.h>

#include "Rmoleculeset.h"

#include <moleculeutils.h>

//--------
//external
//--------

void gramSubtree_self(SEXP s, double lambda, int depthMax, bool filterTotters, bool branchFlag, bool branchKernelUntilN, bool silentMode);

void gramSubtree_test(SEXP s, double lambda, int depthMax, bool filterTotters, bool branchFlag, bool branchKernelUntilN, bool silentMode);

void initialize_tuples(int nmax);

//--------
//internal
//--------

static vector< vector< vector< vector<int> > > > tuples;

class Nextatom {
	// a class to describe the list of neighbors of a given atoms that share the same bond and atom type
public:
	int b;  // the bond type
	string a; // the atom type
	vector<int> i; // the list of atom indices
};
// We also need a basic ordering for this class in order to sort and have fast matching
bool operator<(const Nextatom& a, const Nextatom& b) {
	if (a.a < b.a)
		return true;
	if (a.a > b.a)
		return false;
	return (a.b < b.b);
}

void initialize_extended(MoleculeSet* aSet, vector< vector< vector< Nextatom> > > * edges);

double subTreeKernel(	Molecule* mol1,
			Molecule* mol2,
			vector< vector< Nextatom> >* edges1,
			vector< vector< Nextatom> >* edges2,
			int depthMax,
			double lambda,
			bool noTotters,
			bool untilFlag,
			bool branchFlag);



#endif


