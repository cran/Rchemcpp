/****************************************************************************************
					spectrumhelper.h 
					----------------
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



#ifndef SPECTRUMHELPER_H
#define SPECTRUMHELPER_H

#include <Rcpp.h>
#include <Rinternals.h>

#include "Rmoleculeset.h"

//---------
// external
//---------

void gramSpectrum_self(SEXP s, int depthMax, int kernelType, double kernelParam, bool onlyN, bool silentMode);
void gramSpectrum_test(SEXP s, int depthMax, int kernelType, double kernelParam, bool onlyN, bool silentMode);

//---------
// internal
//---------

typedef struct{
	Atom* endAtom;
        Atom* previousAtom;
	double weight;
}path;

typedef struct{
  vector<path> path_pointers;
  string molName;
  int molInd;
}pathsInMol;

void gramComputeSpectrum_self( MoleculeSet*, int, int, int, double, vector< pathsInMol >*, vector<string>*, vector<int>*, bool, bool silentMode = false);
void gramComputeSpectrum_test( MoleculeSet*, MoleculeSet*, int, int, int, double, vector< pathsInMol >*, vector< pathsInMol >*, vector<string>*, vector<int>*, bool, bool silentMode = false);

void init_path_pointers( MoleculeSet*, vector< pathsInMol >*, string, int);
void updatePaths( MoleculeSet*, string, int, vector< pathsInMol >*, vector< pathsInMol >*, int, int);

void updateGram_self( MoleculeSet*, vector< pathsInMol >*, int, double, int);
void updateGram_test( MoleculeSet*, MoleculeSet*, vector< pathsInMol >*, vector< pathsInMol >*, int, double, int);
void updateSelfKernel(MoleculeSet* , vector< pathsInMol >*, int, double, int);



#endif


