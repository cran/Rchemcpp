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



#ifndef SPECTRUM3DHELPER_H
#define SPECTRUM3DHELPER_H

#include <Rcpp.h>
#include <Rinternals.h>

#include "Rmoleculeset.h"

//---------
// external
//---------

void gramSpectrum3D_self(SEXP s, int depthMax, int kernelType, int nBins, double distMin, double distMax, bool silentMode);
void gramSpectrum3D_test(SEXP s, int depthMax, int kernelType, int nBins, double distMin, double distMax, bool silentMode);

//---------
// internal
//---------

typedef struct{
	Atom* endAtom;
        Atom* startAtom;
        Atom* previousAtom;
	double weight;
}path3D;

typedef struct{
  vector<path3D> path_pointers;
  string molName;
  int molInd;
}pathsInMol3D;


void gramComputeSpectrum3D_self( MoleculeSet*, int, int, int, vector< pathsInMol3D >*, vector<string>*, vector<int>*, bool silentMode = false);
void gramComputeSpectrum3D_test( MoleculeSet*, MoleculeSet*, int, int, int, vector< pathsInMol3D >*, vector< pathsInMol3D >*, vector<string>*, vector<int>*, bool silentMode = false);

void init_path_pointers3D( MoleculeSet*, vector< pathsInMol3D >*, string);
void updatePaths3D( MoleculeSet*, string, int, vector< pathsInMol3D >*, vector< pathsInMol3D >*, int, int);

void updateGram3D_self( MoleculeSet*, vector< pathsInMol3D >*, int);
void updateGram3D_test( MoleculeSet*, MoleculeSet*, vector< pathsInMol3D >*, vector< pathsInMol3D >*, int);
void updateSelfKernel3D(MoleculeSet* , vector< pathsInMol3D >*, int);

#endif


