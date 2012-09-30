/****************************************************************************************
					  cerror.cpp 
					-------------
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



#include "cerror.h"

CError::CError( int anError, string aComment ) {
	why = anError;
	comment = aComment;
}

CError::~CError(){
}

void CError::describe(){
	switch(why){
		case BADFILE:
		  cerr << "BAD FILE error:"; break;
		case MISSINGDESCRIPTOR:
			cerr << "MISSING DESCRIPTOR error:"; break;
		case NOTENOUGHATOMSINMOLECULE:
			cerr << "NOT ENOUGH ATOMS IN MOLECULE error: "; break;
		case FILENOTFOUND:
			cerr << "FILE NOT FOUND error: "; break;
		case ERRORNA:
			cerr << "EMPTY VALUE ERROR: "; break;
		case UNKNOWNDATATYPE:
			cerr << "UNKNOWN DATA TYPE ERROR: "; break;
		case ERRORATOMNOTFOUND:
			cerr << "ATOM NOT FOUND: "; break;
		case VALUENOTALLOWED:
			cerr << "VALUE NOT ALLOWED: "; break;
		case NOTIMPLEMENTED:
			cerr << "FEATURE NOT YET IMPLEMENTED: "; break;
		case NOTFOUND:
			cerr << "NOT FOUND: "; break;
		case COULDNOTOPENFILE:
			cerr << "COULD NOT OPEN FILE: "; break;
		case DUPLICATEENTRIES:
			cerr << "DUPLICATE ENTRIES: "; break;
		case DEPRECATED:
			cerr << "DEPRECATED: "; break;
		case THREADCREATIONERROR:
			cerr << "THREADCREATIONERROR: "; break;
		case EOFERROR:
			cerr << "EOFERROR: "; break;
		case MISSINGDATA:
			cerr << "MISSING DATA: "; break;
		case NOTCALCULATED:
			cerr << "NOTCALCULATED: "; break;
		case BADVALUE:
			cerr << "BADVALUE: "; break;
		case BONDALREADYEXISTS:
			cerr << "BONDALREADYEXISTS: "; break;
		case IOERROR:
			cerr << "IOERROR: "; break;
		case NOSTRUCTURE:
			cerr << "NOSTRUCTURE: "; break;
		case ATOMNOTFOUND:
			cerr << "ATOMNOTFOUND: "; break;
		case SSSRNOTDETECTED:
			cerr << "SSSRNOTDETECTED: "; break;
		case MISSINGRING:
			cerr << "MISSINGRING: "; break;
		case ATOMNOTINRING:
			cerr << "ATOMNOTINRING: "; break;
		case ALLNEIGHBOURSVISITED:
			cerr << "ALLNEIGHBOURSVISITED: "; break;
		case BONDNOTFOUND:
			cerr << "BONDNOTFOUND: "; break;
		case MOLECULENOTFOUND:
			cerr << "MOLECULENOTFOUND: "; break;

		default:
			cerr << "ERROR type:" << getType(); break;
	}
	cerr << " " << getComment() << endl;
}
