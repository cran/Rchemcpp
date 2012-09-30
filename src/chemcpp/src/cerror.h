
/****************************************************************************************
					  cerror.h 
					------------
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


#ifndef CERROR_H
#define CERROR_H

#include <string>
#include <iostream>
using std::string;
using std::endl;
using std::cout;
using std::cerr;

// Error codes
#define BADFILE 1
#define MISSINGDESCRIPTOR 2
#define ATOMALREADYEXISTS 3
#define NOTENOUGHATOMSINMOLECULE 4
#define FILENOTFOUND 5
#define ERRORNA 6
#define UNKNOWNDATATYPE 7
#define ERRORATOMNOTFOUND 8
#define VALUENOTALLOWED 9
#define NOTIMPLEMENTED 10
#define NOTFOUND 11
#define COULDNOTOPENFILE 12
#define DUPLICATEENTRIES 13
#define DEPRECATED 14
#define THREADCREATIONERROR 15
#define EOFERROR 16
#define MISSINGDATA 17
#define NOTCALCULATED 18
#define BADVALUE 19
#define BONDALREADYEXISTS 20
#define IOERROR 21
#define NOSTRUCTURE 22
#define ATOMNOTFOUND 23
#define SSSRNOTDETECTED 24
#define MISSINGRING 25
#define ATOMNOTINRING 26
#define ALLNEIGHBOURSVISITED 27
#define BONDNOTFOUND 28
#define MOLECULENOTFOUND 29


/** CError class thrown when errors occur in chemcpp

		@author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
		@version 0.3
		@date 17 Jan 2004

		CLASS NAME: 	CError

		FOR:					SNSF SPONSORED PROJECT

		PURPOSE:			implements error message thrown by the classes
									in the chemcpp project.

		When an error occurs an instance of this class is created and thrown.
		The class has two attributes: a type (integer why) and a comment (string comment).
		These two attributes can be retrieved using the getType() and getComment() functions.
		An error message corresponding to the error can also be printed on cerr using the
		describe() function.

		The following error types are defined:

		- BADFILE 1
		- MISSINGDESCRIPTOR 2
		- ATOMALREADYEXISTS 3
		- NOTENOUGHATOMSINMOLECULE 4
		- FILENOTFOUND 5
		- ERRORNA 6
		- UNKNOWNDATATYPE 7
		- ERRORATOMNOTFOUND 8
		- VALUENOTALLOWED 9
		- NOTIMPLEMENTED 10
		- NOTFOUND 11
		- COULDNOTOPENFILE 12
		- DUPLICATEENTRIES 13
		- DEPRECATED 14
		- THREADCREATIONERROR 15
		- EOFERROR 16
		- MISSINGDATA 17
		- NOTCALCULATED 18
		- BADVALUE 19
		- BONDALREADYEXISTS 20
		- IOERROR 21
		- NOSTRUCTURE 22
		- ATOMNOTFOUND 23
		- SSSRNOTDETECTED 24
		- MISSINGRING 25
		- ATOMNOTINRING 26
		- ALLNEIGHBOURSVISITED 27
		- BONDNOTFOUND 28
		- MOLECULENOTFOUND 29



	  */

class CError {

/**
		\example cerror_example.cpp
*/

public:

	/** class constructor.
	*/
	CError( int anError, string aComment );

	/** class destructor.
	*/  
	~CError();

	/** returns the error type.
	*/
	int getType() { return why; }
	
	/** returns the comment associated with the error.
	*/
	string getComment() { return comment; }

	/** prints a description of the error to stderr.
	*/
	void describe();

private:
	/** error type.
	*/
	int why;
	
	/** error comment.
	*/
	string comment;

};

#endif
