/****************************************************************************************
					  jlpioutils.h 
					----------------
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

#ifndef JLPIOUTILS_H
#define JLPIOUTILS_H

#include <iostream>
#include <string>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


using namespace std;

/**Provides utilities for io operations
  *@author Jean-Luc Perret
  */
class JLPIOUtils {
public: 
	/** class constructor.
	*/
	JLPIOUtils();

	/** class desctructor.
	*/
	~JLPIOUtils();

	/** list files in a directory.
  		usage:
			vector\<string\> fileVector;
			readDirectory ( "/home/anydir", fileVector);
			start and end allow to limit the number of file returned:
			limiting from start to end (included) starting count from 0
			the defauld (-1) indicates no limits.
	*/
  	static void readDirectory( const string directoryLocation, vector< string > *result, string extension = "", long start = -1, long end = -1 );

	/** concatenates the entries of a vector of float into a string, separated by a given separator.
	*/
	static string vectorToString(vector<float>* a, string separator = "\t");
	
	/** concatenates the entries of a vector of double into a string, separated by a given separator.
	*/
	static string vectorToString(vector<double>* a, string separator = "\t");
};

#endif
