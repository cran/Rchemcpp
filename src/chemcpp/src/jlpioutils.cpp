/****************************************************************************************
					  jlpioutils.cpp 
					------------------
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



//#define DEBUG 1

#include "jlpioutils.h"
#include "cerror.h"

#include <sstream>

JLPIOUtils::JLPIOUtils(){
}
JLPIOUtils::~JLPIOUtils(){
}


string JLPIOUtils::vectorToString( vector<float>* a, string separator ){
	stringstream result;

	vector<float>::iterator i;

	bool first = true;
	for( i = a->begin(); i!= a->end(); i++ ){
		if( first ){
			result << (*i);
			first = false;
		}else{
			result << separator << (*i);
		}
	}
	return ( result.str() );
}

string JLPIOUtils::vectorToString( vector<double>* a, string separator ){
	stringstream result;

	vector<double>::iterator i;

	bool first = true;
	for( i = a->begin(); i!= a->end(); i++ ){
		if( first ){
			result << (*i);
			first = false;
		}else{
			result << separator << (*i);
		}
	}
	return ( result.str() );
}


/** list files in a directory */
void JLPIOUtils::readDirectory( const string directoryLocation, vector< string > *result, string extension, long start, long end ){

    #ifdef DEBUG
    cout << "JLPIOUtils::readDirectory " << start << " " << end << endl;
    #endif

	DIR* pDir = opendir ( directoryLocation.c_str() );

	if ( !pDir ){
//		return false;
	}

	dirent* pEntry;
	struct stat buf;

	long i = 0;
	while ( ( pEntry = readdir ( pDir ) ) ) {
		stat( pEntry->d_name, &buf );

		string aFileName = pEntry->d_name;
		#ifdef DEBUG
			cout << "found " << aFileName << endl;
		#endif

			//#ifdef SGI
			//if( S_ISDIR( buf.st_mode ) ){   //for SGI, if it is not a directory
			//#else
		     if( S_ISDIR( buf.st_mode ) ){   //for LINUX, if it is not a directory
			//#endif
		      #ifdef DEBUG
		        cout << " -> it is not a directory" << endl;
		      #endif
			if( aFileName != "." && aFileName != ".." ){

				if( extension == "" ){
					string sFound = directoryLocation + "/" + aFileName;
					result->push_back ( sFound );
				}else{

					if( aFileName.length() > extension.length()-1 ) {
						#ifdef DEBUG
							cout << "comparing " << aFileName.substr( aFileName.length()-3,3 ).c_str() << " and " << extension.c_str() << endl;
						#endif
						if( aFileName.substr( aFileName.length()-3,3 ) == extension ){
							#ifdef DEBUG
	  						cout << "ADDED " << directoryLocation + "/" + aFileName << endl;
							#endif
							string sFound = directoryLocation + "/" + aFileName;
							//cout << i << " " << (i >= start || start < 0 ) << " " << ( i <= end || end < 0 ) << endl;
							if( (i >= start || start < 0 ) && ( i <= end || end < 0 )){
							    result->push_back ( sFound );
							}
							i++;
						}

					}
				}
			}
		     #ifdef DEBUG
		     }else{
		     cout << " -> it is a directory" << endl; 
		     #endif
		  }

        }

	closedir( pDir );
}


