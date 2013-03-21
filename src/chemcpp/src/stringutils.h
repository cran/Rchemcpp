/****************************************************************************************
					  stringutils.h 
					-----------------
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


#ifndef STRINGUTILS_H
#define STRINGUTILS_H


#include <string>
#include <vector>
#include <sstream>

#include <algorithm>

typedef unsigned int uint;

using namespace std;


/**Static functions to be used to process strings

	  @author Dr Jean-Luc Perret (luc@kuicr.kyoto-u.ac.jp), Kyoto University, Japan
		@version 0.3
		@date 17 Jan 2004

  */

class StringUtils {
public:
	/** splits a string into a vector<string> using a string delimiter.
	*/
	static int Split( const string input, const string delimiter, vector<string>& results );

	/** merges words in a container< string >, inserting separator.
	*/
	static string mergeWords(	vector<string>& words, const string separator );

	/** transforms a string into an int using a stringstream.
	*/
	static int toInt(	const string input );

  	/** transforms a string into a float using a stringstream.
	*/
	static float toFloat(	const string input );

	/** transforms a int into an string using a stringstream.
	*/
	static string toString(	const int input );

	/** transforms a float into an string using a stringstream.
	*/
	static string toString(	const float input );

	/** removes all spaces in a string.
		\todo StringUtils::rmSpace can be improved for more efficiency !!!
	*/
	static string rmSpace( const string input );

	/** removes all tail spaces in a string.
	*/
	static string rmTailSpace( const string input );

	/** "chomp" the string, i.e., removes the last character if it is 'end of line'.
	*/
	static string chomp(string inString);

	/** replaces slash by underscore in the string.
	*/
	static string slashToUnderscore( const string input );

	/** returns the content of a field starting at position start and of length maxlength in string aString.
	*/
	static string field(string aString, int start, int maxlength);

  	/** completes aString with aChar so that the total length is exactly equal to aLength. If aString is longer, it gets cropped. 
	*/
  	static string fill( string aString, int aLength, string aChar = " " );

  	/** prepends aChar to aString so that the total length is exactly equal to aLength. If aString is longer, it gets cropped.
	*/
  	static string preFill( string aString, int aLength, string aChar = " " );
  
	/** prepends aChar to anInt so that the total length is exactly equal to aLength. If aString is longer, it gets cropped.
	*/
	static string preFill( int anInt, int aLength, string aChar = " " );
	
	/** prepends aChar to aFloat so that the total length is exactly equal to aLength. If aString is longer, it gets cropped.
	*/
	static string preFill( float aFloat, int aLength, string aChar = " " );
 
	/** returns the directory portion of a full path.
	*/
  	static string getPath( string aLocation );
  	
	/** returns the extension of a file given a full path.
	*/
	static string getExtension( string aLocation );
	
	/** returns the path without extension.
	*/
	static string getNoExtension( string aLocation );
	
	/** returns the filename portion of a full path (including extension).
	*/	
	static string getFileName( string aLocation );
  	
	/** returns aLength characters at the right end of aString.
	*/
  	static string right( string aString, uint aLength );

	/** converts aString to Upper case.
	*/	
	static string toUpper( string aString );

	/** converts aString to Lower case.
	*/	
	static string toLower( string aString );

	/** returns the string starting at the 1st non space character of aString.
	*/  
  	static string getFirstNonSpace( string aString );

	//StringUtils();
	//~StringUtils();
};

#endif
