/****************************************************************************************
					  stringutils.cpp 
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

#include "stringutils.h"
#include <iostream>

//StringUtils::StringUtils(){
//}
//StringUtils::~StringUtils(){
//}

int StringUtils::toInt(	const string input ){
	std::stringstream oss ( input );
	int rv;
	oss >> rv;
	return rv;
}

float StringUtils::toFloat(	const string input ){
	std::stringstream oss ( input );
	float rv;
	oss >> rv;
	return rv;
}

string StringUtils::toString(	const int input ){
	std::stringstream oss;
	string rv;
	oss << input;
	oss >> rv;
	return rv;
}
string StringUtils::toString(	const float input ){
	std::stringstream oss;
	string rv;
	oss << input;
	oss >> rv;
	return rv;
}


int StringUtils::Split(const string input, const string delimiter, vector<string>& results){
	if (delimiter.empty()) {
        results.push_back(input);
        return results.size();
    }
    string::const_iterator substart = input.begin(), subend;
    bool keep_empty = true;
    while (true) {
        subend = search(substart, input.end(), delimiter.begin(), delimiter.end());
        string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            results.push_back(temp);
        }
        if (subend == input.end()) {
            break;
        }
        substart = subend + delimiter.size();
    }
    return results.size();
}

string StringUtils::mergeWords(	vector<string>& words, const string separator ){
	stringstream result;
	vector< string >::iterator i;
	bool first = true;
	for( i = words.begin(); i != words.end(); i++ ){
		if( first == true ){
			result << (*i);
			first = false;
		}else{
			result << separator << (*i);
		}
	}
	return( result.str() );
}


string StringUtils::rmSpace( const string input ){
	stringstream output;
	string letter;
	for( uint i=0; i<input.length(); i++){
		letter = input[i];
		if( letter != " " ){
			output << letter;
		}
	}
	return( output.str() );
}

string StringUtils::slashToUnderscore( const string input ){
	stringstream output;
	string letter;
	for( uint i=0; i<input.length(); i++){
		letter = input[i];
		if( letter == "/" ){
			output << "_";
		}else{
			output << letter;
		}
	}
	return( output.str() );
}

string StringUtils::rmTailSpace( const string aninput ){

	string input = chomp(aninput);

	string output;
	string letter;
	int status = 0;
	//cout << "R:" <<input<<":R";
	for( int i=input.length(); i>=0; i--){
		//cout << input[i] << "#";
		letter = input[i];
		//cout << dec << letter << "#" << endl;
		if(status<2){
			if( letter != " "){
				output = letter + output;
				status++;  // finished removing spaces
				//cout << ".";
			}else{
				//cout << "#";
			}
		}else{
  		 output = letter + output;
			 //cout << ".";
		}
	}
	//cout << ":" << output << ";" << endl;
	return( output );
}

string StringUtils::chomp(string inString) {
	while ((inString.substr(inString.length(),1) == "\n") || (inString.substr(inString.length(),1) == "\r") ) {
		inString = inString.substr(1,inString.length()-1);
	}
	return inString;
}


string StringUtils::field( string aString, int start, int length ){
	int sl = aString.length();
	if(sl >= start + length){
		return(aString.substr(start, length));
	}else{
		if(start < sl){
			return(aString.substr(start, sl-start));
		}else{
			return("");
		}
	}
}

/** append aChar to aString so that the total length is exactly equal to aLength. If aString is longer, it gets cropped. */
string StringUtils::fill( string aString, int aLength, string aChar ){
	int sl = aString.length();

	if( sl == aLength){
		return( aString );
	}else if ( sl > aLength ){
		return( aString.substr( 0, aLength ) );
	}else{
		string cs = "";
		for( int i = 0; i < aLength - sl; i++ ){
			cs = cs + aChar;
		}
		return( aString + cs );
	}
}

/** prepend aChar to aString so that the total length is exactly equal to aLength. If aString is longer, it gets cropped. */
string StringUtils::preFill( string aString, int aLength, string aChar ){

	//cout << "preFill: " << aString << endl;

	int sl = aString.length();

	if( sl == aLength){
		return( aString );
	}else if ( sl > aLength ){
		return( aString.substr( 0, aLength ) );
	}else{
		string cs = "";
		for( int i = 0; i < aLength - sl; i++ ){
			cs = cs + aChar;
		}
		return( cs + aString );
	}
}

string StringUtils::preFill( int anInt, int aLength, string aChar ){
	return( StringUtils::preFill( StringUtils::toString( anInt ), aLength, aChar ) );
}

string StringUtils::preFill( float aFloat, int aLength, string aChar ){
	//cout << "Float prefill" << aFloat << " " << StringUtils::toString( aFloat ) << endl;
	return( StringUtils::preFill( StringUtils::toString( aFloat ), aLength, aChar ) );
}

/** returns the directory portion of a full path */
string StringUtils::getPath( string aLocation ){
	vector<string> words;
	StringUtils::Split( aLocation, "/", words );

	std::stringstream res;

	if( words.size() > 1){
		vector<string>::iterator itw = words.begin();
		for(; itw!= words.end()-1; itw++ ){
			res << (*itw) << "/";
		}
	}else{
		return( "./" );
	}
	return( res.str() );
}

/** returns the directory portion of a full path */
string StringUtils::getFileName( string aLocation ){

	//#ifdef DEBUG
		//cout << "StringUtils::getFileName: getting file name from " << aLocation << endl;
	//#endif


	vector<string> words;

	// replace all // by /

	// replace all // with /
	/*bool end = false;
	while( end == false ){
		StringUtils::Split( aLocation, "//", words );
		if( words.size() > 1 ){
			aLocation = mergeWords( words , "/" );
			#ifdef DEBUG
				cout << " -> " << aLocation << endl;
			#endif
			words.clear();
		}else{
			#ifdef DEBUG
				cout << "end" << endl;
			#endif
			end = true;
		}
	}
	words.clear();*/


	#ifdef DEBUG
		cout << "now splitting with delimiter / " << endl;
		//aLocation = "/scratch/fhfjdk";
		cout << aLocation << endl;
	#endif
	StringUtils::Split( aLocation, "/", words );
	#ifdef DEBUG
		cout << "splitting done: " << words.size() << " words: " << words[0] << endl;
	#endif


	std::stringstream res;

	if( words.size() > 1){
		vector<string>::iterator itw = words.end()-1;
		//res << (*itw) << "/";
		res << (*itw);
	}else{
		//cout << " ---> " << aLocation << endl;
		return( aLocation );
	}
	//cout << " ---> " << res.str() << endl;
	return( res.str() );
}

/** returns the directory portion of a full path */
string StringUtils::getExtension( string aLocation ){
	string aFileName = getFileName( aLocation );

	vector<string> words;
	StringUtils::Split( aLocation, ".", words );

	std::stringstream res;

	if( words.size() > 1){
		vector<string>::iterator itw = words.end()-1;
		res << (*itw);
	}else{
		return( "" );
	}
	return( res.str() );
}

/** returns the path without extension
*/
string StringUtils::getNoExtension( string aLocation ){
	string aFileName = getFileName( aLocation );

	#ifdef DEBUG
		cout << "StringUtils::getNoExtension: removing extension from " << aFileName << endl;
	#endif

	vector<string> words;
	StringUtils::Split( aLocation, ".", words );

	std::stringstream res;

	if( words.size() > 1){
		vector<string>::iterator itw;
		for( itw = words.begin(); itw != words.end()-1; itw++ ){
			res << (*itw) << ".";
		}
	}else{
		return( aLocation );
	}
	return( getPath( aLocation ) + res.str().substr( 0, res.str().size()-1 ) );
}


/** returns aLength characters at the right end of aString */
string StringUtils::right( string aString, uint aLength ){
	if( aString.length() >= aLength){
		return( aString.substr( aString.length() - aLength, aLength ) );
	}else{
		return( aString );
	}
}

#include <ctype.h>

string StringUtils::toUpper( string aString ){
	for( unsigned i=0; i<aString.size(); i++ ) {
    aString[i] = toupper( aString[i] );
  }
	return( aString );
}

string StringUtils::toLower( string aString ){
	for( unsigned i=0; i<aString.size(); i++ ) {
    aString[i] = tolower( aString[i] );
  }
	return( aString );
}



/*#include <locale>
#include <iterator>    // for back_inserter
#include <string>
#include <algorithm>
//#include <cctype>      // old <ctype.h>
#include <ctype.h>

struct ToUpper
   {
       ToUpper(std::locale const& l) : loc(l) {;}
       //char operator() (char c) const  { return std::toupper(c,loc); }
			 char operator() (char c) const  { return toupper(c); }
   private:
       std::locale const& loc;
   };

struct ToLower
   {
       ToLower(std::locale const& l) : loc(l) {;}
       //char operator() (char c) const  { return std::tolower(c,loc); }
       char operator() (char c) const  { return tolower(c); }
   private:
       std::locale const& loc;
   };

  */


string StringUtils::getFirstNonSpace( string aString )
{
	for( long i = 0; i < aString.size() ; i++ ){
		string c = aString.substr(i,1);
		if( c != " " ){
			return( c );
		}
	}
}
