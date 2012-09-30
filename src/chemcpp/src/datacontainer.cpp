
/****************************************************************************************
					  datacontainer.cpp
					--------------------
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

#include "datacontainer.h"
#include "stringutils.h"


//#define DEBUG 1

DataContainer::DataContainer(){
	#ifdef DEBUG
		cout << "CREATING KIND DESCRIPTOR MAPS" << endl;
	#endif
	kindStringDescriptors = new map< string, Descriptor< string >* >;
	kindIntDescriptors = new map< string, Descriptor< int >* >;
	kindFloatDescriptors = new map< string, Descriptor< float >* >;
  #ifdef DEBUG
		cout << "string map has " << kindStringDescriptors->size() << " elements"  << endl;
	#endif

	//cout << "ELEMENT ####" << endl;

	flagElement = true;  // this datacontainer is a new datacontainer. Useful for destruction

	//cout << "Int map at :" << kindIntDescriptors << endl;

}
DataContainer::~DataContainer(){
	#ifdef DEBUG
		cout << "DataContainer::~DataContainer()" << endl;
	#endif
	deleteAllDescriptors();
	#ifdef DEBUG
	  cout << "shall I delete kind descriptors? ";
	#endif
	if( flagElement == true ){
		#ifdef DEBUG
		  cout << "yes" << endl;
		#endif
		deleteAllKindDescriptors();
		
	}
  	#ifdef DEBUG
		else{
			cout << "no" << endl;
		}
	#endif
	#ifdef DEBUG
		cout << "DataContainer::~DataContainer() DONE" << endl;
	#endif
}

DataContainer::DataContainer( DataContainer& aDataContainer ){
	#ifdef DEBUG
		cout << "DataContainer::DataContainer( DataContainer& aDataContainer )" << endl;
	#endif
	//cout << "DataContainer copy constructor " << endl;
	// set kindDescritors pointer to the same as the source

  	kindStringDescriptors = aDataContainer.kindStringDescriptors;
	kindIntDescriptors = aDataContainer.kindIntDescriptors;
	kindFloatDescriptors = aDataContainer.kindFloatDescriptors;

	//cout << kindStringDescriptors << endl;

	// copy all non kind descriptors
	#ifdef DEBUG
		cout << "copying all non kind descriptors" << endl;
	#endif

	map<const string, Descriptor<int>* >::iterator iti;
	int i = 0;
	//cout << "  int" << endl;
	//cout << (aDataContainer).intDescriptors.size() << endl;
	for( iti = (aDataContainer).intDescriptors.begin(); iti != (aDataContainer).intDescriptors.end(); iti++ ){
			//cout << i << endl;
			i++;

			int aValue = 0;
			bool empty = false;

			try{
				aValue = (*iti).second->getValue(true);
			}catch( CError e ){
				if( e.getType() == ERRORNA ){
					empty = true;
					continue;
				}else{
					e.describe();
					throw( e );
				}
			}
			Descriptor<int>* toto =	addIntDescriptor( (*iti).second->getLabel() , aValue, (*iti).second->getUnit(), (*iti).second->getComment() );
			if( empty == true ){
				  toto->setEmpty();
			}

	}
	#ifdef DEBUG
		cout << "  int DONE" << endl;
	#endif

	#ifdef DEBUG
		cout << "  float" << endl;
		cout << (aDataContainer).floatDescriptors.size() << endl;
	#endif
	map<const string, Descriptor<float>* >::iterator itf;
	for( itf = (aDataContainer).floatDescriptors.begin(); itf != (aDataContainer).floatDescriptors.end(); itf++ ){
			//cout << i << endl;
			i++;
			//cout << (*itf).second->getLabel() << endl;
			float aValue = 0.0;
			bool empty = false;

			try{			
				aValue = (*itf).second->getValue(true);
			}catch( CError e ){
				if( e.getType() == ERRORNA ){
					empty = true;
					continue;
				}else{
					e.describe();
					throw( e );
				}
			}
			Descriptor<float>* toto =	addFloatDescriptor( (*itf).second->getLabel() , aValue, (*itf).second->getUnit(), (*itf).second->getComment() );
			if( empty == true ){
				  toto->setEmpty();
			}
	}
	#ifdef DEBUG
		cout << "float done" << endl;
	#endif
	
	map<const string, Descriptor<string>* >::iterator its;
	for( its = aDataContainer.stringDescriptors.begin(); its != aDataContainer.stringDescriptors.end(); its++ ){
			//cout << i << endl;
			i++;
			string aValue = "";
			bool empty = false;

			try{
				aValue = (*its).second->getValue(true);
			}catch( CError e ){
				if( e.getType() == ERRORNA ){
					empty = true;
					continue;
				}else{
					e.describe();
					throw( e );
				}
			}
			Descriptor<string>* toto =	addStringDescriptor( (*its).second->getLabel() , aValue, (*its).second->getUnit(), (*its).second->getComment() );
			if( empty == true ){
				  toto->setEmpty();
			}
	}
	//cout << "copied " << i << " non kind descriptors" << endl;

	#ifdef DEBUG
		cout << "string done" << endl;
	#endif
  

	#ifdef DEBUG
		cout << "copying all non kind descriptors DONE" << endl;
	#endif


	flagElement = false; // this datacontainer is a copy of an existing datacontainer. Useful for destruction
  #ifdef DEBUG
		cout << "DataContainer::DataContainer( DataContainer& aDataContainer ) DONE" << endl;
	#endif
}

/*DataContainer& DataContainer::operator= (const DataContainer& aDataContainer){
	// copy kind descriptors

	// first string descriptors
	map<const string, Descriptor<string>* >::iterator its;
//	for( its = (aDataContainer).kindStringDescriptors.begin(); its != (aDataContainer).kindStringDescriptors.end(); its++ ){
//  	addKindIntDescriptor( (*its).second -> getLabel(), (*its).second );
//	}

	// then integer descriptors

	// then float descriptors

}*/



/** Adds an integer kind descriptor with a label, value, unit and comment */
Descriptor< int >*  DataContainer::addKindIntDescriptor(string aLabel, int aValue, string aUnit, string aComment){
	
	//cout << "setting " << (*kindIntDescriptors)[aLabel] << " to " << atoi(aValue.c_str()) << endl;
	//(*kindIntDescriptors)[aLabel] = new Descriptor<int>( aLabel, atoi(aValue.c_str()), aUnit, aComment );
	//cout << "map now contains " << endl;
	//(*kindIntDescriptors)[aLabel]->describe();
	(*kindIntDescriptors)[aLabel] = new Descriptor<int>( aLabel, aValue, aUnit, aComment );
	return( (*kindIntDescriptors)[aLabel] );
}
Descriptor< int >*  DataContainer::addKindIntDescriptor(Descriptor< int >* aDescriptor){
	(*kindIntDescriptors)[aDescriptor->getLabel()] = aDescriptor;
	return( aDescriptor );
}

/** Adds a float kind descriptor with a label, value, unit and comment */
Descriptor< float >*  DataContainer::addKindFloatDescriptor(string aLabel, float aValue, string aUnit, string aComment){
	(*kindFloatDescriptors)[aLabel] = new Descriptor<float>( aLabel, aValue, aUnit, aComment );
	return( (*kindFloatDescriptors)[aLabel] );
}
Descriptor< float >*  DataContainer::addKindFloatDescriptor(Descriptor< float >* aDescriptor){
	(*kindFloatDescriptors)[aDescriptor->getLabel()] = aDescriptor;
	return( aDescriptor );
}

/** Adds a string kind descriptor with a label, value, unit and comment */
Descriptor< string >* DataContainer::addKindStringDescriptor(string aLabel, string aValue, string aUnit, string aComment){
	(*kindStringDescriptors)[aLabel] = new Descriptor<string>( aLabel, aValue, aUnit, aComment );
	return( (*kindStringDescriptors)[aLabel] );
}
Descriptor< string >*  DataContainer::addKindStringDescriptor(Descriptor< string >* aDescriptor){
	(*kindStringDescriptors)[aDescriptor->getLabel()] = aDescriptor;
	return( aDescriptor );
}




/** Adds an integer descriptor with a label, value, unit and comment */
Descriptor< int >*  DataContainer::addIntDescriptor(string aLabel, int aValue, string aUnit, string aComment){
	Descriptor< int >* d = new Descriptor<int>( aLabel, aValue, aUnit, aComment );
	intDescriptors[aLabel] = d;
	return( d );
}
/** Adds an float descriptor with a label, value, unit and comment */
Descriptor< float >*  DataContainer::addFloatDescriptor(string aLabel, float aValue, string aUnit, string aComment){
	Descriptor< float >* d = new Descriptor<float>( aLabel, aValue, aUnit, aComment );
	floatDescriptors[aLabel] = d;
	return( d );
}
/** Adds an integer descriptor with a label, value, unit and comment */
Descriptor< string >*  DataContainer::addStringDescriptor(string aLabel, string aValue, string aUnit, string aComment){
	//cout << "DataContainer::addStringDescriptor"<<endl;
	//cout << " adding string " << aLabel << ", " << aValue << ", " << aUnit << ", " << aComment << endl;
	Descriptor< string >* d = new Descriptor<string>( aLabel, aValue, aUnit, aComment );
	stringDescriptors[aLabel] = d;
	//cout << "DataContainer::addStringDescriptor DONE"<<endl;
	return( d );
}

/** sets value aValue to an existing descriptor with label aLabel */
Descriptor< int >* DataContainer::setIntDescriptor(string aLabel, int aValue, string aUnit, string aComment, bool addIfMissing, bool silentError ){
	// check if string descriptor exists. If not create it
	Descriptor< int >* d;

	if( hasIntDescriptor( aLabel ) ){
		d = intDescriptors[aLabel];
		d -> setValue(aValue);
	}else{
		if( addIfMissing ){
			d = addIntDescriptor( aLabel, aValue, aUnit, aComment );
		}else{

			CError e( MISSINGDESCRIPTOR, "DataContainer::setIntDescriptor: no descriptor " + aLabel );
			if( !silentError ){
				e.describe();
			}
			throw(e);
		}
	}
	return( d );
}
Descriptor< float >* DataContainer::setFloatDescriptor(string aLabel, float aValue, string aUnit, string aComment, bool addIfMissing, bool silentError ){
	// check if string descriptor exists. If not create it
	Descriptor< float >* d;

	if( hasFloatDescriptor( aLabel ) ){
		d = floatDescriptors[aLabel];
		//cout << "found float descriptor " << d->toString() << endl;
		d -> setValue( aValue );
	}else{
		if( addIfMissing ){
			d = addFloatDescriptor( aLabel, aValue, aUnit, aComment );
		}else{

			CError e( MISSINGDESCRIPTOR, "DataContainer::setFloatDescriptor: no descriptor " + aLabel );
			if( !silentError ){
				e.describe();
			}
			throw(e);
		}
	}
	return ( d );
}
Descriptor< string >* DataContainer::setStringDescriptor(string aLabel, string aValue, string aUnit, string aComment, bool addIfMissing, bool silentError ){
	// check if string descriptor exists. If not create it
	Descriptor< string >* d;
	if( hasStringDescriptor( aLabel ) ){
		d = stringDescriptors[aLabel];
		d -> setValue(aValue);
	}else{
		if( addIfMissing ){
			d = addStringDescriptor( aLabel, aValue, aUnit, aComment );
		}else{

			CError e( MISSINGDESCRIPTOR, "DataContainer::setStringDescriptor: no descriptor " + aLabel );
			if( !silentError ){
				e.describe();
	    }
			throw(e);
		}
	}
	return( d );
}

/** write a only descriptors of the data container to cout (not Kind descriptors) */
void DataContainer::describeShort() throw( CError ) {
	map<const string, Descriptor<string>* >::iterator its;
  map<const string, Descriptor<int>* >::iterator iti;
  map<const string, Descriptor<float>* >::iterator itf;

	// first write string descriptors
	for( its = stringDescriptors.begin(); its != (stringDescriptors).end(); its++ ){
		(*its).second->describeShort();
	}
	// then write integer descriptors
	for( iti = intDescriptors.begin(); iti != (intDescriptors).end(); iti++ ){
		(*iti).second->describeShort();
	}
	// finally write float descriptors
  for( itf = floatDescriptors.begin(); itf != (floatDescriptors).end(); itf++ ){
		(*itf).second->describeShort();
	}
  cout << "-------------------------" << endl; 	
}

/** write a description of the data container to cout */
void DataContainer::describe() throw( CError ) {

	string valueS;

	cout << "KIND DESCRIPTORS:" << endl;
	cout << "string kindDescriptors:" << endl;
	// first write string kindDescriptors
	map<const string, Descriptor<string>* >::iterator its;
	for( its = (*kindStringDescriptors).begin(); its != (*kindStringDescriptors).end(); its++ ){
	  (*its).second->describe();
	}

	cout << "int kindDescriptors:" << endl;
	// then write integer kindDescriptors 
  map<const string, Descriptor<int>* >::iterator iti;
	for( iti = (*kindIntDescriptors).begin(); iti != (*kindIntDescriptors).end(); iti++ ){
		(*iti).second->describe();
	}
	cout << "float kindDescriptors" << endl;
	// finally write float kindDescriptors
	map<const string, Descriptor<float>* >::iterator itf;
  for( itf = (*kindFloatDescriptors).begin(); itf != (*kindFloatDescriptors).end(); itf++ ){
		(*itf).second->describe();
	}

	cout << "DESCRIPTORS:" << endl;
	// first write string descriptors
	for( its = stringDescriptors.begin(); its != (stringDescriptors).end(); its++ ){
		(*its).second->describe();
	}
	// then write integer descriptors
	for( iti = intDescriptors.begin(); iti != (intDescriptors).end(); iti++ ){
		(*iti).second->describe();
	}
	// finally write float descriptors
  for( itf = floatDescriptors.begin(); itf != (floatDescriptors).end(); itf++ ){
		(*itf).second->describe();
	}
  cout << "-------------------------" << endl;

}

Descriptor< string >* DataContainer::getStringDescriptor( string aLabel, bool silentError ) throw( CError ){

	/*bool found;

	// look among kind descriptors first
  map<const string, Descriptor<string>* >::iterator its;
	for( its = (*kindStringDescriptors).begin(); its != (*kindStringDescriptors).end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	}

	// look among descriptors defined for this instance
	for( its = stringDescriptors.begin(); its != stringDescriptors.end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	}

  // if no corresponding descriptors found, throw exception
	CError e( MISSINGDESCRIPTOR, "no descriptor " + aLabel );
	e.describe();
	throw(e);

	*/

	//cout << "looking among " << kindStringDescriptors->size() << " string kind descriptors" << endl;
	map< string, Descriptor<string>* >::iterator answer = kindStringDescriptors->find(aLabel);
	//cout << "and?" << endl;
	if( answer != kindStringDescriptors->end() ){
		//cout << "found" << endl;
		return(answer->second);
	}else{
		//cout << "looking in non kind descriptors" << endl;
		answer = stringDescriptors.find(aLabel);
	  if( answer != stringDescriptors.end() ){
			//cout << "found" << endl;
			return(answer->second);
		}else{
			//cout << "not found" << endl;
			// if no corresponding descriptors found, throw exception
			CError e( MISSINGDESCRIPTOR, "no descriptor " + aLabel );
			if( silentError == false ){
				e.describe();
			}
			throw(e);
		}
	}

}


long DataContainer::getPossibleValuesInIntDescriptor( string aDescriptorName, vector< int >* p ){

	long res = 0;

	if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),8), 8 ) == ".integer"){
		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 8 );
	}

	map< const string, Descriptor<int>* >::iterator iti;
	int i = 0;
	for( iti = intDescriptors.begin(); iti != intDescriptors.end(); iti++ ){
		if( (*iti).second->getLabel() == aDescriptorName ){
			// check if this value was already found, if not add it to the vector of possible values
			int newValue = (*iti).second->getValue();
			vector<int>::iterator vi;
			long nfound = 0;
			for( vi = p->begin(); vi != p->end(); vi++){
				if( (*vi) == newValue ){
					nfound += 1;
					break;
				}
			}
			if( nfound == 0 ){
				p->push_back( newValue );
				res += 1;
			}

		}
	}
	return( res );
}


Descriptor< int >* DataContainer::getIntDescriptor( string aLabel, bool silentError ) throw( CError ){

/*	bool found;

	// look among kind descriptors first
  map<const string, Descriptor<int>* >::iterator its;
	for( its = (*kindIntDescriptors).begin(); its != (*kindIntDescriptors).end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	}

	// look among descriptors defined for this instance
	for( its = intDescriptors.begin(); its != intDescriptors.end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	}

  // if no corresponding descriptors found, throw exception
	CError e( MISSINGDESCRIPTOR, "no descriptor " + aLabel );
	e.describe();
	throw(e);

	*/
	map< string, Descriptor<int>* >::iterator answer = kindIntDescriptors->find(aLabel);
	if( answer != kindIntDescriptors->end() ){
		return(answer->second);
	}else{
		answer = intDescriptors.find(aLabel);
	  if( answer != intDescriptors.end() ){
			return(answer->second);
		}else{
			// if no corresponding descriptors found, throw exception
			CError e( MISSINGDESCRIPTOR, "no descriptor " + aLabel );
			if( silentError == false ){
				e.describe();
			}
			throw(e);
		}
	}

}

Descriptor< float >* DataContainer::getFloatDescriptor( string aLabel, bool silentError ) throw( CError ){

	//bool found;
	map< string, Descriptor<float>* >::iterator answer;

	// look among kind descriptors first
  /*map<const string, Descriptor<float>* >::iterator its;
	for( its = (*kindFloatDescriptors).begin(); its != (*kindFloatDescriptors).end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	}

	// look among descriptors defined for this instance
	for( its = floatDescriptors.begin(); its != floatDescriptors.end(); its++ ){
   	if( (*its).first == aLabel ){
			found = true;
			return ( (*its).second );
		}
	} */
	answer = kindFloatDescriptors->find(aLabel);
	if( answer != kindFloatDescriptors->end() ){
		return(answer->second);
	}else{
		answer = floatDescriptors.find(aLabel);
	if( answer != floatDescriptors.end() ){
			return(answer->second);
		}else{
			// if no corresponding descriptors found, throw exception
			CError e( MISSINGDESCRIPTOR, "no descriptor " + aLabel );
			if( silentError == false ){
				e.describe();
			}
			throw(e);
		}
	}

}


/** deletes the int descriptor aString from intDescriptors (deletes the object and removes pointer from the Hash) */
bool DataContainer::deleteDescriptor( string aLabel, bool found ){

	// first check for the descriptor among the int descriptors
	map<const string, Descriptor<int>* >::iterator iti;
	for( iti = (intDescriptors).begin(); iti != (intDescriptors).end(); iti++ ){
	if( (*iti).first == aLabel ){
			#ifdef DEBUG
				cout << "deleting int descriptor [ " << (*iti).second->toString() << "]" << endl;
			#endif
			delete (*iti).second;
			(intDescriptors).erase(iti);
			found = true;
			return( true );
		}
	}

	 // then check for the descriptor among the float descriptors
	map<const string, Descriptor<float>* >::iterator itf;
	for( itf = (floatDescriptors).begin(); itf != (floatDescriptors).end(); itf++ ){
	if( (*itf).first == aLabel ){
			#ifdef DEBUG
				cout << "deleting float descriptor [ " << (*itf).second->toString() << "]" << endl;
			#endif
			delete (*itf).second;
			(floatDescriptors).erase(itf);
			found = true;
			return( true );
		}
	}

	 // finally check for the descriptor among the string descriptors
	map<const string, Descriptor<string>* >::iterator its;
	for( its = (stringDescriptors).begin(); its != (stringDescriptors).end(); its++ ){
	if( (*its).first == aLabel ){
			#ifdef DEBUG
				cout << "deleting string descriptor [ " << (*its).second->toString() << "]" << endl;
			#endif
			delete (*its).second;
			(stringDescriptors).erase(its);
			found = true;
			return( true );
		}
	}

	return( found );
}

/** deletes all descriptors from ***Descriptors (deletes the object and removes pointer from the Hash) */
void DataContainer::deleteAllDescriptors(){
	// first delete int descriptors
	#ifdef DEBUG
		int i = 0;
	#endif
	map<const string, Descriptor<int>* >::iterator iti;
	for( iti = intDescriptors.begin(); iti != intDescriptors.end(); iti++ ){
			#ifdef DEBUG
				cout << "deleting int descriptor [ " << (*iti).second->toString() << "]" << endl;
				i++;
			#endif
			delete (*iti).second;
   		//(intDescriptors).erase(iti);
	}
	intDescriptors.clear();

	// then delete float descriptors
	map<const string, Descriptor<float>* >::iterator itf;
	for( itf = floatDescriptors.begin(); itf != floatDescriptors.end(); itf++ ){
			#ifdef DEBUG
				cout << "deleting float descriptor [ " << (*itf).second->toString() << "]" << endl;
				i++;
			#endif
			delete (*itf).second;
   		//(floatDescriptors).erase(itf);
	}
	floatDescriptors.clear();

	// finally delete string descriptors
	map<const string, Descriptor<string>* >::iterator its;
	for( its = stringDescriptors.begin(); its != stringDescriptors.end(); its++ ){
			#ifdef DEBUG
				cout << "deleting string descriptor [ " << (*its).second->toString() << "]" << endl;
				i++;
			#endif
			delete (*its).second;
   		//(stringDescriptors).erase(its);
	}
	stringDescriptors.clear();
	#ifdef DEBUG
		cout << " DELETED " << i << " DESCRIPTORS" << endl;
	#endif
}


/** deletes all kindDescriptors from ***Descriptors (deletes the object and removes pointer from the Hash) */
void DataContainer::deleteAllKindDescriptors(){
	// first delete int descriptors
	//cout << "REMOVING ALLKINDDESCRIPTORS ####" << endl;

	#ifdef DEBUG
		cout << "DataContainer::deleteAllKindDescriptors" << endl;
	#endif
	int i = 0;
	map<const string, Descriptor<int>* >::iterator iti;
	for( iti = (kindIntDescriptors)->begin(); iti != (kindIntDescriptors)->end(); iti++ ){
			#ifdef DEBUG
				cout << "deleting kindInt descriptor [ " << (*iti).second->toString() << "]" << endl;
			#endif
			delete (*iti).second;
			i++;
   		//(kindIntDescriptors).erase(iti);
	}
	kindIntDescriptors->clear();

	// then delete float descriptors
	map<const string, Descriptor<float>* >::iterator itf;
	for( itf = (kindFloatDescriptors)->begin(); itf != (kindFloatDescriptors)->end(); itf++ ){
			#ifdef DEBUG
				cout << "deleting kindFloat descriptor [ " << (*itf).second->toString() << "]" << endl;
			#endif
			delete (*itf).second;
			i++;;
   		//(kindFloatDescriptors).erase(itf);
	}
	kindFloatDescriptors->clear();

	// finally delete string descriptors
	map<const string, Descriptor<string>* >::iterator its;
	for( its = (kindStringDescriptors)->begin(); its != (kindStringDescriptors)->end(); its++ ){
			#ifdef DEBUG
				cout << "deleting kindString descriptor [ " << (*its).second->toString() << "]" << endl;
			#endif
			delete (*its).second;
			i++;
   		//(kindStringDescriptors).erase(its);
	}
	kindStringDescriptors->clear();

	delete kindIntDescriptors;
	delete kindFloatDescriptors;
	delete kindStringDescriptors;

	#ifdef DEBUG
		cout << " DELETED " << i << " KIND DESCRIPTORS" << endl;
	#endif
}



/** returns true if the datacontainer has a string descriptor with label aLabel */
bool DataContainer::hasStringDescriptor( string aLabel ){

	if( stringDescriptors.find( aLabel ) == stringDescriptors.end() ){
		return( false );
	}else{
		return( true );
	}

}

/** returns true if the datacontainer has a float descriptor with label aLabel */
bool DataContainer::hasFloatDescriptor( string aLabel ){

	if( floatDescriptors.find( aLabel ) == floatDescriptors.end() ){
		return( false );
	}else{
		return( true );
	}

}

/** returns true if the datacontainer has a int descriptor with label aLabel */
bool DataContainer::hasIntDescriptor( string aLabel ){

	if( intDescriptors.find( aLabel ) == intDescriptors.end() ){
		return( false );
	}else{
		return( true );
	}

}


/** adds a descriptor to the datacontainer. If aName terminates with .string, a string descriptor will be added. If the aName terminates with a .integer, an integer descriptor will be added. If aName terminates with .float, a float descriptor will be added. */
void DataContainer::addUnknownTypeDescriptor( string aDescriptorName, string aValue ){

	if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),8), 8 ) == ".integer"){

			aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 8 );
			setIntDescriptor( aDescriptorName, StringUtils::toInt( aValue ), "", "", true, true );

	}else if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),4), 4 ) == ".int" ){

			aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 4);
			setIntDescriptor( aDescriptorName, StringUtils::toInt( aValue ), "", "", true, true );

	}else if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),6), 6 ) == ".float"){

		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 6 );
		setFloatDescriptor( aDescriptorName, StringUtils::toFloat( aValue ), "", "", true, true );

	}else if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),4), 4 ) == ".flo" ){

		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 4 );
		setFloatDescriptor( aDescriptorName, StringUtils::toFloat( aValue ), "", "", true, true );

	}else	if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),7), 7 ) == ".string"){

		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 7 );
		setStringDescriptor( aDescriptorName, aValue, "", "", true, true );

	}else if( aDescriptorName.substr( aDescriptorName.length() - min( (int) aDescriptorName.length(),4), 4 ) == ".str" ){

		aDescriptorName = aDescriptorName.substr( 0, aDescriptorName.size() - 4 );
		setStringDescriptor( aDescriptorName, aValue, "", "", true, true );

	}else{
		setStringDescriptor( aDescriptorName, aValue, "", "", true, true );

	}
}
