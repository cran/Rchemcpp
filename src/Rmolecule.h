
#ifndef RMOLECULE_H
#define RMOLECULE_H

#include <Rcpp.h>

#include <molecule.h>

class Rmolecule : public Molecule {
public:
	
	Rmolecule(){};
	
	void addAtom(string aSymbol)
	{Molecule::addAtom(aSymbol);}

	void linkAtoms( int firstAtom, int secondAtom, int aBondLabel)
	{Molecule::linkAtoms(firstAtom, secondAtom, aBondLabel);}

	void writeSD( string aFileName )
	{Molecule::writeSD(aFileName);}

	//Extension 1.0.7
	std::vector<string> listStringDescriptors()
	{
		std::vector<string> vec;
		for (map< string, Descriptor< string >* >::iterator i = beginStringDescriptor(); i != endStringDescriptor(); i++)
		{
			if( !(*i).second->isEmpty() )
			{
				vec.push_back((*i).second->getLabel());
			}
			
		}	
		return vec;
	}

	//uncommon for molecules
	/*
	std::vector<string> listFloatDescriptors()
	{
		std::vector<string> vec;
		for (map< string, Descriptor< float >* >::iterator i = beginFloatDescriptor(); i != endFloatDescriptor(); i++)
		{
			if( !(*i).second->isEmpty() )
			{
				vec.push_back((*i).second->getLabel());
			}
			
		}	
		return vec;
	}

	std::vector<string> listIntDescriptors()
	{
		std::vector<string> vec;
		for (map< string, Descriptor< int >* >::iterator i = beginIntDescriptor(); i != endIntDescriptor(); i++)
		{
			if( !(*i).second->isEmpty() )
			{
				vec.push_back((*i).second->getLabel());
			}
			
		}	
		return vec;
	}
	*/
	
	
	string getStringDescriptorValue(string label){ return getStringDescriptor(label)->getValue(); }
	string getStringDescriptorUnit(string label){ return getStringDescriptor(label)->getUnit(); }
	string getStringDescriptorComment(string label){ return getStringDescriptor(label)->getComment(); }

	//uncommon for molecules -> commented out
	//float getFloatDescriptorValue(string label){ return getFloatDescriptor(label)->getValue(); }
	//float getFloatDescriptorUnit(string label){ return getFloatDescriptor(label)->getUnit(); }
	//float getFloatDescriptorComment(string label){ return getFloatDescriptor(label)->getComment(); }
	//int getIntDescriptorValue(string label){ return getIntDescriptor(label)->getValue(); }
	//int getIntDescriptorUnit(string label){ return getIntDescriptor(label)->getUnit(); }
	//int getIntDescriptorComment(string label){ return getIntDescriptor(label)->getComment(); }
	
	
	
	void setStringDescriptor(string aLabel, string aValue, string aUnit, string aComment)
	{
		//Code taken over from datacontainer.cpp (corrected bug: only value was set if existing)
		Descriptor< string >* d;
		if( hasStringDescriptor( aLabel ) ){
			d = stringDescriptors[aLabel];
			d -> setValue(aValue);
			d -> setUnit(aUnit);
			d -> setComment(aComment);
		}else{
			d = addStringDescriptor( aLabel, aValue, aUnit, aComment );
		}
	}
	
	void deleteStringDescriptor(string aLabel)
	{
		//Code taken from datacontainer.cpp (corrected bug: Original function deletes all types of descriptors with that label)
		map<const string, Descriptor<string>* >::iterator its;
		for( its = (stringDescriptors).begin(); its != (stringDescriptors).end(); its++ ){
			if( (*its).first == aLabel ){
				delete (*its).second;
				(stringDescriptors).erase(its);
				//found = true;
				//return( true );
			}
		}
		
	}
	
	//...


//private:

};

#endif

