
#include "Relements.h"

void loadGramAtoms( std::string aFileName ){
	elements.loadGramAtoms( aFileName );
}

RCPP_MODULE(Relements) {
	Rcpp::function( "loadGramAtoms", &loadGramAtoms );
}
