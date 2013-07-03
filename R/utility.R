#' @title createRMolecule
#' 
#' @description Creates an Rchemcpp::Rmolecule from an atom-vector and a bond-matrix 
#' 
#' @usage createRMolecule(atoms, bonds)
#' 
#' @param atoms A vector containing the symbol names of all atoms in the molecule
#' @param bonds A matrix with the same number of rows and columns as the atoms-vector
#' containing the type of bonds between the atoms
#' @return an instance of "molecule"
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' 
#' @examples 
#' m <- createRMolecule(c("C","C"),matrix(c(0,3,3,0),nrow=2))
#' @export


createRMolecule = function(atoms, bonds)
{
	if(!is.vector(atoms)) stop("atoms must be vector")
	if(!is.matrix(bonds)) stop("bonds must be matrix")

	bonds = apply(bonds,1,as.integer)

	if (length(atoms) != nrow(bonds)) stop("matrix must have as many rows as vector")
	if (length(atoms) != ncol(bonds)) stop("matrix must have as many cols as vector has rows")

	if (! isSymmetric(bonds) ) stop("matrix must be symetric")
	if (! all(diag(bonds) == 0)) stop("cannot connect atoms with themselves")	

	if ( (max(bonds) > 3) || (min(bonds) < 0 ) ) stop("Bond types can only range from 1 to 3") 

	m = new(Rchemcpp::Rmolecule)

	for (k in atoms)
	{
		m$addAtom(k)
	}

	for (i in 1:length(atoms)){
		for (j in 1:length(atoms)){
			if ((i < j) && (bonds[i,j] != 0)){
				m$linkAtoms(i-1,j-1,bonds[i,j]);
			}
		}
	}

	return (m)
}



