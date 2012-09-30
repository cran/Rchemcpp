\name{createRMolecule}
\alias{createRMolecule}
\title{Rchemcpp::createRMolecule}
\usage{
  createRMolecule(atoms, bonds)
}
\arguments{
  \item{atoms}{A vector containing the symbol names of all
  atoms in the molecule}

  \item{bonds}{A matrix with the same number of rows and
  columns as the atoms-vector containing the type of bonds
  between the atoms}
}
\value{
  an instance of "molecule"
}
\description{
  This function creates an Rchemcpp::Rmolecule from an atom-vector
  and a bond-matrix.
  The resulting object can be added to an Rmoleculeset afterwards.
}
\author{
  Michael Mahr <rchemcpp@bioinf.jku.at>
}
\examples{
mset = new (Rchemcpp::Rmoleculeset)
mole = createRMolecule(c("H","O","H"),matrix(data =
c(0,1,0,
  1,0,1,
  0,1,0), ncol = 3))
mset$addMoleculeCopy(mole)
}
