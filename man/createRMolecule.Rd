\name{createRMolecule}
\alias{createRMolecule}
\title{createRMolecule}
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
  Creates an \"Rmolecule\" from an atom-vector and a
  bond-matrix
}
\examples{
m <- createRMolecule(c("C","C"),matrix(c(0,3,3,0),nrow=2))
}
\author{
  Michael Mahr <rchemcpp@bioinf.jku.at>
}

