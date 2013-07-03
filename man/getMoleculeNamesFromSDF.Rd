\name{getMoleculeNamesFromSDF}
\alias{getMoleculeNamesFromSDF}
\title{getMoleculeNamesFromSDF - a helper function}
\usage{
  getMoleculeNamesFromSDF(sdfile)
}
\arguments{
  \item{sdfile}{A character containing the name of the SDF
  file.}
}
\value{
  A character vector with one name per molecule.
}
\description{
  This function helps to extract a certain property from an
  SDF file. Usually the molecule class, like
  "active/non-active" or a property of the molecule, like
  "biological activity", is also stored in the SDF file.
  These values often serve as targets for a prediction
  task. This function is a small wrapper that extracts the
  information.
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
moleculeNames <- getMoleculeNamesFromSDF(sdf)
}
\author{
  Guenter Klambauer <rchemcpp@bioinf.jku.at>
}

