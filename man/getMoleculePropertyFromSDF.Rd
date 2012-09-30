\name{getMoleculePropertyFromSDF}
\alias{getMoleculePropertyFromSDF}
\title{Rchemcpp::getMoleculePropertyFromSDF}
\usage{
  getMoleculePropertyFromSDF(sdfile, property)
}
\arguments{
  \item{sdfile}{A character containing the name of the SDF
  file.}

  \item{property}{The name of the slot in the SDF.}
}
\value{
  A character vector with one value per molecule.
}
\description{
  getMoleculePropertyFromSDF - a helper function

  This function helps to extract a certain property from an
  SDF file. Usually the molecule class, like
  "active/non-active" or a property of the molecule, like
  "biological activity", is also stored in the SDF file.
  These values often serve as targets for a prediction
  task. This function is a small wrapper that extracts the
  information.
}
\author{
  Guenter Klambauer <rchemcpp@bioinf.jku.at>
}

