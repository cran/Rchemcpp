\name{readRmoleculeset}
\alias{readRmoleculeset}
\title{Generating an Rmoleculeset from an SDF file}
\usage{
  readRmoleculeset(sdfFileName, detectArom = TRUE,
    bound = 70, type = 2)
}
\arguments{
  \item{sdfFileName}{The name of the SDF file containing
  the molecules.}

  \item{detectArom}{If the molecules in the SDF file have
  no annotated aromatic bonds, the ChemmineR function
  \code{rings} is used for detecting aromaticity. (Default
  = TRUE).}

  \item{bound}{Detection of aromaticity can be time
  consuming if the molecules are large. Detection is only
  done if the number of atoms is below the given number.
  (Default = 70).}

  \item{type}{Experimental parameter to switch between to
  types of the function.}
}
\value{
  An instance of Rmoleculeset.
}
\description{
  This function uses the ChemmineR package to read an SDF
  file and converts it into an \code{Rmoleculeset} that can
  be used as input for the kernel functions \code{sd2gram},
  \code{sd2gramSpectrum}, ..., \code{sd2gram3Dpharma}.
}
\author{
  Guenter Klambauer <rchemcpp@bioinf.jku.at>
}

