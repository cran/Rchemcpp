\name{Rcpp_Rmolecule-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_Rmolecule-class}
\alias{Rmolecule}

\title{Class \code{"Rcpp_Rmolecule"}}
\description{
This class is a Rcpp modules wrapper for the chemcpp c++ class "Molecule".
It allows creating molecules from scratch or manipulating existing ones.
Currently it exposes only a small fraction of functionality of the base class.
Please note that only a part of the original chemcpp class "Molecule"
is exposed until now.
}
\section{Extends}{
chemcpp c++ class "Molecule"
}
\author{
Michael Mahr;
base class written by Jean-Luc Perret and Pierre Mahe
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set = new (Rmoleculeset)
mol = new (Rmolecule)
mol$addAtom("H")
set$addMoleculeCopy(mol)
}
\keyword{classes}
\section{Methods}{
  \describe{
%%    \item{\code{initialize(...)}:}{ ~~ }
%%    \item{\code{finalize()}:}{ ~~ }
    \item{\code{writeSD(...)}:}{ Write molecule to sd file }
    \item{\code{linkAtoms(...)}:}{ Create a bond between two atoms; Atom index is zero-based }
    \item{\code{addAtom(...)}:}{ Add an atom by specifying its character symbol }

    \item{\code{listStringDescriptors(...)}:}{ Return a vector of all string descriptors of the molecule }
    \item{\code{getStringDescriptorValue(...)}:}{ Return the value of one string descriptor }
    \item{\code{getStringDescriptorUnit(...)}:}{ Return the unit of one string descriptor }
    \item{\code{getStringDescriptorComment(...)}:}{ Return the comment of one string descriptor }

    \item{\code{setStringDescriptor(...)}:}{ Create or replace a string descriptor of the molecule by specifying the name, value, unit and comment }
    \item{\code{deleteStringDescriptor(...)}:}{ Delete one string descriptor from the molecule }

  }
}
