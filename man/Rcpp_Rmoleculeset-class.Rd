\name{Rcpp_Rmoleculeset-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_Rmoleculeset-class}
\alias{Rmoleculeset}

\title{Class \code{"Rcpp_Rmoleculeset"}}
\description{
This class is a Rcpp modules wrapper for the chemcpp c++ class "MoleculeSet".
It allows reading molecule-files and computing simple comparison-matrices.
When calling the function "setComparisonSet" however, the argument object is
copied (instead of storing a reference).
Please note that only a part of the original chemcpp class "MoleculeSet"
is exposed until now.
}
\section{Extends}{
chemcpp c++ class "MoleculeSet"
}
\author{
Michael Mahr
base class written by Jean-Luc Perret and Pierre Mahe
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
set = new (Rmoleculeset)
set$addSD(sdf, FALSE)
}
\keyword{classes}
\section{Methods}{
  \describe{
%%    \item{\code{initialize(...)}:}{ ~~ }
%%    \item{\code{finalize()}:}{ ~~ }
    \item{\code{writeSelfKernelList(...)}:}{ ~~ }
    \item{\code{writeGramMatrix(...)}:}{ Write the gram matrix to a file, if one has been computed }
    \item{\code{setMorganLabels(...)}:}{ ~~ }
    \item{\code{setMorganChargesLabels(...)}:}{ ~~ }
    \item{\code{setKashimaKernelParam(...)}:}{ ~~ }
    \item{\code{setComparisonSetSelf(...)}:}{ Set the comparison set to be the set itself; NOTE: this is the preferred way to compare a set with itself, because faster implementations are used for comparison this way }
    \item{\code{setComparisonSetCopy(...)}:}{ Set the comparison set to be a different set of molecules; NOTE: this function copies the object specified as argument }
    \item{\code{readPartialCharges(...)}:}{ ~~ }
    \item{\code{numMolecules(...)}:}{ Returns the number of contained molecules }
    \item{\code{normalizeTanimoto_raw(...)}:}{ ~~ }
    \item{\code{normalizeTanimotoMinMax(...)}:}{ ~~ }
    \item{\code{normalizeTanimoto(...)}:}{ ~~ }
    \item{\code{normalizeGram_raw(...)}:}{ ~~ }
    \item{\code{normalizeGram(...)}:}{ ~~ }
    \item{\code{noTottersTransform(...)}:}{ ~~ }
    \item{\code{initializeSelfKernel(...)}:}{ ~~ }
    \item{\code{initializeGram(...)}:}{ ~~ }
    \item{\code{hideHydrogens(...)}:}{ Hide hydrogen atoms in all contained molecules }
    \item{\code{gramCompute3D(...)}:}{ ~~ }
    \item{\code{gramCompute(...)}:}{ ~~ }
    \item{\code{getGramNormal(...)}:}{ Return the normalized gram matrix, if one has been computed }
    \item{\code{getGram(...)}:}{ Return the gram matrix, if one has been computed }
    \item{\code{getComparisonSet(...)}:}{ Return A POINTER to the comparison set contained in the set; NOTE: this pointer expires when the set is destroyed or a different comparison set is set }
    \item{\code{bondsListing(...)}:}{ Return a list of all bonds which are present in the set }
    \item{\code{atomsLabelsListing(...)}:}{ Return a list of all atom symbols which are present in the set }
    \item{\code{addSD2(...)}:}{ Load a file containing molecules }
    \item{\code{addSD(...)}:}{ Load a file containing molecules }
    \item{\code{addMoleculeCopy(...)}:}{ Add a copy of a molecule object to the set  }
    \item{\code{addKCF2(...)}:}{ Load a file containing molecules }
    \item{\code{addKCF(...)}:}{ Load a file containing molecules }
  }
}
