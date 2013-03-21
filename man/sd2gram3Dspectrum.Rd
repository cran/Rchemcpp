\name{sd2gram3Dspectrum}
\alias{sd2gram3Dspectrum}
\title{sd2gram3Dspectrum - Similarity of molecules by fast
approximations of the pharmacophore kernel}
\usage{
  sd2gram3Dspectrum(sdf, sdf2, chargesFileName = "",
    chargesFileName2 = "",
    kernelType = c("3Pspectrum", "3Pbinary", "3Ptanimoto", "2Pspectrum", "2Pbinary", "2Ptanimoto"),
    depthMax = as.integer(3), nBins = as.integer(20),
    distMin = 0, distMax = 20, flagRemoveH = FALSE,
    morganOrder = as.integer(0), chargesThreshold = 0,
    silentMode = FALSE, returnNormalized = TRUE,
    detectArom = TRUE)
}
\arguments{
  \item{sdf}{File containing the molecules. Must be in MDL
  file format (MOL and SDF files). For more information on
  the file format see
  http://en.wikipedia.org/wiki/Chemical_table_file. Default
  = "missing".}

  \item{sdf2}{A second file containing molecules. Must also
  be in SDF. If specified the molecules of the first file
  will be compared with the molecules of this second file.
  Default = "missing".}

  \item{chargesFileName}{A character with the name of the
  file containing the atom charges. Default = missing.}

  \item{chargesFileName2}{A character with the name of the
  file containing the atom charges. Default = missing.}

  \item{kernelType}{Type of kernel to be used. Possible
  choices are 3-points spectrum kernel ("3Pspectrum"),
  3-points binary kernel ("3Pbinary"), 3-points Tanimoto
  kernel ("3Ptanimoto"), 2-points spectrum kernel
  ("2Pspectrum"), 2-points binary kernel ("2Pbinary"),
  2-points Tanimoto kernel ("2Ptanimoto"). Default =
  "3Pspectrum".}

  \item{depthMax}{The maximal length of the molecular
  fragments. Default = 3.}

  \item{nBins}{number of bins used to discretize the
  inter-atomic lengths. An adequate value for the number of
  bins is between 20 and 30. Default = 20.}

  \item{distMin}{minimum distance for inter-atomic distance
  range. Default = 0.}

  \item{distMax}{maximum distance in angstrom for
  inter-atomic distance range. Default = 20.}

  \item{chargesThreshold}{specifies a threshold above which
  partial charges are considered as positive/negative. By
  default this threshold is zero, and every positive (resp.
  negative) partial charge is seen as a positive (resp.
  negative) charge. However, it might be interesting to
  consider a threshold of 0.2 for example, in which case
  only partial charges greater than 0.2 (resp. smaller than
  -0.2) would be seen as positive (resp. negative). Default
  = 0.}

  \item{flagRemoveH}{A logical that indicates whether
  H-atoms should be removed or not.}

  \item{morganOrder}{The order of the DeMorgan Indices to
  be used. If set to zero no DeMorgan Indices are used. The
  higher the order the more different types of atoms exist
  and consequently the more dissimilar will be the
  molecules.}

  \item{silentMode}{Whether or not the program should print
  progress reports to the standart output.}

  \item{returnNormalized}{A logical specifying whether a
  normalized kernel matrix should be returned. Default =
  TRUE.}

  \item{detectArom}{Whether aromatic rings should be
  detected and aromatic bonds should a special bond type.
  (Default = TRUE).}
}
\value{
  A numeric matrix containing the similarity values between
  the molecules.
}
\description{
  This tool implements the six discrete approximations of
  the pharmacophore kernel presented in "The pharmacophore
  kernel for virtual screening with support vector
  machines" (\cite{Mahe, 2006}).
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gram3Dspectrum(sdf)
}
\author{
  Michael Mahr <rchemcpp@bioinf.jku.at> c++ function
  written by Jean-Luc Perret and Pierre Mahe
}
\references{
  (Mahe, 2006) -- P. Mahe, L. Ralaivola, V. Stoven, and
  J.-P. Vert. The pharmacophore kernel for virtual
  screening with support vector machines. Technical Report,
  HAL:ccsd-00020066, Ecole des Mines de Paris, March 2006.
}

