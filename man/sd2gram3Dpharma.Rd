\name{sd2gram3Dpharma}
\alias{sd2gram3Dpharma}
\title{sd2gram3Dpharma - Similarity of molecules by the exact pharmacophore
kernel.}
\usage{
  sd2gram3Dpharma(sdf, sdf2, chargesFileName = "",
    chargesFileName2 = "",
    edgeKernelType = c("RBF", "triangular"),
    edgeKernelParameter = 1, atomKernelMatrix = "",
    flagRemoveH = FALSE, morganOrder = as.integer(0),
    morganChargesThreshold = 0, silentMode = FALSE,
    returnNormalized = TRUE, detectArom = TRUE)
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

  \item{edgeKernelType}{Options to specify the kernel
  function comparing distances between atoms. Choices are
  "RBF" or "triangular". Default = "RBF".}

  \item{edgeKernelParameter}{Specifies the parameter
  associated to these kernels. Either the bandwith of the
  RBF kernel or the cut-off of the triangular kernel.
  Default = 1.}

  \item{atomKernelMatrix}{A string that sets the similarity
  measure between atoms that should be used. Dfault =
  "missing".}

  \item{flagRemoveH}{A logical that indicates whether
  H-atoms should be removed or not.}

  \item{morganOrder}{The order of the DeMorgan Indices to
  be used. If set to zero no DeMorgan Indices are used. The
  higher the order the more different types of atoms exist
  and consequently the more dissimilar will be the
  molecules.}

  \item{morganChargesThreshold}{specifies a threshold above
  which partial Morgan charges are considered as
  positive/negative.}

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
  This tool implements the (exact version of) pharmacophore
  kernel for 3D structures of molecules (\cite{Mahe,
  2006}).
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gram(sdf)
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

