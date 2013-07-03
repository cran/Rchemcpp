\name{sd2gram}
\alias{sd2gram}
\title{sd2gram - Similarity of molecules by the marginalized kernel and
proposed extensions.}
\usage{
  sd2gram(sdFileName, sdFileName2 = "", stopP = 0.1,
    filterTottering = FALSE, converg = as.integer(1000),
    atomKernelMatrix = "", flagRemoveH = FALSE, morganOrder
    = as.integer(0), fileType = c("sd","genericsd","kcf"),
    silentMode = FALSE, returnNormalized = TRUE,
    moleculeNameProperty = "", moleculeNameProperty2 = "")
}
\arguments{
  \item{sdFileName}{File containing the molecules. Must be
  in MDL file format (MOL and SDF files). For more
  information on the file format see
  http://en.wikipedia.org/wiki/Chemical_table_file. Default
  = "missing".}

  \item{sdFileName2}{A second file containing molecules.
  Must also be in SDF. If specified the molecules of the
  first file will be compared with the molecules of this
  second file. Default = "missing".}

  \item{stopP}{... . Default = 0.1.}

  \item{filterTottering}{A logical specifying whether
  tottering paths should be removed. Default = FALSE.}

  \item{converg}{A numeric value specifying when
  convergence is reached. The algorithm stops when the
  kernel value does not change by more than 1/c, where c is
  the value specified by the converg option. Default =
  1000.}

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

  \item{fileType}{Which filetype was submitted.}

  \item{silentMode}{Whether or not the program should print
  progress reports to the standart output. Default =
  FALSE.}

  \item{returnNormalized}{A logical specifying whether a
  normalized kernel matrix should be returned. Default =
  TRUE.}

  \item{moleculeNameProperty}{A string which specifies the
  name of the property of the molecules in the sdFile from
  which the row names (and column names) are read from.
  Default = "".}

  \item{moleculeNameProperty2}{A string which specifies the
  name of the property of the molecules in the sdFile 2
  from which the column names are read from. Default = "".}
}
\value{
  A numeric matrix containing the similarity values between
  the molecules.
}
\description{
  This tools compute the marginalized kernel
  (\cite{Kashima, 2004}) and its proposed extensions
  (\cite{Mahe, 2005)}.
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gram(sdf, moleculeNameProperty="Compound Name")
}
\author{
  Michael Mahr <rchemcpp@bioinf.jku.at> c++ function
  written by Jean-Luc Perret and Pierre Mahe
}
\references{
  (Kashima, 2004) -- H. Kashima, K. Tsuda, and A. Inokuchi.
  Kernels for graphs. In B. Schoelkopf, K. Tsuda, and J.P.
  Vert, editors, Kernel Methods in Computational Biology,
  pages 155-170. MIT Press, 2004.

  (Mahe, 2005) -- P. Mahe, N. Ueda, T. Akutsu, J.-L.
  Perret, and J.-P. Vert. Graph kernels for molecular
  structure- activity relationship analysis with support
  vector machines. J Chem Inf Model, 45(4):939-51, 2005.
}

