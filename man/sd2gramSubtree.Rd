\name{sd2gramSubtree}
\alias{sd2gramSubtree}
\title{sd2gramSubtree - Similarity of molecules by several graph kernels
based on the count of common subtrees}
\usage{
  sd2gramSubtree(sdf, sdf2,
    kernelType = c("sizebased", "branchingbased"),
    branchKernelUntilN = FALSE, lambda = 1,
    depthMax = as.integer(3), flagRemoveH = FALSE,
    filterTottering = FALSE, morganOrder = as.integer(0),
    silentMode = FALSE, returnNormalized = TRUE,
    detectArom = TRUE)
}
\arguments{
  \item{sdf}{File containing the molecules. Must be in MDL
  file format (MOL and SDF files). For more information on
  the file format see
  http://en.wikipedia.org/wiki/Chemical_table_file. Default
  = "missing"}

  \item{sdf2}{A second file containing molecules. Must also
  be in SDF. If specified the molecules of the first file
  will be compared with the molecules of this second file.
  Default = "missing".}

  \item{kernelType}{Determines whether subtrees of the
  molecule are penalized size-based or branching-based.
  Default = "sizebased".}

  \item{branchKernelUntilN}{Logical whether tree patterns
  of until N should be considered. Default = FALSE.}

  \item{lambda}{Weighted contribution of tree-patterns
  depending on their sizes Default = 1.}

  \item{depthMax}{tree-patterns of depth. Default = 3.}

  \item{flagRemoveH}{A logical that indicates whether
  H-atoms should be removed or not.}

  \item{filterTottering}{A logical that indicates whether
  tottering paths should be removed. Default = FALSE.}

  \item{morganOrder}{The order of the DeMorgan Indices to
  be used. If set to zero no DeMorgan Indices are used. The
  higher the order the more different types of atoms exist
  and consequently the more dissimilar will be the
  molecules.}

  \item{silentMode}{Whether the program should print
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
  This tools computes several graph kernels based on the
  detection of common subtrees: the so-called tree-pattern
  graph kernels, originally introduced in (\cite{Ramon,
  2003}), and revisited in (\cite{Mahe, 2006}).
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
  (Mahe, 2006) -- P. Mahe and J.-P. Vert. Graph kernels
  based on tree patterns for molecules. Technical Report,
  HAL:ccsd-00095488, Ecoles des Mines de Paris, September
  2006. (Ramon, 2003) -- J. Ramon and T. Gaertner.
  Expressivity versus efficiency of graph kernels. In T.
  Washio and L. De Raedt, editors, Proceedings of the First
  International Workshop on Mining Graphs, Trees and
  Sequences, pages 65-74, 2003.
}

