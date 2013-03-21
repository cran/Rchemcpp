\name{sd2gramSpectrum}
\alias{sd2gramSpectrum}
\title{sd2gramSpectrum - Similarity of molecules by walk-based graph kernels}
\usage{
  sd2gramSpectrum(sdf, sdf2,
    kernelType = c("spectrum", "tanimoto", "minmaxTanimoto", "marginalized", "lambda"),
    margKernelEndProbability = 0.1, lambdaKernelLambda = 1,
    depthMax = as.integer(3), onlyDepthMax = FALSE,
    flagRemoveH = FALSE, morganOrder = as.integer(0),
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

  \item{kernelType}{Sets which kernel is to be used.
  Options are "spectrum (Spectrum kernel) , "tanimoto"
  (Tanimoto kernel), "minmaxTanimoto" (MinMax Tanimoto
  kernel), "marginalized (Marginalized kernel
  approximation) and "lambda" (LambdaK kernel). Default =
  "spectrum".}

  \item{margKernelEndProbability}{The ending probability
  for the marginalized kernel. Default = 0.1.}

  \item{lambdaKernelLambda}{The lambda parameter of the
  LambdaK kernel. Default = 1.0.}

  \item{depthMax}{The maximal length of the molecular
  fragments.}

  \item{onlyDepthMax}{Whether fragments up to the given
  length should be used or only fragments of the given
  length.}

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
  This function computes several walk-based graph kernel
  functions based on finite length walks and a fast
  implementation for input sd file(s).
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gramSpectrum(sdf)
}
\author{
  Michael Mahr <rchemcpp@bioinf.jku.at> c++ function
  written by Jean-Luc Perret and Pierre Mahe
}

