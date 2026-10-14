\name{CanonicaliseEcology}
\alias{CanonicaliseEcology}
\title{Fold the phi label symmetry}
\usage{
CanonicaliseEcology(samples, extension)
}
\arguments{
\item{samples}{Matrix or data frame of samples, e.g. \code{result$samples}.}
\item{extension}{The \code{EcologyExtension()} object used for the run.}
}
\value{\code{samples}, mapped to the representative with \code{phi >= 1}.}
\description{
The likelihood and prior are unchanged by \code{phi -> 1/phi} with
\code{theta -> 1 - theta} and \code{z} 1 and 2 swapped; this maps each sample
to the representative with \code{phi >= 1}.
}
