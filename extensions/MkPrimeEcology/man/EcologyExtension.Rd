\name{EcologyExtension}
\alias{EcologyExtension}
\title{Ecology-biased equilibrium model}
\usage{
EcologyExtension(
  data,
  ecology,
  neomorphic = integer(0),
  refEcology = NULL,
  magnitudeMode = c("global", "per_ecology"),
  pi0Prior = c(360, 120),
  thetaPrior = c(2, 2),
  sigmaPhi = 1.5,
  zWeight = 5,
  init = list()
)
}
\arguments{
\item{data}{A \code{phyDat} object or \code{MkPrime::MkPrimeData()} result.}
\item{ecology}{Ecology of each tip, named by tip label; \code{NA} for unknown.}
\item{neomorphic}{Indices of neomorphic characters, as for \code{MkPrime::RunMkPrime()}.}
\item{refEcology}{The reference ecology, in which no character is tilted.}
\item{magnitudeMode}{\code{"global"} for one \code{phi}, or \code{"per_ecology"}.}
\item{pi0Prior, thetaPrior}{Beta shape parameters of the priors on \code{pi0} and \code{theta}.}
\item{sigmaPhi}{Standard deviation of the log-normal prior on \code{phi}.}
\item{zWeight}{Relative frequency of the Gibbs move on \code{z}.}
\item{init}{Named list of initial values: \code{phi}, \code{pi0}, \code{theta}.}
}
\value{An \code{MkPrimeExtension} object for \code{MkPrime::MkPrimeModel(extensions = ...)}.}
\description{
Builds the ecology-biased equilibrium (EBE) model as an MkPrime engine
extension. It scores MkPrime's neomorphic characters: on an edge in ecology
\code{s}, character \code{c}'s equilibrium of presence is tilted by
\code{logit(pi1) += s_z log(phi)} according to its influence category
\code{z[c, s]}, and edge ecology is a mixture weighted by a JC(K)
reconstruction of tip ecology. With every \code{z = 0} it is MkPrime's MkN.
}
