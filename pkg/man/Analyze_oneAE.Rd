\name{Analyze_oneAE}
\alias{Analyze_oneAE}
\title{
Signal detection using via variable selection in logistic regression. 
The Bayesian Information Criterion maximization is assessed using Metropolis-Hastings
algorithm.
}
\description{
Spontaneous adverse event reports have a high potential for detecting adverse drug reactions.
However, due to their dimension, the analysis of such databases requires statistical methods.
The MHTrajectoryR package propose to use a logistic regression whose sparsity is viewed as a model selection
challenge. Since the model space is huge, a Metropolis-Hastings algorithm carries out the model selection
by maximizing the BIC criterion.

}
\usage{
Analyze_oneAE(ae, drug, maxit, alpha, nbinit)
}
\arguments{
  \item{ae}{
The response binary vector which indicates if the adverse event is observed (value 1)
and (value 0) otherwise. It must have the sime length with the number of rows of matrix of drugs
consumption.
}
  \item{drug}{
The matrix of drugs consumption or the matrix of binary covariates. 
Each row corresponds to one individual drugs consumption. Each column corresponds
to one drug.
}
  \item{maxit}{
  Number of iteration of the Metropolis-Hastings algorithm. 
  In other words, the length of one trajectory of the Metropolis-Hastings
  Markov Chain.
}
  \item{alpha}{
The parameter that define the neighbordhood. 
}
  \item{nbinit}{
The number of random initialisations.
}
}
\value{ list of  (best) the best model that maximize the BIC.  (all) all explored models through the trajectory.
(signals) the detected signals. 
}
\references{
Matthieu Marbac, Pascale Tubert-Bitter, Mohammed Sedki: 
Bayesian model selection in logistic regression for the detection
of adverse drug reactions. \url{http://arxiv.org/abs/1505.03366} (accepted for publication
in Biometrical Journal).
}
\author{
Mohammed Sedki and Matthieu Marbac
}
\seealso{
\code{\link{FindSignals}}
}

\keyword{Binary data}
\keyword{Logistic regression} 
\keyword{Metropolis-Hastings algorithm} 
\keyword{Bayesian Information Criterion}
\keyword{Variable selection}
\keyword{Spontaneous reporting}
