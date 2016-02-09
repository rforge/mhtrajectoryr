\name{MHTrajectoryR-package}
\alias{MHTrajectoryR-package}
\alias{MHTrajectoryR}
\docType{package}
\title{Bayesian model selection in logistic regression for the detection
of adverse drug reactions
}
\description{
Spontaneous adverse event reports have a high potential for detecting adverse drug reactions.
However, due to their dimension, the analysis of such databases requires statistical methods.
The MHTrajectoryR package propose to use a logistic regression whose sparsity is viewed as a model selection
challenge. Since the model space is huge, a Metropolis-Hastings algorithm carries out the model selection
by maximizing the BIC criterion.
}
\details{
\tabular{ll}{
Package: \tab MHTrajectoryR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-07\cr
License: \tab GPL (>= 2)\cr
}
The main function is Analyze_oneAE.
}
\author{
Matthieu Marbac and Mohammed Sedki
Maintainer: Mohammed Sedki <mohammed.sedki@u-psud.fr>
}
\references{
Matthieu Marbac, Pascale Tubert-Bitter, Mohammed Sedki: 
Bayesian model selection in logistic regression for the detection
of adverse drug reactions. (http://arxiv.org/abs/1505.03366) (accepted for publication
in Biometrical Journal).
}
\keyword{Binary data, logistic regression, Metropolis-Hastings algorithm, Bayesian Information Crite-
rion, model selection,spontaneous reporting}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
