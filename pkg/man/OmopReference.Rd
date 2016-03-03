\name{OmopReference}
\alias{OmopReference}
\docType{data}
\title{
The OMOP reference set
}
\description{
 The OMOP Ryan et al.(2013) reference set of test cases that contains both positive and negative controls.
}
\usage{data(OmopReference)}
\format{
  A set of 399 negative and positive controls or couples of adverse event and drug.
  \describe{
    \item{\code{AE}}{a character vector  indicating adverse events}
    \item{\code{ATC}}{The Anatomical Therapeutic Chemical (ATC) code}
    \item{\code{statu}}{a numeric vector indicating if positive (1) or negative (0) control}
  }
}
\details{
To evaluate and compare the performances of the competing methods, we use the OMOP 
Ryan et al.(2013) reference set of test cases that contains both positive and negative 
controls. Four adverse events (i.e. d = 4) were studied in this reference set : 
acute myocardial infarction (AMI), acute kidney injury (AKI), acute liver injury (ALI), and 
upper gastro-intestinal bleeding (GIB). There are three-hundred and ninety-nine test cases 
where 165 positive controls and 234 negative controls were identified across the
four adverse events of interest. More details are given by Table 1. Ryan et al. (2013) 
indicate that the majority of positive controls for AKI and GIB were supported by randomized 
clinical trial evidence, while the majority of positive controls for ALI and AMI were only based 
on published case reports. Then a positive control is a pair (drug and adverse event) which is 
supported by randomized clinical trial evidence or based on published case reports whereas a 
negative control is not met by neither of the these methods.
}
\source{
Ryan, P., Schuemie, M., Welebob, E., Duke, J., Valentine, S., and Hartzema, A. (2013). 
Defining a reference set to support methodological research in drug safety. Drug safety, 36(1), 33-47.
}
\references{
Ryan, P., Schuemie, M., Welebob, E., Duke, J., Valentine, S., and Hartzema, A. (2013). 
Defining a reference set to support methodological research in drug safety. Drug safety, 36(1), 33-47.
}
\examples{
data(OmopReference)
}
\keyword{Reference set}
\keyword{OMOP}
\keyword{Positive controls}
\keyword{Negative controls}