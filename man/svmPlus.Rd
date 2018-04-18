\name{svmPlus}
\alias{svmPlus}
\title{ Creates and returns an instance of the class QPSvmPlus. In future, LIBSVM and LibLinear based faster implementaion is planned to be supported.}
\description{
Creates and returns an instance of the class QPSvmPlus. In future, LIBSVM and LibLinear based faster implementaion is planned to be supported.
}
\usage{
svmPlus(cost = 1, gamma = 1, kernel_x = "rbf", degree_x = 3, gamma_x = 0.001, kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = 0.001, tol = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cost}{
  cost of constraints violation
}
  \item{gamma}{
  parameter needed for priviledged information
}
  \item{kernel_x}{
  the kernel used for standard training data
}
  \item{degree_x}{
parameter needed for polynomial kernel for training data
}
  \item{gamma_x}{
parameter needed for rbf kernel for training data
}
  \item{kernel_xstar}{
the kernel used for priviledged information (PI)
}
  \item{degree_xstar}{
parameter needed for polynomial kernel for PI
}
  \item{gamma_xstar}{
parameter needed for rbf kernel for PI
}
  \item{tol}{
tolerance of dual variables
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
instance of the class QPSvmPlus
}
\references{
[1] Vladimir et. al, 2009, 22, pp 544–557.
\url{<https://doi.org/10.1016/j.neunet.2009.06.042>}
[2] Li et. al, 2016.
\url{<https://github.com/okbalefthanded/svmplus_matlab>}
}
\author{
Niharika Gauraha and Ola Spjuth
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cost = 1, gamma = 1, kernel_x = "rbf", degree_x = 3,
    gamma_x = 0.001, kernel_xstar = "rbf", degree_xstar = 3,
    gamma_xstar = 0.001, tol = 1e-05)
{
    return(QPSvmPlus$new(cost = 1, gamma = 1, kernel_x = "rbf",
        degree_x = 3, gamma_x = 0.001, kernel_xstar = "rbf",
        degree_xstar = 3, gamma_xstar = 0.001, tol = 1e-05))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
