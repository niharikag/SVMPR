\name{SVMP}
\alias{SVMP}
\title{ Creates and returns an instance of the class specified in the svm_type.}
\description{
Creates and returns an instance of the class specified
in the svm_type.
In future, the current solver used for quadratic programming (quadprog)  will be replaced by
the equivaent quadprog solver defined in CVXR package. Also, LIBSVM and LIBLINEAR based faster implementaions are planned to be supported.
}
\usage{
SVMP(cost = 1, gamma = 1, kernel_x = "rbf", degree_x = 3,
gamma_x = 0.001, kernel_xstar = "rbf", degree_xstar = 3,
gamma_xstar = 0.001, tol = 1e-05, svm_type = "QP")
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
  \item{svm_type}{
optimization techiniques used: QP, LibSVM, LibLinear etc.
Currently it supports only QP.
}
}
\value{
an instance of the class specified in the svm_type.
Currently it suports only "QP", hence returns instance of the class
QPSvmPlus. The return instance can be used to call  fit, project
and predict methods of the QPSvmPlus.
}

\author{
Niharika Gauraha and Ola Spjuth
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# This example is similar to the example given in the section 3.3 of the article:
# https://doi.org/10.1007/s10472-017-9541-2

#Generate train data
  mean1 = rep(0, 2)
  mean2 = rep(1, 2)
  cov2  = cov1 = .5 * diag(2)
  n = 20
  X1 = mvrnorm(n, mean1, Sigma = cov1)
  X2 = mvrnorm(n, mean2, Sigma = cov2)
  X_train = rbind(X1, X2)
  y_train = matrix(c(rep(1, n), rep(-1, n)), 2*n, 1)
# geberate privileged information data
  X1Star = matrix(0, n, 2)
  X2Star = matrix(0, n, 2)
  for(i in 1:n)
  {
    X1Star[i, 1] = norm(X1[i,] - mean1, type = "2")
    X1Star[i, 2] = norm(X2[i,] - mean2, type = "2")
  }
  for(i in 1:n)
  {
    X2Star[i, 1] = norm(X1[i, ] - mean2, type = "2")
    X2Star[i, 2] = norm(X2[i, ] - mean1, type = "2")
  }
  XStar = rbind(X1Star, X2Star)
# generate test data
  n_test = 10
  X1 = mvrnorm(n_test, mean1, Sigma = cov1)
  X2 = mvrnorm(n_test, mean2, Sigma = cov2)
  X_test = rbind(X1, X2)
  y_test = matrix(c(rep(1, n_test), rep(-1, n_test)), 2*n_test, 1)
# create instance of the class type QP
  qp = SVMP(cost = 100, gamma = .01,
            kernel_x = "rbf", gamma_x = .01,
            kernel_xstar = "rbf", gamma_xstar = .01,
            tol = .00001, svm_type = "QP")
# call the fit function
  qp$fit(X_train, XStar, y_train)
# call the predict function
  y_predict = qp$predict(X_test)
  print(length(y_predict[y_predict == y_test]))
  print("correct classification out of 20")
}

