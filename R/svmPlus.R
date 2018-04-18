#' Implementation of SVM plus for Classification
#'
#' The svmplus package implements svm+ for classification problems.
#'
#' In SVM+, privileged information is used to estimate a linear model
#' of the slack variables.
#' Currently, the pakcage is built upon quadprog method. We plan to extend package to use other
#' optimization algorithms, (i.e. CVXR, LibSVM and LibLinear) for optimization,
#' and also we plan to support for multi-class classification problems.
#'
#' @docType package
#' @name conformalClassification
"_PACKAGE"

#' svmPlus for binary classification problems
#' @param cost cost of constraints violation
#' @param gamma parameter needed for priviledged information
#' @param kernel_x the kernel used for standard training data
#' @param degree_x parameter needed for polynomial kernel for training data
#' @param gamma_x parameter needed for rbf kernel for training data
#' @param kernel_xstar the kernel used for priviledged information (PI)
#' @param degree_xstar parameter needed for polynomial kernel for PI
#' @param gamma_xstar parameter needed for rbf kernel for PI
#' @param tol tolerance of dual variables
#' @return instance of the class QPSvmPlus
#' @export
svmPlus = function(cost = 1, gamma = 1,
                   kernel_x = "rbf", degree_x = 3, gamma_x = .001,
                   kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .001,
                   tol = .00001)
{
  return( QPSvmPlus$new(cost = 1, gamma = 1,
                     kernel_x = "rbf", degree_x = 3, gamma_x = .001,
                     kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .001,
                     tol = .00001) )

}
