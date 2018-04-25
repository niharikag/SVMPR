#' svmPlus for classification problems
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
SVMP = function(cost = 1, gamma = 1,
                   kernel_x = "rbf", degree_x = 3, gamma_x = .001,
                   kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .001,
                   tol = .00001, svm_type = "QP")
{
  if(svm_type != "QP")
    stop("Currently SVMP supports only QP based optimization")

  return( QPSvmPlus$new(cost, gamma,
                     kernel_x, degree_x, gamma_x,
                     kernel_xstar, degree_xstar, gamma_xstar,
                     tol) )

}
