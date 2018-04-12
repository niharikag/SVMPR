#############################################
### Base SVMPlus class
#############################################

.BaseSvmPlus <- setRefClass(
  "BaseSvmPlus",
  fields = list(
    Cost = "numeric",
    Gamma = "numeric",
    Kernel_x = "character",
    Degree_x = "numeric",
    Gamma_x = "numeric",
    Kernel_xstar = "character",
    Degree_xstar = "numeric",
    Gamma_xstar = "numeric",
    TOL = "numeric",
    Support_vectors = "ANY",
    Support_y = "ANY",
    Dual_coef = "ANY",
    Intercept = "numeric"
    #data.X = "ANY",
    #data.XStar = "ANY",
    #data.y = "ANY"
  ),
  methods = list(
    initialize = function(cost, gamma,
                          kernel_x, degree_x, gamma_x,
                          kernel_xstar, degree_xstar, gamma_xstar,
                          tol)
    {
      "This method is called when you create an instance of the class."
      .self$Cost <<- cost
      .self$Gamma <<- gamma
      .self$Kernel_x <<- kernel_x
      .self$Degree_x <<- degree_x
      .self$Gamma_x <<- gamma_x
      .self$Kernel_xstar <<- kernel_x
      .self$Degree_xstar <<- degree_x
      .self$Gamma_xstar <<- gamma_x
      .self$TOL <<- tol
    },
    fit = function(X, XStar, y)
    {
      if(is.null(X) || is.null(XStar) || is.null(y))
        stop("NULL X, XStar or y")
      #if(model$modelType != "Classification" )

      print("fit function")
    },
    predict = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")
      #if(model$modelType != "Classification" )

      print("predict function")
    },
    project = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")
      #if(model$modelType != "Classification" )

      print("project function")
    },
    decision_function = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")
      #if(model$modelType != "Classification" )

      print("decision function")
    }
  )
)

#static methods
linear_kernel <- function(x1, x2, param = None)
{
  return (t(x1) %*% x2)
}


# Polynomial kernel
poly_kernel <- function(x1, x2, param = 2)
{
  return (1 +  (t(x1) %*% x2)) ^ param
}



# Radial basis kernel
rbf_kernel <- function(x1, x2, param = 1)
{
  return (exp(-(norm((x1 - x2), type="2") ^ 2) * param))
}



