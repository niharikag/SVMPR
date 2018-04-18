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
    Dual_coef = "ANY"
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
      .self$Kernel_xstar <<- kernel_xstar
      .self$Degree_xstar <<- degree_xstar
      .self$Gamma_xstar <<- gamma_xstar
      .self$TOL <<- tol
    },
    fit = function(X, XStar, y)
    {
      if(is.null(X) || is.null(XStar) || is.null(y))
        stop("NULL X, XStar or y")

      print("fit function")
    },
    predict = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")

      print("predict function")
    },
    project = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")

      print("project function")
    }
  )
)



