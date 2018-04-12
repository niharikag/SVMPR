library(quadprog)
source("R/BaseSvmPlus.R")
#############################################
### QPSVMPlus class: optimization based
### on quadprog
#############################################

QPSvmPlus <- setRefClass(
  "QPSvmPlus", contains = "BaseSvmPlus",
  fields = list(
    QP = "ANY"
  ),
  methods = list(
    initialize = function(cost = 1, gamma = 1,
                          kernel_x = "rbf", degree_x = 3, gamma_x = .01,
                          kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .01,
                          tol = .00001)
    {
      "This method is called when you create an instance of the class."

      #.BaseSvmPlus$initialize(cost, gamma,
      #                       kernel_x, degree_x, gamma_x,
      #                       kernel_xstar, degree_xstar, gamma_xstar,
      #                       tol)
      .self$Cost <<- cost
      .self$Gamma <<- gamma
      .self$Kernel_x <<- kernel_x
      .self$Degree_x <<- degree_x
      .self$Gamma_x <<- gamma_x
      .self$Kernel_xstar <<- kernel_xstar
      .self$Degree_xstar <<- degree_xstar
      .self$Gamma_xstar <<- gamma_xstar
      .self$TOL <<- tol
      .self$Support_vectors = NULL
      .self$Support_y = NULL
      .self$Dual_coef = NULL
      .self$Intercept = 0
    },
    fit = function(X, XStar, y)
    {
      if(is.null(X) || is.null(XStar) || is.null(y))
        stop("NULL X, XStar or y")
      print("fit function")

      #QP supports: min( − q^T x + 1 / 2 x^T P x ) with the constraints  aMatrix^T x > = bVec

      n_samples = nrow(X)
      n_features = ncol(X)


      if (.self$Kernel_x == "linear")
      {
        kernel_method = linear_kernel
        kernel_param = NULL
      }
      else if (.self$Kernel_x == "poly")
      {
        kernel_method = poly_kernel
        kernel_param = .self$Degree_x
      }
      else
      {
        kernel_method = rbf_kernel
        if (.self$Gamma_x == 'auto')
          .self$Gamma_x = 1 / n_features
        kernel_param = .self$Gamma_x
      }

      if (.self$Kernel_xstar == "linear")
      {
        kernel_method_star = linear_kernel
        kernel_param_star = NULL
      }
      else if (.self$Kernel_xstar == "poly")
      {
        kernel_method_star = poly_kernel
        kernel_param_star = .self$Degree_xstar
      }
      else
      {
        kernel_method_star = rbf_kernel
        if (.self$Gamma_xstar == 'auto')
          .self$Gamma_xstar = 1 / XStar.shape[1]
        kernel_param_star = .self$Gamma_xstar
      }
      # compute the matrix K and KStar (n_samples X n_samples) using kernel function
      K = matrix(0, n_samples, n_samples)
      KStar = matrix(0, n_samples, n_samples)

      for (i in 1:n_samples)
        for (j in 1:n_samples)
        {
          K[i, j] = kernel_method(X[i, ], X[j, ], kernel_param)
          KStar[i, j] = kernel_method_star(XStar[i, ], XStar[j, ], kernel_param_star)
        }
      print("K")
      print(K)
      print("K Star")
      print(KStar)

      P1 = cbind( (K * (y %*% t(y)) ) + (KStar / .self$Gamma),
                           (KStar / .self$Gamma) )
      P2 = cbind( (KStar / .self$Gamma), (KStar / .self$Gamma) )
      P = rbind(P1, P2)
      print("matrix P")
      print(P)

      A = rbind( matrix(1, 1, 2 * n_samples),
                          cbind( t(y), matrix(0, 1, n_samples)) )
      b = matrix(c(n_samples * .self$Cost, 0), 2, 1)
      G = diag(1, 2 * n_samples)
      #print(G)
      h = matrix(0, 1, 2 * n_samples)
      Q =  (- (.self$Cost / (2 * .self$Gamma))) * colSums(cbind((KStar + t(KStar)), (KStar + t(KStar)))) - (
         rbind(matrix(1, n_samples, 1), matrix(0, n_samples, 1)) )
      q = Q

      print(" matrixq")
      print(q)

      aMat = cbind(t(A), G)
      bVec = rbind( b, t(h))

      print(" matrix aMat")
      print(aMat)

      print(" matrix bVex")
      print(bVec)
      #print(q)
      #print(G)
      #print(h)
      #print(A)
      #print(b)
      retValue = solve.QP(P, -q, Amat = aMat, bvec = bVec, meq = 2)
      print(retValue$solution)
      alpha = retValue$solution[1:n_samples]
      beta = retValue$solution[(n_samples+1):(2*n_samples)]
      print("alpha")
      print(alpha)
      print("beta")
      print(beta)
      #compute b_star first
      wxstar = (1 / .self$Gamma) * (KStar %*% (alpha + beta - .self$Cost))
      print("wxstar")
      print(wxstar)
      wxstar_idx = beta > .self$TOL
      print("wxstar index")
      print(wxstar_idx)

      if (all(wxstar_idx) == FALSE )  # no beta > tol
        b_star = max(-wxstar)
      else
        b_star = mean(-wxstar)

      print("b_star")
      print(b_star)

      sv = alpha > .self$TOL # tolerance
      sv_x = X[sv, ]
      sv_y = y[sv]

      wx = t(K * (alpha * y))
      print("matrix K * (alpha * y)")
      print(wx)
      print("y*(1 - wxstar - b_star)")
      print(y*(1 - wxstar - b_star))

      temp = t(y * (1 - wxstar - b_star))
      wx = -wx + matrix(rep(temp, n_samples),
                        ncol =  ncol(temp) , byrow = TRUE )
      print("matrix wx")
      print(wx)

      if(all(sv) == FALSE ) # no alpha > tol
      {
        lb = max(wx[y == 1])
        ub = min(wx[y == -1])
        b = (lb + ub) / 2
      }
      else
        b = mean(wx[sv])

      print("bias b")
      print(b)

      .self$Support_vectors = sv_x
      .self$Support_y = sv_y
      .self$Dual_coef = alpha[sv]
      .self$Intercept = b
      print(.self$Dual_coef)
    },
    project = function(X)
    {
      print("project function")
      if(is.null(X))
        stop("passed NULL X")

      if (.self$Kernel_x == "linear")
      {
        kernel_method = linear_kernel
        kernel_param = NULL
      }
      else if (.self$Kernel_x == "poly")
      {
        kernel_method = poly_kernel
        kernel_param = .self$Degree_x
      }
      else
      {
        kernel_method = rbf_kernel
        kernel_param = .self$Gamma_x
      }

      y_predict = matrix(0, nrow(X))
      for (i in 1:nrow(X))
      {
        s = 0
        for (j in 1:length(.self$Dual_coef))
        {
          a = .self$Dual_coef[j]
          sv_y = .self$Support_y[j]
          sv = .self$Support_vectors[j, ]
          s = s + a * sv_y * kernel_method(X[i, ], sv, kernel_param)
        }
        y_predict[i] = s
      }

      print(y_predict)
      return(y_predict)
    },
    predict = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")

      print("predict function")
      return(sign(.self$project(X) + .self$Intercept))
    },
    decision_function = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")
      print("decision function")

      return(.self$project(X) + .self$Intercept)
    }
  )
)


