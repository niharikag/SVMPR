##################################################################
### QPSVMPlus class: uses solve.QP from 'quadprog' package.
### In future, LIBSVM and LIBLINEAR based faster
### implementaions are planned to be supported.
#################################################################

QPSvmPlus <- setRefClass(
  "QPSvmPlus", contains = "BaseSvmPlus",
  fields = list(
    Model = "ANY",
    X_train = "ANY",
    N_class = "numeric"
  ),
  methods = list(
    initialize = function(cost = 1, gamma = 1,
                          kernel_x = "rbf", degree_x = 3, gamma_x = .01,
                          kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .01,
                          tol = .00001)
    {
      "Costructor methods is called while creating an instance of the class."

      callSuper(cost, gamma,
                kernel_x, degree_x, gamma_x,
                kernel_xstar, degree_xstar, gamma_xstar,
                tol)
      .self$Model = NULL
      .self$X_train = NULL
      .self$N_class = 1
    },
    fit = function(X, XStar, y)
    {
      if(is.null(X) || is.null(XStar) || is.null(y))
        stop("NULL X, XStar or y")
      print("fit function")

      #QP supports: min( − q^T x + 1 / 2 x^T P x ) with the constraints  aMatrix^T x > = bVec

      .self$X_train = X
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
          .self$Gamma_xstar = 1 / ncol(X)
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

      cls_labels = unique(y)
      n_class = length(cls_labels)

      if(n_class == 2)
      {
        n_class = 1
        cls_label = c(1, -1)
      }
      .self$N_class = n_class
      .self$Model = list()

      for(i in 1:n_class)
      {
          y_temp = y

          y[y != i] = -1
          y[y == i] = 1

          P1 = cbind( (K * (y %*% t(y)) ) + (KStar / .self$Gamma),
                               (KStar / .self$Gamma) )
          P2 = cbind( (KStar / .self$Gamma), (KStar / .self$Gamma) )
          P = rbind(P1, P2)
          P <- nearPD(P)$mat

          A = cbind( matrix(1, 2 * n_samples, 1),
                              rbind(y, matrix(0,  n_samples, 1)) )
          b = matrix(c(n_samples * .self$Cost, 0), 2, 1)
          G = diag(1, 2 * n_samples)
          h = matrix(0, 2 * n_samples, 1)
          Q =  (- (.self$Cost / (2 * .self$Gamma))) * colSums(cbind((KStar + t(KStar)), (KStar + t(KStar)))) - (
             rbind(matrix(1, n_samples, 1), matrix(0, n_samples, 1)) )
          q = Q

          aMat = cbind(A, G)
          bVec = rbind( b, h)

          retValue = solve.QP(P, -q, Amat = aMat, bvec = bVec, meq = 2)
          if(any(is.na(retValue$solution)))
            stop("QP returns NaN")

          alpha = retValue$solution[1:n_samples]
          beta = retValue$solution[(n_samples+1):(2*n_samples)]

          #compute b_star first
          wxstar = (1 / .self$Gamma) * (KStar %*% (alpha + beta - .self$Cost))

          wxstar_idx = beta > .self$TOL

          if (all(wxstar_idx) == FALSE )  # no beta > tol
            b_star = max(-wxstar)
          else
            b_star = mean(-wxstar)

          sv = alpha > .self$TOL # tolerance
          sv_x = X[sv, ]
          sv_y = y[sv]

          wx = t(K %*% (alpha * y))

          temp = t(y * (1 - wxstar - b_star))
          wx = -wx + temp

          if(all(sv) == FALSE ) # no alpha > tol
          {
            lb = max(wx[y == 1])
            ub = min(wx[y == -1])
            b = (lb + ub) / 2
          }
          else
            b = mean(wx[sv])

          m = list()
          support_index = c(1:n_samples)
          m$Support_index = support_index[sv]
          m$Dual_coef = sv_y * alpha[sv]
          m$Intercept = b
          .self$Model = append(.self$Model, list(m))

          y = y_temp
      }
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

      n_test = nrow(X)
      n_samples = nrow(.self$X_train)
      K_test = matrix(0, n_test, n_samples)
      for (i in 1:n_test)
        for (j in 1:n_samples)
        {
          K_test[i, j] = kernel_method(X[i, ], .self$X_train[j, ], kernel_param)
        }

      if(.self$N_class == 1)
      {
        y_project = matrix(0, n_test, 1)
        for (i in 1:n_test)
        {
          s = 0
          for (j in 1:length(.self$Model[[1]]$Dual_coef))
          {
            a = .self$Model[[1]]$Dual_coef[j]
            sv = .self$Model[[1]]$Support_index[j]
            s = s + a * K_test[i, sv]
          }
          y_project[i] = s
        }
      }
      else{
        y_project = matrix(0, n_test, .self$N_class)
        for (cls in 1:.self$N_class) {
          for (i in 1:n_test)
          {
            s = 0
            for (j in 1:length(.self$Model[[cls]]$Dual_coef))
            {
              a = .self$Model[[cls]]$Dual_coef[j]
              sv = .self$Model[[cls]]$Support_index[j]
              s = s + a * K_test[i, sv]
            }
            y_project[i, cls] = s
          }
        }
      }
      return(y_project)
    },
    predict = function(X)
    {
      if(is.null(X))
        stop("passed NULL X")

      print("predict function")
      if(.self$N_class == 1)
      {
        y_predict = sign(.self$project(X) + .self$Model[[1]]$Intercept)
        y_predict[y_predict == 0] = 1
      }
      else
      {
        y_predict = matrix(-1, nrow(X), 1)
        y_project = .self$project(X)
        for (cls in 1:.self$N_class)
        {
          y_temp = y_project[, cls] + .self$Model[[cls]]$Intercept
          y_predict[ y_temp >= 0] = cls
        }
      }
      return(y_predict)
    }
  )
)

