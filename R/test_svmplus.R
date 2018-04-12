  source("R/QPSvmPlus.R")
  qp = QPSvmPlus$new(cost = 1, gamma = 1,
                     kernel_x = "rbf", degree_x = 3, gamma_x = .0019,
                     kernel_xstar = "rbf", degree_xstar = 3, gamma_xstar = .0568,
                     tol = .00001)
  #X = 2*diag(3)
  #X_test = 3*diag(3)

  #X1 = rbind(c(3,1), c(3,-1), c(6,1), c(6,-1))
  #X2 = rbind(c(1,0), c(0,1), c(0,-1), c(-1,0))
  #X = rbind(X1,X2)
  #y = c(1,1,1,1,-1,-1,-1,-1)


  X_train = rbind(c(17, 24,  1,  8, 15),
                      c(23,  5,  7, 14, 16),
                      c(4,  6, 13, 20, 22),
                      c(10, 12, 19, 21,  3),
                      c(11, 18, 25,  2,  9))

  y_train = c(1, 1, -1, 1, -1)

  XStar = diag(5)

  XStar[1, 1] = 1
  XStar[2, 2] = 2
  XStar[3, 3] = 3
  XStar[4, 4] = 4
  XStar[5,5] = 5

  print(XStar)

  qp$fit(X_train, XStar, y_train)
  #qp$predict(X_test)
