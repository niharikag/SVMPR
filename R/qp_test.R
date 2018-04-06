# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(CVXR)
library(quadprog)


Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
solve.QP(Dmat,dvec,Amat,bvec=bvec)

#min( − q^T x + 1 / 2 x^T P x ) with the constraints  aMatrix^T x > = bVec

P = 2* rbind(c(1.0, 0.),c(0., 4.))
q = c(-8.0, -16.)

aMatrix = rbind(c(-1.0, 1., -1., 0.),
                  c(1., 0., 0., -1.))
bVec = c(5.0, 3., 0., 0.)
sol = solve.QP(P, -q, Amat = -aMatrix, bvec = -bVec)
print(sol)
