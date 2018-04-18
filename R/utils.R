# Linear kernel
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

