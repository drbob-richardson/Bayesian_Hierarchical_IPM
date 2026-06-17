# Kernel functionals of a fitted IPM: asymptotic growth rate (dominant
# eigenvalue), stable size distribution (dominant right eigenvector),
# reproductive value (dominant left eigenvector), and the sensitivity /
# elasticity matrices of lambda. Operates on the discrete projection matrix
# A = delta * K.

ipm_kernel_functionals <- function(projection_matrix) {
  A <- projection_matrix
  right <- eigen(A)
  i <- which.max(Re(right$values))
  lambda <- Re(right$values[i])
  w <- Re(right$vectors[, i]); w <- abs(w); w <- w / sum(w)        # stable distribution
  left <- eigen(t(A))
  j <- which.max(Re(left$values))
  v <- Re(left$vectors[, j]); v <- abs(v)
  v <- v / sum(v * w)                                              # scale so <v,w> = 1
  sensitivity <- outer(v, w)                                       # d lambda / d A_ij
  elasticity <- sensitivity * A / lambda                          # proportional sensitivity
  list(
    lambda = lambda,
    stable_distribution = w,
    reproductive_value = v,
    sensitivity = sensitivity,
    elasticity = elasticity
  )
}

# Convenience wrapper that builds the projection matrix from inverse-IPM
# parameters and returns its functionals.
inverse_ipm_functionals <- function(parameter, mesh, kernel_builder = build_centered_ipm) {
  ipm <- kernel_builder(parameter, mesh)
  ipm_kernel_functionals(ipm$delta * ipm$kernel)
}
