# Local diagnostics for weak and confounded parameter directions.

finite_difference_jacobian <- function(
  parameter,
  forward,
  relative_step = 1e-5
) {
  baseline <- forward(parameter)
  stopifnot(is.numeric(baseline), is.numeric(parameter))

  jacobian <- vapply(seq_along(parameter), function(index) {
    step <- relative_step * max(1, abs(parameter[index]))
    perturbed <- parameter
    perturbed[index] <- perturbed[index] + step
    (forward(perturbed) - baseline) / step
  }, numeric(length(baseline)))

  colnames(jacobian) <- names(parameter)
  jacobian
}

sensitivity_svd <- function(jacobian, relative_tolerance = 1e-8) {
  decomposition <- svd(jacobian)
  relative_singular_values <- decomposition$d / max(decomposition$d)
  weak_directions <- decomposition$v[
    , order(relative_singular_values, decreasing = FALSE),
    drop = FALSE
  ]
  rownames(weak_directions) <- colnames(jacobian)

  list(
    singular_values = decomposition$d,
    relative_singular_values = relative_singular_values,
    numerical_rank = sum(relative_singular_values > relative_tolerance),
    condition_number = max(decomposition$d) / min(decomposition$d),
    weak_directions = weak_directions
  )
}
