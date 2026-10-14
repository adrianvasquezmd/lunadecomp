# Reserve temporary columns without changing user-supplied names.
.ld_private_names <- function(occupied, proposed) {
  utils::tail(make.unique(c(occupied, proposed), sep = "_"), length(proposed))
}

# Invert information in rescaled coordinates, then restore coefficient units.
.ld_information_inverse <- function(model_matrix, hessian_weights) {
  column_scale <- apply(abs(model_matrix), 2, max)
  column_scale[column_scale == 0] <- 1
  scaled_matrix <- sweep(model_matrix, 2, column_scale, "/")
  information <- crossprod(scaled_matrix, scaled_matrix * hessian_weights)
  inverse <- tryCatch(
    solve(information),
    error = function(e) MASS::ginv(information)
  )
  inverse / outer(column_scale, column_scale)
}
