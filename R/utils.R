#' Does a one-sided formula contain no terms?
#'
#' `identical(f, ~ 1)` fails because formulas carry their environment, so a
#' `~ 1` typed by the user is not identical to one typed in a function default.
#' This compares the terms instead: `~ 1`, `~1` and `~ 0` all return `TRUE`.
#'
#' @param f A formula.
#' @return Logical.
#' @keywords internal
#' @importFrom stats terms
#' @noRd
is_null_formula <- function(f) {
  length(attr(stats::terms(f), "term.labels")) == 0
}
