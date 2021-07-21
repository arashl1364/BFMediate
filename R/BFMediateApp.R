#' @include utils.R
NULL

#' Run BFMediateApp
#'
#' Launches a Shiny app that shows a demo of what can be done with
#' the \code{BFMediate} package.
#'
#' @export
#'
#'
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   run_BFMediate()
#' }
run_BFMediate <- function() {
  appDir <- system.file("BFMediateApp", package = "BFMediate")
  if (appDir == "") {
    stop("Could not find shiny app directory. Try re-installing `BFMediate`.", call. = FALSE)
  }

  if (requireNamespace("shiny", quietly = TRUE)) {
    shiny::runApp(appDir, display.mode = "normal")
  } else {
    stop("Package 'shiny' is required but not installed on your system.", call. = FALSE)
  }
}
