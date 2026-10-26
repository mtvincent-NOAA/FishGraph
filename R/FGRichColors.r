#' Rich color function to determine color based on the number of categories.
#' Copied from r4SS written by  by Arni Magnusson and Ian Taylor
#'
#' @param n number of categories to output colors
#' @param alpha transparenc level for color
#'
#' @return rich.vector a vector of color values
#'
#' @author Arni Magnusson
#' @author Ian Taylor
#' @author Matt Vincent
#'
#'
#' @examples \donttest{
#' rich.colors.
#' }
#' @export
#'
#'
FGRichColors <- function(n, alpha = 1) {
    x <- seq(0, 1, length = n+1)[-1]
    r <- 1 / (1 + exp(20 - 35 * x))
    g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
    b <- dnorm(x, 0.25, 0.15) / max(dnorm(x, 0.25, 0.15))
    rgb.m <- matrix(c(r, g, b), ncol = 3)
    rich.vector <- apply(rgb.m, 1, function(v) {
        rgb(v[1], v[2], v[3], alpha = alpha)
    })
}
