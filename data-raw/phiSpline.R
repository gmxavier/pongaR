#' @title Computing the phi function spline.
#' @description A utility function that computes the phi function spline used in Hertz's contact model.
#' @return The phi function produced by \code{splinefun}.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname phiSpline
#' @export
phiSpline <- function() {
  ## compute phi function using splinefun
  aux <- seq(from = 0, to = 1, by = 1e-6)
  phiIntegral <- 0*aux
  k <- 1
  for (phi in aux) {
    phiIntegral[k] <- integrate(function(x) (1-x^(5/2))^(-1/2),
                                lower = 0, upper = phi)$value
    k <- k + 1
  }
  y <- c(aux, rev(aux)[-1])
  x <- c(phiIntegral, -rev(phiIntegral)[-1] + 2*rev(phiIntegral)[1])
  phi <- splinefun(x, y, method = "monoH.FC")
  usethis::use_data(phi, internal = TRUE, compress = "xz")
}
