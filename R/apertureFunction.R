#' @title Computing sensor's aperture function.
#' @description A utility function that computes sensor's aperture function.
#' @param x A distance from the center along the sensor's radius, in m.
#' @param gamma_m A parameter from Midlin's plate theory.
#' @param alpha_n A parameter from Midlin's plate theory.
#' @param x0 Sensor position coordinate \code{x0}, in m.
#' @param y0 Sensor position coordinate \code{y0}, in m.
#' @param receiverRadius2 Squared sensor radius, in m2.
#' @return The value of the sensor's aperture function at a distance \code{x} from the center.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname apertureFunction
#' @export
apertureFunction <- function(x, gamma_m, alpha_n, x0, y0, receiverRadius2) {
  ##
  sin(alpha_n*x)*(cos(gamma_m*(y0 - sqrt(receiverRadius2 - (x - x0)^2
                                         + .Machine$double.eps))) -
                    cos(gamma_m*(y0 + sqrt(receiverRadius2 - (x - x0)^2
                                           + .Machine$double.eps))))
}
