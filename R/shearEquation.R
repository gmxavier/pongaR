#' @title Computing shear equation .
#' @description A utility function that computes shear equation (R. D. Mindlin, J. Appl. Mech. 18, 38 (1951)).
#' @param shearCorrectionfactor The material Midlin's plate theory shear correction factor.
#' @param poissonRatio The material Poisson ratio.
#' @return The value of the shear equation used to calculate the material Midlin's
#'         plate theory shear correction factor.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname shearEquation
#' @export
shearEquation <- function(shearCorrectionfactor, poissonRatio) {
  ##
  alpha <- (1 - 2*poissonRatio)/(2 - 2*poissonRatio)
  return(4*sqrt((1 - alpha*shearCorrectionfactor^2)*(1 - shearCorrectionfactor^2)) -
           (2 - shearCorrectionfactor^2)^2)
}
