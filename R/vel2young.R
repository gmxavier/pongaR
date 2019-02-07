#' @title Young modulus and Poisson ratio from longitudinal and transverse velocities and density.
#' @description A utility function that converts longitudinal and transverse velocities and density
#'              to Young modulus and Poisson ratio.
#' @param density The material density in kg/m3, Default: 4460
#' @param longVelocity The material longitudinal velocity in m/s, Default: 6060
#' @param tranVelocity The material transverse velocity in m/s, Default: 3230
#' @return The material Young modulus in Pa and Poisson ratio as a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname vel2young
#' @export
vel2young <- function(density = 4460, longVelocity = 6060, tranVelocity = 3230) {

  ## sanity check
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  if (!is.numeric(longVelocity) | (longVelocity <= 0)) {
    stop("Argument 'longVelocity' must be a positive number.")
  }
  if (!is.numeric(tranVelocity) | (tranVelocity <= 0)) {
    stop("Argument 'tranVelocity' must be a positive number.")
  }

  ## relationships from Disperse User's Manual, p. 181 (version 2.0.20a)
  if (((3*longVelocity^2 - 4*tranVelocity^2) <= 0)) {
    stop("Young's modulus must be a positive number. Check velocity values!")
  }
  youngModulus <- density*tranVelocity^2*(3*longVelocity^2 -
                                            4*tranVelocity^2)/(longVelocity^2 -
                                                                 tranVelocity^2)
  poissonRatio <- (1/2)*(longVelocity^2 - 2*tranVelocity^2)/(longVelocity^2 -
                                                               tranVelocity^2)

  ## make output
  object <- list(youngModulus = youngModulus, poissonRatio = poissonRatio)
  return(object)

}
