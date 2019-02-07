#' @title Lame parameters from longitudinal and transverse velocities and density.
#' @description A utility function that converts longitudinal and transverse velocities and density to Lame parameters.
#' @param density The material density in kg/m3, Default: 4460
#' @param longVelocity The material longitudinal velocity in m/s, Default: 6060
#' @param tranVelocity The material transverse velocity in m/s, Default: 3230
#' @return The material first and second Lame parameters in Pa as a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname vel2lame
#' @export
vel2lame <- function(density = 4460, longVelocity = 6060, tranVelocity = 3230) {

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
  if (((longVelocity^2 - 2*tranVelocity^2) <= 0)) {
    stop("Parameter Lame lambda must be a positive number. Check velocity values!")
  }
  lameLambda <- density*(longVelocity^2 - 2*tranVelocity^2)
  lameMu <- density*tranVelocity^2

  ## make output
  object <- list(lameLambda = lameLambda, lameMu = lameMu)
  return(object)

}
