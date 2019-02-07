#' @title Longitudinal and transverse velocities from Lame parameters.
#' @description A utility function that converts Lame parameters to longitudinal and transverse velocities.
#' @param density The material density in kg/m3, Default: 4460
#' @param lameLambda The material Lame first parameter in Pa, Default: 70725788000
#' @param lameMu The material Lame second parameter in Pa, Default: 46530734000
#' @return The material density in kg/m3, longitudinal and transverse velocities in m/s as
#'         a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname lame2vel
#' @export
lame2vel <- function(density = 4460,
                     lameLambda = 70.725788000e9,
                     lameMu = 46.530734000e9) {

  ## sanity check
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  if (!is.numeric(lameLambda) | (lameLambda <= 0)) {
    stop("Argument 'lameLambda' must be a positive number.")
  }
  if (!is.numeric(lameMu) | (lameMu <= 0)) {
    stop("Argument 'lameMu' must be a positive number.")
  }

  ## relationships from Disperse User's Manual, p. 181 (version 2.0.20a)
  longVelocity <- sqrt((lameLambda + 2*lameMu)/density)
  tranVelocity <- sqrt(lameMu/density)

  ## make output
  object <- list(density = density,
                 longVelocity = longVelocity,
                 tranVelocity = tranVelocity)
  return(object)

}
