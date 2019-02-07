#' @title Longitudinal and transverse velocities from density, Young modulus and Poisson ratio.
#' @description A utility function that converts density, Young modulus and Poisson ratio
#'              to longitudinal and transverse velocities.
#' @param density The material density in kg/m3, Default: 4460
#' @param youngModulus The material Young modulus in Pa, Default: 121127478933
#' @param poissonRatio The material Poisson ratio, Default: 0.3015857
#' @return The material density in kg/m3, longitudinal and transverse velocities in m/s as
#'         a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname young2vel
#' @export
young2vel <- function(density = 4460,
                      youngModulus = 121.127478933e9,
                      poissonRatio = 0.3015857) {

  ## sanity check
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  if (!is.numeric(youngModulus) | (youngModulus <= 0)) {
    stop("Argument 'youngModulus' must be a positive number.")
  }
  if (!is.numeric(poissonRatio) | (poissonRatio <= 0) | (poissonRatio > 1)) {
    stop("Argument 'poissonRatio' must be a positive
         number less than or equal to 1.")
  }

  ## relationships from Disperse User's Manual, p. 181 (version 2.0.20a)
  lameLambda <- youngModulus*poissonRatio/(1 + poissonRatio)/(1 - 2*poissonRatio)
  lameMu <- youngModulus/2/(1 + poissonRatio)

  ## make output
  object <- lame2vel(density = density, lameLambda = lameLambda, lameMu = lameMu)
  return(object)
  }
