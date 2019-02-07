#' @title Generating solid sphere specification.
#' @description A function that generates solid sphere specification.
#' @param diameter PARAM_DESCRIPTION, Default: 0.001
#' @param impactVelocity PARAM_DESCRIPTION, Default: 1
#' @param youngModulus The material Young modulus in Pa, Default: 74643767730
#' @param poissonRatio The material Poisson ratio, Default: 0.2285018
#' @param density The material density in kg/m3, Default: 2480
#' @return The plate specification as a \code{sphere} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeSphere
#' @export
makeSphere <- function(diameter = 0.001,
                       impactVelocity = 1,
                       youngModulus = 74.643767730e9,
                       poissonRatio = 0.2285018,
                       density = 2.48e+3) {
  ## sanity check
  if (!is.numeric(diameter) | (diameter <= 0)) {
    stop("Argument 'diameter' must be a positive number.")
  }
  if (!is.numeric(impactVelocity) | (impactVelocity <= 0)) {
    stop("Argument 'impactVelocity' must be a positive number.")
  }
  if (!is.numeric(youngModulus) | (youngModulus <= 0)) {
    stop("Argument 'youngModulus' must be a positive number.")
  }
  if (!is.numeric(poissonRatio) | (poissonRatio <= 0) | (poissonRatio > 1)) {
    stop("Argument 'poissonRatio' must be a positive
         number less than or equal to 1.")
  }
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  ## compute plate parameters
  hertzKfactor <- (1 - poissonRatio^2)/(pi*youngModulus)
  mass <- (4*pi/3)*(diameter/2)^3*density

  ## compute function output
  object <- list(diameter = diameter,
                 impactVelocity = impactVelocity,
                 youngModulus = youngModulus,
                 poissonRatio = poissonRatio,
                 density = density,
                 hertzKfactor = hertzKfactor,
                 mass = mass)
  class(object) <- "sphere"
  return(object)
}
