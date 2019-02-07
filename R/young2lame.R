#' @title Lame parameters from Young modulus and Poisson ratio.
#' @description A utility function that converts Young modulus and Poisson ratio to Lame parameters.
#' @param youngModulus The material Young modulus in Pa, Default: 121127478933
#' @param poissonRatio The material Poisson ratio, Default: 0.3015857
#' @return The material first and second Lame parameters in Pa as a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname young2lame
#' @export
young2lame <- function(youngModulus = 121.127478933e9,
                       poissonRatio = 0.3015857) {

  ## sanity check
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
  object <- list(lameLambda = lameLambda, lameMu = lameMu)
  return(object)
  }
