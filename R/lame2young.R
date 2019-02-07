#' @title Young modulus and Poisson ratio from Lame parameters.
#' @description A utility function that converts Lame parameters to Young modulus and Poisson ratio.
#' @param lameLambda The material Lame first parameter in Pa, Default: 70725788000
#' @param lameMu The material Lame second parameter in Pa, Default: 46530734000
#' @return The material Young modulus in Pa and the Poisson ratio as
#'         a \code{list} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname lame2young
#' @export
lame2young <- function(lameLambda = 70.725788000e9, lameMu = 46.530734000e9) {

  ## sanity check
  if (!is.numeric(lameLambda) | (lameLambda <= 0)) {
    stop("Argument 'lameLambda' must be a positive number.")
  }
  if (!is.numeric(lameMu) | (lameMu <= 0)) {
    stop("Argument 'lameMu' must be a positive number.")
  }

  ## relationships from Disperse User's Manual, p. 181 (version 2.0.20a)
  youngModulus <- lameMu*(3*lameLambda + 2*lameMu)/(lameLambda + lameMu)
  poissonRatio <- lameLambda/2/(lameLambda + lameMu)

  ## make output
  object <- list(youngModulus = youngModulus, poissonRatio = poissonRatio)
  return(object)

}
