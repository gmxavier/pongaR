#' @title Generating plate specification.
#' @description A function that generates plate specification.
#' @param length Plate length in m, Default: 0.508
#' @param width Plate width in m, Default: 0.381
#' @param materialIsotropic Logical, if \code{TRUE}, the material is isotropic, Default: TRUE
#' @param thickness Plate thickness in m, Default: ifelse(materialIsotropic, 0.005, 0.00254)
#' @param youngModulus The material Young modulus in Pa, Default: 121127478933
#' @param poissonRatio The material Poisson ratio, Default: 0.3015857
#' @param density The material density in kg/m3, Default: ifelse(materialIsotropic, 4460, 1550)
#' @param flexuralRigidity11 The material flexural rigidity 11 in Pa m3, Default: 198.7
#' @param flexuralRigidity22 The material flexural rigidity 22 in Pa m3, Default: 13.23
#' @param flexuralRigidity12 The material flexural rigidity 12 in Pa m3, Default: 3.98
#' @param flexuralRigidity66 The material flexural rigidity 66 in Pa m3, Default: 8.15
#' @param transverseShearmodulus44 The material transverse shear modulus 44 in MPa m, Default: 7940000
#' @param transverseShearmodulus55 The material transverse shear modulus 55 in MPa m, Default: 12640000
#' @param lossFactor The material loss factor, Default: 0
#' @param shearCorrectionfactor The material shear correction factor, Default: NA
#' @return The plate specification as a \code{plate} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makePlate
#' @export
makePlate <- function(length = 0.508,
                      width = 0.381,
                      materialIsotropic = TRUE,
                      thickness = ifelse(materialIsotropic, 5e-3, 2.54e-3),
                      youngModulus = 121.127478933e9,
                      poissonRatio = 0.3015857,
                      density = ifelse(materialIsotropic, 4.46e+3, 1.55e+3),
                      flexuralRigidity11 = 198.7,
                      flexuralRigidity22 = 13.23,
                      flexuralRigidity12 = 3.98,
                      flexuralRigidity66 = 8.15,
                      transverseShearmodulus44 = 7.94e6,
                      transverseShearmodulus55 = 12.64e6,
                      lossFactor = 0,
                      shearCorrectionfactor = NA) {
  ## sanity check
  if (!is.numeric(length) | (length <= 0)) {
    stop("Argument 'length' must be a positive number.")
  }
  if (!is.numeric(width) | (width <= 0)) {
    stop("Argument 'width' must be a positive number.")
  }
  if (materialIsotropic != TRUE & materialIsotropic != FALSE) {
    stop("Argument 'materialIsotropic' must be TRUE or FALSE.")
  }
  if (!is.numeric(thickness) | (thickness <= 0)) {
    stop("Argument 'thickness' must be a positive number.")
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
  if (!is.numeric(flexuralRigidity11) | (flexuralRigidity11 <= 0)) {
    stop("Argument 'flexuralRigidity11' must be a positive number.")
  }
  if (!is.numeric(flexuralRigidity22) | (flexuralRigidity22 <= 0)) {
    stop("Argument 'flexuralRigidity22' must be a positive number.")
  }
  if (!is.numeric(flexuralRigidity12) | (flexuralRigidity12 <= 0)) {
    stop("Argument 'flexuralRigidity12' must be a positive number.")
  }
  if (!is.numeric(flexuralRigidity66) | (flexuralRigidity66 <= 0)) {
    stop("Argument 'flexuralRigidity66' must be a positive number.")
  }
  if (!is.numeric(transverseShearmodulus44) | (transverseShearmodulus44 <= 0)) {
    stop("Argument 'transverseShearmodulus44' must be a positive number.")
  }
  if (!is.numeric(transverseShearmodulus55) | (transverseShearmodulus55 <= 0)) {
    stop("Argument 'transverseShearmodulus55' must be a positive number.")
  }
  if (!is.numeric(lossFactor) | (lossFactor < 0) | (lossFactor > 1)) {
    stop("Argument 'lossFactor' must be a non negative
         number less than or equal to 1.")
  }
  if (!is.na(shearCorrectionfactor)) {
    if (!is.numeric(shearCorrectionfactor) | (shearCorrectionfactor <= 0)) {
      stop("Argument 'shearCorrectionfactor' must be a positive number.")
    }
  }

  ## compute shear correction factor
  if (is.na(shearCorrectionfactor)) {
    ## compute shear correction factor (R. D. Mindlin, J. Appl. Mech. 18, 38 (1951))
    shearCorrectionfactor <- uniroot(f = shearEquation,
                                     interval = c(1e-3, 1),
                                     poissonRatio = poissonRatio)$root
  }

  ## compute other plate parameters
  flexuralRigidity <- youngModulus*thickness^3/(12*(1 - poissonRatio^2))
  shearModulus <- youngModulus/(2*(1 + poissonRatio))
  rotaryInertia <- density*thickness^3/12
  if (materialIsotropic == TRUE) {
    transverseShearmodulus44 <- shearCorrectionfactor^2*shearModulus*thickness
    transverseShearmodulus55 <- transverseShearmodulus44
    flexuralRigidity11 <- flexuralRigidity
    flexuralRigidity12 <- poissonRatio*flexuralRigidity
    flexuralRigidity22 <- flexuralRigidity
    flexuralRigidity66 <- ((1 - poissonRatio)/2)*flexuralRigidity
  }
  hertzKfactor <- (1 - poissonRatio^2)/(pi*youngModulus)
  mass <- density*length*width*thickness
  impedance <- 8*(flexuralRigidity*density*thickness)^(1/2)
  longVelocity <- sqrt(youngModulus*(1 - poissonRatio)/(1 + poissonRatio)/(1 - 2*poissonRatio)/density)
  acousticImpedance <- density*longVelocity

  ## compute function output
  object <- list(length = length,
                 width = width,
                 thickness = thickness,
                 materialIsotropic = materialIsotropic,
                 youngModulus = youngModulus,
                 poissonRatio = poissonRatio,
                 density = density,
                 lossFactor = lossFactor,
                 shearCorrectionfactor = shearCorrectionfactor,
                 flexuralRigidity11 = flexuralRigidity11,
                 flexuralRigidity12 = flexuralRigidity12,
                 flexuralRigidity22 = flexuralRigidity22,
                 flexuralRigidity66 = flexuralRigidity66,
                 shearModulus = shearModulus,
                 rotaryInertia = rotaryInertia,
                 transverseShearmodulus44 = transverseShearmodulus44,
                 transverseShearmodulus55 = transverseShearmodulus55,
                 hertzKfactor = hertzKfactor,
                 mass = mass,
                 impedance = impedance,
                 acousticImpedance = acousticImpedance)
  class(object) <- "plate"
  return(object)
}
