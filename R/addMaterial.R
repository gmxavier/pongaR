#' @title Adding new material data to material database.
#' @description A utility function that adds new material data to \code{pongaR} material database.
#' @param name The material name.
#' @param reference The material data reference.
#' @param materialIsotropic Logical, \code{TRUE} means material is isotropic.
#' @param youngModulus The material Young modulus in Pa, Default: NA
#' @param poissonRatio The material Poisson ratio, Default: NA
#' @param density The material density in kg/m3
#' @param flexuralRigidity11 The material flexural rigidity 11 in Pa m3, Default: NA
#' @param flexuralRigidity22 The material flexural rigidity 22 in Pa m3, Default: NA
#' @param flexuralRigidity12 The material flexural rigidity 12 in Pa m3, Default: NA
#' @param flexuralRigidity66 The material flexural rigidity 66 in Pa m3, Default: NA
#' @param transverseShearmodulus44 The material transverse shear modulus 44 in MPa m, Default: NA
#' @param transverseShearmodulus55 The material transverse shear modulus 55 in MPa m, Default: NA
#' @return Adds the new material data to the material database embedded into the \code{pongaR} library.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname addMaterial
#' @export
addMaterial <- function(name,
                        reference,
                        materialIsotropic,
                        youngModulus = NA,
                        poissonRatio = NA,
                        density,
                        flexuralRigidity11 = NA,
                        flexuralRigidity22 = NA,
                        flexuralRigidity12 = NA,
                        flexuralRigidity66 = NA,
                        transverseShearmodulus44 = NA,
                        transverseShearmodulus55 = NA) {

  ## sanity check
  if (!is.character(name)) {
    stop("Argument 'name' must be a character string.")
  }
  if (!is.character(reference)) {
    stop("Argument 'reference' must be a character string.")
  }
  if (materialIsotropic != TRUE & materialIsotropic != FALSE) {
    stop("Argument 'materialIsotropic' must be TRUE or FALSE.")
  }
  if (materialIsotropic == TRUE) {
    if (!is.numeric(youngModulus) | (youngModulus <= 0)) {
      stop("Argument 'youngModulus' must be a positive number.")
    }
    if (!is.numeric(poissonRatio) | (poissonRatio <= 0) | (poissonRatio > 1)) {
      stop("Argument 'poissonRatio' must be a positive
           number less than or equal to 1.")
    }
    }
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  if (materialIsotropic == FALSE) {
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
  }

  data("material", envir = environment())
  aux <- data.frame(name = name,
                    reference = reference,
                    materialIsotropic = materialIsotropic,
                    youngModulus = youngModulus,
                    poissonRatio = poissonRatio,
                    density = density,
                    flexuralRigidity11 = flexuralRigidity11,
                    flexuralRigidity22 = flexuralRigidity22,
                    flexuralRigidity12 = flexuralRigidity12,
                    flexuralRigidity66 = flexuralRigidity66,
                    transverseShearmodulus44 = transverseShearmodulus44,
                    transverseShearmodulus55 = transverseShearmodulus55,
                    stringsAsFactors = FALSE)
  material <- rbind(material, aux)
  save(material, file = "data/material.rda")

}
