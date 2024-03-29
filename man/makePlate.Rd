% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePlate.R
\name{makePlate}
\alias{makePlate}
\title{Generating plate specification.}
\usage{
makePlate(length = 0.508, width = 0.381, materialIsotropic = TRUE,
  thickness = ifelse(materialIsotropic, 0.005, 0.00254),
  youngModulus = 121127478933, poissonRatio = 0.3015857,
  density = ifelse(materialIsotropic, 4460, 1550),
  flexuralRigidity11 = 198.7, flexuralRigidity22 = 13.23,
  flexuralRigidity12 = 3.98, flexuralRigidity66 = 8.15,
  transverseShearmodulus44 = 7940000,
  transverseShearmodulus55 = 12640000, lossFactor = 0,
  shearCorrectionfactor = NA)
}
\arguments{
\item{length}{Plate length in m, Default: 0.508}

\item{width}{Plate width in m, Default: 0.381}

\item{materialIsotropic}{Logical, if \code{TRUE}, the material is isotropic, Default: TRUE}

\item{thickness}{Plate thickness in m, Default: ifelse(materialIsotropic, 0.005, 0.00254)}

\item{youngModulus}{The material Young modulus in Pa, Default: 121127478933}

\item{poissonRatio}{The material Poisson ratio, Default: 0.3015857}

\item{density}{The material density in kg/m3, Default: ifelse(materialIsotropic, 4460, 1550)}

\item{flexuralRigidity11}{The material flexural rigidity 11 in Pa m3, Default: 198.7}

\item{flexuralRigidity22}{The material flexural rigidity 22 in Pa m3, Default: 13.23}

\item{flexuralRigidity12}{The material flexural rigidity 12 in Pa m3, Default: 3.98}

\item{flexuralRigidity66}{The material flexural rigidity 66 in Pa m3, Default: 8.15}

\item{transverseShearmodulus44}{The material transverse shear modulus 44 in MPa m, Default: 7940000}

\item{transverseShearmodulus55}{The material transverse shear modulus 55 in MPa m, Default: 12640000}

\item{lossFactor}{The material loss factor, Default: 0}

\item{shearCorrectionfactor}{The material shear correction factor, Default: NA}
}
\value{
The plate specification as a \code{plate} object.
}
\description{
A function that generates plate specification.
}
