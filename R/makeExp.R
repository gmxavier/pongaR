#' @title Generating experiment data.
#' @description A function that generates experiment data to produce several acoustic emission signal
#'              waveforms.
#' @param impactFrequency Impact frequency, in impacts/s, Default: 2000
#' @param expDuration Experiment duration, in s, Default: 0.001
#' @param sphereMaterial Solid sphere material name, Default: 'glassMcLaskey'
#' @param diameterDist Distribution of the solid sphere diameter, in m, Default: 'rnorm(n, mean = 1e-3, sd = 10e-6)'
#' @param velocityDist Distribution of the impact velocity, in m/s, Default: 'rnorm(n, mean = 1, sd = 1e-3)'
#' @param plateMaterial Plate material name, Default: 'titaniumDisperse'
#' @param plateGeometry Plate geometry, in m, Default: c(0.508, 0.381, 0.005)
#' @param lossFactor Material loss factor, Default: 0
#' @param shearCorrectionfactor Material shear correction factor, Default: NA
#' @param plateTheory Plate theory model name, Default: 'lamb'
#' @param excitabilityData Excitability data name, Default: 'titaniumEF'
#' @param modes Maximum number of modes \code{c(m,n)}, Default: NULL
#' @param sensorDiameter Sensor diameter (in m), Default: 0
#' @param sensorPosition Sensor position coordinates \code{c(x,y)} in m, Default: c(round(0.508/2, 3), round(0.381/2, 3))
#' @param sourceType Source type, Default: 'impulse'
#' @param sourceModel The source model, Default: 'hertz'
#' @param riseTime The plb rise time, in s, Default: 1.5e-06
#' @param forceAmplitude The plb and tone burst force amplitude, in N, Default: 1
#' @param centerFreq The tone burst center frequency, in Hz, Default: 2e+05
#' @param cyclesNum The tone burst number of cycles, Default: 5
#' @param sensorName The sensor dataset name, Default: 'r15iuc'
#' @param filterRange The band-pass filter range, in Hz, Default: c(50000, 1e+06)
#' @param amplifierGain The amplfier gain in dB, Default: 80
#' @param xunit Character indicating which unit is used in argument \code{duration}.
#'              If \code{xunit = "time"}, the unit is time in s,
#'              otherwise the number of samples, Default: 'time'
#' @param samp.rate Sampling rate of the \code{Wave}, in Hz, Default: 1e+07
#' @param bit Resolution of the \code{Wave} and rescaling unit, Default: 32. This may be
#'            1 for rescaling to numeric values in [-1,1],
#'            8 (i.e. 8-bit) for rescaling to integers in [0, 254],
#'            16 (i.e. 16-bit) for rescaling to integers in [-32767, 32767],
#'            24 (i.e. 24-bit) for rescaling to integers in [-8388607, 8388607],
#'            32 (i.e. 32-bit) for rescaling either to integers in [-2147483647, 2147483647],
#'            64 (i.e. 64-bit) for rescaling to numeric values in [-1, 1] (FLOAT_IEEE Wave format), and
#'            0 for not rescaling at all.
#' @return The experimental data as a \code{exp} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeExp
#' @export
makeExp <- function(impactFrequency = 2e3,
                    expDuration = 1e-3,
                    sphereMaterial = 'glassMcLaskey',
                    diameterDist = 'rnorm(n, mean = 1e-3, sd = 10e-6)',
                    velocityDist = 'rnorm(n, mean = 1, sd = 1e-3)',
                    plateMaterial = 'titaniumDisperse',
                    plateGeometry = c(0.508, 0.381, 5e-3),
                    lossFactor = 0,
                    shearCorrectionfactor = NA,
                    plateTheory = "lamb",
                    excitabilityData = "titaniumEF",
                    modes = NULL,
                    sensorDiameter = 0,
                    sensorPosition = c(round(0.508/2, 3), round(0.381/2, 3)),
                    sourceType = "impulse",
                    sourceModel = "hertz",
                    riseTime = 1.5e-6,
                    forceAmplitude = 1,
                    centerFreq = 2e+05,
                    cyclesNum = 5,
                    sensorName = "r15iuc",
                    filterRange = c(50e3, 1000e3),
                    amplifierGain = 80,
                    xunit = "time",
                    samp.rate = 10e6,
                    bit = 32) {

  ## laod material data
  data('material', envir = environment())

  ## make plate
  plate <- makePlate(length = plateGeometry[1],
                     width = plateGeometry[2],
                     thickness = plateGeometry[3],
                     youngModulus = material$youngModulus[material$name == plateMaterial],
                     poissonRatio = material$poissonRatio[material$name == plateMaterial],
                     density = material$density[material$name == plateMaterial],
                     lossFactor = lossFactor,
                     shearCorrectionfactor = shearCorrectionfactor)

  ## make sphere
  sphere <- makeSphere(diameter = 0.001,
                       impactVelocity = 1,
                       youngModulus = material$youngModulus[material$name == sphereMaterial],
                       poissonRatio = material$poissonRatio[material$name == sphereMaterial],
                       density = material$density[material$name == sphereMaterial])

  ## make medium
  mediumData <- makeMedium(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = c(0, 0),
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = expDuration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit)

  ## make source
  sourceData <- makeSource(sourceModel = sourceModel,
                           plate = plate,
                           sphere = sphere,
                           riseTime = riseTime,
                           forceAmplitude = forceAmplitude,
                           centerFreq = centerFreq,
                           cyclesNum = cyclesNum,
                           duration = expDuration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit)

  ## make readout
  readoutData <- makeReadout(sensorName = sensorName,
                             filterRange = filterRange,
                             amplifierGain = amplifierGain)

  ## make output (don't change the list order)
  object <- list(mediumData = mediumData,
                 sourceData = sourceData,
                 readoutData = readoutData,
                 impactFrequency = impactFrequency,
                 expDuration = expDuration,
                 samp.rate = samp.rate,
                 diameterDist = diameterDist,
                 velocityDist = velocityDist,
                 sensorPosition = sensorPosition)
  class(object) <- "exp"
  return(object)
}
