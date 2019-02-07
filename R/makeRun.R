#' @title Generating run data.
#' @description A function that generates run data to produce several acoustic emission signal
#'              waveforms defined by \code{expData}.
#' @param expData An object of class \code{exp} produced by \code{makeExp}, Default: makeExp()
#' @return The run data as a \code{run} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeRun
#' @export
makeRun <- function(expData = makeExp()) {

  ## sanity check
  if (class(expData) != "exp") {
    stop("Argument 'expData' must be an object of class exp")
  }

  ## make experiment matrix
  impactFrequency <- expData$impactFrequency
  expDuration<- expData$expDuration
  sensorPosition <- expData$sensorPosition
  sourceN <- impactFrequency*expDuration
  sourceXY <- c(sample(seq(0, sensorPosition[1], 1e-3), sourceN, replace = TRUE),
                sample(seq(0, sensorPosition[2], 1e-3), sourceN, replace = TRUE))
  sourceXY <- matrix(data = sourceXY, ncol = 2)
  colnames(sourceXY) <- c('sourceX', 'sourceY')
  diameterDist <- expData$diameterDist
  velocityDist <- expData$velocityDist
  diameter <- eval(parse(text = paste("n<-", sourceN, ";", diameterDist, sep = "")))
  diameter <- round(diameter, digits = 6)
  impactVelocity <- eval(parse(text = paste("n<-", sourceN, ";", velocityDist, sep = "")))
  expMatrix <- cbind(sourceXY, diameter, impactVelocity)

  ## define readout_volt calculation function
  fun <- function(x, par) {

    requireNamespace(tuneR)
    requireNamespace(parallel)

    ## define plate
    plate <- par$mediumData$plate

    ## redefine sphere diameter and impact velocity
    sphere <- par$sourceData$sphere
    sphere$diameter <- x[3]
    sphere$impactVelocity <- x[4]

    ## calculate source data
    sourceData <- makeSource(sourceModel = par$sourceData$sourceModel,
                             plate = plate,
                             sphere = sphere,
                             riseTime = par$sourceData$riseTime,
                             forceAmplitude = par$sourceData$forceAmplitude,
                             centerFreq = par$sourceData$centerFreq,
                             cyclesNum = par$sourceData$cyclesNum,
                             duration = par$sourceData$duration,
                             xunit = par$sourceData$xunit,
                             samp.rate = par$sourceData$samp.rate,
                             bit = par$sourceData$bit)

    ## medium data calculation
    par$mediumData$sourcePosition <- c(x[1], x[2])
    mediumData <- makeMedium(plate = plate,
                             excitabilityData = par$mediumData$excitabilityData,
                             sensorDiameter = par$mediumData$sensorDiameter,
                             sensorPosition = par$mediumData$sensorPosition,
                             sourcePosition = par$mediumData$sourcePosition,
                             sourceType = par$mediumData$sourceType,
                             modes = par$mediumData$modes,
                             plateTheory = par$mediumData$plateTheory,
                             duration = par$mediumData$duration,
                             xunit = par$mediumData$xunit,
                             samp.rate = par$mediumData$samp.rate,
                             bit = par$mediumData$bit)
    # mediumData <- makeMedium(plate = plate,
    #                          excitabilityData = par$mediumData$excitabilityData,
    #                          sensorDiameter = par$mediumData$sensorDiameter,
    #                          sensorPosition = par$mediumData$sensorPosition,
    #                          sourcePosition = c(x[1], x[2]),
    #                          sourceType = par$mediumData$sourceType,
    #                          modes = par$mediumData$modes,
    #                          plateTheory = par$mediumData$plateTheory,
    #                          duration = par$mediumData$duration,
    #                          xunit = par$mediumData$xunit,
    #                          samp.rate = par$mediumData$samp.rate,
    #                          bit = par$mediumData$bit)

    ## readout data definition
    readoutData <- par$readoutData

    ## signal data calculation
    signalData <- makeSignal(sourceData = sourceData,
                             mediumData = mediumData,
                             readoutData = readoutData)

    ## make output
    object <- signalData$readout_volt@left
    return(object)
  }

  ## define experiment matrix and parameters
  x <- expMatrix
  dimnames(x) <- NULL
  par <- expData[1:3]

  ## start cluster
  cl <- makeCluster(getOption("cl.cores", (detectCores() - 1)))

  ## export needed functions to cluster (not necessary with AEA package)
  clusterExport(cl, c('makePlate',
                      'makeSphere',
                      'makeSource',
                      'phiSpline',
                      'makeMedium',
                      'apertureFunction',
                      'makeReadout',
                      'makeSignal',
                      'vel2lame',
                      "BesselJ",
                      "BesselY"))

  ## calculate readout_volt for each impact
  aux <- parApply(cl, x, 1, fun, par = par)

  ## stop cluster
  stopCluster(cl)

  ## make output
  object <- list(expData = expData,
                 expMatrix = expMatrix,
                 out = aux)
  class(object) <- "run"
  return(object)
}
