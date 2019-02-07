#' @title Generating source impulse response.
#' @description A function that generates source impulse response.
#' @param sourceModel The source model, Default: c("toneburst", "plb", "hertz", "hunter", "reed", "akay")
#' @param plate The plate specification produced by \code{makePlate}, Default: NA
#' @param sphere The solid sphere specification produced by \code{makeSphere}, Default: NA
#' @param riseTime The plb rise time, in s, Default: 1.5e-06
#' @param forceAmplitude The plb and tone burst force amplitude, in N, Default: 1
#' @param centerFreq The tone burst center frequency, in Hz, Default: 2e+05
#' @param cyclesNum The tone burst number of cycles, Default: 5
#' @param duration Duration of the \code{Wave} in \code{xunit}, Default: 0.0001024
#' @param xunit Character indicating which unit is used in argument \code{duration}.
#'              If \code{xunit = "time"}, the unit is time in s,
#'              otherwise the number of samples, Default: c("time", "samples")
#' @param samp.rate Sampling rate of the \code{Wave}, in Hz, Default: 1e+07
#' @param bit Resolution of the \code{Wave} and rescaling unit, Default: 32. This may be
#'            1 for rescaling to numeric values in [-1,1],
#'            8 (i.e. 8-bit) for rescaling to integers in [0, 254],
#'            16 (i.e. 16-bit) for rescaling to integers in [-32767, 32767],
#'            24 (i.e. 24-bit) for rescaling to integers in [-8388607, 8388607],
#'            32 (i.e. 32-bit) for rescaling either to integers in [-2147483647, 2147483647],
#'            64 (i.e. 64-bit) for rescaling to numeric values in [-1, 1] (FLOAT_IEEE Wave format), and
#'            0 for not rescaling at all.
#' @return The source impulse response waveform as a \code{source} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeSource
#' @export
makeSource <- function(sourceModel = c("toneburst", "plb", "hertz", "hunter", "reed", "akay"),
                       plate = NA,
                       sphere = NA,
                       riseTime = 1.5e-6,
                       forceAmplitude = 1,
                       centerFreq = 0.2e6,
                       cyclesNum = 5,
                       duration = 102.4e-6,
                       xunit = c("time", "samples"),
                       samp.rate = 10e6,
                       bit = 32) {

  phi <- NULL

  ## sanity check
  sourceModel <- match.arg(sourceModel)
  if (sourceModel != "hertz" &
      sourceModel != "hunter" &
      sourceModel != "reed" &
      sourceModel != "akay" &
      sourceModel != "plb" &
      sourceModel != "toneburst") {
    stop("Argument 'sourceModel' must be 'hertz', 'hunter', 'reed', 'akay', 'plb' or ''toneburst.")
  }
  if ((!is.na(plate)) && (class(plate) != "plate")) {
    stop("Argument 'plate' must be an object of class plate.")
  }
  if ((!is.na(sphere)) && (class(sphere) != "sphere")) {
    stop("Argument 'sphere' must be an object of class sphere.")
  }
  if (!is.numeric(riseTime) | (sourceModel == "plb" & riseTime <= 0)) {
    stop("Argument 'riseTime' must be a positive number.")
  }
  if (!is.numeric(forceAmplitude)) {
    stop("Argument 'forceAmplitude' must be a number.")
  }
  if ((abs(round(cyclesNum) - cyclesNum) > 0) | (cyclesNum < 1)) {
    stop("Argument 'cycles' must be a positive integer greater than 1.")
  }
  if (!is.numeric(duration) | (duration <= 0)) {
    stop("Argument 'duration' must be a positive number.")
  }
  xunit <- match.arg(xunit)
  if (xunit != "samples" & xunit != "time") {
    stop("Argument 'xunit' must be 'samples' or 'time'.")
  }
  if (!is.numeric(samp.rate) | (samp.rate <= 0)) {
    stop("Argument 'samp.rate' must be a positive number.")
  }

  if ((sourceModel == "akay")) {
    if ( (is.na(plate)) || (is.na(sphere)) ) {
      stop("Arguments 'plate' and 'sphere' must be specified for this 'sourceModel'.")
    }
    ## compute parameters
    hertzKfactor <- plate$hertzKfactor + sphere$hertzKfactor
    k <- (4/3/pi)*sqrt(sphere$diameter/2)*hertzKfactor^(-1)
    alpham <- (5*sphere$impactVelocity^2*sphere$mass/(4*k))^(2/5)
    contactDuration <- 2.9432*alpham/sphere$impactVelocity
    inelasticityParameter <- 3.218*sphere$mass/(contactDuration*plate$impedance)
    forceAmplitudeHertz <- k*alpham^(3/2)
    forceAmplitude <- forceAmplitudeHertz*exp(-inelasticityParameter)
    ## compute force pulse
    n <- 2
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    force_wf <- ifelse(time <= contactDuration,
                       forceAmplitude*(sin(pi*time/contactDuration)^n), 0)
    force_wf <- Wave(left = force_wf, samp.rate = samp.rate, bit = bit)
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = alpham,
                   contactDuration = contactDuration,
                   riseTime = 0,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
  if ((sourceModel == "hunter")) {
    if ( (is.na(plate)) || (is.na(sphere)) ) {
      stop("Arguments 'plate' and 'sphere' must be specified for this 'sourceModel'.")
    }
    ## compute parameters
    ## Hunter's paper eq. 29 (adapted for the definition of hertzKfactor)
    g <- pi*(plate$hertzKfactor + sphere$hertzKfactor)
    ## Hunter's paper eq. 33
    Z0 <- ((15/16)*g*(plate$mass*sphere$mass)/(plate$mass + sphere$mass))^(2/5)
    Z0 <- Z0*(sphere$diameter/2)^(-1/5)*(sphere$impactVelocity)^(4/5)
    ## hunter's paper eq. 39 (considering the exact value of integral in eq. 34)
    omega0 <- pi/(4*beta(2/5,1/2)/5)*(sphere$impactVelocity/Z0)
    ## Hunter's paper eq. 40
    forceAmplitude <- sphere$mass*Z0*(omega0)^2
    ## compute force pulse
    n <- 1
    contactDuration <- pi/omega0
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    force_wf <- ifelse(time <= contactDuration,
                       forceAmplitude*(sin(pi*time/contactDuration)^n), 0)
    force_wf <- Wave(left = force_wf, samp.rate = samp.rate, bit = bit)
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = Z0,
                   contactDuration = contactDuration,
                   riseTime = 0,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
  if ((sourceModel == "reed")) {
    if ( (is.na(plate)) || (is.na(sphere)) ) {
      stop("Arguments 'plate' and 'sphere' must be specified for this 'sourceModel'.")
    }
    ## compute parameters
    hertzKfactor <- plate$hertzKfactor + sphere$hertzKfactor
    ## Reed's paper K (adapted for the definition of hertzKfactor)
    K <- (4/3/pi)*hertzKfactor^(-1)
    alpha0 <- (5*sphere$mass*sphere$impactVelocity^2/(4*(sphere$diameter/2)^(1/2)*K))^(2/5)
    contactDuration <- (4*beta(2/5,1/2)/5)*(alpha0/sphere$impactVelocity)
    forceAmplitude <- (sphere$diameter/2)^(1/2)*K*(alpha0)^(3/2)
    ## compute force pulse
    n <- 3/2
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    force_wf <- ifelse(time <= contactDuration,
                       forceAmplitude*(sin(pi*time/contactDuration)^n), 0)
    force_wf <- Wave(left = force_wf, samp.rate = samp.rate, bit = bit)
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = alpha0,
                   contactDuration = contactDuration,
                   riseTime = 0,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
  if ((sourceModel == "hertz")) {
    if ( (anyNA(plate)) | (anyNA(sphere)) ) {
      stop("Arguments 'plate' and 'sphere' must be specified for this 'sourceModel'.")
    }
    ## load phi function approximation
    phi <- pongaR:::phi
    ## compute parameters
    hertzKfactor <- plate$hertzKfactor + sphere$hertzKfactor
    ## Reed's paper K (adapted for the definition of hertzKfactor)
    K <- (4/3/pi)*hertzKfactor^(-1)
    alpha0 <- (5*sphere$mass*sphere$impactVelocity^2/(4*(sphere$diameter/2)^(1/2)*K))^(2/5)
    contactDuration <- (4*beta(2/5,1/2)/5)*(alpha0/sphere$impactVelocity)
    forceAmplitude <- (sphere$diameter/2)^(1/2)*K*(alpha0)^(3/2)
    ## compute force pulse
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    force_wf <- ifelse(time <= contactDuration,
                       forceAmplitude*phi(time*sphere$impactVelocity/alpha0), 0)
    force_wf <- Wave(left = force_wf, samp.rate = samp.rate, bit = bit)
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = alpha0,
                   contactDuration = contactDuration,
                   riseTime = 0,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
  if ((sourceModel == "plb")) {
    alpha0 <- 0
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    force_wf <- ifelse(time <= riseTime,
                       (0.5 - 0.5*cos(pi*time/riseTime)), 1)
    force_wf <- Wave(left = force_wf, samp.rate = samp.rate, bit = bit)*forceAmplitude
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = 0,
                   contactDuration = 0,
                   riseTime = riseTime,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
  if ((sourceModel == "toneburst")) {
    alpha0 <- 0
    if (xunit == "samples") {
      time <- seq(from = 0,
                  to = (duration - 1),
                  by = 1)
      contactDuration <- ceiling(contactDuration*samp.rate)
      force_wf <- rep(0, duration)
    } else {
      time <- seq(from = 0,
                  to = (duration*samp.rate - 1),
                  by = 1)/samp.rate
      force_wf <- rep(0, duration*samp.rate)
    }
    input <- sin(2*pi*centerFreq*time)
    input <- ifelse(time <= cyclesNum/centerFreq,
                       input, 0)
    timeLength <- duration*samp.rate
    windowLength <- NROW(seq(0, cyclesNum/centerFreq, 1/samp.rate))
    input <- input*c(ftwindow(windowLength, wn = "hanning"),
                     rep(0, (timeLength - windowLength)))
    force_wf <- Wave(left = input, samp.rate = samp.rate, bit = bit)*forceAmplitude
    object <- list(sphere = sphere,
                   sourceModel = sourceModel,
                   alpha0 = 0,
                   contactDuration = 0,
                   riseTime = riseTime,
                   forceAmplitude = forceAmplitude,
                   centerFreq = centerFreq,
                   cyclesNum = cyclesNum,
                   duration = duration,
                   xunit = xunit,
                   samp.rate = samp.rate,
                   bit = bit,
                   waveform = force_wf)
    class(object) <- "source"
    return(object)
  }
}
