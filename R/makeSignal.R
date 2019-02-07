#' @title Generating acoustic emission signal due solid particle impact onto a plate.
#' @description A function that generates acoustic emission signal due solid particle impact onto a plate.
#' @param runData An object of class \code{run} produced by \code{makeRun}, Default: NULL
#' @param sourceData An object of class \code{source} produced by \code{makeSource}, Default: NULL
#' @param mediumData An object of class \code{medium} produced by \code{makeMedium}, Default: NULL
#' @param readoutData An object of class \code{readout} produced by \code{makeReadout}
#' @param plot If \code{TRUE}, plots the readout impulse response, Default: FALSE
#' @return The acoustic emission signal waveform as a \code{signal} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeSignal
#' @export
makeSignal <- function(runData = NULL,
                       sourceData = NULL,
                       mediumData = NULL,
                       readoutData = NULL,
                       plot = FALSE) {

  if (!is.null(runData)) {
    ## sanity check
    if (class(runData) != "run") {
      stop("Argument 'runData' must be an object of class run")
    }
    if (!is.logical(plot)) {
      stop("Argument 'plot' must be an object of type logical")
    }

    ## define readout_volt calculation function
    fun <- function(x, par) {

      top <- rep(0, (x[1] - 1))
      middle <- x[-1]
      object <- c(top, middle)
      object <- object[1:par]
      return(object)

    }

    ## define impact instants and add to out matrix
    par1 <- floor(runData$expData$expDuration*runData$expData$samp.rate)
    # instant <- c(2, sample(2:par1, size = (nrow(runData$expMatrix) - 1)))
    instant <- c(2, sample(2:(par1/2), size = (nrow(runData$expMatrix) - 1)))
    x <- rbind(instant, runData$out)

    ## start cluster
    cl <- makeCluster(getOption("cl.cores", (detectCores() - 1)))

    ## calculate readout_volt for each impact
    aux <- parApply(cl, x, 2, fun, par = par1)

    ## stop cluster
    stopCluster(cl)

    ## calculate readout_volt for all impacts
    readout_volt <- rowSums(aux)
    readout_volt <- Wave(left = readout_volt,
                         samp.rate = runData$expData$sourceData$samp.rate,
                         bit = runData$expData$sourceData$bit)

    ## plot signal if plot is TRUE
    if (plot) {
      wf <- readout_volt@left
      dt <- 1/runData$expData$sourceData$waveform@samp.rate
      timeArray <- (0:(length(wf)-1))*dt
      timeArray <- timeArray * 1e6 #scale to microseconds
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = "Readout output [V]")
      grid()
    }

    ## make output
    return(readout_volt)

  }
  else {
    ## sanity check
    if (class(sourceData) != "source") {
      stop("Argument 'sourceData' must be an object of class source")
    }
    if (class(mediumData) != "medium") {
      stop("Argument 'mediumData' must be an object of class medium")
    }
    if (class(readoutData) != "readout") {
      stop("Argument 'readoutData' must be an object of class readout")
    }
    if (!is.logical(plot)) {
      stop("Argument 'plot' must be an object of type logical")
    }

    ## calculate medium displacement and velocity due source
    medium_displ <- convolve(sourceData$waveform@left,
                             rev(mediumData$waveform@left),
                             type = "open")
    medium_displ <- (1/sourceData$waveform@samp.rate)*medium_displ
    medium_displ <- medium_displ[1:length(sourceData$waveform@left)]
    medium_veloc <- c(0, diff(medium_displ)*sourceData$waveform@samp.rate)

    ## calculate readout impulse response by filter design
    if (readoutData$samp.rate < sourceData$waveform@samp.rate) {
        fn <- sourceData$waveform@samp.rate/2
        f <- c(readoutData$readoutFR$freq, max(readoutData$readoutFR$freq), fn)/fn
        m <- c(readoutData$readoutFR$amp, 0, 0)

    } else {
        f <- readoutData$readoutFR$freq/max(readoutData$readoutFR$freq)
        m <- readoutData$readoutFR$amp
    }
    N <- 2000
    filt <- fir2(N, f, m)
    # ## plots check
    # fh <- freqz(filt)
    # readout_ir <- impz(filt = as.Arma(filt))
    # plot((readout_ir$t/10e6)*1e6, readout_ir$x, type = "l")
    # op <- par(mfrow = c(1, 2))
    # plot(f, m, type = "b", ylab = "magnitude", xlab = "Frequency")
    # lines(fh$f / pi, abs(fh$h), col = "blue")
    # # plot in dB:
    # plot(f, 20*log10(m + 1e-5), type = "b", ylab = "dB", xlab = "Frequency")
    # lines(fh$f / pi, 20*log10(abs(fh$h)), col = "blue")
    # par(op)
    # pad <- 5000
    # m1 <- Mod(fft(c(readout_ir$x, rep(0, pad))))[1:((N+pad)/5/2)]
    # f1 <- seq(0, 1, length.out = length(m1))
    # plot(f1, 20*log10(m1), type = "l", col = "blue")
    # lines(readoutData$readoutFR$freq/max(readoutData$readoutFR$freq), 20*log10(readoutData$readoutFR$amp), col = "red")
    # grid()
    # spec(readout_ir$x, 10e6, dB="max0", flim = c(0,1e3), alim = c(-50,0))

    ## calculate sensor output voltage due plate velocity (velocity sensor)
    if (names(readoutData$readoutFR)[2] == "amp_V/ubar") {
        medium_press <- medium_veloc*mediumData$plate$acousticImpedance*10 # convert m/s to ubar
        readout_volt <- fftfilt(b = filt, x = medium_press)
    }
    if (names(readoutData$readoutFR)[2] == "amp_V/(m/s)") {
        readout_volt <- fftfilt(b = filt, x = medium_veloc)
    }
    ## calculate sensor output voltage due plate displacement (displacement sensor)
    if (names(readoutData$readoutFR)[2] == "amp_V/m") {
        readout_volt <- fftfilt(b = filt, x = medium_displ)
    }

    ## change to waveform
    medium_displ <- Wave(left = medium_displ,
                         samp.rate = sourceData$waveform@samp.rate,
                         bit = 32)
    medium_veloc <- Wave(left = medium_veloc,
                         samp.rate = sourceData$waveform@samp.rate,
                         bit = 32)
    readout_volt <- Wave(left = readout_volt,
                         samp.rate = sourceData$waveform@samp.rate,
                         bit = 32)

    ## plot signals if plot is TRUE
    if (plot) {
      op <- par(mfrow = c(3, 1))
      wf <- sourceData$waveform@left
      dt <- 1/sourceData$waveform@samp.rate
      timeArray <- (0:(length(wf)-1))*dt
      timeArray <- timeArray * 1e6 #scale to microseconds
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = "Force [N]")
      wf <- mediumData$waveform@left
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = expression(paste("Plate i.r. [m] ", N^-1)))
      wf <- medium_displ@left
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = "Displacement [m]")
      wf <- medium_veloc@left
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = expression(paste("Velocity [m ", s^-1, "]")))
      wf <- readout_volt@left
      plot(timeArray,
           wf,
           type='l',
           col='black',
           xlab = expression(paste("Time [", mu, "s]")),
           ylab = "Readout output [V]")
      par(op)
    }

    ## make output
    object <- list(source_waveform = sourceData$waveform,
                   medium_ir = mediumData$waveform,
                   medium_displ = medium_displ,
                   medium_veloc = medium_veloc,
                   readout_volt = readout_volt)
    class(object) <- "signal"
    return(object)
  }

}
