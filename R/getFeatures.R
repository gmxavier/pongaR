#' @title Extracting acoustic emission features.
#' @description Extracts several acoustic emission features from a waveform as defined in
#'              Sause's PhD thesis table 2.2, p. 50.
#' @param wave The waveform as a \code{Wave} object, Default: makeSource()$waveform
#' @param xlim The x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads
#'             to a ‘reversed axis’., Default: c(0, duration(wave))
#' @param samp.rate Sampling rate of the \code{Wave} in Hz, Default: 1e+07
#' @param preAmpGain Pre-amplifier gain in dB, Default: 40
#' @param preTriggerSamples Number of pre-trigger samples, Default: 0
#' @param threshold Thershold in V, Default: 0.01
#' @param maxDuration Max hit duration in s, Default: 0.05
#' @param maxHitNum Maximum number of hits extracted, Default: 5000
#' @param HDT Hit definition time in s, Default: 1e-04
#' @param HLT Hit lockout time in s, Default: 1e-05
#' @param PDT Peak definition time in s, Default: 5e-05
#' @param powerLimits Ranges of the four partial power features,
#'                    Default: rbind(c(0, 250000), c(250000, 5e+05), c(5e+05, 750000), c(750000, 1e+06))
#' @param plot If \code{TRUE}, plots the waveform with indications of the threshold
#'             and some features, Default: FALSE
#' @return A \code{list} object with the extracted acoustic emission features for each
#'         detected hit in the waveform. It also includes the start, peak and end indexes of the hits.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname getFeatures
#' @export
getFeatures <- function(wave = makeSource()$waveform,
                        xlim = c(0, duration(wave)),
                        samp.rate = 10e6,
                        preAmpGain = 40,
                        preTriggerSamples = 0,
                        threshold = 10e-3,
                        maxDuration = 50e-3,
                        maxHitNum = 5000,
                        HDT = 100e-6,
                        HLT = 10e-6,
                        PDT = 50e-6,
                        powerLimits = rbind(c(0, 250e3),
                                            c(250e3, 500e3),
                                            c(500e3, 750e3),
                                            c(750e3, 1000e3)),
                        plot = FALSE) {

  ## sanity check
  if (wave@samp.rate != samp.rate) {
    stop("Argument 'wave' must have samp.rate equal to 'samp.rate'.")
  }
  if ((!is.numeric(preAmpGain)) & (preAmpGain < 0)) {
    stop("Argument 'preAmpGain' must be numeric and non negative.")
  }
  if ((!is.numeric(preTriggerSamples)) & (preTriggerSamples < 0)) {
    stop("Argument 'preTriggerSamples' must be numeric and non negative")
  }
  if ((!is.numeric(threshold)) & (threshold < 0)) {
    stop("Argument 'threshold' must be numeric and non negative.")
  }
  if ((!is.numeric(maxDuration)) & (maxDuration < 0)) {
    stop("Argument 'maxDuration' must be numeric and non negative")
  }
  if ((!is.numeric(maxHitNum)) & (maxHitNum < 0)) {
    stop("Argument 'maxHitNum' must be numeric and non negative")
  }
  if ((!is.numeric(HDT)) & (HDT < 0)) {
    stop("Argument 'HDT' must be numeric and non negative")
  }
  if ((!is.numeric(HLT)) & (HLT < 0)) {
    stop("Argument 'HLT' must be numeric and non negative")
  }
  if ((!is.numeric(PDT)) & (PDT < 0)) {
    stop("Argument 'PDT' must be numeric and non negative")
  }
  if ((!is.numeric(powerLimits)) & all(powerLimits < 0)) {
    stop("Argument 'powerLimits' must be numeric and non negative")
  }
  if (!is.logical(plot)) {
    stop("Argument 'plot' must be an object of type logical")
  }

  # wave <- 1e-6*(c(bob6, rep(0, 500)) + c(rep(0, 500), 2*bob6)) # [V]
  # wave <- 1e-6*bob6 # [V]
  # f <- 10e6 # [Hz]
  # maxDuration <- 50e-3 # [s]
  # threshold <- 2.5e-6 # [V] 10e-6
  # maxHitNum <- 500 # [-]
  # PDT <- 200e-6 # [s]
  # HDT <- 200e-6 # [s]
  # HLT <- 15e-3 # [s]
  # preTriggerSamples <- 0 # [-]
  # preAmpGain <- 40 # [dB]
  # powerLimits <- rbind(c(0,250e3), c(250e3,500e3), c(500e3,750e3), c(750e3,1000e3))
  # plot = 1
  # RMSTime <- 10e-3 # between 10 and 1000 ms [s] (Sause's PhD thesis table 2.2, p. 50)

  # make input
  # input <- inputw(wave = wave, f = samp.rate)
  wave <- matrix(wave@left, ncol = 1) # wave <- input$w
  f <- samp.rate # f <- input$f
  # rm(input)

  # set thresholds
  th <- threshold
  th1 <- 100*th/max(abs(wave))

  # make waveform data
  wf <- wave[, 1]
  time <- (seq(1, length(wf)) - preTriggerSamples)/f
  wfLength <- length(wf)

  state <- 1
  hits <- matrix(data = 0, nrow = maxHitNum, ncol = 22)
  colnames(hits)<-c("firstIdx", "peakIdx", "lastIdx", "amplitude", "counts",
                    "duration", "riseTime", "RMS", "ASL", "averageFreq",
                    "reverbFreq", "initialFreq", "riseAngle", "decayAngle",
                    "absoluteEnergy", "peakFreq", "centroidFreq", "weightedFreq",
                    "pp1Freq", "pp2Freq", "pp3Freq", "pp4Freq")
  hit.counter <- 1
  i <- 1
  pb <- txtProgressBar(style = 3)

  # check hit presence
  if (max(abs(wave)) <= threshold) {
    hits <- hits[1, ]
    # # make plot
    # if (plot) {
    #     # trace plot
    #     plot(time*1e6,
    #          wf,
    #          type = 'l',
    #          col = 'black',
    #          xlab = expression(paste("Time [", mu, "s]")),
    #          ylab = "Amplitude [V]",
    #          ylim = c(-0.9*th, 1.1*th))
    #     grid()
    #     # threshold lines
    #     abline(h = th,
    #            lty = "dashed",
    #            col = "blue")
    #     mtext(text = expression("U"["thr"]),
    #           side = 4,
    #           at = th,
    #           las = 1,
    #           col = "blue")
    #     }
  } else {
    # begin calculation (based on Muravin's code)
    while (i <= wfLength) {
      setTxtProgressBar(pb, i/wfLength)
      if (
        ((wf[i] > 0) & (wf[i] > th) & (state == 1))
        # ||
        # ((wf[i] < 0) && (wf[i] < -th) && (state  == 1))
      ) {
        peak.value <- wf[i]
        peak.value.idx <- i
        first.th.cross.start.time <- time[i]
        first.th.cross.idx <- i
        peak.start.time <- time[i]
        hit.start.time <- time[i]
        state <- 2
      }
      if (i != 1) {
        if (
          ((wf[i] > 0) && (state == 2) && (wf[i] < th) && (wf[i - 1] > th))
          # ||
          # ((wf[i] < 0) && (state == 2) && (wf[i] < -th) && (wf[i - 1] > -th))
        ) {
          hit.start.time <- time[i]
        }
      } else {
        if (
          ((wf[i] > 0) && (state == 2) && (wf[i] < th))
          # ||
          # ((wf[i] < 0) && (state == 2) && (wf[i] < -th))
        ) {
          hit.start.time <- time[i]
        }
      }
      if ((state == 2)
          &&
          ((time[i] - peak.start.time) > PDT)
      ) {
        peak.start.time <- 1e+50
      }
      if (
        ((wf[i] > 0) && (state == 2) && (wf[i] > abs(peak.value))
         &&
         ((time[i] - peak.start.time) < PDT))
        # ||
        # ((wf[i] < 0) && (state == 2) && (wf[i] < -abs(peak.value))
        #  &&
        #  ((time[i] - peak.start.time) < PDT))
      ) {
        peak.value <- wf[i]
        peak.value.idx <- i
        peak.start.time <- time[i]
      }
      if (i != 1) {
        if (((wf[i] > 0) && (state == 2) && (wf[i] > th) && (wf[i - 1] < th))
            ||
            ((wf[i] > 0) && (state == 2) && (wf[i] < th) && (wf[i - 1] > th))
            # ||
            # ((wf[i] < 0) && (state == 2) && (wf[i] > -th) && (wf[i - 1] < -th))
            # ||
            # ((wf[i] < 0) && (state == 2) && (wf[i] < -th) && (wf[i - 1] > -th))
        ) {
          last.th.cross.idx <- i
        }
      } else {
        if (((wf[i] > 0) && (state == 2) && (wf[i] > th))
            ||
            ((wf[i] > 0) && (state == 2) && (wf[i] < th))
            # ||
            # ((wf[i] < 0) && (state == 2) && (wf[i] > -th))
            # ||
            # ((wf[i] < 0) && (state == 2) && (wf[i] < -th))
        ) {
          last.th.cross.idx <- i
        }
      }
      if ( ( state == 2 )
           &&
           ( ( (time[i] - hit.start.time) > HDT)
             || ( (time[i] - first.th.cross.start.time) > maxDuration ) ||
             ( i == length(wf) ) ) )  {

        if ( i == length(wf) ) {
          last.th.cross.idx <- i
        }

        wfHit = wf[first.th.cross.idx:last.th.cross.idx]

        ## basic indexes
        # hit start index
        hits[hit.counter, 1] <- first.th.cross.idx
        # hit max index
        hits[hit.counter, 2] <- peak.value.idx
        # hit stop index
        hits[hit.counter, 3] <- last.th.cross.idx

        ## basic signal parameters (Sause's PhD thesis table 2.1, p. 49)
        # t0 - time of first threshold crossing (arrival time) [s]
        t0 <- time[first.th.cross.idx]
        # NAE - number of signals threshold crossings (counts) [-]
        NAE <- sum(rle(wfHit >= threshold)$values)
        # tAE - time between first and last threshold crossing of signal (duration) [s]
        tAE <- time[last.th.cross.idx] - time[first.th.cross.idx]
        # Umax - maximum signal voltage [V]
        Umax <- wf[peak.value.idx]
        # tpeak - time of maximum signal voltage [s]
        tpeak <- time[peak.value.idx]
        # fpeak - frequency of maximum signal contribution [Hz]
        wfspec <- spec(wave = wfHit, f = f, norm = F, plot = F)
        fpeak <- wfspec[(wfspec[, 2] == max(wfspec[, 2])), 1]*1e3
        # Npeak - number of threshold crossings between t0 and tpeak (counts to peak) [-]
        wf1 <- wf[first.th.cross.idx:peak.value.idx]
        Npeak <- sum(rle(wf1 >= threshold)$values)
        rm(wf1)

        ## acoustic emission signal features (Sause's PhD thesis table 2.2, p. 50)
        # amplitude (log amplitude) [dB ref 1 uV]
        hits[hit.counter, 4] <- 20*log10(abs(Umax)/1e-6) - preAmpGain
        # counts [-]
        hits[hit.counter, 5] <- NAE
        # duration [s]
        hits[hit.counter, 6] <- tAE
        # rise time [s]
        hits[hit.counter, 7] <- tpeak - t0
        # root mean square (RMS) [V]
        hits[hit.counter, 8] <- sqrt(sum((1/f)*wfHit^2)/tAE)
        # average signal level (ASL) [dB]
        wf1 <- abs(wfHit)
        wf1 <- 20*log10(wf1/threshold) - preAmpGain
        hits[hit.counter, 9] <- sum((1/f)*wf1)/tAE
        # average frequency [Hz]
        hits[hit.counter, 10] <- NAE/tAE
        # reverberation frequency [Hz]
        hits[hit.counter, 11] <- (NAE - Npeak)/(t0 + tAE - tpeak)
        # initial frequency [Hz]
        hits[hit.counter, 12] <- Npeak/(tpeak - t0)
        # rise angle [degree] (adapted from original)
        hits[hit.counter, 13] <- atan(abs(Umax)/((tpeak - t0)/tAE))*180/pi
        # decay angle [degree] (adapted from original)
        hits[hit.counter, 14] <- atan(abs(Umax)/((tAE - (tpeak - t0))/tAE))*180/pi
        # absolute energy [aJ]
        hits[hit.counter, 15] <- 1e18*sum((1/f)*wfHit^2)/10e3
        # peak frequency [Hz]
        hits[hit.counter, 16] <- fpeak
        # frequency centroid [Hz]
        hits[hit.counter, 17] <- sum(wfspec[, 1]*wfspec[, 2])/sum(wfspec[, 2])*1e3
        # weighted peak-frequency [Hz]
        hits[hit.counter, 18] <- sqrt(hits[hit.counter, 16]*hits[hit.counter, 17])
        # partial power 1 [%]
        wfpspec <- spec(wave = wfHit, f = f, norm = F, PSD = T, plot = F)
        aux <- sum(wfpspec[, 2])
        hits[hit.counter, 19] <- 100*sum(wfpspec[(wfpspec[, 1] >= powerLimits[1, 1]*1e-3 &
                                                    wfpspec[, 1] <= powerLimits[1, 2]*1e-3), 2])/aux
        # partial power 2 [%]
        hits[hit.counter, 20] <- 100*sum(wfpspec[(wfpspec[, 1] >= powerLimits[2, 1]*1e-3 &
                                                    wfpspec[, 1] <= powerLimits[2, 2]*1e-3), 2])/aux
        # partial power 3 [%]
        hits[hit.counter, 21] <- 100*sum(wfpspec[(wfpspec[, 1] >= powerLimits[3, 1]*1e-3 &
                                                    wfpspec[, 1] <= powerLimits[3, 2]*1e-3), 2])/aux
        # partial power 4 [%]
        hits[hit.counter, 22] <- 100*sum(wfpspec[(wfpspec[, 1] >= powerLimits[4, 1]*1e-3 &
                                                    wfpspec[, 1] <= powerLimits[4, 2]*1e-3), 2])/aux

        hit.counter <- hit.counter + 1
        state <- 3
        lock.time.begin <- time[i]
        # print(wf[first.th.cross.idx])
        # readline(prompt="Press [enter] to continue")
      }
      if (state == 3) {
        if ((time[i] - lock.time.begin) < HLT) {
          state <- 3
        } else {
          state <- 1
        }
      }
      i <- i + 1
    }
    if (hit.counter != 1) {
      hits <- hits[1:(hit.counter - 1), ]
    }
    ## make plot
    if (plot) {
      # trace plot
      plot(time*1e6,
           wf,
           type = 'l',
           col = 'black',
           xlab = expression(paste("Time [", mu, "s]")),
           xlim = xlim*1e6,
           ylab = "Amplitude [V]",
           cex.lab = 1.5,
           cex.axis = 1.5,
           las = 1)
      grid()
      usr <- par("usr")
      # threshold lines
      abline(h = th,
             lty = "dashed",
             col = "black")
      mtext(text = expression("U"["thr"]),
            side = 4,
            at = th,
            las = 1,
            col = "black",
            cex = 1.5)
      # amplitude line
      abline(h = Umax,
             lty = "dashed",
             col = "black")
      mtext(text = expression("U"["peak"]),
            side = 4,
            at = Umax,
            las = 1,
            col = "black",
            cex = 1.5)
      # arrival time line
      abline(v = t0*1e6,
             lty = "dashed",
             col = "black")
      mtext(text = expression("t"[0]),
            side = 3,
            las = 0,
            at = t0*1e6,
            col = "black",
            cex = 1.5)
      # peak time line
      abline(v = tpeak*1e6,
             lty = "dashed",
             col = "black")
      mtext(text = expression("t"["peak"]),
            side = 3,
            las = 0,
            at = tpeak*1e6,
            col = "black",
            cex = 1.5)
      # duration time line and interval
      abline(v = (t0 + tAE)*1e6,
             lty = "dashed",
             col = "black")
      mtext(text = expression(paste("t"[0]," + t"["AE"])),
            side = 3,
            las = 0,
            at = (t0 + tAE)*1e6,
            col = "black",
            cex = 1.5)
    }
  }
  close(pb)

  ## make output
  object <- list(AEfeatures = hits)
  return(object)
}
