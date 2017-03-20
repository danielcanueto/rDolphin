

signal_fitting = function(parS, Xdata,multiplicities,roof_effect,Ydata,freq) {
  
  i = as.numeric(parS[seq(1, length(parS) - 4, 5)])
  p = as.numeric(parS[seq(2, length(parS) - 3, 5)])
  w = as.numeric(parS[seq(3, length(parS) - 2, 5)])*0.5/freq
  g = as.numeric(parS[seq(4, length(parS) - 1, 5)])
  j = as.numeric(parS[seq(5, length(parS) - 0, 5)])/freq
  signals_parameters=rbind(i,p,w,g,j)
  fitted_signals = matrix(0, dim(signals_parameters)[2], length(Xdata))

  # multiplicities = as.numeric(parS[seq(6, length(parS) - 1, 7)])
  # roof_effect = as.numeric(parS[seq(7, length(parS) - 0, 7)])
  NumSignals = length(parS) / 5


  for (s in seq_along(multiplicities)) {
    if (roof_effect[s] > 0) {
      if (multiplicities[s] == 1)   {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 2) {
        fitted_signals[s, ] = peakpvoigt(c(
          signals_parameters[1, s]/(1 - roof_effect[s]),
          (signals_parameters[2, s] - signals_parameters[5, s]/2),
          signals_parameters[3, s],
          signals_parameters[4, s]
        ),
          Xdata) + peakpvoigt(
            c(
              signals_parameters[1, s],
              (signals_parameters[2, s] + signals_parameters[5, s]/2),
              signals_parameters[3, s],
              signals_parameters[4, s]
            ),
            Xdata
          )
      } else if (multiplicities[s] == 3) {
        y=1/(2 + roof_effect[s])
        x= 1-y
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s] *x,
            (signals_parameters[2, s] - signals_parameters[5, s]),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s] *y,
            (signals_parameters[2, s] + signals_parameters[5, s]),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 4) {
        #     fitted_signals[s, ] = peakpvoigt(
        #       c(
        #         signals_parameters[1, s] / 3,
        #         (signals_parameters[2, s] - 3 * signals_parameters[5, s]),
        #         signals_parameters[3, s],
        #         signals_parameters[4, s]
        #       ),
        #       Xdata
        #     ) + peakpvoigt(c(
        #       signals_parameters[1, s],
        #       (signals_parameters[2, s] - signals_parameters[5, s]),
        #       signals_parameters[3, s],
        #       signals_parameters[4, s]
        #     ),
        #     Xdata) + peakpvoigt(c(
        #       signals_parameters[1, s],
        #       (signals_parameters[2, s] + signals_parameters[5, s]),
        #       signals_parameters[3, s],
        #       signals_parameters[4, s]
        #     ),
        #     Xdata) + peakpvoigt(
        #       c(
        #         signals_parameters[1, s] / 3,
        #         (signals_parameters[2, s] + 3 * signals_parameters[5, s]),
        #         signals_parameters[3, s],
        #         signals_parameters[4, s]
        #       ),
        #       Xdata
        #     )
      }
    } else if (roof_effect[s] == 0) {
      if (multiplicities[s] == 0) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 1) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 2) {
        fitted_signals[s, ] = peakpvoigt(c(
          signals_parameters[1, s],
          (signals_parameters[2, s] - signals_parameters[5, s]/2),
          signals_parameters[3, s],
          signals_parameters[4, s]
        ),
          Xdata) + peakpvoigt(c(
            signals_parameters[1, s],
            (signals_parameters[2, s] + signals_parameters[5, s]/2),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
            Xdata)
      } else if (multiplicities[s] == 3) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s] / 2,
            (signals_parameters[2, s] - signals_parameters[5, s]),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s] / 2,
            (signals_parameters[2, s] + signals_parameters[5, s]),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 4) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s] / 3,
            (signals_parameters[2, s] - 3 * signals_parameters[5, s]/2),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(c(
          signals_parameters[1, s] ,
          (signals_parameters[2, s] - signals_parameters[5, s]/2) ,
          signals_parameters[3, s],
          signals_parameters[4, s]
        ),
        Xdata) + peakpvoigt(c(
          signals_parameters[1, s],
          (signals_parameters[2, s] + signals_parameters[5, s]/2),
          signals_parameters[3, s],
          signals_parameters[4, s]
        ),
        Xdata) + peakpvoigt(
          c(
            signals_parameters[1, s] / 3,
            (signals_parameters[2, s] + 3 * signals_parameters[5, s]/2) ,
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      }
    } else if (roof_effect[s] < 0) {
      if (multiplicities[s] == 1) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 2) {
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s],
            (signals_parameters[2, s] - signals_parameters[5, s]/2),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(c(
          signals_parameters[1, s] *
            (1 - roof_effect[s]) ,
          (signals_parameters[2, s] + signals_parameters[5, s]/2),
          signals_parameters[3, s],
          signals_parameters[4, s]
        ),
          Xdata)
      } else if (multiplicities[s] == 3) {
        y=1/(2 + roof_effect[s])
        x= 1-y
        fitted_signals[s, ] = peakpvoigt(
          c(
            signals_parameters[1, s]*x,
            (signals_parameters[2, s] - signals_parameters[5, s]),
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s],
            signals_parameters[2, s],
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        ) + peakpvoigt(
          c(
            signals_parameters[1, s]*y,
            (signals_parameters[2, s] + signals_parameters[5, s]) ,
            signals_parameters[3, s],
            signals_parameters[4, s]
          ),
          Xdata
        )
      } else if (multiplicities[s] == 4) {
        # fitted_signals[s, ] = peakpvoigt(
        #   c(
        #     signals_parameters[1, s] / 3 ,
        #     (signals_parameters[2, s] - 3 * signals_parameters[5, s]),
        #     signals_parameters[3, s],
        #     signals_parameters[4, s]
        #   ),
        #   Xdata
        # ) + peakpvoigt(c(
        #   signals_parameters[1, s],
        #   (signals_parameters[2, s] - signals_parameters[5, s]),
        #   signals_parameters[3, s],
        #   signals_parameters[4, s]
        # ),
        # Xdata) + peakpvoigt(c(
        #   signals_parameters[1, s],
        #   (signals_parameters[2, s] + signals_parameters[5, s]),
        #   signals_parameters[3, s],
        #   signals_parameters[4, s]
        # ),
        # Xdata) + peakpvoigt(
        #   c(
        #     signals_parameters[1, s] / 3,
        #     (signals_parameters[2, s] + 3 * signals_parameters[5, s]),
        #     signals_parameters[3, s],
        #     signals_parameters[4, s]
        #   ),
        #   Xdata
        # )
      }
    }
  }


  return(fitted_signals)


}
