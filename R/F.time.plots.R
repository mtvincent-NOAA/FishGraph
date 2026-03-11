#' Plots of fishing mortality rate over time
#'
#' The routine \code{F.time.plots} provides time plots of fishing mortality rate
#' \emph{F} over time. Plots are of absolute estimated \emph{F} and \emph{F}
#' relative to various reference points. Through the argument \code{F.references},
#' arbitrary reference points can be included in plots, one by one or together.

#'
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param start.drop Number of years at the start of the data to be omitted from
#' plots, as when a model includes an initialization period.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param legend.pos A text string compatible with the \code{legend} function of \code{R}.
#' Defines the position of the legend (ex. "bottomright", "bottom", etc.)
#' @param F.references A list of character-string vectors to specify reference
#' points to appear on \emph{F} plots.
#' @param F.additional A vector of character-string names specifying additional F metrics to be plotted.
#' An extension of *.ratio indicates that the metric is a relative measure.
#'
#' @return Graphics
#'
#' @author M Prager
#' @author E Williams
#' @author K Shertzer
#' @author R Cheshire
#' @author K Purcell
#'
#' @examples \donttest{
#' F.time.plots(gag)
#' }
#' @export
F.time.plots <-
function(x, DataName = deparse(substitute(x)), draft = TRUE,
   start.drop = 0, graphics.type = NULL, use.color = TRUE, legend.pos = "topleft",
   F.references = NULL, F.additional=NULL)
###########################################################################################
#  ARGUMENTS:
#  x - an R list with output from the assessment models
#     The list x must have a component x$t.series that is a data frame
#     The data frame must have components
#        x$t.series$year (integer) with years of analysis
#        x$t.series$F.* (real) with fishing mortality rates
#     The list x must have a component x$parms that is a list
#     The list x$parms must have a (vector) component x$parms$SSBmsy that has
#        a single value of spawning-stock biomass at MSY.
#  DataName - a string representation an identifier for the data (run) in use.
#  draft - TRUE if figures are to have main titles
#  graphics.type - a character vector with graphics-file types
#  start.drop - integer, number of years to drop from start of file
#  use.color - TRUE of graphs are in color
###########################################################################################
{  ### Check for needed data components
   if (! ("t.series" %in% names(x)))
   {  Errmsg = (paste("Component ", deparse(substitute(x)), "$t.series not found.",
         sep = ""))
      warning(Errmsg, immediate. = TRUE)
      return(invisible(-1))
   }
   # Find column labeled 'year':
   yrcol <- grep("year",names(x$t.series))
   if (length(yrcol) != 1)
   {  warning("No year column found in x$t.series", immediate. = TRUE)
      return(invisible(-1))
   }
   # Get indices of columns beginning with "F":
   Fcols <-  grep("^F\\.",names(x$t.series))
   if (length(Fcols) == 0)
   {  warning("No F columns found in F.plots!", immediate. = TRUE)
      return(invisible(-1))
   }
   ### Make local copies of data:
   yrndx <- (start.drop + 1):nrow(x$t.series)
   year <- x$t.series$year[yrndx]
   Fdata <- x$t.series[yrndx, Fcols, drop = FALSE]
   # Extract F.full if present:
   Ffull.found <- FALSE
   Ffullcol <- grep("F.full", names(Fdata))
   if (length(Ffullcol) > 0)
   {  Ffull <- Fdata[,Ffullcol]
      Fdata <- Fdata[,-Ffullcol, drop = FALSE]  # Remove it from main F data frame
      Ffull.found <- TRUE
   }
   # Extract F.Fmsy if present:
   {  FFmsy.found = FALSE}
   FFmsycol <- grep("F.Fmsy$", names(Fdata))
   if (length(FFmsycol) > 0)
   {  FFmsy <- Fdata[,FFmsycol]
      Fdata <- Fdata[,-FFmsycol, drop = FALSE]
      FFmsy.found = TRUE
   }

   # Extract F.additional if present:
   if (length(F.additional>0))
   {
     Fadd=matrix(0, nrow=length(yrndx), ncol=length(F.additional))
     for (i in 1:length(F.additional)) {
       Faddname <- F.additional[i]
       Faddcol <- which(names(Fdata)==Faddname)
       if (length(Faddcol) > 0)
       {  Fadd[,i] <- Fdata[,Faddcol]
          Fdata <- Fdata[,-Faddcol, drop = FALSE]
       } else {Errmsg = (paste("User-defined component ", Faddname, " not found in $t.series.",
                               sep = ""))
               warning(Errmsg, immediate. = TRUE)
               return(invisible(-1))

       }
      }
    }

   nplots <- ncol(Fdata)             # number of individual series to plot

   ### If writing graphics files, make sure there is a directory for them:
   if (length(graphics.type > 0))
   {  write.graphs <- TRUE
      GraphicsDirName <- paste(DataName, "-figs/F", sep="")
   } else
   { write.graphs <- FALSE
   }
   ### Set graphics parameters, constants, data structures:
   plot.options = FGGetOptions()
   savepar <- FGSetPar(draft)
   PlotTitle <- ""
   colvec <- FGGetPal(nplots, use.color)
   if(use.color) Y1Col <- plot.options$color$clr.line
   else Y1Col <- plot.options$bw$clr.line
   fgcex = 1.0
   Fmat <- matrix(NA, nplots, length(year))    # matrix to save series F's
   bar.names <- rep("", nplots)              # vector of bar names

   ### Make plots of individual F series:
   lab.x <- "Year"
   lab.y <- "Fishing mortality rate"
   for (iplot in 1:nplots)
   {  fishery   <- unlist(strsplit(names(Fdata[iplot]), ".", fixed = TRUE))
      fishery   <- paste(fishery[-1], collapse=".")
      if(draft)
      { PlotTitle <- FGMakeTitle(paste("Fishery:", fishery), DataName)
      }
      ymin <- 0.0
      ymax <- max(Fdata[,iplot], 0.01, na.rm = TRUE)  # (0.01) --> very small Fs only

      FGTimePlot(year, Fdata[,iplot], lab.x = lab.x, lab.y = lab.y,
         use.color = use.color, ylim = c(ymin, ymax), main = PlotTitle,
         FGtype = "1line", Y1Col = Y1Col)
      if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
         GraphName = paste("F.", fishery, sep=""), graphics.type)
   }  # end for (iplot...)

   ### Make plot of F.full
   if (Ffull.found)
   {  if(draft) PlotTitle <- FGMakeTitle("Full F", DataName)
      ymin = 0.0
      ymax <- max(Ffull, 0.01, na.rm = TRUE)
      FGTimePlot(year, Ffull, lab.x = lab.x, lab.y = lab.y, use.color = use.color,
         ylim = c(ymin, ymax), main = PlotTitle, FGtype = "1line", Y1Col = Y1Col)
      if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
         GraphName = "F.full", graphics.type)
   }
   ### Make plot of F.Fmsy:
   if (FFmsy.found)
   {  if(draft) PlotTitle <- FGMakeTitle("F/Fmsy", DataName)
      lab.y <- expression(bold(italic(F) / italic(F)[MSY]))
      ymin = 0.0
      ymax <- max(FFmsy, 0.01, na.rm = TRUE)

      FGTimePlot(year, FFmsy, lab.x = lab.x, lab.y = lab.y, use.color = use.color,
         ylim = c(ymin, ymax), main = PlotTitle, href = 1.0, hrefnames = NULL,
         FGtype = "1line", Y1Col = Y1Col)
      if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
         GraphName = "F.Fmsy", graphics.type)
    }

   ### Make plots of Ffull with user-specified reference lines:
   if (Ffull.found)
   {  lab.x <- "Year"
      lab.y <- "Fishing mortality rate (full)"
      if (draft) PlotTitle <- FGMakeTitle("Full F", DataName)
      if (is.list(F.references))
      {  nrefs = length(F.references)
         for (iplot in 1:nrefs)
         {  hrefnames <- F.references[[iplot]]
            # Get a vector of the indices of the references in x$parms:
            hrefindex <- NULL
            for (k in (F.references[[iplot]]))
            {  hrefindex <- c(hrefindex, which(names(x$parms) == k))
            }
            href <- unlist(x$parms[hrefindex])
            ymax <- max(Ffull, 0.01, href, na.rm = TRUE)

            FGTimePlot(year, Ffull, lab.x = lab.x, lab.y = lab.y, use.color = use.color,
               ylim = c(0, ymax), main = PlotTitle, FGtype = "1line", Y1Col = Y1Col,
               hrefnames = hrefnames, href = href, legend.pos = legend.pos)

            if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
               GraphName = paste("F.full.ref", sprintf("%03i",iplot), sep=""),
               graphics.type)
         }
      }
   }

   ### Make plots of user-specified F metrics:
   if (length(F.additional>0))
   {  lab.x <- "Year"

      for (i.Fadd in 1:length(F.additional)) {
          Faddname <- F.additional[i.Fadd]

          addname.dum=unlist(strsplit(F.additional[i.Fadd],".", fixed=TRUE))
          if (tail(addname.dum,1)=="ratio") {ratio=TRUE} else {ratio=FALSE}
          addname=paste(addname.dum[-1],collapse=".")

          if (ratio) {
            lab.y <- paste("Relative fishing mortality rate (",addname,")",sep="")
            if (draft) PlotTitle <- FGMakeTitle(Faddname, DataName)
            ymax <- max(Fadd[,i.Fadd], 0.01, na.rm = TRUE)
            FGTimePlot(year, Fadd[,i.Fadd], lab.x = lab.x, lab.y = lab.y, use.color = use.color,
                       ylim = c(0, ymax), main = PlotTitle, FGtype = "1line", Y1Col = Y1Col,
                       href=1.0, hrefnames=NULL, legend.pos = legend.pos)
            if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                                         GraphName = Faddname, graphics.type)


          } else {

            lab.y <- paste("Fishing mortality rate (",addname,")",sep="")
            if (draft) PlotTitle <- FGMakeTitle(Faddname, DataName)
            ymax <- max(Fadd[,i.Fadd], 0.01, na.rm = TRUE)
            FGTimePlot(year, Fadd[,i.Fadd], lab.x = lab.x, lab.y = lab.y, use.color = use.color,
                       ylim = c(0, ymax), main = PlotTitle, FGtype = "1line", Y1Col = Y1Col,
                       legend.pos = legend.pos)
            if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                                         GraphName = Faddname, graphics.type)


            if (is.list(F.references))
            {  nrefs = length(F.references)
               for (iplot in 1:nrefs)
               {  hrefnames <- F.references[[iplot]]
                  # Get a vector of the indices of the references in x$parms:
                  hrefindex <- NULL
                  for (k in (F.references[[iplot]]))
                    {hrefindex <- c(hrefindex, which(names(x$parms) == k))}
                  href <- unlist(x$parms[hrefindex])
                  ymax <- max(Fadd[,i.Fadd], 0.01, href, na.rm = TRUE)

                  FGTimePlot(year, Fadd[,i.Fadd], lab.x = lab.x, lab.y = lab.y, use.color = use.color,
                             ylim = c(0, ymax), main = PlotTitle, FGtype = "1line", Y1Col = Y1Col,
                             hrefnames = hrefnames, href = href, legend.pos = legend.pos)

                  if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                                               GraphName = paste(Faddname,".ref", sprintf("%03i",iplot), sep=""),
                                               graphics.type)
               }
            }
          }

      } #end for i.Fadd
   } #end if F.additional



   ### Make stacked barplot of F's:
   if (nplots > 1)
   {  Fmat <- t(as.matrix(Fdata))
      bar.names <- sub("F\\.", "", rownames(Fmat))

      if(draft) PlotTitle <- FGMakeTitle("F by fishery", DataName)
      # Convert NAs to zeros for plotting
      Fmat <- ifelse(is.na(Fmat), 0, Fmat)

      barplot(Fmat, beside = FALSE, axis.lty = 1, col = colvec,
         main = PlotTitle, ylab = "Fishing mortality rate", xlab = "Year", las = 1)
      # Note: legend colors & names are reversed to correspond with plot:
      legend(legend.pos, inset = c(.05,.1), rev(bar.names), fill = rev(colvec),
         bg = "white")
      abline(h = 0)
      if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
         GraphName = "F.stacked", graphics.type)
    }

   par(savepar)
   return(invisible(NULL))
}  # End of function definition

