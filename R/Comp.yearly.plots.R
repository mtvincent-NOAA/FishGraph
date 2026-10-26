#' Plots of age- and length-composition fits by year
#'
#' The routine \code{Comp.yearly.plots} generates plots of age- and length-composition
#' fits by year and data series. Optionally, the sample size N, effective sample
#' size Neff, angular deviation can be printed on the plot surface. The background is
#' shaded if the model does not fit to the composition (Neff < 0).  Plots of the
#' ratio Neff/N over time may be made. Plots are available in two formats: one
#' plot per page, or a compact format with 21 plots per page.
#'
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param units a character string (e.g. \code{"cm"}) for labeling the X-axis of
#' length-composition plots.
#' @param print.n When \code{TRUE}, \emph{N} is printed in each plot.
#' @param print.neff When \code{TRUE}, \emph{N_eff} (rounded) is printed in each plot.
#' @param plot.missing When \code{TRUE}, plots out all comps. Otherwise comps where
#' \code{is.na(N)} are not printed.
#' @param plot.pred when \code{TRUE}, plots the fit to the composition data, when \code{FALSE} only plots the observed compositions
#' @param plot.neff When \code{TRUE}, a timeplot of \emph{N_eff} is made for each data
#' series
#' @param compute.neff When \code{TRUE}, \strong{FishGraph} computes \emph{N_eff}
#' internally, rather than using values stored by the user.  The computations are
#' based on the sample sizes and observed and estimated proportions.
#' @param print.year When \code{TRUE}, the year is printed in each plot.
#' @param compact When \code{TRUE}, plots are arranged by default in 5x3 matrix on each page.
#' @param uniform When \code{TRUE}, all years of the same data series are scaled
#' the same.
#' @param connect.obsd When \code{TRUE}, observed points are connected by lines.
#'
#' @return Graphics
#'
#' @author M Prager
#' @author E Williams
#' @author K Shertzer
#' @author R Cheshire
#' @author K Purcell
#'
#'
#' @examples \donttest{
#' Comp.yearly.plots(gag)
#' }
#' @export
Comp.yearly.plots <-
function(x, DataName = deparse(substitute(x)), draft = TRUE,
   graphics.type = NULL, use.color = TRUE, units = x$info$units.length,
   print.n = !compact, print.neff = !compact, plot.missing = FALSE, plot.pred = TRUE,
   plot.neff = FALSE, compute.neff = FALSE, print.year = compact, compact = TRUE,
   uniform = TRUE, connect.obsd = FALSE)
#-------------------------------------------------------------------------------
#  ARGUMENTS:
#  x - an R list with output from the assessment models
#     The list x must have a component x$comp.mats that is a list of matrices
#     The list must have matrices in pairs in order
#        x$comp.mats$xxx.ob observed age- or length-comps as proportions
#        x$comp.mats$xxx.pr predicted age- or length-comps as proportions
#        Each of these matrices should have the years as column names and
#        the ageclasses or length bins as row names
#     To plot sample sizes, the list x must have a component x$t.series
#        that is a data frame.  For each xxx.ob matrix (above), there should
#        be a variable xxx.n in x$t.series with the sample sizes.  There should
#        also be a variable x$t.series.year with years.  These may be a superset
#        of the years in each matrix.
#  DataName - a string representation an identifier for the data (run) in use.
#  draft - if TRUE the figure has a main title
#  graphics.type - a character vector with graphics-file types
#  use.color - TRUE of graphs are in color
#  units - a character string for labeling the units of length bins.
#  print.n - TRUE to print actual sample size on each plot.
#  print.neff - TRUE to print effective sample size on each plot.
#  plot.neff - TRUE to plot effective sample size by year for each series.
#  compute.neff - TRUE to compute eff. sample size internally.
#  print.year - TRUE to print year on each plot.
#  compact - TRUE to arrange many plots on a page.
#  uniform - TRUE to use same scale for all plots for the same fleet.
#  connect.obsd - TRUE to connect observed points with a line.
#-------------------------------------------------------------------------------
{
   ### Check for data presence:
   Errstring = ("No composition data found.  Terminating Comp.yearly.plots");
   if (! ("comp.mats" %in% names(x))) stop(Errstring)

   ### Make local copies of list of matrices:
   cm <- x$comp.mats
   compnames <- names(cm)

   ### Set graphics parameters for compact or non-compact mode:
   plot.options = FGGetOptions()
   savepar <- FGSetPar(draft)
   if (compact)
   {  par(mar = c(3.5, 4, 0.85, 0.7), cex = 0.9, cex.main = 0.9, cex.lab = 0.85,
         cex.axis = 0.75, mfcol = plot.options$graphics$fglayout, lab = c(5, 2, 4),
         mgp=c(1.8, 0.75, 0), las = 0)
      draft <- FALSE  ### Plots have no main titles in compact mode
   }    else
   {  par(las = 1)
   }
   PlotsPerPage <- prod(par("mfrow"))
   PlotsWaiting <- 0
   PlotsPage <- 1
   ### Set colors:
   if (use.color)
   {  col.obsd <- plot.options$color$clr.obsd
      col.pred <- plot.options$color$clr.pred
      Y1Col  <- plot.options$color$clr.line
   }    else
   {  col.obsd <- plot.options$bw$clr.obsd
      col.pred <- plot.options$bw$clr.pred
      Y1Col  <- plot.options$bw$clr.line
   }
   ### If writing graphics files, make sure there is a directory for them:
   if (length(graphics.type > 0))
   {  write.graphs <- TRUE
      GraphicsDirName <- paste(DataName, "-figs/compyr", sep="")
   }   else
   {  write.graphs <- FALSE
   }
   ### Make plots:
   nindexes <- length(cm) %/% 2              # integer division operator
   for (jindex in 1:nindexes)                ### By data series
   {  prop.obsd <- cm[[jindex*2 - 1]]
      prop.pred <- cm[[jindex*2]]
      ymax <- min(max(prop.obsd,prop.pred, na.rm=TRUE),
                  1.3*quantile(c(prop.obsd, prop.pred), 0.98, na.rm=TRUE))
      years     <- as.numeric(rownames(prop.obsd))
      nyears    <- length(years)
      bins      <- as.numeric(colnames(prop.obsd))
      fileroot  <- FGTrimName(names(cm)[jindex*2], removePrefix = 0, removeSuffix = 1)
      ptype     <- substr(fileroot, 1, 1)        # ptype = "l" or "a" for length or age
      fishery   <- FGTrimName(fileroot, removePrefix = 1, removeSuffix = 0)
      titleroot <- paste( "Data: ", DataName, "    Sample: ", fishery, sep = "")
      ### Get Nobs for this matrix from the x$t.series dataframe:
      ## if (print.n || print.neff || plot.neff)
      ## {
          yearvec = rownames(prop.obsd)
         strtmp = paste("x$t.series$", fileroot, ".n[x$t.series$year %in% yearvec]", sep="")
         Nobs = eval(parse(text = strtmp))
      ## } else {Nobs=rep(NA,nyears)}
      ### Get Neff for this matrix from the x$t.series dataframe IF Neff isn't computed internally.
      ### otherwise, set up a blank vector for use when it IS computed:
      ## if (print.neff || plot.neff)
      ##    {
             if (compute.neff)
                { Neff <- rep(NA_integer_, nyears) }
             else
                { strtmp <- paste("x$t.series$", fileroot,
                    ".neff[x$t.series$year %in% yearvec]", sep="")
                  Neff <- eval(parse(text = strtmp))
                  if (is.null(Neff) || any(is.na(Neff)))
                    {   warning("Comp.yearly.plots: comp.neff=FALSE but no data or NA in\n",
                            paste("   x$t.series$", fileroot, ".neff", sep = ""), call. = FALSE)
                        Neff <- rep(NA_integer_, nyears)
                    } else (Neff<-round(Neff,digits=1))

                }
         ## }
      ### Set up an empty vector for angular deviation:
      #if (print.angle) angle <- rep(NA, nyears)

      ### Make a label plot for the fishery (index)
      if(compact)
         {  ### But if it would be in the last frame, save the graph & move to new page:
         if (PlotsWaiting == PlotsPerPage - 1)
         {  frame()
            if (write.graphs)
            {  gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
               FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                  graphics.type)
            }  # end if (write.graph)
            PlotsWaiting = 0
            PlotsPage = PlotsPage + 1
         }  # end if(PlotsWaiting = PlotsPerPage - 1)

         plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
         text(0.5, 0.5, bquote(phantom(.) %down% phantom(.) ~~ .(fileroot)
            ~~ phantom(.) %down% phantom(.)), adj=c(0.5, 0.5))
         box()
         PlotsWaiting = PlotsWaiting + 1
         }  # end (if compact)

     for(iyear in 1:nyears)                ### By year
     {   thisyear = years[iyear]
         ### Set main title:
         if (draft)
         {  mtitle <- paste(titleroot, "   Year:", thisyear) }
         else
         {  mtitle <- "" }
         ## Set x-axis title:
         if(ptype == "l")
            {  title.x <- FGMakeLabel("Length bin", units) }
         else
         {  title.x <- "Age class" }

         if (compute.neff)
         {    num <- den <- 0.0
             den <- sum( (prop.obsd[iyear,] - prop.pred[iyear,])^2)
             num <- sum ( prop.pred[iyear,] * (1 - prop.pred[iyear,]))
             Neff[iyear] <- round(num/den, digits=1)
         }

         ## Plot axes and notations:
         if (! uniform) ymax <- max(c(prop.obsd[iyear,], prop.pred[iyear,]))
         #if (length(bins) > 4) {
         #   xlim <- range(bins)
         #   }
         #else { if (length(bins) == 1) xlim <- c(bins[1] - 2, bins[1] + 2) else
         #       xlim <- c(bins[1], bins[1] + (bins[2]-bins[1]) * 4)}
         if ((!is.na(Nobs[iyear]) & Nobs[iyear]>0) & (!is.na(Neff[iyear]) & Neff[iyear]>0) | plot.missing ) { #Only print out if not na sample size or if plot.missing is TRUE
         plot(bins, prop.pred[iyear,], ylim = c(0, ymax),
            xlab = title.x, ylab = "Proportion", main = mtitle, type = "n")
         grid(col = "lightgray", lty = 1)
         ### Plot data:
         if(plot.pred) lines (bins, prop.pred[iyear,], type = "l", lty = 1, col = col.pred, lwd = 2)
         points(bins, prop.obsd[iyear,], type = "p", pch = 21, cex = 1.25, col = col.obsd,
            lwd = 2)
         if (connect.obsd) points(bins, prop.obsd[iyear,], type = "l", col = col.obsd,
            lwd = 1)
         PlotsWaiting <- PlotsWaiting + 1
         ### Add optional notations:
         mypar = par()
         # Set user scale to plot size in inches:
         par(usr = c(0, mypar$fin[1], 0, mypar$fin[2]))
         xloc <- 0.975 * mypar$fin[1]
         yloc <- 0.980 * mypar$fin[2]
         # Set linespacing to 170% of typical text size in inches:
         vspace <- strheight("M/gh3", units="inches") * 1.7
#          if (print.angle)
#             {  ### Add notation of angular deviation
#                angle[iyear] <- (vecAngle(prop.obsd[iyear,],prop.pred[iyear,], degrees = TRUE))
#                anntext <- bquote(Deviation == .(round(angle[iyear], digits = 1)) *degree)
#                text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
#                yloc <- yloc - vspace
#             }
         if (print.n)
            {  ### Add notation of number observed
               anntext <- bquote( italic(N) == .(Nobs[iyear]))
               text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
               yloc <- yloc - vspace
            }
         if (print.neff || plot.neff)
            {  ### Add effective N and ratio Neff/N to plot

               if (print.neff)
               {  anntext <- bquote(Effective ~~ italic(N) == .(Neff[iyear]))
                  text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
                  yloc <- yloc - vspace
                  #anntext <- bquote(italic(N)[eff] / italic(N) ==
                  #   .(round(Neff[iyear]/Nobs[iyear], digits=2)))
                  #text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
                                        #yloc <- yloc - vspace
                  if (Neff[iyear]<=0|is.na(Neff[iyear])==TRUE)
               {
                   polygon(c(0,0,mypar$fin[1],mypar$fin[1]),c(0,mypar$fin[2],mypar$fin[2],0),col=rgb(.211, .211, .211, .2))
               }

               }
            }
         ## background shading for years not fitted
         if (is.na(Nobs[iyear])|Nobs[iyear]==0)
         {
             polygon(c(0,0,mypar$fin[1],mypar$fin[1]),c(0,mypar$fin[2],mypar$fin[2],0),col=rgb(.211, .211, .211, .2))
         }
         if (print.year)
            {  anntext <- thisyear
               text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
               yloc <- yloc - vspace
            }

         ### Write plot to file(s)
         if (write.graphs)
            {  if ( !compact || (PlotsWaiting == PlotsPerPage) )
               {  if (compact)
                        gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
                     else
                        gname  = paste(fileroot, ".", thisyear, sep = "")
                  FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                     graphics.type)
                  PlotsWaiting <- 0
                  PlotsPage = PlotsPage + 1
                  }
            }
         }                              #end if statement for sample size
         }  # end for iyear in 1:nyears

         if (plot.neff && nyears > 1)
         {  ### Plot relative effective sample size over time
            lab.y <- expression(italic(N)[eff] / italic(N))
            max.y <- max(Neff / (Nobs+1e-6), 1.0,na.rm=TRUE)
            if (!draft) mtitle <- ""
            else mtitle <- FGMakeTitle(fileroot, DataName)
            # Prevent plotting fractional years:
            if (length(years) == 1)
                {   xlim <- years[1] + c(-2, 2)
                }
            else
                { if (range(years)[2]-range(years)[1] > 4)
                    {   xlim <- range(years)
                    }
                    else
                    {   xlim <- years[1] + c(0, 4)
                    }
                }
            FGTimePlot(years, Neff/(Nobs+1e-6),  lab.x = "Year", lab.y = lab.y,
               href = NULL, hrefnames = NULL, Y1Col = Y1Col, ylim = c(0, max.y),
               xlim = xlim, FGtype = "stick", main = mtitle)
            PlotsWaiting <- PlotsWaiting + 1
            ### Write plot to file(s)
            if (write.graphs)
            {  if ( !compact || (PlotsWaiting == PlotsPerPage) )
               {  if (compact)
                  {  gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
                  }
                  else
                  {   gname  = paste(fileroot, ".nratio", sep = "")
                  }
                  FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                     graphics.type)
                  PlotsWaiting <- 0
                  PlotsPage = PlotsPage + 1
               }
            }     # end if (write.graphs)
         }     # end if (plot.neff && ...)

     } # end for jindex in 1:nindexes

            ### Write any remaining plots to file(s)
            if (write.graphs && compact && PlotsWaiting)
               {  gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
                  FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                     graphics.type)
               }  # end if (write.graphs && compact)

   par(savepar)     # reset graphics device
   return(invisible(0))
} # end function definition

