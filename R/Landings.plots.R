#' Landings and discards trajectories
#'
#' The function \code{Landings.plots} provides time trajectories of landings and
#' discards by fishery. Optional arguments specify whether data represent
#' observed and predicted landings or simply independent series of landings.
#' The same holds for discards.
#'
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}.
#' @param from.zero When \code{TRUE}, the Y-axis of each plot starts at zero.
#' @param L.units Character vector containing units of measure associate with
#' each landings series found in \code{x$t.series}.
#' @param D.units Same as the preceding, but for discards.
#' @param start.drop Number of years at the start of the data to be omitted from
#' plots, as when a model includes an initialization period.
#' @param L.obs.pre When \code{TRUE}, landings data are interpreted as pairs of
#' columns (observed and predicted). When \code{FALSE}, each column represents
#' observed landings.
#' @param D.obs.pre Same as the preceding, but for discards.
#' @param notconfidential a dataframe with names the same as landings and entries are booleans where true is confidential data. Will need to be input.
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
#' Landings.plots(gag)
#' }
#' @export
Landings.plots <-
function(x, DataName = deparse(substitute(x)), draft = TRUE,
    graphics.type = NULL, use.color = TRUE, from.zero = TRUE,
    L.units = x$info$units.landings, D.units = x$info$units.discards,
    start.drop = 0, L.obs.pre = TRUE, D.obs.pre = TRUE,confidential=NULL)
{
    ### Check for required data:
    if (! ("t.series" %in% names(x)))
    {  Errmsg <- paste("Data frame t.series not found in", deparse(substitute(x)))
      warning(Errmsg, immediate. = TRUE)
      return(invisible(-1))
    }
    if (! "year" %in% names(x$t.series))
    {  Errmsg <- paste("Column 'year' not found in ", deparse(substitute(x)),
         "$t.series", sep = "")
      warning(Errmsg, immediate. = TRUE)
      return(invisible(-1))
    }

    ### Local copies of data:
    year <- x$t.series$year
    # Get columns beginning with "L."
    landings <-  x$t.series[, grep("^L\\.", names(x$t.series))]
    discards <-  x$t.series[, grep("^D\\.", names(x$t.series))]
    if(start.drop > 0)
    {   year <- year[-(1:start.drop)]
        landings <- landings[-(1:start.drop),]
        discards <- discards[-(1:start.drop),]
    }
    nlandings <- ncol(landings)
    ndiscards <- ncol(discards)
    ## Check confidential is same length as landings and discards
    if(dim(landings)[1]!=dim(confidential)[1])
    {   warning("number of years in Landings do not match number of years in confidential.", immediate. = TRUE)
        return(invisible(-1))
    }

    ### Does data set include any landings?
    if (nlandings == 0)
    {   warning("No landings data (x$t.series.L*) found.", immediate. = TRUE)
        return(invisible(-1))
    }
    ### if (L.obs.pre) there should be an even number of landings columns
    if ( L.obs.pre && (nlandings %% 2) != 0 )
    {   warning("Odd number of Landings columns", immediate. = TRUE)
        return(invisible(-1))
    }
    if ( D.obs.pre && (ndiscards %% 2) != 0 )
    {   warning("Odd number of discards columns", immediate. = TRUE)
        return(invisible(-1))
    }
    #  Find the indices of the years to be plotted -- same for all plots.
    #  Start: the first year with a non-missing value of any index.
    #  Stop: the last year in the data.
    temp <- apply(is.na(cbind(landings, discards)), 1, all)  # Vector w/ TRUE if all cols NA
    first  <- which.min(temp)                 # First year with any data not NA
    ss     <- first:nrow(landings)
    landings <- landings[ss,]
    discards <- discards[ss,]
    confidential <- confidential[ss,]
    year <- year[ss]
    if (L.obs.pre) nlandings <- nlandings / 2
    if (D.obs.pre) ndiscards <- ndiscards / 2

    # Expand units if only one value given
    if (length(L.units) == 1) L.units <- rep(L.units, nlandings)
    if (length(D.units) == 1) D.units <- rep(D.units, ndiscards)

    # Set graphics parameters
    savepar <- FGSetPar(draft)
    PlotTitle <- ""

    if (! is.null(graphics.type))
    {   write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/L", sep="")
    }
    else
    {   write.graphs <- FALSE }

    #===========================================================================
    ### Plots of observed and predicted landings:
    for (iplot in 1:nlandings)
    {   if (L.obs.pre)
        {   FGtype = "linepoint"
            y1 <- landings[, 2 * iplot]      # predicted L
            y2 <- landings[, 2 * iplot - 1]  # observed L
            series = FGTrimName(names(landings)[2 * iplot], removePrefix = 0,
                                removeSuffix = 1)
            if (series %in% names(confidential)) {
                y2[confidential[,series]] <- NA  # confidential data don't plot
            }
        }
        else
        {   FGtype = "1line"
            y1 = landings[, iplot]
            y2 = NULL
            series = names(landings)[iplot]
            if (series %in% names(confidential)) {
                y1[confidential[,series]] <- NA  # confidential data don't plot
            }
        }
        if (from.zero) yrange = c(0, max(y1, y2, na.rm = TRUE))
        else yrange <- range(y1, y2, na.rm = TRUE)

        if (draft) PlotTitle <-FGMakeTitle(paste("Fishery: ", series),DataName)

        FGTimePlot(x = year, y = y1, y2 = y2, lab.x = "Year",
            lab.y = FGMakeLabel("Landings", L.units[iplot]), ylim = yrange,
            FGtype = FGtype, main = PlotTitle, use.color = use.color)
        if (write.graphs)
        {   FGSavePlot(GraphicsDirName, DataName, GraphName =
                series, graphics.type)
        }
   } # End for iplot
    # ===================================================================================
    # Plots of discards:

    if (ndiscards > 0)
    {   for (iplot in 1:ndiscards)
        {   if (D.obs.pre)
            {   FGtype = "linepoint"
                y1 <- discards[, 2 * iplot]        # predicted D
                y2 <- discards[, 2 * iplot - 1]    # observed D
                series = FGTrimName(names(discards)[2 * iplot], 0, 1)
                if (series %in% names(confidential)) {
                    y2[confidential[,series]] <- NA  # confidential data don't plot
                }
            }
            else
            {   FGtype = "1line"
                y1 = discards[, iplot]
                y2 = NULL
                series = names(discards)[iplot]
                if (series %in% names(confidential)) {
                    y1[confidential[,series]] <- NA  # confidential data don't plot
                }
            }
            if (from.zero)
            {   yrange = c(0, max(y1, y2, na.rm = TRUE))
            }
            else
            {   yrange <- range(y1, y2, na.rm = TRUE)
            }

            if (draft) PlotTitle <- FGMakeTitle(paste("Fishery: ", series), DataName)

            FGTimePlot(x = year, y = y1, y2 = y2, lab.x = "Year",
                lab.y = FGMakeLabel("Discards", D.units[iplot]), ylim=yrange,
                FGtype = FGtype, main = PlotTitle, use.color = use.color)
            if (write.graphs)
            {   FGSavePlot(GraphicsDirName, DataName, GraphName = series, graphics.type)
            }
        }   # end for loop
     }  # end if (ndiscards > 0)

    par(savepar)
    return(invisible(NULL))
}  # End of function

