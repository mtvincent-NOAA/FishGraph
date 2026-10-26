#' Comparison plots of fits to time series from 2 assessment models
#'
#' Plots of observed and predicted time series from 2 assessment models
#'
#' @param x an R list with output from the assessment model
#' @param y an R list with output from another assesssment model
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param ModelNames a vector of string names for the names of models to be used
#'  in the legend of the plot. Defaults to arguments \code{c('x','y')}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param connect.obsd When \code{TRUE} a line connecting observed points is plotted
#' @param err.bar When \code{TRUE} error bars around the observation points are plotted
#'
#' @return Graphics
#'
#' @author M Vincent
#' @author N Klibanski
#'

#'
#' @examples \donttest{
#' Compare.tseries.plots(gag1,gag2)
#' }
#' @export Compare.tseries.plots
#'

########################################################################################
## R function to plot observed and predicted catch and indices from fish stock assessment models.
##
## Assumes a dataframe object called "t.series" is a component of "x" and "y"
##
## M.T. Vincent and N. Klibanski, October 2026
########################################################################################
Compare.tseries.plots <- function(x,y, DataName = deparse(substitute(x)),
                                  ModelNames = c('x','y'), draft = TRUE,
                                  graphics.type = NULL, use.color= TRUE,
                                  connect.obsd = FALSE, from.zero = TRUE,
                                  err.bar = TRUE)
######################################################################################
    ## ARGUMENTS
    ## x - an R list output from the assessment models
    ##    The list x must have a component x$t.series that is a data frame
    ##    The components in x should match with those in y
    ## y - an R list output from the assessment models
    ##    The list y must have a component y$t.series that is a data frame
    ##    The components in y should match with those in x
    ## DataName - a string representation an identifier for the titles
    ## ModelNames - a vector of strings used in legend to distinguish models
    ##    if nothing present defaults to c('x','y')
    ## draft - TRUE if figures are to have main titles
    ## graphics.type - a character vector with graphics-file types, e.g., c("png","eps")
    ## use.color - TRUE if graphs of index fits are in color
    ## connect.obsd - TRUE to connect the observed points
    ## from.zero - TRUE to start the Y-axis of each plot from zero
    ## err.bar - TRUE to plot the error bars around the observation if exists

##########################################################################################
{
    ### Check for x$t.series data frame:
    if (! ("t.series" %in% names(x))) {
        ErrText = (paste("Component ", deparse(substitute(x)),
            "$t.series not found.", sep = ""))
        warning(ErrText, immediate. = TRUE)
        return(invisible(-1))
    }
    ### Check for y$t.series data frame:
    if (! ("t.series" %in% names(y))) {
        ErrText = (paste("Component ", deparse(substitute(y)),
            "$t.series not found.", sep = ""))
        warning(ErrText, immediate. = TRUE)
        return(invisible(-1))
    }

    ## Assume units of the two models are the same. If needs to be added as an arguement this can be added later, just ask
    L.units = x$info$units.landings
    D.units = x$info$units.discards

    ### Local copy of data:
    xts <- x$t.series
    ### Local copy of data:
    yts <- y$t.series

    ### Check for other needed data:
    if (! ("year" %in% names(x$t.series))) {
        warning("Year variable not found in Index.plots!", immediate. = TRUE)
        return(invisible(-1))
    }
    ### Check for other needed data:
    if (! ("year" %in% names(y$t.series))) {
        warning("Year variable not found in Index.plots!", immediate. = TRUE)
        return(invisible(-1))
    }

    ### Check that x and y share the sane
    if (nrow(xts)!=nrow(yts)) {
        warning("The length of t.series for models x and y are not idenitcal!",
                immediate. = TRUE)
        return(invisible(-1))
    }

    t.series.names <- unique(c(names(xts),names(yts)))
    t.series.year <- sort(unique(xts$year,yts$year))

    ## determine names of t.series variables associated with fits
    t.series.names.ob <- t.series.names[grepl(t.series.names,pattern=".ob")]
    t.series.names.pr <- t.series.names[grepl(t.series.names,pattern=".pr")]

    nplots <- length(t.series.names.ob) #Number of distinct indices

    ## Set graphics parameters
    savepar <- FGSetPar(draft)
    plot.options <- FGGetOptions()
    ## Set the colors
    if (use.color) {
        col.obsd <- plot.options$color$clr.obsd
        col.pred <- plot.options$color$clr.line
        col.pred2  <- plot.options$color$clr.line2
        col.obsd2  <- plot.options$color$clr.obsd2
        err.bar.col <- plot.options$color$clr.obsd
        err.bar.col2 <- plot.options$color$clr.obsd2
    } else {
        col.obsd <- plot.options$bw$clr.obsd
        col.pred <- plot.options$bw$clr.line
        col.pred2  <- plot.options$bw$clr.line2
        col.obsd2  <- plot.options$bw$clr.obsd2
        err.bar.col <- plot.options$bw$clr.obsd
        err.bar.col2 <- plot.options$bw$clr.obsd2
    }


    ### If writing graphics files, make sure there is a directory for them:
    if (! is.null(graphics.type)){
        write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/compare", sep="")
        FGCheckGraphDir(GraphicsDirName)
    } else {write.graphs <- FALSE }

    if (connect.obsd) FGtype <- "2line" else FGtype <- "linepoint"

    for (iplot in t.series.names.ob) {
        ob.x <- xts[,iplot]
        pr.name <- gsub(".ob",".pr",iplot)
        pr.x <- xts[,pr.name]
        ob.y <- yts[,iplot]
        pr.y <- yts[,pr.name]

        tmp <- is.na(cbind(ob.x,ob.y)) #matrix of TRUE (missing or FALSE (not missing)
        temp <- apply(tmp, 1, all)     #Vector with TRUE if all cols are TRUE (missing)
        first <- which.min(temp)       #First year with any column not NA
        yrndx <- first:nrow(xts)
        nyear <- length(yrndx)
        ob.x <- ob.x[yrndx]
        ob.y <- ob.y[yrndx]
        pr.x <- pr.x[yrndx]
        pr.y <- pr.y[yrndx]
        years <- xts$year[yrndx]

        ## Get name of quantity plotted and make plot title
        QuantType <-switch(substr(iplot, 1, 1),
                           U = "Index",
                           L = "Fishery Landings",
                           D = "Fishery Discards")
        QuantName <- FGTrimName(iplot, removePrefix = 1, removeSuffix = 1)

        if (draft) PlotTitle <- FGMakeTitle(paste(QuantType,QuantName, sep = ": "), DataName) else PlotTitle <-""

        ob.maxrange <- c(ob.x,ob.y)
        if (err.bar) {
            CVName <- paste0("cv.",FGTrimName(iplot, removeSuffix = 1))
            CVxcol <- which(names(xts) == CVName)
            if (length(CVxcol) != 0) {
                xSD <- xts[yrndx,CVxcol] * ob.x
                ob.x.lower <- ob.x - 2*xSD
                ob.x.upper <- ob.x + 2*xSD
            } else {ob.x.lower <- ob.x.upper <- NULL}
            CVycol <- which(names(yts) == CVName)
            if (length(CVycol) != 0) {
                ySD <- yts[yrndx,CVycol] * ob.y
                ob.y.lower <- ob.y - 2*ySD
                ob.y.upper <- ob.y + 2*ySD
            } else {ob.y.lower <- ob.y.upper <- NULL}
            ob.maxrange <- c(ob.x.upper,ob.y.upper,ob.x.lower,ob.y.lower)
        }

        ## Scale Y-axis:
        if (from.zero) {yrange = c(min(0,ob.maxrange, pr.x,pr.y, na.rm = TRUE),max(ob.maxrange, pr.x,pr.y, na.rm = TRUE))
        } else {yrange <- range(ob.maxrange, pr.x, pr.y, na.rm = TRUE) }
        ## Y-axis label
        YLab <-switch(substr(iplot, 1, 1),
                      U = "Relative abundance (CPUE)",
                      L = FGMakeLabel("Landings", L.units),
                      D = FGMakeLabel("Discards", D.units))

        ## Plot of observed and predicted for the two models
        FGTimePlot(x = years, y = pr.x, y2 = ob.x, lab.x = 'Year',
            lab.y = YLab, href = NULL, hrefnames = NULL,
            use.color = use.color, FGtype = FGtype, Y1Col=col.pred,
            Y2Col=col.obsd, ylim = yrange, main = PlotTitle)
        if (length(ob.x.upper)!=0){
          for (ibars in 1:length(years)){
            lines(c(years[ibars],years[ibars]),c(ob.x.lower[ibars], ob.x.upper[ibars]), col=err.bar.col)
          }
        }

        ## Add in y observed and predicted
        points(x = years, y = pr.y, type = 'o', pch = 16, col = col.pred2, cex = 1.0,lty=2,lwd=2)
        if (FGtype == "2line") lines(years, ob.y, lwd = 2, type = "b", pch = 1,
                                     col = col.obsd2, cex = 1.5)
        if ((substr(FGtype, 1, 9) == "linepoint") && (!is.null(ob.y))){
            points(years, ob.y, lwd = 2, pch = 1, col = col.obsd2, cex = 1.5)}
        if (length(ob.y.upper)!=0){
          for (ibars in 1:length(years)){
            lines(c(years[ibars],years[ibars]),c(ob.y.lower[ibars], ob.y.upper[ibars]), col=err.bar.col2, lty=3)
          }
        }
        Legends=paste(rep(ModelNames,each=2),rep(c('ob','pr'),2))
        colvec = c(col.obsd, col.pred, col.obsd2, col.pred2)
        ## Add in legend
        if (err.bar) {
            if (connect.obsd) {
                legend('topleft',,Legends,col=colvec,lwd=2,lty=c(1,1,1,2),pch =c(1,16,1,16),bty='n',pt.cex=c(1.5,1,1.5,1))
            } else { legend('topleft',,Legends,col=colvec,lwd=c(1,2,1,2),lty=c(1,1,3,2),pch =c(1,16,1,16),bty='n',pt.cex=c(1.5,1,1.5,1)) }
        } else{
            if (connect.obsd) {
                legend('topleft',,Legends,col=colvec,lwd=2,lty=c(1,1,1,2),pch =c(1,16,1,16),bty='n',pt.cex=c(1.5,1,1.5,1))
            } else { legend('topleft',,Legends,col=colvec,lwd=2,lty=c(0,1,0,2),pch =c(1,16,1,16),bty='n',pt.cex=c(1.5,1,1.5,1))}
        }

        if (write.graphs)
        {   FGSavePlot(GraphicsDirName, DataName, GraphName =
                paste(QuantType,QuantName, sep = " "), graphics.type)
        }
    } #End iplot
}
