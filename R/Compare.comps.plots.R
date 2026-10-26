#' Plots comparing fits to age and length composition between two models
#'
#' The functions \code{Compare.comps.plots} generates plots of the observed
#' age and length composition data and the predicted compositions  for two models.
#' The function checks that the data are between the two models are identical.
#' If data are not identical for the observed compositions of the same name between
#' the two models the function will not plot that data set and will give a warning.
#' The plots generated include annual plots and weighted aggregate plots
#'
#' @param x an R list with output from the assessment models.
#' @param y an R list with output from the assessment models with names the same as x
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param ModelNames a vector of string names for the names of models to be used
#' in the legend of the plot. Defaults to arguments \code{c('x','y')}.
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
#' @param plot.neff When \code{TRUE}, a timeplot of \emph{N_eff} is made for each data
#' series
#' @param print.year When \code{TRUE}, the year is printed in each plot.
#' @param compute.neff When \code{TRUE}, \strong{FishGraph} computes \emph{N_eff}
#' internally, rather than using values stored by the user.  The computations are
#' based on the sample sizes and observed and estimated proportions.
#' @param compact When \code{TRUE}, plots are arranged by default in 5x3 matrix on each page.
#' @param uniform When \code{TRUE}, all years of the same data series are scaled
#' the same.
#' @param connect.obsd When \code{TRUE}, observed points are connected by lines.
#' @param plot.pooled When \code{TRUE}, plots pooled weighted compositions
#' @param Nyrs when \code{TRUE}, plots years of pooled weighted compositions
#'
#' @return Graphics
#'
#' @author M Vincent
#' @author M Prager
#' @author E Williams
#' @author K Shertzer
#' @author R Cheshire
#' @author K Purcell
#'
#'
#' @examples \donttest{
#' Compare.comps.plots(gag,gag)
#' }
#' @export Compare.comps.plots
#'

Compare.comps.plots <- function(x, y, DataName = deparse(substitute(x)),
                                ModelNames = c('x','y'), draft = TRUE, graphics.type = NULL,
                                use.color = TRUE, units = x$info$units.length,
                                print.n = !compact, print.neff = !compact,
                                plot.missing = FALSE, plot.neff = FALSE, print.year = TRUE,
                                compute.neff = FALSE, compact = TRUE, uniform = TRUE,
                                connect.obsd = FALSE, plot.pooled = TRUE, Nyrs = FALSE)
{
    ##-------------------------------------------------------------------------------
    ## Check for data presence:
    Errstring = ("No composition data found.  Terminating Comp.yearly.plots");
    if (! ("comp.mats" %in% names(x))) stop(Errstring)
    if (! ("comp.mats" %in% names(y))) stop(Errstring)

    ## Make local copies of list of matrices:
    cmx <- x$comp.mats
    cmy <- y$comp.mats
    comp.names <- unique(names(cmx),names(cmy))
    obs.names <- comp.names[grep("[.]ob",comp.names)]

    ## Set graphics parameters for compact or non-compact mode:
    plot.options = FGGetOptions()
    savepar <- FGSetPar(draft)
    if (compact) {
        par(mar = c(3.5, 4, 0.85, 0.7), cex = 0.9, cex.main = 0.9, cex.lab = 0.85,
            cex.axis = 0.75, mfcol = plot.options$graphics$fglayout, lab = c(5, 2, 4),
            mgp=c(1.8, 0.75, 0), las = 0)
        draft1 <- FALSE  ### Plots have no main titles in compact mode
    } else {
        par(las = 1)
        draft1 <- draft
    }
    PlotsPerPage <- prod(par("mfrow"))
    PlotsWaiting <- 0
    PlotsPage <- 1
    ## Set colors:
    if (use.color) {
        col.obsd <- plot.options$color$clr.obsd
        col.pred  <- plot.options$color$clr.line
        ## col.obsd2  <- plot.options$color$clr.obsd2
        col.pred2  <- plot.options$color$clr.line2
    } else {
        col.obsd <- plot.options$bw$clr.obsd
        col.pred <- plot.options$bw$clr.line
        ## col.obsd2  <- plot.options$bw$clr.obsd2
        col.pred2  <- plot.options$bw$clr.line2
    }
    ## If writing graphics files, make sure there is a directory for them:
    if (length(graphics.type > 0)) {
        write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/compare", sep="")
        FGCheckGraphDir(GraphicsDirName)
    } else {  write.graphs <- FALSE}

    ## Make plots:
    for (obname in obs.names)  {               ### By data series
        ## Get the corresponding predicted name
        prname <-gsub("[.]ob", ".pr", obname)
        ## pull out the data by name of source from both models
        prop.obsd <- cmx[[obname]]
        prop.pred <- cmx[[prname]]
        prop.obsd2 <- cmy[[obname]]
        prop.pred2 <- cmy[[prname]]

        ## Check that both observed data sources are not null and if not identical skip to next data source
        if (!is.null(prop.obsd) & !is.null(prop.obsd2)) {
            if (!identical(prop.obsd,prop.obsd2)) { next }}

        ymax <- min(max(prop.obsd,prop.pred,prop.pred2, na.rm=TRUE),
                    1.3*quantile(c(prop.obsd, prop.pred, prop.pred2), 0.98, na.rm=TRUE))
        years     <- as.numeric(rownames(prop.obsd))
        nyears    <- length(years)
        bins      <- as.numeric(colnames(prop.obsd))
        fileroot  <- FGTrimName(obname, removePrefix = 0, removeSuffix = 1)
        ptype     <- substr(fileroot, 1, 1)        # ptype = "l" or "a" for length or age
        fishery   <- FGTrimName(fileroot, removePrefix = 1, removeSuffix = 0)
        titleroot <- paste( "Data: ", DataName, "    Sample: ", fishery, sep = "")
        ## Get Nobs for this matrix from the x$t.series dataframe:
        ## if (print.n || print.neff || plot.neff) {
            yearvec = rownames(prop.obsd)
            strtmp = paste("x$t.series$", fileroot, ".n[x$t.series$year %in% yearvec]", sep="")
            Nobs = eval(parse(text = strtmp))
            yearvec2 = rownames(prop.obsd2)
            strtmp = paste("y$t.series$", fileroot, ".n[y$t.series$year %in% yearvec2]", sep="")
            Nobs2 = eval(parse(text = strtmp))

        ## } else {Nobs=rep(NA,nyears)}
        ## Get Neff for this matrix from the x$t.series dataframe IF Neff isn't computed internally.
        ## otherwise, set up a blank vector for use when it IS computed:
        ## if (print.neff || plot.neff) {
            if (compute.neff) {
                Neff <- rep(NA_integer_, nyears)
                Neff2 <- rep(NA_integer_, nyears)
            } else {
                yearvec = rownames(prop.obsd)
                strtmp <- paste("x$t.series$", fileroot,
                                ".neff[x$t.series$year %in% yearvec]", sep="")
                Neff <- eval(parse(text = strtmp))
                if (is.null(Neff) || any(is.na(Neff))) {
                    warning("Comp.yearly.plots: comp.neff=FALSE but no data or NA in\n",
                            paste("   x$t.series$", fileroot, ".neff", sep = ""), call. = FALSE)
                    Neff <- rep(NA_integer_, nyears)
                } else {
                    Neff<-round(Neff,digits=1)
                    Neff[Neff < 0] <- 0
                }
                ## Now for the second Neff
                yearvec2 = rownames(prop.obsd2)
                strtmp2 <- paste("y$t.series$", fileroot,
                                ".neff[y$t.series$year %in% yearvec2]", sep="")
                Neff2 <- eval(parse(text = strtmp2))
                if (is.null(Neff2) || any(is.na(Neff2))) {
                    warning("Comp.yearly.plots: comp.neff=FALSE but no data or NA in\n",
                            paste("   y$t.series$", fileroot, ".neff", sep = ""), call. = FALSE)
                    Neff2 <- rep(NA_integer_, nyears)
                } else {
                    Neff2 <- round(Neff2, digits=1)
                    Neff2[Neff2 < 0] <- 0
                }
            }
        ## }

        ## Make a label plot for the fishery (index)
        if(compact) {
            ## But if it would be in the last frame, save the graph & move to new page:
            if (PlotsWaiting == PlotsPerPage - 1) {
                frame()
                if (write.graphs) {
                    gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
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

        for(iyear in 1:nyears) {               ### By year
            thisyear = years[iyear]
            ## Set main title:
            if (draft1) {
                mtitle <- paste(titleroot, "   Year:", thisyear)
            } else { mtitle <- "" }
            ## Set x-axis title:
            if(ptype == "l") {
                title.x <- FGMakeLabel("Length bin", units)
            } else { title.x <- "Age class" }

            if (compute.neff) {
                num <- den <- 0.0
                den <- sum( (prop.obsd[iyear,] - prop.pred[iyear,])^2)
                num <- sum ( prop.pred[iyear,] * (1 - prop.pred[iyear,]))
                Neff[iyear] <- round(num/den, digits=1)
            }

            ## Plot axes and notations:
            if (! uniform) ymax <- max(c(prop.obsd[iyear,], prop.pred[iyear,]))

            if ((!is.na(Nobs[iyear]) & Nobs[iyear]>0) & (!is.na(Neff[iyear]) & Neff[iyear]>0) | plot.missing ) {
                ## Only print out if not na sample size or if plot.missing is TRUE
                plot(bins, prop.pred[iyear,], ylim = c(0, ymax),
                     xlab = title.x, ylab = "Proportion", main = mtitle, type = "n")
                grid(col = "lightgray", lty = 1)
                ## Plot data:
                lines (bins, prop.pred[iyear,], type = "l", lty = 1, col = col.pred, lwd = 2)
                lines (bins, prop.pred2[iyear,], type = "l", lty = 2, col = col.pred2, lwd = 2)
                points(bins, prop.obsd[iyear,], type = "p", pch = 21, cex = 1.25, col = col.obsd, lwd = 2)
                if (is.null(prop.obsd)) { points(bins, prop.obsd2[iyear,], type = "p", pch = 21, cex = 1.25, col = col.obsd, lwd = 2) }
                if (connect.obsd) points(bins, prop.obsd[iyear,], type = "l", col = col.obsd, lwd = 1)
                PlotsWaiting <- PlotsWaiting + 1
                ## Add optional notations:
                mypar = par()
                ## Set user scale to plot size in inches:
                par(usr = c(0, mypar$fin[1], 0, mypar$fin[2]))
                xloc <- 0.975 * mypar$fin[1]
                yloc <- 0.980 * mypar$fin[2]
                ## Set linespacing to 170% of typical text size in inches:
                vspace <- strheight("M/gh3", units="inches") * 1.7

                if (print.n) {
                    ## Add notation of number observed
                    anntext <- bquote( italic(N)~1 == .(Nobs[iyear])~~ italic(N)~2 == .(Nobs2[iyear]))
                    text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
                    yloc <- yloc - vspace
                }
                if (print.neff || plot.neff) {
                    ## Add effective N and ratio Neff/N to plot
                    if (compute.neff) {
                        num <- den <- 0.0
                        den <- sum( (prop.obsd[iyear,] - prop.pred[iyear,])^2)
                        num <- sum ( prop.pred[iyear,] * (1 - prop.pred[iyear,]))
                        Neff[iyear] <- round(num/den, digits=1)
                        ## Now for second Neff
                        num <- den <- 0.0
                        den <- sum( (prop.obsd2[iyear,] - prop.pred2[iyear,])^2)
                        num <- sum ( prop.pred2[iyear,] * (1 - prop.pred2[iyear,]))
                        Neff2[iyear] <- round(num/den, digits=1)
                    }
                    if (print.neff) {
                        anntext <- bquote(Effective ~~ italic(N)~1 == .(Neff[iyear])~~ italic(N)~2 == .(Neff2[iyear]))
                        text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
                        yloc <- yloc - vspace
                        if (Neff[iyear]<=0|is.na(Neff[iyear])|Neff2[iyear]<=0|is.na(Neff2[iyear])) {
                            polygon(c(0,0,mypar$fin[1],mypar$fin[1]),c(0,mypar$fin[2],mypar$fin[2],0),col=rgb(.211, .211, .211, .2))
                        }
                    }
                }
                ## background shading for years not fitted
                if (is.na(Nobs[iyear])|Nobs[iyear]==0|is.na(Nobs2[iyear])|Nobs2[iyear]==0) {
                    polygon(c(0,0,mypar$fin[1],mypar$fin[1]),c(0,mypar$fin[2],mypar$fin[2],0),col=rgb(.211, .211, .211, .2))
                }
                if (print.year) {
                    anntext <- thisyear
                    text(xloc, yloc, anntext, adj = c(1, 1), cex = 0.9)
                    yloc <- yloc - vspace
                }

                ## Write plot to file(s)
                if (write.graphs) {
                    if ( !compact || (PlotsWaiting == PlotsPerPage) ) {
                        if (compact) {
                            gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
                        } else { gname  = paste(fileroot, ".", thisyear, sep = "")}
                        ## Add legend to plot before saving
                        legend(xloc,yloc,c("Observed", ModelNames),col=c(col.obsd,col.pred,col.pred2),pch=c(21,NA,NA),lty=c(NA,1,2),lwd=2, adj=c(1,1))
                        FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                                   graphics.type)
                        PlotsWaiting <- 0
                        PlotsPage = PlotsPage + 1
                    }
                } else if (PlotsWaiting == PlotsPerPage) {
                    ## Add legend if last plot
                    legend(xloc,yloc,c("Observed", ModelNames),col=c(col.obsd,col.pred,col.pred2),pch=c(21,NA,NA),lty=c(NA,1,2),lwd=2, adj=c(1,1))
                }
            }                              #end if statement for sample size
        }  # end for iyear in 1:nyears

        if (plot.neff && nyears > 1) {
            ## Plot relative effective sample size over time
            lab.y <- expression(italic(N)[eff] / italic(N))
            max.y <- max(Neff / (Nobs+1e-6), 1.0, Neff2 / (Nobs2+1e-6), na.rm=TRUE)
            if (!draft1) { mtitle <- ""
            } else mtitle <- FGMakeTitle(fileroot, DataName)
            ## Prevent plotting fractional years:
            if (length(years) == 1){
                xlim <- years[1] + c(-2, 2)
            } else {
                if (range(years)[2]-range(years)[1] > 4) {
                    xlim <- range(years)
                } else {
                    xlim <- years[1] + c(0, 4)
                }
            }
            FGTimePlot(years, Neff/(Nobs+1e-6),  lab.x = "Year", lab.y = lab.y,
                       href = NULL, hrefnames = NULL, Y1Col = col.pred, ylim = c(0, max.y),
                       xlim = xlim, FGtype = "stick", main = mtitle)
            lines(years, Neff2/(Nobs2+1e-6), lwd = 2, type = "h", pch = pch, col = col.pred2, cex = 1.0, lty=2)
            lines(years, Neff2/(Nobs2+1e-6), lwd = 2, type = "p", pch = 21, col = col.pred2, cex = 1.0)
            PlotsWaiting <- PlotsWaiting + 1
            ## Write plot to file(s)
            if (write.graphs) {
                if ( !compact || (PlotsWaiting == PlotsPerPage) ) {
                    if (compact) {
                        gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
                    } else { gname  = paste(fileroot, ".nratio", sep = "")  }
                    legend('topright', ,ModelNames,col=c(col.pred,col.pred2),pch=c(16,21),lty=c(1,2),lwd=2)
                    FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                               graphics.type)
                    PlotsWaiting <- 0
                    PlotsPage = PlotsPage + 1
                }
            }  else{   # end if (write.graphs)
                if (PlotsWaiting == PlotsPerPage) {
                    ## Add legend if last plot
                    legend('topright', ,ModelNames,col=c(col.pred,col.pred2),pch=c(16,21),lty=c(1,2),lwd=2)
                }
            }
        }     # end if (plot.neff && ...)
    } # end for obname in obs.name

    ## Add legend if last plot
    legend(xloc,yloc,c("Observed", ModelNames),col=c(col.obsd,col.pred,col.pred2),pch=c(21,NA,NA),lty=c(NA,1,2),lwd=2, adj=c(1,1))

    ## Write any remaining plots to file(s)
    if (write.graphs && compact && PlotsWaiting)
    {  gname = paste("comp.page", sprintf("%02.0f", PlotsPage), sep = "")
        FGSavePlot(GraphicsDirName, DataName, GraphName = gname,
                   graphics.type)
    }  # end if (write.graphs && compact)

    par(savepar)     # reset graphics device

    #######################################################################################
    ##---Cumulative fit plots--------------------------------------------------------------

    selx=x$sel.age
    sely=y$sel.age
    if (plot.pooled){
        for (obname in obs.names) {
            ## Get the corresponding predicted name
            prname <-gsub("[.]ob", ".pr", obname)
            par(oma=c(2,2.5,0,0))
            if (draft) {par(mar = c(3, 4, 3, 1 ))
            } else{ par(mar = c(3, 1, 1, .5)) }

            m1x.all <- cmx[[obname]]         # Matrix of observed
            m2x.all <- cmx[[prname]]           # matrix of predicted
            m1y.all <- cmy[[obname]]         # Matrix of observed
            m2y.all <- cmy[[prname]]           # matrix of predicted

            ## Extract column names for plots
            binnames=colnames(m1x.all)
            if (!identical(binnames,colnames(m1y.all))) {
                warning(paste0("Column names for the two models do not match for ",obname,
                               ". Can not make plots, skipping for now."))
                next
            }

            ## Get various string representations of data series:
            ##  gfileroot is used in the name for the graphics file(s):
            gfileroot <- FGTrimName(obname, removePrefix = 0, removeSuffix = 1)
            if (gfileroot != FGTrimName(prname, removePrefix = 0, removeSuffix = 1)) {
                warning(paste("There is not an observed and predicted matrix for",gfileroot, "Skipping this plot for now"),sep=" ")
                next #skip to next name
            }
            groot=FGTrimName(obname, removePrefix = 1, removeSuffix = 1)
            gfilename <- paste("pooled.",gfileroot, sep="")
            ## test if a matrix or a vector

            if(!any(paste(c("sel.m","sel.v"),groot,sep=".") %in% names(selx))) {
                if(!any(paste(c("sel.m","sel.v"),groot,sep=".") %in% names(sely))) {
                    warning(paste0("Neither sel.m.",groot," nor sel.v.",groot," were found in the selectivity information in the rdats provided. Skipping this plot for now.") )
                    next
                }
            }


            if (paste("sel.m",groot,sep=".") %in% names(selx)){
                sel.blockx <- selx[[paste("sel.m",groot,sep=".")]]
                sel.blockx <- sel.blockx[rownames(sel.blockx)%in%rownames(m1x.all),,drop=FALSE]
                block.textx <- rownames(unique(as.data.frame(sel.blockx)))
            } else if (paste("sel.v",groot,sep=".") %in% names(selx)){
                block.textx = rownames(m1x.all)[1]
            }

            if (paste("sel.m",groot,sep=".") %in% names(sely)){
                sel.blocky <- sely[[paste("sel.m",groot,sep=".")]]
                sel.blocky <- sel.blocky[rownames(sel.blocky)%in%rownames(m1y.all),,drop=FALSE]
                block.texty <- rownames(unique(as.data.frame(sel.blocky)))
            } else if (paste("sel.v",groot,sep=".") %in% names(sely)){
                block.texty = rownames(m1y.all)[1]
            }

            yr.blockx <- as.numeric(block.textx)
            yr.blocky <- as.numeric(block.texty)

            ## Check if identical. If not skip for now. Might change this in the future
            if (! identical(yr.blockx,yr.blocky)) {
                warning(paste0("The selectivity time blocks for ", obname," were not identical. Skipping for now. If you have ideas on how to do this or need these plots contact Matt."))
                next
            }
            ##add last year to vector of years
            yr.block=c(yr.blockx,tail(x$t.series$year,n=1)+1)
            ## titleroot is used as part of the plot title:
            titleroot <- paste("Fishery: ", gfileroot)

            neffname<-paste(gfileroot,".neff", sep="")
            nname<-paste(gfileroot,".n", sep="")
            nseriesx<-nseriesy <-NULL
            if (neffname%in%names(x$t.series)) {
                nseriesx<-x$t.series[,names(x$t.series)==neffname]
            } else if(nname%in%names(x$t.series)) {
                nseriesx<-x$t.series[,names(x$t.series)==nname]
            }
            if (neffname%in%names(y$t.series)) {
                nseriesy<-y$t.series[,names(y$t.series)==neffname]
            } else if(nname%in%names(y$t.series)) {
                nseriesy<-y$t.series[,names(y$t.series)==nname]
            }

            if (is.null(nseriesx) & is.null(nseriesy)) {
                warning("There is no effective N or sample size in t.series x or y for the ",neffname," or ",nname,". Going to the next plot.")
                next
            }
            yrs.includex<-x$t.series$year[nseriesx>0 & !is.na(nseriesx)]
            yrs.includey<-y$t.series$year[nseriesy>0 & !is.na(nseriesy)]
            ## Figure out how many plots will actually be plotted
            subplotsx=length(which(sapply(1:length(block.textx),function(i) length(which(rownames(m1x.all[rownames(m1x.all)%in%yr.block[i]:(yr.block[i+1]-1),,drop=FALSE])%in%yrs.includex)))>0))
            subplotsy=length(which(sapply(1:length(block.texty),function(i) length(which(rownames(m1y.all[rownames(m1y.all)%in%yr.block[i]:(yr.block[i+1]-1),,drop=FALSE])%in%yrs.includey)))>0))
            subplots = max(subplotsx,subplotsy)

            ## Set the format based on the number of plots
            if(subplots<4){par(mfcol=c(subplots,1));cexLab=1#/.83
            } else if (subplots==4){par(mfrow=c(2,2));cexLab=1/.83
            } else if (subplots<=8){par(mfrow=c(ceiling(subplots/2),2));cexLab=1/.66
            } else if (subplots>9){par(mfrow=c(ceiling(subplots/3),3));cexLab=1/.66}
            if (subplots==3)cexLab=1/.66

            ##subset to years in each selectivity block (last block interval is defined differently because last value of yr.block is
            ## end of interval instead of beginning of next interval)
            for(iyrs in 1:(length(yr.blockx))){

                m1x.sub<-m1x.all[rownames(m1x.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]
                m2x.sub<-m2x.all[rownames(m2x.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]
                m1y.sub<-m1y.all[rownames(m1y.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]
                m2y.sub<-m2y.all[rownames(m2y.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]

                ##Exclude those yrs that don't make the min sample size requirement, designated by n<0
                if (dim(m1x.sub)[1]>1) {
                    m1x<-m1x.sub[rownames(m1x.sub)%in%yrs.includex,,drop=FALSE]
                    m2x<-m2x.sub[rownames(m2x.sub)%in%yrs.includex,,drop=FALSE]
                    wtx=x$t.series[[neffname]][!is.na(x$t.series[[neffname]])&x$t.series[[neffname]]>0&x$t.series$year%in%as.numeric(rownames(m1x))]
                    ## wt=wt[is.na(wt)==FALSE] # I think this is redundant
                } else {
                    m1x<-m1x.sub
                    m2x<-m2x.sub
                }
                if (dim(m1y.sub)[1]>1) {
                    m1y<-m1y.sub[rownames(m1y.sub)%in%yrs.includey,,drop=FALSE]
                    m2y<-m2y.sub[rownames(m2y.sub)%in%yrs.includey,,drop=FALSE]
                    wty=y$t.series[[neffname]][!is.na(y$t.series[[neffname]])&y$t.series[[neffname]]>0&y$t.series$year%in%as.numeric(rownames(m1y))]
                    ## wt=wt[is.na(wt)==FALSE] # I think this is redundant
                } else {
                    m1y<-m1y.sub
                    m2y<-m2y.sub
                }

                if (dim(m1x.sub)[1]>1 | dim(m1y.sub)[1]>1) {
                    poolylabtxt="wgt by effective N"
                } else { poolylabtxt="unweighted"}

                ## Set Y-axis title according to data type:
                if(substr(gfileroot, 1, 1)  == "l"){title.x <- FGMakeLabel("Length bin", units)
                }else{title.x <- "Age class"}

                ## allow plotting of vector instead of data frame object
                if (length(m1x)==0 & length(m1y)==0) {
                    warning(paste0("You do not have any samples to plot for ", gfileroot," in year block ",yr.block[iyrs],". Skipping this for now but you should probably check this in BAM!"))
                    next
                }

                if(dim(m1x)[1]>1){
                    ob.mx=apply(m1x,2,weighted.mean,w=wtx)
                    pr.mx=apply(m2x,2,weighted.mean,w=wtx)
                    x.vecx=as.numeric(dimnames(m1x)[[2]])
                    Poolyrsx=dim(m1x)[1]
                } else {
                    ob.mx=m1x
                    pr.mx=m2x
                    x.vecx=as.numeric(dimnames(m1x)[[2]])
                    Poolyrsx=1
                }

                if(dim(m1y)[1]>1){
                    ob.my=apply(m1y,2,weighted.mean,w=wty)
                    pr.my=apply(m2y,2,weighted.mean,w=wty)
                    x.vecy=as.numeric(dimnames(m1y)[[2]])
                    Poolyrsy=dim(m1y)[1]
                } else {
                    ob.my=m1y
                    pr.my=m2y
                    x.vecy=as.numeric(dimnames(m1y)[[2]])
                    Poolyrsy=1
                }

                ymax=max(c(ob.mx,pr.mx,ob.my,pr.my))

                if(!is.nan(ob.mx[1])){
                    plot(x=x.vecx,y=ob.mx,col=NULL,ylim=c(0,ymax),xlab="",ylab="",las=1,cex.axis=cexLab)
                    if (draft) title(main=titleroot)
                    grid(col = "lightgray", lty = 1)
                    points(x=x.vecx,y=ob.mx,col=col.obsd,pch=21, cex = 1.25, lwd=2)
                    lines(x=x.vecx,y=pr.mx,col=col.pred,lwd=2,type="l")
                    points(x=x.vecx,y=pr.mx,col=col.pred,lwd=2, pch=16)
                    lines(x=x.vecy,y=pr.my,col=col.pred2,lwd=2,type="l",lty=2)
                    points(x=x.vecy,y=pr.my,col=col.pred2,lwd=2, pch=17)
                    YrRange=paste(yr.block[iyrs],yr.block[iyrs+1]-1,sep="-")
                    legend(x="topright",bg="white",legend=c("Observed",ModelNames),
                           col = c(col.obsd,col.pred,col.pred2), lty = c(0,1,2),
                           lwd=2, pch = c(21, 16, 17),cex=cexLab,xjust=0.5,bty='n',
                           title=ifelse(Nyrs,paste(YrRange,paste0("Nyrs = ",Poolyrsx),sep="\n"),YrRange))
                    if (draft){mtext(paste("Pooled",poolylabtxt,sep="-"),side=2,line=4)}
                } else if (!is.nan(ob.my[1])){
                    plot(x=x.vecy,y=ob.my,col=NULL,ylim=c(0,ymax),xlab="",ylab="",las=1,cex.axis=cexLab)
                    if (draft) title(main=titleroot)
                    grid(col = "lightgray", lty = 1)
                    points(x=x.vecy,y=ob.my,col=col.obsd,pch=21, cex = 1.25, lwd=2)
                    lines(x=x.vecy,y=pr.my,col=col.pred2,lwd=2,type="l",lty=2)
                    points(x=x.vecy,y=pr.my,col=col.pred2,lwd=2, pch=17)
                    YrRange=paste(yr.block[iyrs],yr.block[iyrs+1]-1,sep="-")
                    legend(x="topright",bg="white",legend=c("Observed",ModelNames),
                           col = c(col.obsd,col.pred,col.pred2), lty = c(0,1,2),
                           lwd=2, pch = c(21, 16, 17),cex=cexLab,xjust=0.5,bty='n',
                           title=ifelse(Nyrs,paste(YrRange,paste0("Nyrs = ",Poolyrsy),sep="\n"),YrRange))
                } #End if ob.m is not nan
            } #End iyrs loop
            if (!draft) {mtext("Proportion",side=2,outer=TRUE,line=1,cex=cexLab)}
            mtext(title.x,side=1,outer=TRUE,line=0)
            if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = gfilename,
                                         graphics.type)

        } #End if plot pooled
    }     #End obsname loop
    return(invisible(0))
} # end function definition

