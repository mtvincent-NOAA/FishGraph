#' Abundance Index plots from stock assessment models
#'
#' Plots observed and predicted and diagnostic plots for index development
#'
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param connect.obsd When \code{TRUE} a line connecting observed points is plotted
#' @param from.zero When \code{TRUE}, the Y-axis of each plot (except recruitment deviations) starts at zero.
#' @param two.panel when \code{TRUE}, the observed-predicted plot and residual time plot are drawn as two-panels.  When \code{FALSE} they are drawn independently.
#' @param log.resid When \code{TRUE} residuals are computed as \eqn{R = log(U \diagup\hat{U})}.
#' @param err.bar When \code{TRUE}, error bars indicating plus/minus two standard errors of the observed index are plotted
#' @param resid.analysis When \code{TRUE}, additional diagnostic plots of residuals are created
#' @param resid.tests Runs statistical tests on residuals. Set to "runs" for runs test only or "multiple" for additional tests. Defaults to argument \code{NULL} for none.
#' @param runs.box Boolean for if you want to print the runs test 3 times SD on the log residual plot for the 2 box
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
#' Index.plots(gag)
#' }
#' @export
#######################################################################################
#  R function to plot abundance indices from fish stock assessment models.
#  Plots observed and predicted, and residual diagnostics.
#
#  Assumes a dataframe object called "t.series" as a component of "x"
#     and that CPUE data columns begin with "U"
#  Assumes that the first, third, etc., U* columns in x$t.series are observed values
#     and that the following U* columns are predicted values.
#
# M.H. Prager and E. H. Williams, December 2005
# Revised and expanded by A. Stevens, March 2006
# Major revision by M. H. Prager, November, 2006
# Fixed error: pairs plot now made only w/ > 1 series, MHP, 26 Aug 2009
# Revised for R package September 2015 by K. Shertzer
# Revised to add runs test and residual colors 2022 M. Vincent
#######################################################################################
Index.plots <- function(x, DataName = deparse(substitute(x)), draft = TRUE,
   graphics.type = NULL, use.color = TRUE, connect.obsd = FALSE, from.zero = TRUE,
   two.panel = TRUE, log.resid = TRUE, err.bar = TRUE, resid.analysis = TRUE, resid.tests=NULL, runs.box = TRUE)
#######################################################################################
#  ARGUMENTS:
#  x - an R list with output from the assessment models
#     The list x must have a component x$t.series that is a data frame
#     The data frame must have components
#        x$t.series$year with (integer) years of the analysis
#        x$t.series$U.*.ob with observed index values
#        x$t.series$U.*.pr with predicted index values
#        x$t.series$U.*.cv with predicted index values
#        IN ALL CASES, THE U.*.ob MUST PRECEDE THE CORRESPONDING U.*.pr.
#  DataName - a string representation an identifier for the data (run) in use.
#  draft - TRUE if the figures are to have main titles
#  graphics.type - a character vector with graphics-file types, e.g., c("png", "eps")
#  use.color - TRUE if graphs of index fits are in color
#  connect.obsd - TRUE to connect the observed points
#  from.zero - TRUE to start the Y-axis of each plot from zero
#  two.panel - TRUE to plot observed and predicted in upper panel, and residuals in lower panel
#  log.resid - TRUE to plot log(resid) rather than resid / mean(abs(resid))
#  err.bar - TRUE to plot error bars (plus/minus two SEs) on observed indices; requires
#             component x$t.series$U.*.cv in the data frame
#  resid.analysis - TRUE for additional diagnostic plots and tests on residuals, plots use color regardless
#                   of use.color setting
#  resid.tests - tests on residuals. "runs" for runs test only, "multiple" for multiple tests including runs test,
#                default is none.

#######################################################################################
{
    ### Check for x$t.series data frame:
    if (! ("t.series" %in% names(x)))
    {   ErrText = (paste("Component ", deparse(substitute(x)),
            "$t.series not found.", sep = ""))
        warning(ErrText, immediate. = TRUE)
        return(invisible(-1))
    }
    ### Local copy of data:
    ts <- x$t.series

    ### Check for other needed data:
    if (! ("year" %in% names(x$t.series)))
    {   warning("Year variable not found in Index.plots!", immediate. = TRUE)
        return(invisible(-1))
    }
    # Get indices of columns beginning with "U" (columns with CPUE):
    Ucols <-  grep("^U", names(ts))
    if (length(Ucols) == 0)
    {   warning("No data columns found that start with'U'.", immediate. = TRUE)
        return(invisible(-1))
    }
    # Use modulo operator to see if # of columns is odd:
    if (( length(Ucols) %% 2) != 0 )
    {   warning("Odd number of index columns found in Index.plots!", immediate. = TRUE)
        return(invisible(-1))
    }
    # Get indices of columns beginning with "cv.U" (columns with CPUE):
    if (err.bar){
      CVcols <-  grep("^cv.U", names(ts))
      if (length(CVcols) == 0)
      {   warning("Error bars desired (err.bar=TRUE), but no data columns found that start with 'cv'", immediate. = TRUE)
          return(invisible(-1))
      }
    }
    nplots <- length(Ucols) / 2                   # Number of distinct indices
    ### Set graphics parameters
    savepar <- FGSetPar(draft)
    plot.options = FGGetOptions()
    ifelse (use.color, err.bar.col <- plot.options$color$clr.obsd, err.bar.col <- plot.options$bw$clr.obsd)

    PlotTitle <- ""
    # Define matrix to split screen later:
    smatrix <- rbind(c(0.0, 1.0, 0.40, 1.0), c(0.0, 1.0, 0.0, 0.40))
    colvec <- FGGetPal(nplots, use.color)
    ### If writing graphics files, make sure there is a directory for them:
    if (! is.null(graphics.type))
    {   write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/index", sep="")
        FGCheckGraphDir(GraphicsDirName)
    } else {write.graphs <- FALSE }

    if (connect.obsd) FGtype <- "2line" else FGtype <- "linepoint"

    #   Find the indices of the years to be plotted -- same for all plots.
    #   Start: the first year with a non-missing value of any index.
    #   Stop: the last year in the t.series data frame.
    temp   <- is.na(ts[Ucols])          # Array of TRUE (missing) or FALSE (not missing)
    temp   <- apply(temp, 1, all)       # Vector with TRUE if all cols are TRUE (missing)
    first  <- which.min(temp)           # First year with any column not NA
    yrndx  <- first:nrow(ts)            # These years (rows) of the d.f. will be plotted
    nyear  <- length(yrndx)
    ts     <- ts[yrndx,]                # Drop years that won't be used

    resids <- matrix(NA, nplots, nyear) # empty matrix to save residuals
    bar.names <- rep("", nplots)        # empty vector of bar names
#----------------------------------------------------------------------------------
# Make the obsd-pred plots:
#----------------------------------------------------------------------------------
    if (two.panel) localpar <- par(no.readonly = TRUE)
    for (iplot in 1:nplots)
    {   # Get indices for plotted columns:
        i2 <- Ucols[2*iplot]
        i1 <- i2 - 1
        #Get CPUE and replace missing -99999 with NA
        cpue.obs <- replace(ts[, Ucols[2*iplot - 1]],ts[, Ucols[2*iplot - 1]]==-99999,NA)
        cpue.pre <- replace(ts[, Ucols[2*iplot]],ts[, Ucols[2*iplot]]==-99999,NA)

        # Save residuals for barplot
        if (log.resid) {resids[iplot,] <- log((cpue.obs + 1e-20) / (cpue.pre + 1e-20))
        } else {resids[iplot,] <- cpue.obs - cpue.pre}

        # Get name of index and make plot title:
        IndexName  <- FGTrimName(names(ts[i1]), removePrefix = 1, removeSuffix = 1)
        bar.names[iplot] <- IndexName

        if (draft) PlotTitle <- FGMakeTitle(paste("Index:", IndexName), DataName) else PlotTitle <- ""

        cpue.obs.maxrange <- cpue.obs
        if (err.bar){
          CVName <- paste("cv.U.",IndexName,sep="")
          CVcol <- which(names(ts) == CVName)
          if (length(CVcol)==0)
          {warning(paste("Error bars desired (err.bar=TRUE), but no CV found for index series ", IndexName, sep=""), immediate. = TRUE)
           return(invisible(-1))
          }
          cpue.sd <- ts[, CVcol] * cpue.obs
          cpue.obs.upper <- cpue.obs + 2*cpue.sd
          cpue.obs.lower <- cpue.obs - 2*cpue.sd
          cpue.obs.maxrange <- cpue.obs.upper
        }

        # Scale Y-axis:
        if (from.zero) {yrange = c(0, max(cpue.obs.maxrange, cpue.pre, na.rm = TRUE))
        } else {yrange <- range(cpue.obs.maxrange, cpue.pre, na.rm = TRUE) }

        # Plot of obs and predicted:
        lab.x = "Year"
        if (two.panel)
        {   split.screen(smatrix)
            screen(1)
            par(mar = localpar$mar + c(-3, 0, 0, 0), cex = 1, cex.axis = 0.85,
                cex.main = 1)
            lab.x = ""
        }
        FGTimePlot(x = ts$year, y = cpue.pre, y2 = cpue.obs, lab.x = lab.x,
            lab.y = "Relative abundance (CPUE)", href = NULL, hrefnames = NULL,
            use.color = use.color, FGtype = FGtype,
            ylim = yrange, main = PlotTitle)
        if (err.bar){
          for (ibars in 1:length(ts$year)){
            lines(c(ts$year[ibars],ts$year[ibars]),c(cpue.obs.lower[ibars], cpue.obs.upper[ibars]), col=err.bar.col)
          }
        }

        ### Write plot to file
        if (!two.panel && write.graphs)
            {   FGSavePlot(GraphicsDirName, DataName,
                GraphName = paste("U", IndexName, sep="."), graphics.type)
            }

        # Plot of residuals:
        if (two.panel)
        {   screen(2)
            PlotTitle = ""
            if (draft) par(mar = localpar$mar + c(0, 0, -1, 0)) else par(mar = localpar$mar)
            par(cex = 1, cex.main = 1, cex.axis = 0.85, lab = c(5, 4, 7))
        }

        # Scaled residuals if not in log space:
        if (log.resid)
        {   scaled.resids <- resids[iplot,]
            lab.y <- "Log residual"
        } else {
            meany <-mean(abs(resids[iplot,]), na.rm = TRUE)
            scaled.resids <- (resids[iplot,] + 1e-20) / (meany + 1e-20)
            lab.y <- "Scaled residual"
        }

        maxy <- max(abs(scaled.resids), na.rm = TRUE)

        if (runs.box){
            runstest=runs.sig3(scaled.resids)
            maxy <- max(c(abs(scaled.resids),runstest$sig3lim), na.rm = TRUE)
        }

        ylim = maxy * 1.05 * c(-1, 1)

        FGTimePlot(x = ts$year, y = scaled.resids, lab.x = "Year",
            lab.y = lab.y, href = NULL, hrefnames = NULL, use.color = use.color,
            FGtype = "stick", ylim = ylim,  main = PlotTitle)
        abline(h = 0)

        if (runs.box){
            lims=runstest$sig3lim
            cols =  ifelse(runstest$p.runs<0.01,rgb(1,0,0,0.5),ifelse(runstest$p.runs<0.05,rgb(1,0.647,0,0.5),rgb(0,1,0,0.5)))
            rect(min(ts$year-1),lims[1],max(ts$year+1),lims[2],col=cols,border=cols)
            points(ts$year[scaled.resids<lims[1]|scaled.resids>lims[2]],scaled.resids[scaled.resids<lims[1]|scaled.resids>lims[2]],pch=16,col='red',cex=1)
        }

        ### Write plot to file(s)
        if (two.panel) GraphName <- paste("U", IndexName, sep=".") else GraphName = paste("U", IndexName, "resid", sep=".")
        if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName, graphics.type)
        if (two.panel) close.screen(all.screens = TRUE)

    } # end for iplot
    if (two.panel) par(localpar)
#----------------------------------------------------------------------------------
# Make the barplot of residuals:
#----------------------------------------------------------------------------------
#     if (draw.barplot)
#     {   if (log.resid)
#         {   ylab = "Log residual"
#         }
#         else
#         {   # Standardize residuals:  subtract mean, divide by sd.
#             meanresids <- mean(resids, na.rm = TRUE)
#             sdresids <- sd(as.vector(resids), na.rm = TRUE)
#             resids <- (resids - meanresids)
#             resids <- resids/sdresids
#             ylab = "Standardized residual"
#         }
#         plotresids <- resids                         # Copy for plotting
#         plotresids[is.na(resids)] <- 0.0             # Change NA to zero in copy only
#
#         # Set up plot
#         par(mfrow = c(2,1))  # Divide plot region for plots of + and - residuals
#         if(draft) PlotTitle <- FGMakeTitle("Residuals in abundance indices.", DataName)
#         else PlotTitle <- ""
#
#         # Plot positive residuals:
#         xresids <- ifelse(plotresids > 0.0, plotresids, 0.0)  #Get positive resids only
#         # Get height of tallest bar, to place legend:
#         ymax <- max(apply(xresids, 2, sum, na.rm = TRUE))
#         par(mar=c(1, 4.5, 2, 0) + 0.2)
#         barplot(xresids, beside = FALSE, names.arg = ts$year, cex.names = 0.9,
#             axis.lty = 1, col = colvec, ylab = "", main = "", las = 1)
#         legend(x = 1, y = ymax, xjust = 0, yjust = 1, fill = colvec, legend = bar.names)
#         mtext(PlotTitle, side = 3, line = 1, adj = 0.5, font = 2)
#
#         ### Plot negative residuals:
#         xresids <- ifelse(plotresids < 0.0, plotresids, 0.0)
#         par(mar=c(1, 4.5, 1.3, 0) + 0.2)
#         mp <- barplot(xresids, beside = FALSE, names.arg = ts$year, axis.lty = 1,
#             col = colvec, ylab = "", main = "", xlab = "",las = 1, axisnames = FALSE)
#         axis(side = 3, tick = TRUE, at = mp, labels = rep("", nyear), tcl = -0.25)
#         par(mfrow = c(1,1))
#         mtext(ylab, side = 2, line = 3.5, adj = 0.5, font = 2)
#
#         if (write.graphs) { FGSavePlot(GraphicsDirName, DataName,
#             GraphName = "U.resids", graphics.type) }
#     }  # End if (barplot)

#----------------------------------------------------------------------------------------------------------------------------
#Log-residual analysis plots (EHW added on 9-10-14)
#---------------------------------------------------------------------------------------------------------------------------
if (resid.analysis) {
  #--start residual analysis plots------------------------------
  #---- Residual analysis page 1--------------------------------------------------

  for (iplot in 1:nplots)
  {
    # Get indices for plotted columns:
    i2 <- Ucols[2*iplot]
    i1 <- i2 - 1

    cpue.obs <- ts[, Ucols[2*iplot - 1]]
    cpue.pre <- ts[, Ucols[2*iplot]]

    IndexName  <- FGTrimName(names(ts[i1]), removePrefix = 1, removeSuffix = 1)

    ob=log(cpue.obs)
    pr=log(cpue.pre)
    res=pr-ob
    res.vals=(!is.na(res))
    res=res[res.vals]
    yr=ts$year[res.vals]
    n.res=length(res)

    if (draft) PlotTitle <- FGMakeTitle(paste("Residual Analysis (1 of 2), Index:", IndexName), DataName) else PlotTitle <- ""

    GraphName=paste("U.",IndexName,".resid1", sep="")

    layout(rbind(c(1,2),c(3,4)))
    par(mai=c(1,1,1,0.25))

   if(log.resid==TRUE) {
      resid.label="Log-Residuals (Ob-Pr)"
      } else {resid.label="Residuals (Ob-Pr)"}


    plot(yr,res,pch=1,col="darkblue",xlab="Year",ylab=resid.label,cex=1.5)
    abline(h=0,lwd=2,lty=2)
    grid()
    points(yr,res,pch=1,col="darkblue",cex=1.5)
    lines(lowess(res~yr),col="red")

    mtext(PlotTitle,side=3,line=-2,outer=TRUE,cex=1.5)

    plot(pr[res.vals],res,pch=1,col="darkblue",xlab="Predicted Value",ylab=resid.label,cex=1.5)
    abline(h=0,lwd=2,lty=2)
    grid()
    points(pr[res.vals],res,pch=1,col="darkblue",cex=1.5)
    lines(lowess(res~pr[res.vals]),col="red")

    qqnorm(res,pch=16)
    qqline(res,col=2)

    acf(res,main="Autocorrelation Function")
   ### Write plot to file(s)
   if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName, graphics.type)

   #---- Residual analysis page 1--------------------------------------------------

    if (draft) PlotTitle <- FGMakeTitle(paste("Residual Analysis (2 of 2), Index:", IndexName), DataName) else PlotTitle <- ""
    GraphName=paste("U.",IndexName,".resid2", sep="")
    ##Breusch-Pagan Test for heteroskedasticity
    bptest.out=bptest(yr~res)
    #Harrison-McCabe test for heteroskedasticity
    hmctest.out=hmctest(yr~res)

    ##Breusch-Godfrey test for higher-order serial correlation
    bgtest.out=bgtest(yr~res,order=1)
    ##Durbin-Watson test for autocorrelation of disturbances
    dwtest.out=dwtest(yr~res)

    #from package nortest
    #Lilliefors (Kolmogorov-Smirnov) test for normality
   if (n.res>4) {lillie.test.out=lillie.test(res)}
   if (n.res<=4)
   {print(paste("Index.plots: Lilliefors (Kolmogorov-Smirnov) test not performed on index ", IndexName,
                "; test requires N > 4.", sep=""))}

    #Anderson-Darling test for normality
    if (n.res>7) {ad.test.out<-ad.test(res)}
    if (n.res<=7)
      {print(paste("Index.plots: Anderson-Darling test not performed on index ", IndexName,
      "; test requires N > 7.", sep=""))}
    #Pearson chi-square test for normality
    pearson.test.out=pearson.test(res)

    #Shapiro-Wilk test for normality
    shapiro.test.out=shapiro.test(res)

    #Phillips-Perron test for the null hypothesis that x has a unit root
    pp.test.out=pp.test(res)

    #---second page of plots-----------------------------------------------
    layout(rbind(c(1,2),c(3,3)))
    par(oma=c(0,0,1,0))
    ##using package lawstat
    #runs test
    runs.test.v2.out=runs.test.v2(res,plot.it=T)

    #mtext(index.title,side=3,line=-1,outer=TRUE,cex=1.5)
    mtext(PlotTitle,side=3,line=-1,outer=TRUE,cex=1.5)

    max.res=max(abs(res))
    range=seq(from=-max.res,to=max.res,length.out=100)
    res.hist=hist(res,plot=F)
    multiplier=res.hist$counts/res.hist$density
    plot(res.hist,xlab="Residual (ob-pr)",main="Histogram of residuals w/ normal curve")
    curve(dnorm(x,mean=0,sd=sd(res))* multiplier[1],col="blue",lwd=2, add=TRUE, yaxt="n")

    if (!is.null(resid.tests)) {
    par(mai=c(0.25,0.25,0.25,0.25))
    plot(1,1,ylim=c(-15,0),xlim=c(0,3),xlab="",ylab="",bty="n",col.axis="transparent",xaxt="n",yaxt="n",col=0)

    text(x=1.75,y=0,labels="Traffic light",cex=1,font=2)
    text(x=2.25,y=0,labels="P-value",cex=1,font=2)

    if (resid.tests=="runs") {
      text(x=1.5,y=-1,labels="Runs test:",adj=c(1,0.5),cex=1,font=2)
      if(runs.test.v2.out$p.value<0.05){sign.col="red"
        sign.pos=1.85}
      if(runs.test.v2.out$p.value>=0.05&runs.test.v2.out$p.value<=0.1){sign.col="yellow"
       sign.pos=1.75}
      if(runs.test.v2.out$p.value>0.1){sign.col="green3"
        sign.pos=1.65}
      points(x=sign.pos,y=-1,pch=16,col=sign.col,cex=2.25)

      points(x=c(1.65,1.75,1.85),y=c(-1,-1,-1),cex=2.25,col=1)

      text(x=2.25,y=-1,labels=round(runs.test.v2.out$p.value,4),cex=1,font=1)

    } else if (resid.tests=="multiple"){
      text(x=1.5,y=-1,labels="Runs test:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-2,labels="Breusch-Pagan test for heteroskedasticity:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-3,labels="Harrison-McCabe test for heteroskedasticity:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-4,labels="Breusch-Godfrey test for higher-order serial correlation:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-5,labels="Durbin-Watson test for autocorrelation of disturbances:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-6,labels="Lilliefors (Kolmogorov-Smirnov) test for normality:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-7,labels="Anderson-Darling test for normality:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-8,labels="Pearson chi-square test for normality:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-9,labels="Shapiro-Wilk test for normality:",adj=c(1,0.5),cex=1,font=2)
      text(x=1.5,y=-10,labels="Phillips-Perron test for null hypothesis x has a unit root:",adj=c(1,0.5),cex=1,font=2)


    if(runs.test.v2.out$p.value<0.05){sign.col="red"
                                      sign.pos=1.85}
    if(runs.test.v2.out$p.value>=0.05&runs.test.v2.out$p.value<=0.1){sign.col="yellow"
                                                                sign.pos=1.75}
    if(runs.test.v2.out$p.value>0.1){sign.col="green3"
                                        sign.pos=1.65}
    points(x=sign.pos,y=-1,pch=16,col=sign.col,cex=2.25)

    if(bptest.out$p.value<0.05){sign.col="red"
                                sign.pos=1.85}
    if(bptest.out$p.value>=0.05&bptest.out$p.value<=0.1){sign.col="yellow"
                                                         sign.pos=1.75}
    if(bptest.out$p.value>0.1){sign.col="green3"
                                  sign.pos=1.65}
    points(x=sign.pos,y=-2,pch=16,col=sign.col,cex=2.25)

    if(hmctest.out$p.value<0.05){sign.col="red"
                                 sign.pos=1.85}
    if(hmctest.out$p.value>=0.05&hmctest.out$p.value<=0.1){sign.col="yellow"
                                                          sign.pos=1.75}
    if(hmctest.out$p.value>0.1){sign.col="green3"
                                   sign.pos=1.65}
    points(x=sign.pos,y=-3,pch=16,col=sign.col,cex=2.25)

    if(bgtest.out$p.value<0.05){sign.col="green3"
                                sign.pos=1.65}
    if(bgtest.out$p.value>=0.05&bgtest.out$p.value<=0.1){sign.col="yellow"
                                                         sign.pos=1.75}
    if(bgtest.out$p.value>0.1){sign.col="red"
                                  sign.pos=1.85}
    points(x=sign.pos,y=-4,pch=16,col=sign.col,cex=2.25)

    if(dwtest.out$p.value<0.05){sign.col="green3"
                                sign.pos=1.65}
    if(dwtest.out$p.value>=0.05&dwtest.out$p.value<=0.1){sign.col="yellow"
                                                         sign.pos=1.75}
    if(dwtest.out$p.value>0.1){sign.col="red"
                                  sign.pos=1.85}
    points(x=sign.pos,y=-5,pch=16,col=sign.col,cex=2.25)

   if (n.res > 4) {
      if(lillie.test.out$p.value<0.05){sign.col="red"
                                       sign.pos=1.85}
      if(lillie.test.out$p.value>=0.05&lillie.test.out$p.value<=0.1){sign.col="yellow"
                                                                sign.pos=1.75}
      if(lillie.test.out$p.value>0.1){sign.col="green3"
                                         sign.pos=1.65}
      points(x=sign.pos,y=-6,pch=16,col=sign.col,cex=2.25)
   }
    if (n.res > 7) {
      if(ad.test.out$p.value<0.05){sign.col="red"
                                   sign.pos=1.85}
      if(ad.test.out$p.value>=0.05&ad.test.out$p.value<=0.1){sign.col="yellow"
                                                            sign.pos=1.75}
      if(ad.test.out$p.value>0.1){sign.col="green3"
                                     sign.pos=1.65}
      points(x=sign.pos,y=-7,pch=16,col=sign.col,cex=2.25)
    }
    if(pearson.test.out$p.value<0.05){sign.col="red"
                                      sign.pos=1.85}
    if(pearson.test.out$p.value>=0.05&pearson.test.out$p.value<=0.1){sign.col="yellow"
                                                               sign.pos=1.75}
    if(pearson.test.out$p.value>0.1){sign.col="green3"
                                        sign.pos=1.65}
    points(x=sign.pos,y=-8,pch=16,col=sign.col,cex=2.25)

    if(shapiro.test.out$p.value<0.05){sign.col="red"
                                      sign.pos=1.85}
    if(shapiro.test.out$p.value>=0.05&shapiro.test.out$p.value<=0.1){sign.col="yellow"
                                                               sign.pos=1.75}
    if(shapiro.test.out$p.value>0.1){sign.col="green3"
                                        sign.pos=1.65}
    points(x=sign.pos,y=-9,pch=16,col=sign.col,cex=2.25)

    if(pp.test.out$p.value<0.05){sign.col="green3"
                                 sign.pos=1.65}
    if(pp.test.out$p.value>=0.05&pp.test.out$p.value<=0.1){sign.col="yellow"
                                                          sign.pos=1.75}
    if(pp.test.out$p.value>0.1){sign.col="red"
                                   sign.pos=1.85}
    points(x=sign.pos,y=-10,pch=16,col=sign.col,cex=2.25)

      points(x=c(1.65,1.75,1.85),y=c(-1,-1,-1),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-2,-2,-2),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-3,-3,-3),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-4,-4,-4),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-5,-5,-5),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-6,-6,-6),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-7,-7,-7),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-8,-8,-8),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-9,-9,-9),cex=2.25,col=1)
      points(x=c(1.65,1.75,1.85),y=c(-10,-10,-10),cex=2.25,col=1)


      text(x=2.25,y=-1,labels=round(runs.test.v2.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-2,labels=round(bptest.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-3,labels=round(hmctest.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-4,labels=round(bgtest.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-5,labels=round(dwtest.out$p.value,4),cex=1,font=1)
      if (n.res > 4) {text(x=2.25,y=-6,labels=round(lillie.test.out$p.value,4),cex=1,font=1)
      } else {text(x=2.25,y=-6,labels="NA",cex=1,font=1)}
      if (n.res > 7) {text(x=2.25,y=-7,labels=round(ad.test.out$p.value,4),cex=1,font=1)
      } else {text(x=2.25,y=-7,labels="NA",cex=1,font=1)}
      text(x=2.25,y=-8,labels=round(pearson.test.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-9,labels=round(shapiro.test.out$p.value,4),cex=1,font=1)
      text(x=2.25,y=-10,labels=round(pp.test.out$p.value,4),cex=1,font=1)
    } else {(print("Warning: Index.plots argument resid.tests not specified correctly. See help (?Index.plots)."))
    }# end resid.tests=runs or multiple

    } #end resids.test

    ### Write plot to file(s)
    if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName, graphics.type)
  }

} #end resid.analysis
#----------------------------------------------------------------------------------
# Make the pairs plot of residuals:
#----------------------------------------------------------------------------------
# if (draw.pairs && nplots > 1)
# {   ts <- ts[,Ucols]    # Only CPUE series
#     # Only odd (obs) CPUE series
#     ts <- ts[, seq(from = 1, to = ncol(ts) - 1, by = 2)]
#     tsnames <- names(ts)
#     # Take ".ob" from end of names:
#     for (i in 1:ncol(ts))
#     {    tsnames[i] <- FGTrimName(tsnames[i], removeSuffix = 1)
#     }
#     names(ts) <- tsnames
#     clr <- "black"
#     if (use.color) clr <- FGoptions$color$clr.line
#     PlotTitle = ""
#     if (draft) PlotTitle <- FGMakeTitle("Comparison of abundance indices", DataName)
#     pairs(ts, main = PlotTitle, lwd = 2, cex.labels = 1.6, cex = 1.5, col = clr,
#           cex.axis = 1.2)
#     if (write.graphs) { FGSavePlot(GraphicsDirName, DataName,
#                                    GraphName = "U.pairs", graphics.type) }
# }

#----------------------------------------------------------------------------------
    # Add names to the matrix of residuals before returning it to caller
    rownames(resids) <- bar.names
    colnames(resids) <- ts$year
     # Restore graphics parameters and return:
     par(savepar)
    return(invisible(t(resids)))
}   # END OF FUNCTION
#=================================================================================

