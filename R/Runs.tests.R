runs.test.v2=function (y, plot.it = FALSE, alternative = c("two.sided", "positive.correlated", "negative.correlated"))
{
    alternative <- match.arg(alternative)
    DNAME = deparse(substitute(y))
    y <- na.omit(y)
    med <- median(y, na.rm = TRUE)
    for (k in 2:length(y)) {
        if ((y[k] == med) & (y[k - 1] < med)) {
            y[k] = y[k - 1]
        }
        else if ((y[k] == med) & (y[k - 1] > med)) {
            y[k] = y[k - 1]
        }
    }
    q <- rep(0.05, length(y))
    p <- rep(-0.05, length(y))
    d <- y
    q[I(d < med) | I(d == med)] <- NA
    p[I(d >= med)] <- NA
    if (plot.it) {
        plot(q, type = "p", pch = "A", col = "red", ylim = c(-0.5, 0.5),
             xlim = c(1, length(y)), xlab = "", ylab = "",main="Runs Test")
        points(p, pch = "B", col = "blue")
        abline(h = 0)
    }
    m <- length(na.omit(q))
    n <- length(na.omit(p))
    R <- 1
    s <- sign(y - med)
    for (k in 1:(length(y) - 1)) {
        if (s[k] != s[k + 1]) {
            R <- R + 1
        }
    }
    E <- 1 + 2 * n * m/(n + m)
    s2 <- (2 * n * m * (2 * n * m - n - m))/((n + m)^2 * (n +
                                                          m - 1))
    statistic <- (R - E)/sqrt(s2)
    if (alternative == "positive.correlated") {
        p.value = pnorm(statistic)
        METHOD = "Runs Test - Positive Correlated"
    }
    else if (alternative == "negative.correlated") {
        p.value = 1 - pnorm(statistic)
        METHOD = "Runs Test - Negative Correlated"
    }
    else {
        p.value = 2 * min(pnorm(statistic), 1 - pnorm(statistic))
        alternative = "two.sided"
        METHOD = "Runs Test - Two sided"
    }
    STATISTIC = statistic
    names(STATISTIC) = "Standardized Runs Statistic"
    structure(list(statistic = STATISTIC, p.value = p.value,
                   method = METHOD, data.name = DNAME), class = "htest")
}


#' Function to do runs.test and calculate 3 sigma limits
#' runs test is conducted with library(randtests)
#' @param x residuals from CPUE fits
#' @param type only c("resid","observations")
#' @param mixing c("less","greater","two.sided"). Default less is checking for postive autocorrelation only
#' @return runs test p value and 3 x sigma limits
#'
#' @author Henning Winker (JRC-EC) and Laurence Kell (Sea++) Edited by Matt Vincent (NOAA SEFSC)
runs.sig3 <- function(x,type="resid",mixing="two.sided") {
  if(type=="resid"){
      mu = 0}else{mu = mean(x, na.rm = TRUE)}
  alternative= c("two.sided", "positive.correlated", "negative.correlated")[which(c("two.sided", "positive.correlated", "negative.correlated")%in%mixing)]
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    # Make the runs test non-parametric
    runstest = runs.test.v2(x,FALSE,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }

  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}


#' Simple plot function for runs test that plots the test for a single index or mean length/age time series. This function needs to be given the residuals precalculated for a single object with names Yr and residuals
#'
#' @param year a vector with the year of the residual
#' @param scaled.resids a vector of the
#' @param lab.y a string of the y axis label
#' @param PlotTitle a string of the plot title
#'
#' @return Graphics
#'
#' @author M T Vincent (many ideas from Henning Winker and Laurence Kell)
#'
#' @export

FGRunsPlot <- function(year,scaled.resids,lab.y,PlotTitle){

    maxy <- max(abs(scaled.resids), na.rm = TRUE)
    runstest=runs.sig3(scaled.resids)
    maxy <- max(c(abs(scaled.resids),runstest$sig3lim), na.rm = TRUE)

    ylim = maxy * 1.05 * c(-1, 1)

    FGTimePlot(x = year, y = scaled.resids, lab.x = "Year",
               lab.y = lab.y, href = NULL, hrefnames = NULL, use.color = FALSE,
               FGtype = "stick", ylim = ylim,  main = PlotTitle)
    abline(h = 0)

    lims=runstest$sig3lim
    cols =  ifelse(runstest$p.runs<0.05,rgb(1,0,0,0.5),ifelse(runstest$p.runs<0.1,rgb(1,1,0,0.5),rgb(0,1,0,0.5)))
    rect(min(year-1),lims[1],max(year+1),lims[2],col=cols,border=cols)
    points(year[scaled.resids<lims[1]|scaled.resids>lims[2]],scaled.resids[scaled.resids<lims[1]|scaled.resids>lims[2]],pch=16,col='red',cex=1)
}

#' Function that takes the output from an assessment and produces plots of the runs test for the age and length composition.
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param log.resid
#'
#'
#'@export
ALRuns.plots <- function(x, DataName = deparse(substitute(x)), draft = TRUE, graphics.type = NULL, log.resid = TRUE){

    Errstring = ("No composition data found.  Terminating ALRuns.plots");
    if (! ("comp.mats" %in% names(x))) stop(Errstring)

    savepar <- FGSetPar(draft)

    ### Make local copy of needed data components
    cm <- x$comp.mats
    ts=x$t.series
    ## Is number of columns odd?  This is a problem -- they should be in pairs!
    if ( (length(cm) %% 2) != 0 ) stop("Odd number of matrices found in ALRuns.plots!\n")

    ### Set graphics parameters
    plot.options = FGGetOptions()
    savepar <- FGSetPar(draft)

    ### If writing graphics files, make sure there is a directory for them:
    if (! is.null(graphics.type)){  write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/comp", sep="")
        FGCheckGraphDir(GraphicsDirName)
    }else{  write.graphs <- FALSE
    }

    lcomp=cm[grep("lcomp",names(cm))]
    lcomp=lcomp[sapply(lcomp,dim)[1,]!=1]
    ## Is number of columns odd?  This is a problem -- they should be in pairs!
    if ( (length(lcomp) %% 2) != 0 ) stop("Odd number of length compositions found in ALRuns.plots!\n")
    nlplot=length(lcomp)/2
    if (nlplot>0){               #Check to make sure there are length comps present.
        ## If length comps present figure out how many there are
        if(nlplot<4){par(mfcol=c(nlplot,1))
        } else if (nlplot>3&nlplot<7){par(mfcol=c(ceiling(nlplot/2),2))
        } else if (nlplot>6){par(mfcol=c(ceiling(nlplot/3),3))}
        par(cex = 1, cex.main = 1, cex.lab = 0.85, cex.axis = 0.75, mar=c(3.5,4,0.85,0.7), oma=c(0,0,2,0), las=1, lab = c(5, 2, 4),  mgp=c(1.8, 0.75, 0))

        for (iplot in 1:nlplot){
            lenObs=lcomp[[2*iplot-1]]
            lenPred=lcomp[[2*iplot]]
            ## Check to see if it is a matrix. If it's just a vector then it is a pooled comp and not worth plotting
            years=as.numeric(row.names(lenObs))
            lenbins=as.numeric(colnames(lenObs))
            meanlenObs=rep(NA,length(years))
            meanlenPred=rep(NA,length(years))
            for(i in 1:length(years)){
                ## Calculate the mean observed and predicted length
                meanlenObs[i]=sum(lenObs[i,]*lenbins)
                meanlenPred[i]=sum(lenPred[i,]*lenbins)
            }
            if(log.resid){
                resid= log(meanlenObs)-log(meanlenPred)
            } else {
                resid= meanlenObs-meanlenPred
                resid = (resid+1e-20)/(mean(abs(resid),na.rm=TRUE)+1e-20)
            }
            lab.y=ifelse(log.resid,"Log residual","Scaled residual")
            FGRunsPlot(years,resid,lab.y,FGTrimName(names(lcomp[2*iplot-1]),1,1))
        }

        if (draft) title("Runs Tests for Length Compostions",line=.5,outer=TRUE)
        if (write.graphs) FGSavePlot(GraphicsDirName, DataName, "length.runs.test", graphics.type)
    }

    ## Make plot for age composition means
    acomp=cm[grep("acomp",names(cm))]
    acomp=acomp[sapply(acomp,dim)[1,]!=1]
    ## Is number of columns odd?  This is a problem -- they should be in pairs!
    if ( (length(acomp) %% 2) != 0 ) stop("Odd number of age compositions found in ALRuns.plots!\n")
    naplot=length(acomp)/2
    if (naplot>0){
         ## If age compositions are present then set up the layout accordingly
        if(naplot<4){par(mfcol=c(naplot,1))
        } else if (naplot>3&naplot<7){par(mfcol=c(ceiling(naplot/2),2))
        } else if (naplot>6){par(mfcol=c(ceiling(naplot/3),3))}
        par(cex = 1, cex.main = 1, cex.lab = 0.85, cex.axis = 0.75, mar=c(3.5,4,.85,.7), oma=c(0,0,2,0), las=1, lab = c(5, 2, 4),  mgp=c(1.8, 0.75, 0))
        for (iplot in 1:naplot){
            ageObs=acomp[[2*iplot-1]]
            agePred=acomp[[2*iplot]]

            years=as.numeric(row.names(ageObs))
            ages=as.numeric(colnames(ageObs))
            meanageObs=rep(NA,length(years))
            meanagePred=rep(NA,length(years))
            for(i in 1:length(years)){
                ## Calculate the mean observed and predicted length
                meanageObs[i]=sum(ageObs[i,]*ages)
                meanagePred[i]=sum(agePred[i,]*ages)
            }
            if(log.resid){
                resid= log(meanageObs)-log(meanagePred)
            } else {resid= meanageObs-meanagePred}
            lab.y=ifelse(log.resid,"Log residual","Scaled residual")
            FGRunsPlot(years,resid,lab.y,FGTrimName(names(acomp[2*iplot-1]),1,1))
        }

        if (draft) title("Runs Tests for Age Compostions",line=.5,outer=TRUE)
        if (write.graphs) FGSavePlot(GraphicsDirName, DataName, "age.runs.test", graphics.type)
    }
    par(savepar)
    return(invisible(NULL))
}

