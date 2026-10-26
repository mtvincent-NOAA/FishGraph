#' Plots of age- and length-composition fits
#'
#' The function \code{Comp.plots} generates bubble plots of residuals of age-
#' and length-composition fits for the entire time frame of the assessment.
#' Bubble areas are scaled to the largest residual within each plot. Optionally,
#' a small inset plot displays the correlation or angular
#' deviation between observed and predicted values each year. In addition,
#' \code{Comp.plots} generates plots of predicted and observed mean compositions
#' pooled across years and weighted by sample size.
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
#' @param resid.type type of residuals to plot for in each of the diagnostics, options are "OSA","standard","pearson", and "all". Input can be a single value or vector of values
#' @param corr when \code{FALSE}, angular deviation is displayed in the inset plot; otherwise,
#' correlations.
#' @param c.min lower bound on the y-axis range of the inset plot, applies only when corr is \code{TRUE}
#' @param max.bub cex value for maximum bubble size, default is 1.5.
#' @param b.plot create boxplots of composition residusls, default is FALSE.
#' @param resid.plots create qqplot and scatterplot of residuals, default is FALSE.
#' @param Nyrs boolean for whether the number of pooled years is to be added to legend. Default is \code{FALSE}
#'
#'
#' @return Graphics
#'
#' @author M Prager
#' @author E Williams
#' @author K Shertzer
#' @author R Cheshire
#' @author K Purcell
#' @author M Vincent
#'
#' @examples \donttest{
#' Comp.plots(gag)
#' }
#' @export
###########################################################################################
                                        #  R function to make bubble plots of age- or length-composition matrices
                                        #  from fish stock-assessment models.
                                        #  Part of FishGraph.
                                        #  REQUIRES source("fgsupport.r") before running this.
                                        #  E. H. Williams, December 2002
                                        #  Revised by M. H. Prager, December, 2005
                                        #  Revised by A. Stephens, May, 2006
                                        #  Major revision by M. H. Prager, November 2006
                                        #  Revised by R. Cheshire, November 2012
                                        #  Revised by K. Shertzer, September 2015
                                        #  Revised by R. Cheshire, January 2019 (add boxplots, weighting of pooled plots if 1) neff in t.series
                                        #                                         or 2. 'n' in t.series, otherwise no weights,
                                        #                                         pooled comps split by selectivity blocks)
#######################################################################################
Comp.plots <- function(x, DataName = deparse(substitute(x)), draft = TRUE,
                       graphics.type = NULL, use.color = TRUE, units = x$info$units.length,
                       resid.type = "OSA", corr = TRUE, c.min = 0.25, max.bub = 1.5, b.plot = FALSE, resid.plots=FALSE, Nyrs=FALSE)
#######################################################################################
                                        #  ARGUMENTS:
                                        #  x - an R list with output from the assessment models
                                        #     The list x must have a component x$comp.mats that is a list of matrices
                                        #     These matrices must be found in pairs.  First, xxxx.ob, then xxxx.pr.
                                        #  DataName - a string representation an identifier for the data (run) in use.
                                        #  draft - TRUE if the figure has a main title.
                                        #  graphics.type - a character vector with graphics-file types
                                        #  use.color - TRUE of graphs are in color
                                        #  units - character string with length units for Y-axis of length-comp plots
                                        #  p.resid - when FALSE, bubbles areas scaled to the largest residuals, otherwise scaled to the
                                        #          Pearson residuals.
                                        #  corr - TRUE plots correlation coefficients in bottom pane,  FALSE plots
                                        #        angular deviation
                                        #  c.min - lower y axis value for correlation plots, default is 0.25.  All correlations
                                        #          below this value are plotted at the minimum as a different symbol and color.
                                        #  max.bub - cex value for maximum bubble size, default is 8.0
#######################################################################################
{  Errstring = ("No composition data found.  Terminating Comp.plots");
    if (! ("comp.mats" %in% names(x))) stop(Errstring)

### Make local copy of needed data components
    cm <- x$comp.mats
    ts=x$t.series
    sel=x$sel.age
    ## Replace all with all 3 types of

    residsType=unlist(ifelse(tolower(resid.type)=="all",list(c("osa","standard","pearson")),tolower(resid.type)))
    if(!all(residsType%in%c("osa","standard","pearson")))stop('resid.type needs to be one of the following "osa","standard","pearson", or "all".')
    ## Is effective sample size present if resid.type=OSA
    if(any(residsType%in%c("pearson","osa"))){
        if(length(grep('neff',names(ts)))==0) stop("resid.type of OSA or Pearson requires effective sample size 'neff' components in x$t.series.")
    }

    if(any(residsType%in%"osa") &(!require(compResidual,quietly=TRUE))) stop("You do not have the required package installed to plot One Step Ahead residuals find out how at https://github.com/fishfollower/compResidual")

    ## Is number of columns odd?  This is a problem -- they should be in pairs!
    if ( (length(cm) %% 2) != 0 ) stop("Odd number of matrices found in Comp.plots!\n")

### Set graphics parameters
    plot.options = FGGetOptions()
    savepar <- FGSetPar(draft)

### If writing graphics files, make sure there is a directory for them:
    if (! is.null(graphics.type)){  write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/comp", sep="")
        FGCheckGraphDir(GraphicsDirName)
    }else{  write.graphs <- FALSE
    }

### Set the dimensions for splitting the graphics screen into two plotting regions:
    smatrix <- rbind(c(0.0, 1.0, 0.25, 1.0), c(0.0, 1.0, 0.0, 0.25))
    nplots <- length(cm) %/% 2               # integer division operator
    ## Added this as a short cut to control colors without having to reference outside controls.

    if (use.color){
        clr.pos <- plot.options$color$clr.pos    # color for positive residuals
        clr.neg <- plot.options$color$clr.neg    # color for negative residuals
        clr.ang <- plot.options$color$clr.ang    # color for angle plot
        col.obsd <- plot.options$color$clr.obsd
        col.pred <- plot.options$color$clr.pred
        Y1Col  <- plot.options$color$clr.line
    }  else{
        clr.pos <- plot.options$bw$clr.pos    # color for positive residuals
        clr.neg <- plot.options$bw$clr.neg    # color for negative residuals
        clr.ang <- plot.options$bw$clr.ang    # color for angle plot
        col.obsd <- plot.options$bw$clr.obsd
        col.pred <- plot.options$bw$clr.pred
        Y1Col  <- plot.options$bw$clr.line
    }

    ## Set up the matrices to store the residual average length and age over time between observed and predicted
    nageComp <- length(grep('acomp[.].*[.]ob',names(cm)))
    namesAgeComp <- unname(sapply(names(cm)[grep('acomp[.].*[.]ob',names(cm))], function (x) FishGraph:::FGTrimName(x, removePrefix = 1, removeSuffix = 1)))
    AgeResids <- matrix(NA, length(ts$year),nageComp,dimnames= list(ts$year, namesAgeComp))

    nlenComp <- length(grep('lcomp[.].*[.]ob',names(cm)))
    namesLenComp <- unname(sapply(names(cm)[grep('lcomp[.].*[.]ob',names(cm))], function (x) FishGraph:::FGTrimName(x, removePrefix = 1, removeSuffix = 1)))
    LenResids <- matrix(NA, length(ts$year),nlenComp,dimnames= list(ts$year, namesLenComp))

    for (iplot in 1:nplots) {

        ##---Cumulative fit plots-------------------------------------------------------------------
        par(oma=c(2,2.5,0,0))
        if (draft){ par(mar = c(3, 5, 3, 1 ))
        }else par(mar = c(3, 3, 1, 1))
        ## check to make sure that ob comes first
        if(length(grep("ob",names(cm)[iplot * 2-1]))==1){
            m1.all <- cm[[iplot*2-1]]         # Matrix of observed
            m2.all <- cm[[iplot*2]]           # matrix of predicted
        } else{
            m1.all <- cm[[iplot*2]]         # Matrix of observed
            m2.all <- cm[[iplot*2-1]]       # matrix of predicted
        }
        ## Extract column names for plots
        binnames=colnames(m1.all)

        ## Get various string representations of data series:
        ##  gfileroot is used in the name for the graphics file(s):
        gfileroot <- FGTrimName(names(cm)[iplot * 2], removePrefix = 0, removeSuffix = 1)
        if (gfileroot != FGTrimName(names(cm)[iplot * 2-1], removePrefix = 0, removeSuffix = 1)) {stop(paste("There is not an observed and predicted matrix for",gfileroot),sep=" ")}
        groot=FGTrimName(names(cm)[iplot * 2], removePrefix = 1, removeSuffix = 1)
        gfilename <- paste("pooled.",gfileroot, sep="")
        ## test if a matrix or a vector
        if(!any(paste(c("sel.m","sel.v"),groot,sep=".") %in% names(sel))) {
            warning(paste0("Neither sel.m.",groot," nor sel.v.",groot," were found in the selectivity information in the rdat provided. Skipping this plot for now.") )
            next
        }

        if (paste("sel.m",groot,sep=".") %in% names(sel)){
            sel.block <- sel[[paste("sel.m",groot,sep=".")]]
            sel.block <- sel.block[rownames(sel.block)%in%rownames(m1.all),,drop=FALSE]
            block.text <- rownames(unique(as.data.frame(sel.block)))
        } else if (paste("sel.v",groot,sep=".") %in% names(sel)){
            block.text = rownames(m1.all)[1]
        }
        yr.block <- as.numeric(block.text)
        ##add last year to vector of #if(length(yr.block)>1){
        yr.block=c(yr.block,last(ts$year)+1) #max(as.numeric(rownames(as.data.frame(m1.all))))+1)
        ## titleroot is used as part of the plot title:
        titleroot <- paste("Fishery: ", gfileroot)

        neffname<-paste(gfileroot,".neff", sep="")
        nname<-paste(gfileroot,".n", sep="")
        if (neffname%in%names(x$t.series)) {
            nseries<-x$t.series[,names(x$t.series)==neffname]
        } else if(nname%in%names(x$t.series)) {
            nseries<-x$t.series[,names(x$t.series)==nname]
        } else {
            warning("There is no effective N or sample size in t.series for the ",neffname," or ",nname,". Going to the next plot.")
            next
        }
        yrs.include<-x$t.series$year[nseries>0 & !is.na(nseries)]
        ## Figure out how many plots will actually be plotted
        subplots=length(which(sapply(1:length(block.text),function(i) length(which(rownames(m1.all[rownames(m1.all)%in%yr.block[i]:(yr.block[i+1]-1),,drop=FALSE])%in%yrs.include)))>0))

        ## Set the format based on the numnber of plots
        if(subplots<4){par(mfcol=c(subplots,1));cexLab=1#/.83
        } else if (subplots==4){par(mfrow=c(2,2));cexLab=1/.83
        } else if (subplots<=8){par(mfrow=c(ceiling(subplots/2),2));cexLab=1/.66
        } else if (subplots>9){par(mfrow=c(ceiling(subplots/3),3));cexLab=1/.66}
                                        #if (subplots==1)cexLab=1
        if (subplots==3)cexLab=1/.66
        ##subset to years in each selectivity block (last block interval is defined differntly because last value of yr.block is
        ## end of interval instead of beginning of next interval)
        for(iyrs in 1:(length(yr.block)-1)){

            m1.sub<-m1.all[rownames(m1.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]
            m2.sub<-m2.all[rownames(m2.all)%in%yr.block[iyrs]:(yr.block[iyrs+1]-1),,drop=FALSE]

            ##Exclude those yrs that don't make the min sample size requirement, designated by n<0
            if (dim(m1.sub)[1]>1) {
                m1<-m1.sub[rownames(m1.sub)%in%yrs.include,,drop=FALSE]
                m2<-m2.sub[rownames(m2.sub)%in%yrs.include,,drop=FALSE]
                poolylabtxt="wgt by effective N"
                wt=x$t.series[[neffname]][!is.na(x$t.series[[neffname]])&x$t.series[[neffname]]>0&x$t.series$year%in%as.numeric(rownames(m1))]
                ## wt=wt[is.na(wt)==FALSE] # I think this is redundant
            } else {
                m1<-m1.sub
                m2<-m2.sub
                poolylabtxt="unweighted"
            }

            ## Set Y-axis title according to data type:
            if(substr(gfileroot, 1, 1)  == "l"){title.x <- FGMakeLabel("Length bin", units)
            }else{title.x <- "Age class"}

            ## allow plotting of vector instead of data frame object
            if (length(m1)==0) {
                warning(paste0("You do not have any samples to plot for ", gfileroot," in year block ",yr.block[iyrs],". Skipping this for now but you should probably check this in BAM!"))
                next
            }

            if(dim(m1)[1]>1){
                ob.m=apply(m1,2,weighted.mean,w=wt)
                pr.m=apply(m2,2,weighted.mean,w=wt)
                x.vec=as.numeric(dimnames(m1)[[2]])
                Poolyrs=dim(m1)[1]
            } else {
                ob.m=m1
                pr.m=m2
                x.vec=as.numeric(dimnames(m1)[[2]])
                Poolyrs=1
            }

            ymax=max(c(ob.m,pr.m))

            if(!is.nan(ob.m[1])){
                                        #plot(x=x.vec,y=ob.m,type="b",col=1,lwd=2,cex=2,ylim=c(0,ymax),xlab=title.y,ylab="Mean Composition",main=titleroot)
                ## par(las=FGSetLas(c(ob.m, pr.m)))
                plot(x=x.vec,y=ob.m,col=NULL,ylim=c(0,ymax),xlab="",ylab="",las=1,cex.axis=cexLab)
                if (draft) title(main=titleroot)
                ##polygon(x=c(min(x.vec),x.vec,max(x.vec)),y=c(0,ob.m,0),border="lightgray",col="lightgray")
                grid(col = "lightgray", lty = 1)
                points(x=x.vec,y=ob.m,col=col.obsd,pch=21, cex = 1.25, lwd=2)
                lines(x=x.vec,y=pr.m,col=Y1Col,lwd=2,type="l"); points(x=x.vec,y=pr.m,col=Y1Col,lwd=2, pch=16)
                YrRange=paste(yr.block[iyrs],yr.block[iyrs+1]-1,sep="-")
                legend(x="topright",bg="white",legend=ifelse(Nyrs,paste(YrRange,paste0("Nyrs = ",Poolyrs),sep="\n"),YrRange),cex=cexLab,xjust=0.5,bty='n')
                ## if (Nyrs) legend(x="topleft",bg="white",legend=c(paste0("Nyrs = ",Poolyrs)),cex=cexLab,xjust=1,bty='n')

                if (draft){mtext(paste("Pooled",poolylabtxt,sep="-"),side=2,line=4)}
            }}
        if (!draft) {mtext("Proportion",side=2,outer=TRUE,line=1,cex=cexLab)}
        mtext(title.x,side=1,outer=TRUE,line=0)
        if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = gfilename,
                                     graphics.type)

        for (RT in residsType){
            ##--------------bubble plots--------------------------------------------------------------
            ## Set Y-axis title according to data type:
            ylabresidtxt= switch(RT,"osa"="OSA Residuals","pearson"="Pearson Residuals","standard"="Deviance Residuals")

            ## Extract the used years
            ## if (nname%in%names(ts)) {
            m1<-m1.all[rownames(m1.all)%in%yrs.include,,drop=FALSE]
            m2<-m2.all[rownames(m2.all)%in%yrs.include,,drop=FALSE]
            ## } else {
            ##     m1<-m1.all
            ##     m2<-m2.all
            ## }

            ## Set Y-axis title according to data type:
            if(substr(gfileroot, 1, 1)  == "l"){ title.y <- FGMakeLabel(paste(ylabresidtxt,"Length bin",sep=" - "), units)
            } else {title.y <-  paste(ylabresidtxt,"Age class",sep=" - ")}


            ## Get size and color of the bubbles:
            if(RT=="pearson"){
                wt=ts[names(ts)==paste(gfileroot,"neff",sep=".")] #get sample sized associated with each year
                wt=wt[wt>0&is.na(wt)==FALSE]
                z1=(m1-m2)/sqrt(m2*(1-m2)/wt)
            } else if (RT=="standard"){z1 <- m1 - m2
            } else if (RT=="osa"){

                ## Need to get the dirichlet multinomial variance scalar to calculate the the OSA residuals
                dmroot=paste0("log_dm_",groot,"_",substr(gfileroot,1,2))
                if (length(grep(dmroot,names(x$parms),ignore.case=TRUE))>0){ dmvar=x$parms[[grep(dmroot,names(x$parms),ignore.case=TRUE)]]
                } else { dmvar=x$parm.cons[8,grep(dmroot,names(x$parm.cons),ignore.case=TRUE)]
                }
                if (length(dmvar)==0) { warning("Neither parms nor parm.cons contains ",dmroot," so you cannot calculate the OSA residuals. Skipping this for now.")
                    break
                }
                if(nname%in%names(x$t.series)) {
                    ns =x$t.series[,names(x$t.series)==nname]
                    ns=ns[!is.na(ns) & ns>0]
                } else { warning(paste0("OSA needs to have the sample size to be calculated for ",gfileroot));break}

                if(!is.null(dim(m1))){
                    if(length(ns)!=dim(m1)[1]) {
                        warning(paste("The vector of sample sizes does not match the dimensions of the matrix of compositions for",gfileroot,". Skipping this for now.",sep=" "))
                        break
                    }
                    obs=round(sweep(m1,1,ns,"*"))
                    alpha=sweep(m2,1,ns*exp(dmvar),"*")+0.00001
                    ## Calculate residuals using compResidual package
                    z1=t(compResidual::resDirM(t(obs),t(alpha)))
                } else {
                    if(length(ns)!=length(m1)) {
                        warning(paste("The vector of sample sizes does not match the vector length of compositions for",gfileroot,". Skipping this for now.",sep=" "))
                        break
                    }
                    obs=round(m1*ns)
                    alpha=m2*ns*exp(dmvar)+0.00001
                    z1=t(compResidual::resDirM((obs),(alpha)))
                }
                if (any(is.infinite(z1))){
                    z1[is.infinite(z1)]=NA
                    warning(paste0("There were infinite values in the OSA for ",gfileroot,". This was replaces with an NA to alow for plotting. This plot should probably not be used."))
                }
            }

            ## Get coordinates of bubbles:
            if(is.matrix(z1)){
                if ("cres"%in%class(z1)){ #class of OSA residuals
                    irn <- as.integer(rownames(m1))                        # year names
                    x1 <- rep(irn, ncol(m1)-1)                              # year names
                    y1 <- rep(as.numeric(head(binnames,-1)), each=nrow(m1))    # age- or length-class names

                } else { ##("matrix"%in%class(z1)) {
                    irn <- as.integer(rownames(z1))                        # year names
                    x1 <- as.integer(rep(irn, ncol(z1)))                   # year names
                    y1 <- rep(as.numeric(colnames(z1)), each=nrow(z1))    # age- or length-class names
                }
            } else {
                irn <- as.integer(yrs.include[!is.na(yrs.include)])                        # year names
                x1 <- as.integer(rep(irn, length(z1)))                   # year names
                y1 <- sort(as.numeric(names(z1)))    # age- or length-class names
            }

            z2 <- max.bub*(sqrt(abs(z1))/sqrt(max(abs(z1), na.rm = TRUE)))  ###plots area of bubble
            colvec <- ifelse(z1 < 0.0, clr.neg, clr.pos)

            ## Compute angular deviation for each year using the formula:
            ## angle - arccos (dotprod(a,b)/(sqrt(dotprod(a,a) * dotprot(b,b)))
            if (is.matrix(m1)){
                m12dp <- apply(m1*m2, 1, sum)   # This returns a vector of length nyear
                m11dp <- apply(m1*m1, 1, sum)   # Ditto
                m22dp <- apply(m2*m2, 1, sum)   # Ditto
            } else {
                m12dp <- m1%*%m2                      #Dot product works if vector
                m11dp <- m1%*%m1
                m22dp <- m2%*%m2
            }
            angdev = acos(m12dp/sqrt(m11dp*m22dp)) * 180.0 / pi   # Converts to degrees
            if(is.matrix(m1)){
                fit.metric=rep(0,dim(m1)[1])
                for (i in 1:length(fit.metric)){
                    fit.metric[i]=cor(m1[i,],m2[i,],method='pearson')
                }
            } else {
                fit.metric=cor(m1,m2,method='pearson')
            }


            ## Note: Set x limits so axis range is >= 4
            ## This prevents plotting fractional years (e.g. 1995.5)
            xmin = irn[1]-1
            xmax = max(irn,na.rm=TRUE)+1
            layout(matrix(c(1,2,3),3,1),heights=c(1,4,1))
            par(cex = 1, cex.main = 1, cex.axis = 0.85)
        {  if (draft) par(mar = c(0, 4, 3, 1 ))
           else par(mar = c(0, 4, 1, 1))
        }
            bub.scale=c(max.bub,max.bub*0.75,max.bub*0.5,max.bub*0.25,max.bub*0.1)
            legx=c(1,2,3,4,5)
            legy=rep(1,length(bub.scale))
            find.bub.label=function(x){x*max(abs(z1),na.rm=TRUE)/max.bub}
            bub.label=round(find.bub.label(bub.scale),2)

            plot(legx,legy,pch=19,cex=bub.scale,ylim=c(0.8,1.2),xlim=c(0,7),xaxt="n",
                 yaxt="n",xlab="",ylab="",bty="n")
            text(1:5,rep(0.85,6),bub.label)
            if (draft) {
                if (use.color){
                    title(main = FGMakeTitle(paste(titleroot, "    Orange: underestimate"),
                                             DataName))
                } else {
                    title(main = FGMakeTitle(paste(titleroot, "    Light: underestimate"),
                                             DataName))
                }
            }   #end if draft

            par(mar=c(0,4,0,1))
            ## Draw the main (bubble) plot:
            plot(x1, y1, xlab = "", ylab = title.y, type = "n", las = 1, xaxt = "n",
                 xlim = c(xmin, xmax))
            grid(col = "lightgray", lty = 1)
            points(x1, y1, cex = z2, col = colvec, pch = 19)

            par(mar=c(3,4, 0 , 1 ))
            if(corr==FALSE){
                plot(irn, angdev, xlab = "Year", xaxt = "n", ylab = "Ang.dev.",
                     ylim = c(0, 90), axes = FALSE, frame.plot = TRUE, type = "n",
                     xlim = c(xmin,xmax))
                axis(side = 2, at = seq(from = 0, to = 90, by = 10), las = 1,
                     tck = 1, col = "gray75", lty = "dotted", labels = FALSE)
                axis(side = 2, at = c(20), tck = 1, col = "gray55", labels = FALSE)
                axis(side = 2, at = c(0, 30, 60, 90), las = 1, labels = TRUE, cex.axis = 0.75 )
                box()
                points(irn, angdev, pch = 18, col = clr.ang, cex = 1.3)
            }else
            {
                plot(yrs.include[is.na(yrs.include)==FALSE],fit.metric,ylim=c(c.min,1.1),type="n",pch=16,ylab="Corr.",las=1, cex.lab=0.9,xlim = c(xmin,xmax))
                grid(col = "lightgray", lty = 1)
                fit.col=rep(0,length(fit.metric))
                fit.pch=rep(0,length(fit.metric))
                for (i in 1:length(fit.metric)){

                    if(fit.metric[i]<c.min) {fit.pch[i]=13; fit.col[i]="red"; fit.metric[i]=c.min}
                    else{fit.pch[i]=18; fit.col[i]="black"; fit.metric[i]=fit.metric[i]}}
                points(yrs.include[is.na(yrs.include)==FALSE],fit.metric,pch=fit.pch,col=fit.col)
            }
            if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = paste(RT,gfileroot,sep='.'),
                                         graphics.type)
                                        #--------------boxplots----------------------------------------------------------
            par(savepar)
            par(mar=c(4,4,2,2))
            par(las=0)
            par(oma=c(2,2,0,0))
            if(b.plot){

                gfilename <- paste0(RT,".boxplot.",gfileroot)
### Set Y-axis title according to data type:
                if(substr(gfileroot, 1, 1)  == "l") {title.x <- FGMakeLabel("Length bin", units)
                } else {title.x <- "Age class"}

                par(cex = 1, cex.main = 1, cex.axis = 0.85)
                if (draft) { par(mar = c(4, 5, 3, 1 ))
                } else par(mar = c(4, 4, 1, 1))

                par(mfrow=c(2,1))
                boxplot(z1,xlab=title.x,ylab="",names=colnames(z1),show.names=TRUE)
                abline(h=0,lty=3)
                if (draft)
                { title(main = FGMakeTitle(titleroot,DataName))
                }   #end if draft
                boxplot(t(z1),xlab='Year',main="",ylab="",names=rownames(z1),show.names=TRUE)
                abline(h=0,lty=3)
                mtext(ylabresidtxt,side=2,outer=TRUE,line=0)
                if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = gfilename,
                                             graphics.type)
            }                                   # if b.plot=TRUE
###################################Residual plots of qq and scatter plot###################

            if (resid.plots){
                par(savepar)
                par(mfrow=c(2,1),mar=c(4,4,3,1),oma=c(0,0,0,0))
                col <- rep(1:nrow(z1), ncol(z1))
                ## Make QQ plot of residuals
                qqnorm(z1,col=col,main="")
                if(draft) title(main=FGMakeTitle(paste("QQ for ",titleroot),DataName))
                abline(0,1)
                legend("topleft",col=col,legend=binnames,pch=1,bty='n',cex=0.8,ncol=2)
                ## Make residual plots
                plot(as.vector(z1), col=col, ylab="Residuals", xlab="")
                if(draft) title(main=FGMakeTitle(titleroot,DataName))
                if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = paste(RT,"QQScatter",gfileroot,sep="."),
                                             graphics.type)

            }                               #If adding residual plots
        }                                   #Loop over residual types

#################################Joint residual plots of mean age  and length###########

        ## Calculate the mean length/age by year for the observed and predicted
        BinVals= as.numeric(dimnames(m1.all)[[2]])
        meanVals1= c(m1.all%*%BinVals)
        meanVals2= c(m2.all%*%BinVals)

        ## Calculate the residual mean length/age by year and store in the correct matrix
        if(substr(gfileroot, 1, 1)  == "l"){ #Length
            LenResids[ts$year%in%as.numeric(dimnames(m1.all)[[1]]),dimnames(LenResids)[[2]]==groot]=log(meanVals1)-log(meanVals2)
        }else{## Age
            AgeResids[ts$year%in%as.numeric(dimnames(m1.all)[[1]]),dimnames(AgeResids)[[2]]==groot]=log(meanVals1)-log(meanVals2)
        }

    }                                   #Loop over data matrices

    par(savepar)
    par(mar=c(4.5,5,.5,.5))
    ## ------------------------------------------------------------------------------
    ## Make plot of joint residuals mean length by year: MTV Add 9/2026
    ## ------------------------------------------------------------------------------
    palette(FGRichColors(ncol(LenResids)))

    FGTimePlot(x = ts$year, y = LenResids[,1], lab.y = 'Mean Length Residuals', Y1Col = 1, FGtype = 'stick', ylim = c(min(LenResids,na.rm=TRUE) , max(LenResids,na.rm=TRUE)) )
    abline(h=0, lty=2)
    loessDF= data.frame(y=c(LenResids),x=rep(ts$year,ncol(LenResids)))
    bploc=ts$year[apply(LenResids,1,function(x) !all(is.na(x)))]
    boxplot(y~x,data=loessDF,xlab='',ylab='',add=TRUE,at=bploc,xaxt='n',yaxt='n')
    for (i in 1:ncol(LenResids)) {
        lines(ts$year, LenResids[,i], lwd = 2, type = 'h', pch = 16, col = i, cex = 1.0)
        ## points(ts$year, resids[i,], lwd = 2, type = 'p', pch = 16, col = i, cex = 1.0)
        points(ts$year, LenResids[,i], type = 'p', pch = 21, bg = i, cex = 1.0)
    }

    loessModel <- loess(y~x,data=loessDF)
    lines(ts$year, predict(loessModel, newdata = data.frame(x = ts$year)), col = 'black', lwd = 2)

    RMSE=sqrt( sum(LenResids^2,na.rm=TRUE)/ length(which(!is.na(LenResids))))
    text(max(ts$year),max(LenResids,na.rm=TRUE),labels=paste0("RMSE = ",round(RMSE,2)),adj= c(1,.5))
    legend('bottomleft',, c(dimnames(LenResids)[[2]],'Loess'), col = c(1:ncol(LenResids),'black'), lty=c(rep(NA,ncol(LenResids)),1), pch = c(rep(16,ncol(LenResids)),NA),lwd=2, bty='n')
    GraphName = "MeanLen.JointResid"
    if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName, graphics.type)

    ## ------------------------------------------------------------------------------
    ## Make plot of joint residuals mean age by year: MTV Add 9/2026
    ## ------------------------------------------------------------------------------
    palette(FGRichColors(ncol(AgeResids)))

    FGTimePlot(x = ts$year, y = AgeResids[,1], lab.y = 'Mean Age Residuals', Y1Col = 1, FGtype = 'stick', ylim = c(min(AgeResids,na.rm=TRUE) , max(AgeResids,na.rm=TRUE)) )
    abline(h=0, lty=2)
    loessDF= data.frame(y=c(AgeResids),x=rep(ts$year,ncol(AgeResids)))
    bploc=ts$year[apply(AgeResids,1,function(x) !all(is.na(x)))]
    boxplot(y~x,data=loessDF,xlab='',ylab='',add=TRUE,at=bploc,xaxt='n',yaxt='n')
    for (i in 1:ncol(AgeResids)) {
        lines(ts$year, AgeResids[,i], lwd = 2, type = 'h', pch = 16, col = i, cex = 1.0)
        ## points(ts$year, resids[i,], lwd = 2, type = 'p', pch = 16, col = i, cex = 1.0)
        points(ts$year, AgeResids[,i], type = 'p', pch = 21, bg = i, cex = 1.0)
    }

    loessModel <- loess(y~x,data=loessDF)
    lines(ts$year, predict(loessModel, newdata = data.frame(x = ts$year)), col = 'black', lwd = 2)

    RMSE=sqrt( sum(AgeResids^2,na.rm=TRUE)/ length(which(!is.na(AgeResids))))
    text(max(ts$year),max(AgeResids,na.rm=TRUE),labels=paste0("RMSE = ",round(RMSE,2)),adj= c(1,.5))
    legend('bottomleft',, c(dimnames(AgeResids)[[2]],'Loess'), col = c(1:ncol(AgeResids),'black'), lty=c(rep(NA,ncol(AgeResids)),1), pch = c(rep(16,ncol(AgeResids)),NA),lwd=2, bty='n')
    GraphName = "MeanAge.JointResid"
    if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName, graphics.type)

    par(savepar)
    palette('R4')
    return(invisible(NULL))
}
                                        # End function definition
###########################################################################################
