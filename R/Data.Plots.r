#' Plot of available data from the output of an assessment
#' The function \code{Data.Plots} generates a timeseries plot of the different data sources. Optionally, the size of the points for the age, length, and weight composition can be proportion to the sample size of input data
#'
#' @param x an R list with output from the assessment model.
#' @param DataName string used in plot titles. Defaults to argument \code{x}.
#' #' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param scaleSize a boolean whether to scale size of point for age and length composition based on the sample size. Default is FALSE
#' @param cexval cex multiplier of squares to make continuous bar if desired
#' @param fill a boolean whether to fill the point, this is useful for diagnosing how much overlap occurs in the points, Default if TRUE
#' @param sampleName a character string of the length and age composition name to search for in the tseries matrix. Do not include the '.'. Default is 'neff'
#'
#' @return Graphics
#'
#' @author MT Vincent
#'
#' @examples \donttest{
#' Data.Plots(gag)
#' }
#' @export
#'
###########################################################################################
## R function to make plot of data input into assessment
## Part of FishGraph


Data.Plots  <-  function(x, DataName = deparse(substitute(x)), draft = TRUE,
                       graphics.type = NULL, use.color = TRUE, scaleSize=FALSE,cexval=1,fill=TRUE,sampleName="neff")
{
###########################################################################################


### Check for x$t.series data frame:
    if (! ("t.series" %in% names(x)))
    {   ErrText = (paste("Component ", deparse(substitute(x)),
                         "$t.series not found.", sep = ""))
                         warning(ErrText, immediate. = TRUE)
                         return(invisible(-1))
    }
    ts = x$t.series

    ## Replace time series values with NA for less than 0 values
    ts[ts<0]=NA

    ### Check for other needed data:
    if (! ("year" %in% names(x$t.series)))
    {   warning("Year variable not found in Index.plots!", immediate. = TRUE)
        return(invisible(-1))
    }

    ## Get the years of available catch by fishery
    years = ts$year

    ## midyr = mean(years)
    minyr=min(years)-(min(years)%%5)
    maxyr=max(years)+ifelse((max(years)%%5)==0,0,5-(max(years)%%5))
    midyr=.5
    xyears=(years-minyr)/(maxyr-minyr)


    ## Get years of available CPUE by fishery
    CPUEcols=grep("^U.*ob",names(ts))
    CPUEs=ts[,CPUEcols]
    tmp=names(ts)[CPUEcols]
    ## remove the first 2 and last 3 characters from string to get name of fleet
    CPUENames=substr(substr(tmp,3,nchar(tmp)),1,nchar(tmp)-5)

    ## Landings
    Landingcols=grep("^L.*ob",names(ts))
    Landing=ts[,Landingcols]
    tmp=names(ts)[Landingcols]
    LandingNames=substr(substr(tmp,3,nchar(tmp)),1,nchar(tmp)-5)

    ## Discards
    Discardcols=grep("^D.*ob",names(ts))
    Discard=ts[,Discardcols]
    tmp=names(ts)[Discardcols]
    DiscardNames=substr(substr(tmp,3,nchar(tmp)),1,nchar(tmp)-5)

    ## Length compositions
    Lengthcols=grep(paste0("^lcomp.*",sampleName),names(ts))
    LengthN=ts[,Lengthcols]
    tmp=names(ts)[Lengthcols]
    LengthNames=substr(substr(tmp,7,nchar(tmp)),1,nchar(tmp)-7-nchar(sampleName))

    ## Age compositions
    Agecols=grep(paste0("^acomp.*",sampleName),names(ts))
    AgeN=ts[,Agecols]
    tmp=names(ts)[Agecols]
    AgeNames=substr(substr(tmp,7,nchar(tmp)),1,nchar(tmp)-7-nchar(sampleName))

    ## Close kin mark recapture
    CKMRcol=grep('s1',names(x$ckmr_data))
    if (length(CKMRcol) > 0){
        CKMRYears=(unique(x$ckmr_data[,CKMRcol])-minyr)/(maxyr-minyr)
        CKMRN=aggregate(x=x$ckmr_data$ckmr.npairs,by=list(x$ckmr_data$s1),FUN=sum)
        names(CKMRN) <- c("Year","Comparisons")
    }



    ## Check which data sources are used in this assessment
    PlotCPUE=ifelse(length(CPUENames)==0,FALSE,TRUE)
    PlotLanding=ifelse(length(LandingNames)==0,FALSE,TRUE)
    PlotDiscard=ifelse(length(DiscardNames)==0,FALSE,TRUE)
    PlotLength=ifelse(length(LengthNames)==0,FALSE,TRUE)
    PlotAge=ifelse(length(AgeNames)==0,FALSE,TRUE)
    PlotCKMR=ifelse(length(CKMRcol)==0,FALSE,TRUE)

    ## Calculate what the height of the y axis must be to plot all of these
    ## Gives 2 spaces for title and 1 for rows of data
    pt=1
    ttl=2
    ycount =  length(CPUENames)*pt + length(LandingNames)*pt + length(DiscardNames)*pt + length(LengthNames)*pt + length(AgeNames)*pt + PlotCPUE*ttl + PlotLanding*ttl + PlotDiscard*ttl + PlotLength*ttl + PlotAge*ttl + PlotCKMR*(ttl + pt)
    spacer=1/(ycount-2)

    AllNames=c(LandingNames,DiscardNames,CPUENames,LengthNames,AgeNames)

    if (use.color){
        colvec <- rainbow(length(unique(AllNames)), s = 0.9, v = 0.65, start = 0.1, end = 1)
        if (length(unique(AllNames)) == 1) colvec = "gray55"
    } else {
        colvec  <- gray.colors(length(unique(AllNames)), start = 0.2, end = 0.8)
    }
    if(PlotCKMR) {
        AllNames=c(AllNames,"CKMR")
        colvec=c(colvec,"grey50")
    }

    fishcols=data.frame(names=unique(AllNames),col=colvec)

    ## Make the plot with the right dimensions for the plot

    savepar = FGSetPar(draft)

    PlotTitle=ifelse(draft, paste(DataName,"Data Availability"),"")


    if (draft){ par(mar=c(2.5,1+max(nchar(AllNames)),3,3+max(nchar(AllNames))))
    } else{ par(mar=c(2.5,1+max(nchar(AllNames)),.5,3+max(nchar(AllNames))))}


    plot(NA,xlim=c(0,1),ylim=c(-1,0),main=PlotTitle,axes=FALSE,frame.plot=TRUE,ylab='',xlab='Year',cex=1.5)

    axis(1,at=seq(0,1,length.out=length(seq(minyr,maxyr,5))),labels=seq(minyr,maxyr,5),cex.axis=1.1)
    counter=spacer
    savecount=NULL
    ## Put in Landings title and data available points
    if (PlotLanding){
        counter=counter-spacer
        text(midyr,counter,"Landings",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        if (length(LandingNames)==1){
            points(xyears[!is.na(Landing)],rep(counter,length(xyears[!is.na(Landing)])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==LandingNames,2],"transparent"),cex=cexval,col=fishcols[fishcols$names==LandingNames,2])
            savecount=c(savecount,counter)
            counter=counter-spacer
        } else {
            for (i in 1:length(LandingNames)) {
                points(xyears[!is.na(Landing[,i])],rep(counter,length(xyears[!is.na(Landing[,i])])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==LandingNames[i],2],"transparent"),cex=cexval,col=fishcols[fishcols$names==LandingNames[i],2])
                savecount=c(savecount,counter)
                counter=counter-spacer
            }                           #For Landing length
        }                           #if length ==1
    }                 #End if plot Landing

    ## Put in Discard title and data available point
    if (PlotDiscard){
        counter=counter-spacer
        text(midyr,counter,"Discards",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        if (length(DiscardNames)==1){
            points(xyears[!is.na(Discard)],rep(counter,length(xyears[!is.na(Discard)])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==DiscardNames,2],"transparent"),cex=cexval,col=fishcols[fishcols$names==DiscardNames,2])
            savecount=c(savecount,counter)
            counter=counter-spacer
        } else {
            for (i in 1:length(DiscardNames)) {
                points(xyears[!is.na(Discard[,i])],rep(counter,length(xyears[!is.na(Discard[,i])])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==DiscardNames[i],2],"transparent"),cex=cexval,col=fishcols[fishcols$names==DiscardNames[i],2])
                savecount=c(savecount,counter)
                counter=counter-spacer
            }                           #For Discard length
        }                           #if length ==1
    }                 #End if plot Discard

    ## Put in CPUE title and available data
    if (PlotCPUE){
        counter=counter-spacer
        text(midyr,counter,"Abundance Indices",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        if (length(CPUENames)==1){
            points(xyears[!is.na(CPUEs)],rep(counter,length(xyears[!is.na(CPUEs)])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==CPUENames,2],"transparent"),cex=cexval,col=fishcols[fishcols$names==CPUENames,2])
            savecount=c(savecount,counter)
            counter=counter-spacer
        } else {
            for (i in 1:length(CPUENames)) {
                points(xyears[!is.na(CPUEs[,i])],rep(counter,length(xyears[!is.na(CPUEs[,i])])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==CPUENames[i],2],"transparent"),cex=cexval,col=fishcols[fishcols$names==CPUENames[i],2])
                savecount=c(savecount,counter)
                counter=counter-spacer
            }                               #For CPUE length
        }                               #if length==1
    }                 #End if plot CPUE

    ## Put in Length comp data availability with option for size of points
    if (PlotLength){
        counter=counter-spacer
        text(midyr,counter,"Length Composition",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        if (length(LengthNames)==1){
            ifelse(scaleSize,cexmult <- (LengthN[!is.na(LengthN)&LengthN>0]/max(LengthN,na.rm=TRUE))*cexval ,cexmult <- cexval)
            points(xyears[!is.na(LengthN)&LengthN>0],rep(counter,length(xyears[!is.na(LengthN)&LengthN>0])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==LengthNames,2],"transparent"),cex=cexmult,col=fishcols[fishcols$names==LengthNames,2])
            savecount=c(savecount,counter)
            counter=counter-spacer
        } else {
            for (i in 1:length(LengthNames)) {
                ifelse(scaleSize,cexmult <- (LengthN[(!is.na(LengthN[,i])&LengthN[,i]>0),i]/max(LengthN,na.rm=TRUE))*cexval ,cexmult <- cexval)
                points(xyears[(!is.na(LengthN[,i])&LengthN[,i]>0)],rep(counter,length(xyears[(!is.na(LengthN[,i])&LengthN[,i]>0)])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==LengthNames[i],2],"transparent"),cex=cexmult,col=fishcols[fishcols$names==LengthNames[i],2])
                savecount=c(savecount,counter)
                counter=counter-spacer
            }                               #For Length length
        }                               #if length==1
    }                 #End if plot Length

    ## Put in Age comp data availability with option for size of points
    if (PlotAge){
        counter=counter-spacer
        text(midyr,counter,"Age Composition",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        if (length(AgeNames)==1){
            ifelse(scaleSize,cexmult <- (AgeN[!is.na(AgeN)]/max(AgeN,na.rm=TRUE))*cexval ,cexmult <- cexval)
            points(xyears[!is.na(AgeN)],rep(counter,length(xyears[!is.na(AgeN)])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==AgeNames,2],"transparent"),cex=cexmult,col=fishcols[fishcols$names==AgeNames,2])
            savecount=c(savecount,counter)
            counter=counter-spacer
        } else {
            for (i in 1:length(AgeNames)) {
                ifelse(scaleSize,cexmult <- (AgeN[!is.na(AgeN[,i]),i]/max(AgeN,na.rm=TRUE))*cexval ,cexmult <- cexval)
                points(xyears[!is.na(AgeN[,i])],rep(counter,length(xyears[!is.na(AgeN[,i])])),pch=22,lwd=3,bg=ifelse(fill,fishcols[fishcols$names==AgeNames[i],2],"transparent"),cex=cexmult,col=fishcols[fishcols$names==AgeNames[i],2])
                savecount=c(savecount,counter)
                counter=counter-spacer
            }                               #For Age length
        }                               #if length==1
    }                 #End if plot Age

    ## Put in CKMR sampling year if available
    if (PlotCKMR){
        counter=counter-spacer
        text(midyr,counter,"Close Kin Mark Recapture Samples",cex=1.25,adj=c(0.5,0))
        counter=counter-spacer
        ifelse(scaleSize,cexmult <- (CKMRN$Comparisons/max(CKMRN$Comparisons,na.rm=TRUE))*cexval ,cexmult <- cexval)
        points(CKMRYears,rep(counter,length(CKMRYears)),pch=22,lwd=3,bg=ifelse(fill,'grey50',"transparent"),cex=cexmult,col='grey50')
        savecount=c(savecount,counter)
        counter=counter-spacer
    }

    axis(2,at=savecount,labels=AllNames,las=2,cex.axis=1.1)
    legend(1.05, -0.5, bty='n', xpd=NA, legend=fishcols[,1], pch=22, pt.bg=fishcols[,2],col=fishcols[,2],cex=1.1,yjust=0.5,lwd=3,lty=NA)

    ### If writing graphics files, make sure there is a directory for them:
    if (! is.null(graphics.type)){
        write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs", sep="")
        FGCheckGraphDir(GraphicsDirName)
    } else {  write.graphs <- FALSE }


    if (write.graphs) FGSavePlot(GraphicsDirName, DataName, GraphName = "data_availability", graphics.type)
         # Restore graphics parameters and return:
     par(savepar)
    return(invisible(0))

}
                                        #End Function

