#' Time plots of model catch, landings, and discards by fishery
#'
#' The function \code{CLD.total.plots} generates barplots of estimated catch, landings,
#' and discards (in numbers and weight) over time. Each bar is subdivided by fishery.
#'
#' @param x an R list with output from the assessment models.
#' @param DataName string used in plot titles.  Defaults to argument \code{x}.
#' @param draft modifies plots for use in a report.  When \code{FALSE} main titles
#' are omitted.
#' @param graphics.type a vector of graphics file types to which graphics are saved.
#' When \code{NULL}, no plots are saved.
#' @param use.color plots are made in grayscale when \code{FALSE}
#' @param first.year An integer (e.g. 1960) indicating the first year for plots.
#' The subsetting is done via row names of the matrices being plotted.
#' @param units.CLD.n A text string (e.g. \code{"million fish"}) used in labeling
#' Y-axis of plots of numbers caught, landed, or discarded.
#' @param units.CLD.w A text string (e.g. \code{"t"}) used in labeling the Y-axis
#' of plots of biomass caught, landed, or discarded.
#' @param CLD.n.references A list of three character-string names that specify reference
#' points in numbers, in order of catch, landings, and discards, to be included on numbers plots.
#' Input NULL for none, e.g., CLD.w.references=list(NULL,"msy.num","Dmsy.num").
#' @param CLD.w.references A list of three character-string names that specify reference
#' points in biomass, in order of catch, landings, and discards, to be included on biomass plots.
#' Input NULL for none, e.g., CLD.w.references=list(NULL,"msy", NULL).
#' @param plot.proportion When \code{TRUE}, additional plots are made, scaled as
#' proportions.
#'
#' @return None
#'
#' @author M Prager
#' @author E Williams
#' @author K Shertzer
#' @author R Cheshire
#' @author K Purcell
#'
#'
#' @examples \donttest{
#' CLD.total.plots(gag)
#' }
#'
#' @export
CLD.total.plots <-
function(x, DataName=deparse(substitute(x)),
    draft = TRUE, graphics.type = NULL, use.color = TRUE,
    first.year = 0, units.CLD.n = x$info$units.numbers,
    units.CLD.w = x$info$units.biomass,
    CLD.n.references = NULL, CLD.w.references = NULL,
    plot.proportion = TRUE)

{   ### Check for needed data component:
    if (! ("CLD.est.mats" %in% names(x)))
    {   Errmsg <- paste("List 'CLD.est.mats' not found in data object:",
            deparse(substitute(x)))
        warning(Errmsg, immediate. = TRUE)
        return(invisible(-1))
    }

    h.ref=ifelse (is.list(CLD.n.references), TRUE, FALSE)
    if (h.ref){
      ref.parms=CLD.n.references %in% names(x$parms)
      for (i in 1:3) {
        if ((!is.null(CLD.n.references[[i]])) && (!ref.parms[i])) {
          Errstring = (paste("CLD.n.reference element ", CLD.n.references[[i]], " not found.", sep = ""))
          warning(Errstring, immediate. = TRUE)
          return(invisible(-1))
        }
      }
    }
    h.ref=ifelse (is.list(CLD.w.references), TRUE, FALSE)
    if (h.ref){
      ref.parms=CLD.w.references %in% names(x$parms)
      for (i in 1:3) {
        if ((!is.null(CLD.w.references[[i]])) && (!ref.parms[i])) {
          Errstring = (paste("CLD.w.reference element ", CLD.w.references[[i]], " not found.", sep = ""))
          warning(Errstring, immediate. = TRUE)
          return(invisible(-1))
        }
      }
    }

    ### Set up plotting-related stuff:
    savepar <- FGSetPar(draft)
    PlotTitle <- ""
    if (! is.null(graphics.type))
    {   write.graphs <- TRUE
        GraphicsDirName <- paste(DataName, "-figs/CLD", sep="")
    } else {write.graphs <- FALSE}
#---------------------------------------------------------------------------------------
    mats <- x$CLD.est.mats
    # Set up information list for labeling plots
    specs = list(
        data.types = c("C", "L", "D"),
        data.units = c("n", "w"),
        txt.types = c("Catch", "Landings", "Dead discards"),
        txt.units = c("numbers", "weight"),
        #uom = c(x$info$units.numbers, x$info$units.biomass)
        uom = c(units.CLD.n, units.CLD.w)
    )
    n.t = length(specs$data.types)
    n.u = length(specs$data.units)
    # Go through the six combinations of types and units:
    for (tt in 1:n.t)
    {   for (uu in 1:n.u)
        {   # Construct prefix, e.g., "Lw" or "Cn":
            prefix <- paste(specs$data.types[tt], specs$data.units[uu], sep="")
            # Construct grep-search string, e.g., "^Lw" or "^Cn":
            search <- paste("^", prefix, sep = "")
            # Construct description, e.g. "Catch in numbers":
            desc <- paste(specs$txt.types[tt], "in",  specs$txt.units[uu])
            # Construct name of possible totals matrix, e.g., "Ln.total":
            tot.name = paste(prefix, "total", sep=".")
            # Get indices of columns beginning with the prefix:
            indexlist <-  grep(search, names(mats))
            # Remove the index that is the totals matrix:
            indexlist = indexlist[names(mats)[indexlist] != tot.name]
            # indexlist now has the index of matrices of one of {Cn, Cw, Ln, or ...}
            n.mats <- length(indexlist)
            if (n.mats > 0)
            {   # Find the first and last years of all matrices, store in yr.range:
                yr.range <- NULL
                for (i in indexlist)
                {   yr.range <- range(yr.range,
                        range(as.integer(rownames(mats[[i]]))))
                }
                n.yr = diff(yr.range) + 1
                # Create matrix to hold annual totals:
                totals  <- matrix(0, n.mats, n.yr)
                colnames(totals) <- seq(from = yr.range[1], to = yr.range[2])
                # Create vector to hold names for legend:
                clr.names <- rep("", n.mats)
                # Accumulate totals and store names for legend:
                n.clr <- 0
                for (i in 1:length(indexlist))
                {   ix = indexlist[i]
                    thismat <- mats[[ix]]
                    firstyr <- min(as.integer(rownames(thismat)))
                    offset <- firstyr - yr.range[1]
                    matlen = nrow(thismat)
                    totals[i, (1:matlen) + offset] <- rowSums(thismat)
                    n.clr <- n.clr + 1
                    # Store the fishery name into clr.names:
                    clr.names[n.clr] <- FGTrimName(names(mats)[ix], removePrefix = 1)
                } # End for i
                rownames(totals) <- clr.names
                # Remove rows before argument first.year
                totals <- totals[, as.integer(colnames(totals)) >= first.year, drop = FALSE]
                if(draft) PlotTitle <- FGMakeTitle(paste(desc, "by fishery"), DataName)

                hrefval=NULL
                hrefstring=NULL
                if (is.list(CLD.n.references) && (prefix == "Cn")){
                    hrefstring <- CLD.n.references[[1]]
                    hrefindex <- which(names(x$parms) == CLD.n.references[[1]])
                    hrefval <- unlist(x$parms[hrefindex])
                }
                if (is.list(CLD.n.references) && (prefix == "Ln")){
                  hrefstring <- CLD.n.references[[2]]
                  hrefindex <- which(names(x$parms) == CLD.n.references[[2]])
                  hrefval <- unlist(x$parms[hrefindex])
                }
                if (is.list(CLD.n.references) && (prefix == "Dn")){
                  hrefstring <- CLD.n.references[[3]]
                  hrefindex <- which(names(x$parms) == CLD.n.references[[3]])
                  hrefval <- unlist(x$parms[hrefindex])
                }
                if (is.list(CLD.w.references) && (prefix == "Cw")){
                  hrefstring <- CLD.w.references[[1]]
                  hrefindex <- which(names(x$parms) == CLD.w.references[[1]])
                  hrefval <- unlist(x$parms[hrefindex])
                }
                if (is.list(CLD.w.references) && (prefix == "Lw")){
                  hrefstring <- CLD.w.references[[2]]
                  hrefindex <- which(names(x$parms) == CLD.w.references[[2]])
                  hrefval <- unlist(x$parms[hrefindex])
                }
                if (is.list(CLD.w.references) && (prefix == "Dw")){
                  hrefstring <- CLD.w.references[[3]]
                  hrefindex <- which(names(x$parms) == CLD.w.references[[3]])
                  hrefval <- unlist(x$parms[hrefindex])
                }

                FGBarplot(t(totals), lab.y = FGMakeLabel(desc, specs$uom[uu]),
                    use.color = use.color, PlotTitle = PlotTitle,
                    leg.title = "Fishery", href=hrefval, hrefname=hrefstring)

                if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                    GraphName = paste("CLD", prefix, sep = "."), graphics.type)
                if (plot.proportion)
                {
                FGBarplot(t(totals), lab.y = "Proportion",
                    use.color = use.color, PlotTitle = PlotTitle,
                    leg.title = "Fishery", proportion = TRUE)
                if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                    GraphName = paste("CLD", prefix, "prop", sep = "."), graphics.type)
                } # end if plot proportion

            }  # End if n.mats > 0 ...
        }  # End for uu
    } # End for tt

    ## Make an optional plot that plots Landings and discards if they are in the same units
    for (uu in 1:n.u)
    {   # Construct prefix, e.g., "Lw" or "Cn":
        prefix <- paste0(c("L","D"), specs$data.units[uu])
                                        # Construct grep-search string, e.g., "^Lw" or "^Cn":
        search <- paste("^", prefix, sep = "")
                                        # Construct description, e.g. "Catch in numbers":
        desc <- paste("Landings and Discards in",  specs$txt.units[uu])
                                        # Construct name of possible totals matrix, e.g., "Ln.total":
        ## tot.name = paste(prefix, "total", sep=".")
                                        # Get indices of columns beginning with the prefix:
        indexlist <-  c(grep(search[1], names(mats)),grep(search[2], names(mats)))
                                        # Remove the index that is the totals matrix:
        indexlist = indexlist[-grep("total",names(mats)[indexlist])]
                                        # indexlist now has the index of matrices of one of {Cn, Cw, Ln, or ...}
        n.mats <- length(indexlist)
        if (n.mats > 0)
        {   # Find the first and last years of all matrices, store in yr.range:
            yr.range <- NULL
            for (i in indexlist){
               yr.range <- range(yr.range,
                                  range(as.integer(rownames(mats[[i]]))))
            }
            n.yr = diff(yr.range) + 1
                                        # Create matrix to hold annual totals:
            totals  <- matrix(0, n.mats, n.yr)
            colnames(totals) <- seq(from = yr.range[1], to = yr.range[2])
            ## Create vector to hold names for legend:
            clr.names <- rep("", n.mats)
            ## Accumulate totals and store names for legend:
            n.clr <- 0
            for (i in 1:length(indexlist)){
               ix = indexlist[i]
                thismat <- mats[[ix]]
                firstyr <- min(as.integer(rownames(thismat)))
                offset <- firstyr - yr.range[1]
                matlen = nrow(thismat)
                totals[i, (1:matlen) + offset] <- rowSums(thismat)
                n.clr <- n.clr + 1
               ## Store the fishery name into clr.names:
                clr.names[n.clr] <- paste(substr(names(mats)[ix],1,1),FGTrimName(names(mats)[ix], removePrefix = 1),sep='.')
            } # End for i
            rownames(totals) <- clr.names
            ## Remove rows before argument first.year
            totals <- totals[, as.integer(colnames(totals)) >= first.year, drop = FALSE]
            if(draft) PlotTitle <- FGMakeTitle(paste(desc, "by fishery"), DataName)

            hrefval=NULL
            hrefstring=NULL
            if (uu == 1){
                if (!is.null(CLD.n.references[[1]])){
                    hrefstring <- CLD.n.references[[1]]
                    hrefindex <- which(names(x$parms) == CLD.n.references[[1]])
                    hrefval <- unlist(x$parms[hrefindex])
                } else if (!is.null(CLD.n.references[[2]]) & !is.null(CLD.n.references[[3]])){
                    Lrefstring <- CLD.n.references[[2]]
                    Lrefindex <- which(names(x$parms) == CLD.n.references[[2]])
                    Lrefval <- unlist(x$parms[Lrefindex])
                    Drefstring <- CLD.n.references[[3]]
                    Drefindex <- which(names(x$parms) == CLD.n.references[[3]])
                    Drefval <- unlist(x$parms[Drefindex])
                    hrefval <- Lrefval + Drefval
                    hrefstring <- paste("L",Drefstring,sep=" + ")
                } else if (!is.null(CLD.n.references[[2]])){
                    hrefstring <- CLD.n.references[[2]]
                    hrefindex <- which(names(x$parms) == CLD.n.references[[2]])
                    hrefval <- unlist(x$parms[hrefindex])
                }
            }
            if (uu == 2){
                if (!is.null(CLD.w.references[[1]])){
                    hrefstring <- CLD.w.references[[1]]
                    hrefindex <- which(names(x$parms) == CLD.w.references[[1]])
                    hrefval <- unlist(x$parms[hrefindex])
                } else if (!is.null(CLD.w.references[[2]]) & !is.null(CLD.w.references[[3]])){
                    Lrefstring <- CLD.w.references[[2]]
                    Lrefindex <- which(names(x$parms) == CLD.w.references[[2]])
                    Lrefval <- unlist(x$parms[Lrefindex])
                    Drefstring <- CLD.w.references[[3]]
                    Drefindex <- which(names(x$parms) == CLD.w.references[[3]])
                    Drefval <- unlist(x$parms[Drefindex])
                    hrefval <- Lrefval + Drefval
                    hrefstring <- paste("L",Drefstring,sep=" + ")
                } else if (!is.null(CLD.w.references[[2]])){
                    hrefstring <- CLD.w.references[[2]]
                    hrefindex <- which(names(x$parms) == CLD.w.references[[2]])
                    hrefval <- unlist(x$parms[hrefindex])
                }
            }

            FGBarplot(t(totals), lab.y = FGMakeLabel(desc, specs$uom[uu]),
                      use.color = use.color, PlotTitle = PlotTitle,
                      leg.title = "Fishery", href=hrefval, hrefname=hrefstring)

            if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                                         GraphName = paste("CLD", paste0("LD",specs$data.units[uu]), sep = "."), graphics.type)
            if (plot.proportion)
            {
                FGBarplot(t(totals), lab.y = "Proportion",
                          use.color = use.color, PlotTitle = PlotTitle,
                          leg.title = "Fishery", proportion = TRUE)
                if (write.graphs) FGSavePlot(GraphicsDirName, DataName,
                                             GraphName = paste("CLD", paste0("LD",specs$data.units[uu]), "prop", sep = "."), graphics.type)
            } # end if plot proportion

        }  # End if n.mats > 0 ...
    }




    par(savepar)
}  # End function definition

