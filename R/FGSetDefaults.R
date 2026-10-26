FGSetDefaults <-
function()
{
    fgopts <- list(
        graphics = list(
            fglayout = c(5,3)),
        color = list(
            clr.grid = "lightgray",         # for all grids
            lty.grid = "solid",             # ditto
            clr.pos = "orange",     # for bubble plots
            clr.neg = "royalblue3",         # for bubble plots
            clr.ang = "royalblue4",         # for lower panel, bubble plots
            clr.line = "royalblue4",
            clr.obsd = "darkorchid4",
            clr.pred = "black",
            clr.lightline = "wheat3",
            clr.line2 = "#1B9E77",
            clr.obsd2 = "#D95F02"),
        bw = list(
            clr.grid = "lightgray",
            lty.grid = "solid",
            clr.pos = gray(0.95),
            clr.neg = gray(0.40),
            clr.ang = "black",
            clr.line = "black",
            clr.obsd = "black",
            clr.pred = "black",
            clr.lightline = "gray60",
            clr.line2 = "gray75",
            clr.obsd2 = "gray35")
        )
    options(FGoptions = fgopts)
    #assign("FGoptions", fgopts, pos = 1)
    return(invisible(NULL))
}   # end function FGSetDefaults

