#' Plot side-by-side comparison of genotypes between rqtl2 and DOQTL methods for a single chromosome
#'
#' @param geno1 Imputed phase-known genotypes
#'     (as produced by \code{\link[qtl2geno]{maxmarg}}) or a list of
#'     three-dimensional arrays (as produced by \code{\link[qtl2geno]{guess_phase}}).
#' @param geno2 Imputed phase-known genotypes
#'       (from DOQTL)
#' @param map1 Marker map for qtl2 method (a list of vectors of marker positions).
#' @param map2 Marker map for doqtl method (a list of vectors of marker positions).
#' @param ind Individual to plot, either a numeric index or an ID.
#' @param chr Selected chromosome to plot; a vector of character strings.
#' @param col Vector of colors for the different genotypes.
#' @param na_col Color for missing segments.
#' @param border Color of outer border around chromosome rectangles.
#' @param shift If TRUE, shift the chromosomes so they all start at 0.
#' @param bgcolor Color for background rectangle
#' @param chrwidth Total width of rectangles for each chromosome, as a
#'        fraction of the distance between them.
#  @param file Output file name, as a string
#' @param ... Additional graphics parameters
#'
#' @export
#' @importFrom graphics plot rect par axis title abline box


plot_compare_x <-
    function(geno1, geno2, map1, map2, ind=NULL, chr=NULL, col=NULL, na_col="white",
             border="black", shift=FALSE, bgcolor="gray90",
             chrwidth=0.5, file = "compare_plot.png", ...)
{

    chr="X"

    #Subset to selected	individuals
    selected_ind <- ind
    if(!is.null(ind)) {
        geno1 <- subset(geno1, ind=selected_ind)
        geno2 <- subset(geno2, ind=selected_ind)
    }

    #Check if the maps are null
    if(is.null(map1)) stop("map1 is NULL")
    if(is.null(map2)) stop("map2 is NULL")

    #Check chromosome width requirement... see parameter description above
    stopifnot(chrwidth > 0 && chrwidth < 1)

    # ignore class of geno objects
    for(k in 1:length(geno1[,20]))
	geno1[k,] <- unclass(geno1[k,])
    	geno2[k,] <- unclass(geno2[k,])

###WARNING: Need to check objects (geno marker chr) for "X" chromosome existence ... will add later


    
    #Check for same number of markers between map and geno
	
    nmar_map1 <- sapply(map1, length)                                                     #list with number of markers per chromosome in map
    for(k in 1:length(geno1[,20])) {
    	nmar_geno1 <- sapply(geno1[k,], ncol)                                                     #list with number of markers per chromosome in geno
    	if(any(nmar_geno1 != nmar_map1))                                                      #stops if mismatch
        	stop("Mismatch between numbers of markers between geno1 and map1 on chr ",
             		paste(names(geno1[k,])[nmar_geno1 != nmar_map1], collapse=", "))
    	for(i in seq_along(geno1[k,])) {
           if(!all(names(map1[[i]]) == colnames(geno1[k,][[i]])))
        	    stop("Mismatch between marker names on chr ", names(geno1[k,])[i])
    	}
    }

    nmar_map2 <- sapply(map2, length)                                                     #list with number of markers per chromosome in map
    for(k in 1:length(geno2[,20])) {
        nmar_geno2 <- sapply(geno2[k,], ncol)                                                     #list with number of markers per chromosome in geno
	if(any(nmar_geno2 != nmar_map2))                                                      #stops if mismatch
                stop("Mismatch between numbers of markers between geno2 and map2 on chr ",
                        paste(names(geno2[k,])[nmar_geno2 != nmar_map2], collapse=", "))
        for(i in seq_along(geno2[k,])) {
           if(!all(names(map2[[i]]) == colnames(geno2[k,][[i]])))
                    stop("Mismatch between marker names on chr ", names(geno2[k,])[i])
        }
    }

    # shift map to start at 0
    if(shift) map1 <- lapply(map1, function(a) a-min(a,na.rm=TRUE))
    if(shift) map2 <- lapply(map2, function(a) a-min(a,na.rm=TRUE))


  plot_compare_internal <-
        function(geno1, geno2, map1, map2,ind, col=NULL, na_col="white",
                 border="black", bgcolor="gray90",
                 chrwidth=0.5,
                 xlab="Individual", ylab="Position (Mbp)",
                 ylim=NULL, main="", las=1, xaxs="i", yaxs="r",
                 mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
                 hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
                 vlines=NULL, vlines.col="gray80", vlines.lwd=1, vlines.lty=1,
                 ...)

  {
	#essentially extends possible arguments in functions...? Check this later
        dots <- list(...)

        #Set y limits for plot
        #Updated to take max from both methods as lim
        #ATTN !!! NEED TO UPDATE.... JUST SELECTED 1st ONE
        if(is.null(ylim)) ylim <- rev(range(unlist(map1), na.rm=TRUE))

        #Set x limits accoring to number of individuals
        nchr <- length(map1)
        nind <- length(geno1[,20])
        xlim <- c(0.5, nind+0.5)

        #Create Base Plot (type=n is has nothing)
        plot(0, 0, type="n", xlab="", ylab="", main=main,
             xaxs=xaxs, yaxs=yaxs,
             xaxt="n", yaxt="n", xlim=xlim, ylim=ylim)

        #Sets a vector of the form c(x1,x2,y1,y2) giving the extremes of the user coordinates of the plotting region
        u <- par("usr")
        if(!is.null(bgcolor))
            rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)

        # include axis labels?
        if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
        if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")

#### REVISIT LATER
#       # add x axis unless par(xaxt="n")
#       if(dots$xaxt != "n") {
#           odd <- seq(1, nind, by=2)
#           axis(side=1, at=odd, names(map1)[odd],                                      #check what is happening here
#                mgp=mgp.x, las=las, tick=FALSE)
#           if(nind > 1) {
#               even <- seq(2, nchr, by=2)
#               axis(side=1, at=even, names(map)[even],                                 #check what is happening here
#                    mgp=mgp.x, las=las, tick=FALSE)
#           }
#	}


	# add y axis unless par(yaxt="n")
        if(dots$yaxt != "n") {
            axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
        }

	# grid lines
        if(!(length(vlines)==1 && is.na(vlines))) {
            if(is.null(vlines)) vlines <- 1:nchr
            abline(v=vlines, col=vlines.col, lwd=vlines.lwd, lty=vlines.lty)
        }

	if(!(length(hlines)==1 && is.na(hlines))) {
            if(is.null(hlines)) hlines <- pretty(ylim)
            abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
        }

	# x and y axis labels
        title(xlab=xlab, mgp=mgp.x)
        title(ylab=ylab, mgp=mgp.y)


        #Sets colors, recycles if more than 8 genotypes (or more than prespecified colors)
        max_geno <- max(max(unlist(geno1[1,]), na.rm=TRUE), max(unlist(geno2[1,]), na.rm=TRUE))
        if(is.null(col)) {
            if(max_geno <= 8) {
                col <- qtl2plot::CCcolors
            }
            else {
                warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
                col <- rep(qtl2plot::CCcolors, max_geno)
            }
	}
	else if(max_geno > length(col)) {
            warning("not enough colors; recycling them")
            col <- rep(col, max_geno)
        }

        #Sets number of chromosomes
        #likely wont be needed anymore
        nchr <- length(map1)
        chr_num <- which(names(map1) == chr)

        #for each individual
        #creates a plot object for each individual
        for(i in 1:nind) {
            #subset to only a particular chromosome
            g1 <- geno1[i,][[chr_num]]
            g2 <- geno2[i,][[chr_num]]

            this_chrwidth <- chrwidth

                rect(i-this_chrwidth/2, min(map1[[chr_num]], na.rm=TRUE),
                     i, max(map1[[chr_num]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map2[[chr_num]], na.rm=TRUE),
                     i, max(map2[[chr_num]], na.rm=TRUE),
                     col=na_col, border=border, lend=1, ljoin=1)

                addgenorect(g1[1,,1], map1[[chr_num]], i-this_chrwidth/2, i,            #g is ind x marker x parent ... for male x chrom we only have mom
                            col=col)
                addgenorect(g2[1,,1], map2[[chr_num]], i+this_chrwidth/2, i,
                            col=col)

                rect(i-this_chrwidth/2, min(map1[[chr_num]], na.rm=TRUE),
                     i, max(map1[[chr_num]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)
                rect(i+this_chrwidth/2, min(map2[[chr_num]], na.rm=TRUE),
                     i, max(map2[[chr_num]], na.rm=TRUE),
                     col=NULL, border=border, lend=1, ljoin=1)

	}


	#Creates box around whole plot
        box()


  } #end internal

  png(paste0("plots/",file))
  plot_compare_internal(geno1, geno2, map1, map2, col=col, na_col=na_col, border=border,
                          bgcolor=bgcolor, chrwidth=chrwidth, file = file,...)
  dev.off()

} #end plot_compare_x


# add rectangles for the genotypes
addgenorect <-
    function(geno, map, x1, x2, col)
{
    intervals <- geno2intervals(geno, map)
    if(is.null(intervals) || nrow(intervals) < 1) return(NULL)

    for(i in 1:nrow(intervals)) {
        rect(x1, intervals[i,1],
             x2, intervals[i,2],
             col=col[intervals[i,3]],
             border=NA, lend=1, ljoin=1)
    }
}

# convert vector of integer genotypes to intervals with common genotypes
# (start, end, genotype)
geno2intervals <-
    function(geno, map)
{
    if(all(is.na(geno))) return(NULL)

    stopifnot(length(geno) == length(map))

    # drop missing values
    map <- map[!is.na(geno)]
    geno <- geno[!is.na(geno)]

    d <- diff(geno)
    xo_int <- which(d != 0)

    data.frame(lo=map[c(1,xo_int+1)],
               hi=map[c(xo_int, length(map))],
               geno=geno[c(xo_int, length(map))])

}


