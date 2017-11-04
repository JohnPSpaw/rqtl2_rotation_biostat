#John Spaw

#Plot RMS difference between rqtl2 and doqtl geno probs at each marker



plot_RMS_x <- 
	function(geno1, geno2, map1, map2, ind=NULL, chr="X", col=NULL, na_col="white",
		border="black", shift=FALSE, bgcolor="gray90", file="rms_plot.pdf", ...)

{
	#Insert warning if not X chromosome


    	#Subset to selected individuals
    	#Will only matter once this function is extended to multiple individuals
	selected_ind <- ind
    	if(!is.null(ind)) {
        	geno1 <- subset(geno1, ind=selected_ind)
        	geno2 <- subset(geno2, ind=selected_ind)
    	}
	nind <- dim(geno1[[20]])[1]

	#Check if the maps are null
    	if(is.null(map1)) stop("map1 is NULL")
    	if(is.null(map2)) stop("map2 is NULL")

    	#ignore class of geno objects
	for(k in 1:nind) {
		geno1[[20]][k,,] <- unclass(geno1[[20]][k,,])
   	 	geno2[[20]][k,,] <- unclass(geno2[[20]][k,,])
	}

    	#Check for same number of markers between map and geno
	#################ALERT: NEED TO ADD

	#Ensure that geno1 and geno2 have the same markers
	#################ALERT: NEED TO ADD

    	# shift map to start at 0
    	if(shift) map1 <- lapply(map1, function(a) a-min(a,na.rm=TRUE))
    	if(shift) map2 <- lapply(map2, function(a) a-min(a,na.rm=TRUE))	


        nmar <- dim(geno1[[20]])[3]
        marker_num <- seq(1:nmar)
        rms <- function(error) { sqrt(mean(error^2))}
	
	rmsd_list <- list()
	
	for(k in 1:nind) {
		rmsd <- rep(NA,nmar)
        	for(j in 1:nmar) {
        	        difference <- geno1[[20]][k,,j] - geno2[[20]][k,,j]
        	        rmsd[j] <- rms(difference)
        	}
		rmsd_list[[k]] <- rmsd
	}

	plot_RMS_internal <- 
		function(geno1,geno2, map1, map2,na_col="white",
                 border="black", bgcolor="gray90",
                 xlab="Marker", ylab="Root Mean Square Diff",
                 ylim=NULL, main="", las=1, xaxs="i", yaxs="r",
                 mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
                 hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
                 vlines=NULL, vlines.col="gray80", vlines.lwd=1, vlines.lty=1,
                 ...)

	{
		#essentially extends possible arguments in functions...? Check this later
        	dots <- list(...)

        	#Set y limits for plot
		ymax <- max(sapply(rmsd_list, max)) 	
        	ylim <- c(0,ymax + 0.2*ymax)
        	xlim <- c(0, nmar)


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

		# add y axis unless par(yaxt="n")
        	if(dots$yaxt != "n") {
        	    axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
        	}

		if(!(length(hlines)==1 && is.na(hlines))) {
        	    if(is.null(hlines)) hlines <- pretty(ylim)
        	    abline(h=hlines, col=hlines.col, lwd=hlines.lwd, lty=hlines.lty)
        	}

		# x and y axis labels
        	title(xlab=xlab, mgp=mgp.x)
        	title(ylab=ylab, mgp=mgp.y)
		
		col <- qtl2plot::CCcolors
		for(k in 1:nind) {
			lines(marker_num,rmsd_list[[k]], lwd=1.5, col=col[k])
		}
	}

  	pdf(paste0("plots/",file))
  	plot_RMS_internal(geno1, geno2, map1, map2, col=col, na_col=na_col, border=border,
                          bgcolor=bgcolor, file = file,...)
	dev.off()


}


