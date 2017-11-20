for(i in inds) {
        print(i)
        pr1 <- get(paste0("pr_rqtl2_",i))
        pr2 <- get(paste0("pr_doqtl_",i))
        ind_name <- paste0("DO-",i)
	filename <- paste0("plots/genocompare_16_",ind_name,".pdf")

	pdf(filename)
		plot_genoprobcomp(probs1=pr1, probs2=pr2, map=pmap_rqtl2, ind=ind_name,chr=16)
	dev.off()
}





