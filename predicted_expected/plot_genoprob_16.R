
inds <- sort(c(22,109,171,186,223,224,233))

for(i in inds) {
        print(i)
        pr1 <- get(paste0("pr_rqtl2_",i))
        pr2 <- get(paste0("pr_doqtl_",i))
        ind_name <- paste0("DO-",i)

        filename1 <- paste0("plots/genoprob_16_",ind_name,"_rqtl2.pdf")
        filename2 <- paste0("plots/genoprob_16_",ind_name,"_doqtl.pdf")

	main1 <- paste0("r/qtl2 ",ind_name)
        main2 <- paste0("DOQTL ",ind_name)

        pdf(filename1)
                plot_genoprob(probs=pr1, map=pmap_rqtl2, ind=ind_name,chr=16, main=main1)
        dev.off()

        pdf(filename2)
                plot_genoprob(probs=pr2, map=pmap_doqtl, ind=ind_name,chr=16, main=main2)
        dev.off()
}


