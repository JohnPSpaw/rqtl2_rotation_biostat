#compare_genoprob(probs1=pr_rqtl2_203, probs2=pr_doqtl_203, cross=attieDO, ind="DO-203",chr=16,minprob=0.9, annotate=TRUE)

for(i in inds) {
	print(i)
	pr1 <- get(paste0("pr_rqtl2_",i))
        pr2 <- get(paste0("pr_doqtl_",i))
	ind_name <- paste0("DO-",i)
	print(compare_genoprob(probs1=pr1, probs2=pr2, cross=attieDO, ind=ind_name,chr=16,minprob=0.9, annotate=TRUE))
}
