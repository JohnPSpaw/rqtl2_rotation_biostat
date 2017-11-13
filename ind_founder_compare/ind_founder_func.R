#absdiff <- readRDS("../abs_diff/computed_data/absdiff.rds")
#attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])

ind_founder_compare <- function(ind, chr, markers){

	#convert index of individual to individual name
	#ind is index with respect to absdiff
	ind_name <- ind_names[ind] #name of ind... such as "DO-101"... primary key between absdiff and attieDO
	ind_attieDO <- which(dimnames(attieDO$cross_info)[[1]] == ind_name) #index in of ind in attieDO

	#raw genotypes (matrix)
	g <-attieDO$geno[[chr]]; fg <-attieDO$founder_geno[[chr]] #individuals x markers ..... 1 = AA, 3 = BB, 2 = AB

	#markers
	markers <- colnames(g)[markers[1]:markers[2]]

	#genotypes 
	gg <- g[ind_attieDO, markers]; gg[gg==0] <- NA
	fgg <- fg[,markers]; fgg[fgg==0] <- NA

	proportion_match <- rep(NA,8)
	for(i in 1:8) {
		proportion_match[i] <- mean(gg == fgg[i,], na.rm=TRUE) 
	}
	names(proportion_match) <- c("A","B","C","D","E","F","G","H")
	print(proportion_match)	



	#load genotype probs
	pr_rqtl2_71 <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_",ind_name,".rds"))
        pr_doqtl_71 <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_",ind_name,".rds"))

	#determine which genotype is predicted in rqtl2 along selected marker range
	
}





ind_founder_compare(ind=10,chr=16,markers=c(577,870))
