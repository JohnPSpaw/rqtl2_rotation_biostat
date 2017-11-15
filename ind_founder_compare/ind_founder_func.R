#absdiff <- readRDS("../abs_diff/computed_data/absdiff.rds")
#attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])

letters <- c("A","B","C","D","E","F","G","H")


ind_founder_compare <- function(ind, chr, markers){

	##### PRINT A VECTOR WITH PROPORTION OF MATCHES BETWEEN IND AND EACH FOUNDER ALONG GIVEN MARKERS ######

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
	print(paste0("Proportion of markers in ",ind," matching founders:"))
	print(round(proportion_match,4))	




	###### PRINT TABLE OF MAX PROB PREDICTED GENOTYPES FOR BOTH RQTL2 AND DOQTL ######

	#load genotype probs
	pr_rqtl2 <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_",ind_name,".rds"))
        pr_doqtl <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_",ind_name,".rds"))
	
	#determine which genotype is predicted in rqtl2 along selected marker range
	writeLines("\n r/qrtl2 predicted geno \n")
	predicted_genos(chr=chr,markers=markers,probs=pr_rqtl2)

        writeLines("\n \n doqtl predicted geno \n")
        predicted_genos(chr=chr,markers=markers,probs=pr_doqtl)
	#writeLines("\n")


	##### PRINT 8x8 MATRIX OF FOUNDER MATCHES ######
	founder_matrix <- matrix(nrow=8,ncol=8)
	for(i in 1:8) {
		for(j in 1:8) {
			founder_matrix[i,j] <- mean(fgg[i,] == fgg[j,], na.rm=TRUE)
		}
	}
	
	rownames(founder_matrix) <- letters
	colnames(founder_matrix) <- letters

	print(round(founder_matrix,3))


}

predicted_genos <- function(chr, markers,probs) { 
	predict <- rep(NA,length(markers))
	print(length(markers))
	max_prob <- rep(NA,length(markers))
	count <- 0
	for(i in markers) {
		#print(i)
		count <- count + 1
                prob_36 <- probs[[chr]][1,,i]
                predict[count] <- names(which(prob_36 == max(prob_36)))
		max_prob[count] <- max(prob_36)
	}
	predict_freq <- table(predict)
	print(predict_freq)
	print(paste0("Prob: ", round(mean(max_prob),3)))
}




#ind_founder_compare(ind=10,chr=16,markers=c(577,870))
