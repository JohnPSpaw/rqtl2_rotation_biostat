#absdiff <- readRDS("../abs_diff/computed_data/absdiff.rds")
#attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")
ind_names <- (dimnames(absdiff[[1]])[1])[[1]]
ind_index <- seq(1,dim(absdiff[[1]])[1])

letters <- c("A","B","C","D","E","F","G","H")
ind_vec_names <- c("Chr", "Ind", "ML", "MR", "A","B","C","D","E","F","G","H","rqtl2_geno","rqtl2_prob","doqtl_geno","doqtl_prob",
                   "fm1","fm1_prop", "fm2","fm2_prop","f3","fm3_prop","fm4","fm4_prop")


auto_ind_founder_compare <- function(ind, chr, markers){

	ind_vec <- c(NA,length(ind_vec_names))	

        ind_vec[1] <- chr
        ind_vec[2] <- ind
        ind_vec[3] <- markers[1]
        ind_vec[4] <- markers[2]

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
                ind_vec[i+4] <- round(mean(gg == fgg[i,], na.rm=TRUE),3)
        }
	

        ###### PRINT TABLE OF MAX PROB PREDICTED GENOTYPES FOR BOTH RQTL2 AND DOQTL ######

        #load genotype probs
        pr_rqtl2 <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_",ind_name,".rds"))
        pr_doqtl <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_",ind_name,".rds"))
	
	print("probs loaded")

        #determine which genotype is predicted in rqtl2 and doqtl along selected marker range
	ind_vec[13] <- predicted_geno(chr,markers,probs=pr_rqtl2)
	print("geno done")
	ind_vec[14] <- predicted_geno_prob(chr,markers,probs=pr_rqtl2)
	print("probs done")
        ind_vec[15] <- predicted_geno(chr,markers,probs=pr_doqtl)
	print("geno 2 done")
        ind_vec[16] <- predicted_geno_prob(chr,markers,probs=pr_doqtl)
	print("probs 2 done")
        ##### PRINT 8x8 MATRIX OF FOUNDER MATCHES ######
        founder_matrix <- matrix(nrow=8,ncol=8)
        for(i in 1:8) {
                for(j in 1:8) {
                        founder_matrix[i,j] <- mean(fgg[i,] == fgg[j,], na.rm=TRUE)
                }
        }

	rownames(founder_matrix) <- letters
        colnames(founder_matrix) <- letters

	rqtl2_geno_str <- strsplit(ind_vec[13],"")
	doqtl_geno_str <- strsplit(ind_vec[15],"")

	ind_vec[17] <- paste0(rqtl2_geno_str[[1]][1],doqtl_geno_str[[1]][1])
	ind_vec[18] <- founder_combo_match(letter1 = rqtl2_geno_str[[1]][1], letter2 = doqtl_geno_str[[1]][1], fmat=founder_matrix)

        ind_vec[19] <- paste0(rqtl2_geno_str[[1]][1],doqtl_geno_str[[1]][2])
	#print(ind_vec)        
	ind_vec[20] <- founder_combo_match(letter1 = rqtl2_geno_str[[1]][1], letter2 = doqtl_geno_str[[1]][2], fmat=founder_matrix)

        ind_vec[21] <- paste0(rqtl2_geno_str[[1]][2],doqtl_geno_str[[1]][1])
        ind_vec[22] <- founder_combo_match(letter1 = rqtl2_geno_str[[1]][2], letter2 = doqtl_geno_str[[1]][1], fmat=founder_matrix)

        ind_vec[23] <- paste0(rqtl2_geno_str[[1]][2],doqtl_geno_str[[1]][2])
        ind_vec[24] <- founder_combo_match(letter1 = rqtl2_geno_str[[1]][2], letter2 = doqtl_geno_str[[1]][2],  fmat=founder_matrix)

	names(ind_vec) <- c("Chr", "Ind", "ML", "MR", "A","B","C","D","E","F","G","H","rqtl2_geno","rqtl2_prob","doqtl_geno","doqtl_prob",
                   "fm1","fm1_prop", "fm2","fm2_prop","f3","fm3_prop","fm4","fm4_prop")

	return(ind_vec)
}


predicted_geno <- function(chr, markers,probs) {
        predict <- rep(NA,length(markers))
        count <- 0
        for(i in markers) {
                count <- count + 1
		print(i)
		#print(length(probs[[chr]][1,,i]))
                prob_36 <- (probs[[chr]][1,,i])
                predict[count] <- names(which(prob_36 == max(prob_36)))
        }
	predict_freq <- table(predict)
	#print(predict_freq)
	return(names(predict_freq[1]))
}

predicted_geno_prob <- function(chr, markers,probs) {
        max_prob <- rep(NA,length(markers))
        count <- 0
        for(i in markers) {
                count <- count + 1
                prob_36 <- (probs[[chr]][1,,i])
                max_prob[count] <- max(prob_36)
        }
	return(round(mean(max_prob),3))
}

founder_combo_match <- function(letter1,letter2,fmat) {
	index1 <- which(letters==letter1)
	index2 <- which(letters==letter2)
	prop <- round(fmat[index1,index2],3)
	return(prop)
}
