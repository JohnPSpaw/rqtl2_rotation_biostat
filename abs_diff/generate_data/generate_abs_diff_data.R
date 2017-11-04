#Second attempt at generating sum abs diff between genotype probs for R/qtl2 and DOQTL
#Based off of Karl Broman's code in /z/Proj/attie/kbroman/AttieDO/DOQTL_v_Rqtl2/calc_prob_diff.R
#
#result: matrix for each chromosome, ind x position, with values in [0,2]
#
#note: ignore's males in X chromosome due to issue in DOQTL source code


#################################################################################
# determine individuals
#################################################################################
#doqtl files and individuals
doqtl_dir <- "/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/"
doqtl_files <- list.files(doqtl_dir, pattern="^probs4qtl2")				#prob files for individuals
doqtl_ind <- as.numeric(sapply(strsplit(doqtl_files, "[-\\.]"), "[", 2))		#vector of individual numbers (483 in total)

#rqtl2 files and individuals
rqtl2_dir <- "/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/"
rqtl2_files <- list.files(rqtl2_dir, pattern="^attieDO_probs_DO-")
rqtl2_ind <- as.numeric(sapply(strsplit(rqtl2_files, "[-\\.]"), "[", 2))

#individuals in common
ind <- rqtl2_ind[rqtl2_ind %in% doqtl_ind] # 483 ind
doind <- paste0("DO-", ind) # as "DO-#"							#vector of string ID's for individualso

#sex of individuals
load("/z/Proj/attie/kbroman/AttieDO/DerivedData/pheno_clin.RData")
sex <- setNames(pheno_clin[,"sex"], pheno_clin[,"mouse"])
sex <- sex[doind]

#################################################################################
#calculate abs diff
#################################################################################
absdiff <- vector("list", 20)
names(absdiff) <- c(1:19,"X")
for(i in seq_along(ind)) {
    message("ind ", i, " of ", length(ind))
    doqtl_prob <- readRDS(paste0(doqtl_dir, "probs4qtl2_", doind[i], ".rds"))
    rqtl2_prob <- readRDS(paste0(rqtl2_dir, "attieDO_probs_", doind[i], ".rds"))

    for(chr in c(1:19, "X")) {
        mar_doqtl <- dimnames(doqtl_prob[[chr]])[[3]]
        mar_rqtl2 <- dimnames(rqtl2_prob[[chr]])[[3]]
        mar <- mar_rqtl2[mar_rqtl2 %in% mar_doqtl]					#common markers

        if(is.null(absdiff[[chr]])) {
            #subset to only females when using X chr
            if(chr=="X") these_ind <- names(sex)[sex=="F"]
            else these_ind <- doind

	    #pre-generate matrix to hold absdiff
            absdiff[[chr]] <- matrix(nrow=length(these_ind), ncol=length(mar))
            dimnames(absdiff[[chr]]) <- list(these_ind, mar)
        }

	if(chr=="X" && sex[doind[i]] == "M") next # skip X chr for males

        absdiff[[chr]][doind[i], ] <-
            colSums(abs(doqtl_prob[[chr]][1,,mar] - rqtl2_prob[[chr]][1,,mar]))
    }
}

saveRDS(absdiff, "../computed_data/absdiff.rds")


#################################################################################
#markers
#################################################################################

# get genetic and physical maps
broman_dir <- "/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/"
pmap <- readRDS(paste0(broman_dir, "attieDO_pmap.rds"))
gmap <- readRDS(paste0(broman_dir, "attieDO_gmap.rds"))

#subset to relevent markers
for(chr in seq_along(pmap)) {
    pmap[[chr]] <- pmap[[chr]][colnames(absdiff[[chr]])]
    gmap[[chr]] <- gmap[[chr]][colnames(absdiff[[chr]])]
}

dir <- "/z/Proj/attie/spaw/abs_diff/computed_data"
saveRDS(pmap, paste0(dir,"/absdiff_pmap.rds"))
saveRDS(gmap, paste0(dir,"/absdiff_gmap.rds"))
