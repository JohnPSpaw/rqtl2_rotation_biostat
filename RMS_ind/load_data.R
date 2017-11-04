#Load data for test≈ing plot_RMS_x.R
#Contains multiple individual for now
#John Spaw

library(qtl2)
	
	#functions for subsetting pr objects to common markers
	geno_common_markers1 <- function(geno1,geno2)
	{
		for(i in seq(along=geno1)) {
			mar1 <- dimnames(geno1[[i]])[[3]]
                        mar2 <-	dimnames(geno2[[i]])[[3]]
			mar <- mar1[mar1 %in% mar2]
			geno1[[i]] <- geno1[[i]][,,mar,drop=FALSE]
		}
		return(geno1)
	}

        geno_common_markers2 <- function(geno1,geno2)
        {
                for(i in seq(along=geno1)) {
                        mar1 <- dimnames(geno1[[i]])[[3]]
                        mar2 <-	dimnames(geno2[[i]])[[3]]
                        mar <- mar1[mar1 %in% mar2]
                        geno2[[i]] <- geno2[[i]][,,mar,drop=FALSE]
                }
                return(geno2)
        }

	#read data files
        pr_rqtl2_71 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-71.rds")
        pr_rqtl2_72 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-72.rds")
        pr_rqtl2_73 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-73.rds")
        pr_rqtl2_74 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-74.rds")
        pr_rqtl2_75 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-75.rds")
        pr_rqtl2_76 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-76.rds")
        pr_rqtl2_77 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-77.rds")
        pr_rqtl2_78 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-78.rds")
        pr_rqtl2_79 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-79.rds")

        pr_doqtl_71 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-71.rds")
        pr_doqtl_72 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-72.rds")
        pr_doqtl_73 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-73.rds")
        pr_doqtl_74 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-74.rds")
        pr_doqtl_75 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-75.rds")
        pr_doqtl_76 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-76.rds")
        pr_doqtl_77 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-77.rds")
        pr_doqtl_78 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-78.rds")
        pr_doqtl_79 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-79.rds")


	#subset to common markers
	pr_rqtl2_71 <- geno_common_markers1(pr_rqtl2_71,pr_doqtl_71)
        pr_rqtl2_72 <- geno_common_markers1(pr_rqtl2_72,pr_doqtl_72)
        pr_rqtl2_73 <- geno_common_markers1(pr_rqtl2_73,pr_doqtl_73)
        pr_rqtl2_74 <- geno_common_markers1(pr_rqtl2_74,pr_doqtl_74)
        pr_rqtl2_75 <- geno_common_markers1(pr_rqtl2_75,pr_doqtl_75)
        pr_rqtl2_76 <- geno_common_markers1(pr_rqtl2_76,pr_doqtl_76)
        pr_rqtl2_77 <- geno_common_markers1(pr_rqtl2_77,pr_doqtl_77)
        pr_rqtl2_78 <- geno_common_markers1(pr_rqtl2_78,pr_doqtl_78)
        pr_rqtl2_79 <- geno_common_markers1(pr_rqtl2_79,pr_doqtl_79)

        pr_doqtl_71 <- geno_common_markers1(pr_doqtl_71,pr_rqtl2_71)
        pr_doqtl_72 <- geno_common_markers1(pr_doqtl_72,pr_rqtl2_72)
        pr_doqtl_73 <- geno_common_markers1(pr_doqtl_73,pr_rqtl2_73)
        pr_doqtl_74 <- geno_common_markers1(pr_doqtl_74,pr_rqtl2_74)
        pr_doqtl_75 <- geno_common_markers1(pr_doqtl_75,pr_rqtl2_75)
        pr_doqtl_76 <- geno_common_markers1(pr_doqtl_76,pr_rqtl2_76)
        pr_doqtl_77 <- geno_common_markers1(pr_doqtl_77,pr_rqtl2_77)
        pr_doqtl_78 <- geno_common_markers1(pr_doqtl_78,pr_rqtl2_78)
        pr_doqtl_79 <- geno_common_markers1(pr_doqtl_79,pr_rqtl2_79)


	#combine into a single list of matrices
        pr_rqtl2 <- rbind(pr_rqtl2_71,pr_rqtl2_72,pr_rqtl2_73,pr_rqtl2_74,pr_rqtl2_75, pr_rqtl2_76,pr_rqtl2_77,pr_rqtl2_78,pr_rqtl2_79)
        pr_doqtl <- rbind(pr_doqtl_71,pr_doqtl_72,pr_doqtl_73,pr_doqtl_74,pr_doqtl_75,pr_doqtl_76,pr_doqtl_77,pr_doqtl_78,pr_doqtl_79)

        #load marker maps
        pmap_rqtl2 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_pmap.rds")
        pmap_doqtl <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/pmap_doqtl.rds")

        #load cross object
        #use same cross object for rqtl2 and doqtl
        attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")

