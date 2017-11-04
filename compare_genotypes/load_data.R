#Load data for test≈ing plot_compare.R
#John Spaw

library(qtl2)

########################################################################
#This block is me being naive.... I need a way to systematically handle individuals
#Will update later

#Individuals 71-75 are all male

        pr_rqtl2_71 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-71.rds")
        pr_rqtl2_72 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-72.rds")
        pr_rqtl2_73 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-73.rds")
        pr_rqtl2_74 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-74.rds")
        pr_rqtl2_75 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-75.rds")
        pr_rqtl2_76 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-76.rds")
        pr_rqtl2_77 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-77.rds")
        pr_rqtl2_78 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-78.rds")
        pr_rqtl2_79 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-79.rds")
	pr_rqtl2 <- rbind(pr_rqtl2_71,pr_rqtl2_72,pr_rqtl2_73,pr_rqtl2_74,pr_rqtl2_75, pr_rqtl2_76,pr_rqtl2_77,pr_rqtl2_78,pr_rqtl2_79)


        pr_doqtl_71 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-71.rds")
        pr_doqtl_72 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-72.rds")
        pr_doqtl_73 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-73.rds")
        pr_doqtl_74 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-74.rds")
        pr_doqtl_75 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-75.rds")
        pr_doqtl_76 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-76.rds")
        pr_doqtl_77 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-77.rds")
        pr_doqtl_78 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-78.rds")
        pr_doqtl_79 <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-79.rds")
	pr_doqtl <- rbind(pr_doqtl_71,pr_doqtl_72,pr_doqtl_73,pr_doqtl_74,pr_doqtl_75,pr_doqtl_76,pr_doqtl_77,pr_doqtl_78,pr_doqtl_79)

        #load marker maps
        pmap_rqtl2 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_pmap.rds")
        pmap_doqtl <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/pmap_doqtl.rds")

        #max marg probs
        m_rqtl2_71 <- maxmarg(pr_rqtl2_71)
        m_rqtl2_72 <- maxmarg(pr_rqtl2_72)
        m_rqtl2_73 <- maxmarg(pr_rqtl2_73)
        m_rqtl2_74 <- maxmarg(pr_rqtl2_74)
        m_rqtl2_75 <- maxmarg(pr_rqtl2_75)
        m_rqtl2_76 <- maxmarg(pr_rqtl2_76)
        m_rqtl2_77 <- maxmarg(pr_rqtl2_77)
        m_rqtl2_78 <- maxmarg(pr_rqtl2_78)
        m_rqtl2_79 <- maxmarg(pr_rqtl2_79)
	m_rqtl2 <- rbind(m_rqtl2_71,m_rqtl2_72,m_rqtl2_73,m_rqtl2_74,m_rqtl2_75,m_rqtl2_76,m_rqtl2_77,m_rqtl2_78,m_rqtl2_79)        

	m_doqtl_71 <- maxmarg(pr_doqtl_71)
        m_doqtl_72 <- maxmarg(pr_doqtl_72)
        m_doqtl_73 <- maxmarg(pr_doqtl_73)
        m_doqtl_74 <- maxmarg(pr_doqtl_74)
        m_doqtl_75 <- maxmarg(pr_doqtl_75)
        m_doqtl_76 <- maxmarg(pr_doqtl_76)
        m_doqtl_77 <- maxmarg(pr_doqtl_77)
        m_doqtl_78 <- maxmarg(pr_doqtl_78)
        m_doqtl_79 <- maxmarg(pr_doqtl_79)
	m_doqtl <- rbind(m_doqtl_71, m_doqtl_72, m_doqtl_73,m_doqtl_74,m_doqtl_75,m_doqtl_76,m_doqtl_77,m_doqtl_78,m_doqtl_79)

        #load cross object
        #use same cross object for rqtl2 and doqtl
        attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")

        #guess phase
        ph_rqtl2_71 <- guess_phase(attieDO, m_rqtl2_71)
        ph_rqtl2_72 <- guess_phase(attieDO, m_rqtl2_72)
        ph_rqtl2_73 <- guess_phase(attieDO, m_rqtl2_73)
        ph_rqtl2_74 <- guess_phase(attieDO, m_rqtl2_74)
        ph_rqtl2_75 <- guess_phase(attieDO, m_rqtl2_75)
        ph_rqtl2_76 <- guess_phase(attieDO, m_rqtl2_76)
        ph_rqtl2_77 <- guess_phase(attieDO, m_rqtl2_77)
        ph_rqtl2_78 <- guess_phase(attieDO, m_rqtl2_78)
        ph_rqtl2_79 <- guess_phase(attieDO, m_rqtl2_79)
	ph_rqtl2 <- rbind(ph_rqtl2_71,ph_rqtl2_72,ph_rqtl2_73,ph_rqtl2_74,ph_rqtl2_75,ph_rqtl2_76,ph_rqtl2_77,ph_rqtl2_78,ph_rqtl2_79)
#        ph_rqtl2 <- rbind(ph_rqtl2_71,ph_rqtl2_72,ph_rqtl2_73,ph_rqtl2_74,ph_rqtl2_75)

	ph_doqtl_71 <- guess_phase(attieDO, m_doqtl_71)
        ph_doqtl_72 <- guess_phase(attieDO, m_doqtl_72)
        ph_doqtl_73 <- guess_phase(attieDO, m_doqtl_73)
        ph_doqtl_74 <- guess_phase(attieDO, m_doqtl_74)
        ph_doqtl_75 <- guess_phase(attieDO, m_doqtl_75)
        ph_doqtl_76 <- guess_phase(attieDO, m_doqtl_76)
        ph_doqtl_77 <- guess_phase(attieDO, m_doqtl_77)
        ph_doqtl_78 <- guess_phase(attieDO, m_doqtl_78)
        ph_doqtl_79 <- guess_phase(attieDO, m_doqtl_79)
	ph_doqtl <- rbind(ph_doqtl_71,ph_doqtl_72,ph_doqtl_73,ph_doqtl_74,ph_doqtl_75,ph_doqtl_76,ph_doqtl_77,ph_doqtl_78,ph_doqtl_79)
#       ph_doqtl <- rbind(ph_doqtl_71,ph_doqtl_72,ph_doqtl_73,ph_doqtl_74,ph_doqtl_75)

########################################################################


