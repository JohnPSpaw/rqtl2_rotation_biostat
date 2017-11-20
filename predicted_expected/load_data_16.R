#Load Cross Data
attieDO <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO.rds")
print("Loaded cross")

#Load rqtl2 and doqtl geno probs
inds <- c(109,146,163,164,170,171,180,185,186,187,195,200,203,204,205,217,223,224,225,226,22,233,254,259,25,264,266)
for(i in inds) {
	var_name_rqtl2 <- paste0("pr_rqtl2_",i)
	var_name_doqtl <- paste0("pr_doqtl_",i)
	
	load_rqtl2 <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_probs_byind/attieDO_probs_DO-",i,".rds"))
	load_doqtl <- readRDS(paste0("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/probs4qtl2_DO-",i,".rds"))        

	assign(var_name_rqtl2, load_rqtl2)
	assign(var_name_doqtl, load_doqtl)
}
print("Loaded Probs")

#Load markers
pmap_rqtl2 <- readRDS("/z/Proj/attie/kbroman/AttieDO/CalcGenoProb/attieDO_pmap.rds")
pmap_doqtl <- readRDS("/z/Proj/attie/kbroman/AttieDO/DerivedData/attie_36state_probs/pmap_doqtl.rds")
print("Loaded markers")

