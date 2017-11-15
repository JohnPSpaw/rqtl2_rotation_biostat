ind_markers <- read.csv("input_data/input_20.csv")
m <- matrix(nrow=nrow(ind_markers), ncol=24)
for(i in 1:nrow(ind_markers)) {
        print(i)
        vec <- auto_ind_founder_compare(ind=ind_markers[i,2],chr=ind_markers[i,1],markers=c(ind_markers[i,3],ind_markers[i,4]))
        m[i,] <- vec
}

names(m) <- c("Chr", "Ind", "ML", "MR", "A","B","C","D","E","F","G","H","rqtl2_geno","rqtl2_prob","doqtl_geno","doqtl_prob",
                   "fm1","fm1_prop", "fm2","fm2_prop","f3","fm3_prop","fm4","fm4_prop")

m <- m[1:nrow(ind_markers),1:24]
write.csv(m,"output_data/explore_20.csv")
