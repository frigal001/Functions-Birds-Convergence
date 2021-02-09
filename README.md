# Deterministic assembly and anthropogenic extinctions drive convergence of island bird communities

Code accompanying the paper: Triantis et al. Deterministic assembly and anthropogenic extinctions drive convergence of island bird communities (Submitted)

The R scripts can be used to run the analyses in the paper. The code Code_CreateArtificialData.R is to create artifical datasets for three hypothetic archipelagos along with their respective species pool. The code Code_AllFunctions.R contains the main functions to calculate and test morphological and phylogenetic convergence among archipelagos.

The code Code_CreateArtificialData.R creates four datasets:

coloAge.Rdata contains the names of the archipelagos, the numbers of colonization events per archipelago and the maximum geological age for each archipelago.

data.community.Rdata contains the data for the three archipelago together. This is a list with (1) the archipelago x species matrix, (2) the traits, (3) the phylogenetic tree and (4) the statu (endemic, non-endemic) and the family for each species. This data are specifically used to calculate the observed morphological and phylogenetic convergence.  

data.archip.Rdata contains a list containing for each archipelgo (1) the phylogenetic tree, (2) the traits, (3) the statu (endemic, non-endemic) and the family for each species, (4) the phylogenetic tree of the species pool, (5) the traits of the species pool (6) the rate from the Brownian motion model for each family contained in the pool and (7) the birth and death rate for each family contained in the pool.

sppoolGlobal.Rdata is the global phylogenetic tree combining all the pools.
