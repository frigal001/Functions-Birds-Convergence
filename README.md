# Deterministic assembly and anthropogenic extinctions drive convergence of island bird communities

Code accompanying the paper: Triantis et al. Deterministic assembly and anthropogenic extinctions drive convergence of island bird communities (Submitted)

The R scripts can be used to run the analyses in the paper. The code Code_CreateArtificialData.R is to create artifical datasets for three hypothetic archipelagos along with their respective species pool. The code Code_AllFunctions.R contains the main functions to calculate and test morphological and phylogenetic convergence among archipelagos.

The code Code_CreateArtificialData.R creates four datasets:

coloAge.Rdata contains the names of the archipelagos, the numbers of colonization events per archipelago and the maximum geological age for each archipelago.

data.community.Rdata contains the data for the three archipelago together. This is a list with (1) the archipelago x species matrix, (2) the traits, (3) the phylogenetic tree and (4) the statu (endemic, non-endemic) and the family for each species. This data are specifically used to calculate the observed morphological and phylogenetic convergence.  

data.archip.Rdata contains a list containing for each archipelgo (1) the phylogenetic tree, (2) the traits, (3) the statu (endemic, non-endemic) and the family for each species, (4) the phylogenetic tree of the species pool, (5) the traits of the species pool (6) the rates from the Brownian motion model for each family contained in the pool and (7) the birth and death rates for each family contained in the pool.

sppoolGlobal.Rdata is the global phylogenetic tree combining all the pools.

The code Code_AllFunctions.R contains the functions:

``` r
fct_phy_colo_spe_nnd_turn
```
to calculate MNTDturn between archipelagos with the reconstructed trees from the simulations
``` r
fct_nnd_turnover
```
to calculate MNTDturn between archipelagos with the reconstructed traits from the simulations
``` r
fct_age
```
to calculate the species age in a given phylogeny


The code has also three core functions:
``` r
fct_cladogenesis
fct_anagenesis
fct_ana_clado
```
these functions generate random traits and tree for each archipelago based on random colonization and speciation while keeping the number of species and number of endemic species constant. Three functions could be probably merged in one to simplify the code. The function needs as inputs the number of species that will use the generate anagenesis and the cladogenetic event(s) and the number of species that will stay unchanged (native).

The function
``` r
prepare_data_simulation
```
is used to prepare and organize the new traits and new trees of each archipelago that will be used to calcualte the simulated MNTDturn

and the function
``` r
Simulation_null_NND
```
is used to simulate N values of random MNTDturn.

The function
``` r
observed_metrics_NND
```
calculates the observed MNTD

and the function 
``` r
ses_function_complete
```
calculate the ses and p-values for a one tailed-test (convergence)


















