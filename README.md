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
these functions generate random traits and trees for each archipelago based on random colonization and speciation while keeping the number of species and number of endemic species constant. These three functions could be probably merged in one to optimize the code. These functions needs as inputs the number of species that will use the generate anagenesis and the cladogenetic event(s) and the number of species that will stay unchanged (native).

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

# Example

# Step 1: upload the datasets created with Code_CreateArtificialData

``` r
load("coloAge.Rdata")
load("data.community.Rdata")
load("data.archip.Rdata")
load("sppoolGlobal.Rdata")
```

# Step 2: we create a list for each archipelago with the N simulated data (traits and trees)
``` r
archipelago.names <- c("A", "B", "C")
List_data <- list()
for (i in 1:length(archipelago.names)){
  List_data[[i]] <- prepare_data_simulation(data.archip, coloAge, names = archipelago.names[i], nsim=100, verbose = T)
}
names(List_data) <- archipelago.names
```
# Step 3: we simulate MNTDturn for traits and trees for all species, endemic and native non-endemic species
``` r
sim.data <- Simulation_null_NND(List_data, sppoolGlobal)
```

# Step 4: we calculate the observed MNTDturn for trait and trees for all species, endemic and native non-endemic species 
``` r
diss_obs <- observed_metrics_NND(data.community)
```
# Step 5: we use the observed and simulated MNTDturn values to calculate the SES and the associated p-values
``` r
res_ses <- ses_function_complete(sim.data, diss_obs)

res_ses$df #the results
```
# Step 6: Plot results
``` r
traits_all <- data.frame(rd = (res_ses$null["nnd.trait",]))
traits_end <- data.frame(rd = (res_ses$null["nnd.trait.end",]))
traits_nat <- data.frame(rd = (res_ses$null["nnd.trait.nat",]))

obs_traits_all <- res_ses$df["nnd.trait",]$obs
obs_traits_end <- res_ses$df["nnd.trait.end",]$obs
obs_traits_nat <- res_ses$df["nnd.trait.nat",]$obs

pv_traits_all <- res_ses$df["nnd.trait",]$pv
pv_traits_end <- res_ses$df["nnd.trait.end",]$pv
pv_traits_nat <- res_ses$df["nnd.trait.nat",]$pv


df.box <- data.frame(sim = c(traits_all$rd, traits_end$rd, traits_nat$rd), 
                     grp = gl(3, 100, labels = c("ALL", "END", "NAT")))
df.box$grp <- factor(df.box$grp, levels = rev(levels(df.box$grp)))
df.points <- data.frame(obs = c(obs_traits_all, obs_traits_end, obs_traits_nat), 
                        grp = c("ALL", "END", "NAT"), limit = rep(0.1, 3), 
                        pvalue = c(paste("P=", pv_traits_all), 
                                   paste("P=", pv_traits_end),
                                   paste("P=", pv_traits_nat)))

ggplot() + 
  geom_violin(data = df.box, mapping = aes(grp, sim), fill= "grey", alpha=0.6, color = NA)  + coord_flip() +
  labs(y="Among archipelago turnover", x = "") + theme_bw() + 
  geom_point(data = df.points, mapping = aes(grp, obs), size = 4, shape = 19, alpha=1) + 
  theme(legend.position = "none") + scale_color_manual(values = c("red")) + 
  geom_label(data = df.points, mapping = aes(grp, limit, label = pvalue), size = 2, 
             hjust = 0, color = "black", fontface = "bold") 
```
















