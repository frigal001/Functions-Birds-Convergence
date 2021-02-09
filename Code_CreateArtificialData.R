library(Rphylopars)
library(mvMORPH)
library(tidyverse)
library(TreeTools)

############################################################################
############################################################################
############################################################################
############################################################################

# create dataset archipelago A, B and C with their respective species pool 
# sppool is the tree of the pool, sparchi the tree of the archipelago, traitpool
# the traits of the pool (5 as in the paper) and traitarchi the traits of the archipelago
# for sppool with combined 2 trees of 50 species to differentiate two different families in the pool
# for each family, we will extract the BM rate. Therefore, in the course of our simulation 
# procedure and for a given archipelago A, if a species S, belonging to the family F, 
# has been randomly selected from the pool to colonize archipelago A, 
# and subsequently to give birth to an endemic clade of four species, we (1) used the 
# geological age of archipelago A, and Lambda, mu estimated for the tree of the family F to 
# create the tree of the endemic clade, and (2) the trait value of S as the ancestral state and sigma2 
# estimated for the tree and traits of family F, to generate trait values for the four endemic 
# species along the BD tree previously generated.

F1 <- pbtree(n=50, scale=100, tip.label = paste("t", 1:50, sep=""))
F2 <- pbtree(n=50, scale=100, tip.label = paste("t", 51:100, sep=""))
F3 <- pbtree(n=50, scale=100, tip.label = paste("t", 101:150, sep=""))
F4 <- pbtree(n=50, scale=100, tip.label = paste("t", 151:200, sep=""))

traitF1 <- replicate(5, rTraitCont(F1, model = "BM", sigma = 0.1, root.value = 0))
traitF2 <- replicate(5, rTraitCont(F2, model = "BM", sigma = 0.1, root.value = 0))
traitF3 <- replicate(5, rTraitCont(F3, model = "BM", sigma = 0.1, root.value = 0))
traitF4 <- replicate(5, rTraitCont(F4, model = "BM", sigma = 0.1, root.value = 0))

## the global tree including all the tree pools involved in the analyses i.e. the 4 families

rootSp <- pbtree(n=4, scale=20)
t1 <- bind.tree(rootSp, F1, 1)
t2 <- bind.tree(t1, F2, 1)
t3 <- bind.tree(t2, F3, 1)
sppoolGlobal <- bind.tree(t3, F4, 1)
plot(sppoolGlobal)

# and the pool got each archipelago A, B and C

sppoolA <- drop.tip(sppoolGlobal, which(!sppoolGlobal$tip.label%in%c(F1$tip.label, F2$tip.label)))
sppoolB <- drop.tip(sppoolGlobal, which(!sppoolGlobal$tip.label%in%c(F1$tip.label, F3$tip.label)))
sppoolC <- drop.tip(sppoolGlobal, which(!sppoolGlobal$tip.label%in%c(F1$tip.label, F4$tip.label)))

FamPoolA <- gl(2, 50, labels = c("F1", "F2"))%>%as.vector
FamPoolB <- gl(2, 50, labels = c("F1", "F3"))%>%as.vector
FamPoolC <- gl(2, 50, labels = c("F1", "F4"))%>%as.vector

traitpoolA <- rbind(traitF1, traitF2)%>%as.data.frame;traitpoolA$Fam = FamPoolA
traitpoolB <- rbind(traitF1, traitF3)%>%as.data.frame;traitpoolB$Fam = FamPoolB
traitpoolC <- rbind(traitF1, traitF4)%>%as.data.frame;traitpoolC$Fam = FamPoolC
colnames(traitpoolA) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5', "Fam")
colnames(traitpoolB) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5', "Fam")
colnames(traitpoolC) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5', "Fam")

rateF1 <- phylopars(data.frame(species = F1$tip.label, traitF1), F1, model = "BM", REML = F)
rateF2 <- phylopars(data.frame(species = F2$tip.label, traitF2), F2, model = "BM", REML = F)
rateF3 <- phylopars(data.frame(species = F3$tip.label, traitF3), F3, model = "BM", REML = F)
rateF4 <- phylopars(data.frame(species = F4$tip.label, traitF4), F4, model = "BM", REML = F)


bdF1 <- bd(birthdeath(F1))
bdF2 <- bd(birthdeath(F2))
bdF3 <- bd(birthdeath(F3))
bdF4 <- bd(birthdeath(F4))


list_rateA <- list("F1" = rateF1, "F2"=rateF2)
list_rateB <- list("F1" = rateF1, "F3"=rateF3)
list_rateC <- list("F1"=rateF1, "F4"=rateF4)

list_bdA <- list("F1" = bdF1, "F2" = bdF2)
list_bdB <- list("F1" = bdF1, "F3" = bdF3)
list_bdC <- list("F1" = bdF1, "F4" = bdF4)

# create data for 3 archipelagos

sparchiA <- pbtree(n=10, scale=10, tip.label = paste("A", 1:10, sep=""))
sparchiB <- pbtree(n=10, scale=10, tip.label = paste("B", 1:10, sep=""))
sparchiC <- pbtree(n=10, scale=10, tip.label = paste("C", 1:10, sep=""))


traitarchiA <- replicate(5, rTraitCont(sparchiA, model = "BM", sigma = 0.1, root.value = 0))
traitarchiB <- replicate(5, rTraitCont(sparchiB, model = "BM", sigma = 0.1, root.value = 0))
traitarchiC <- replicate(5, rTraitCont(sparchiC, model = "BM", sigma = 0.1, root.value = 0))
colnames(traitarchiA) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5')
colnames(traitarchiB) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5')
colnames(traitarchiC) <- c("tr1", 'tr2', "tr3", 'tr4', 'tr5')


traits.archipelago <- rbind(traitarchiA, traitarchiB, traitarchiC)
phy.archipelago <- pbtree(n=30, scale=10, tip.label = rownames(traits.archipelago))
matarchi <- table(gl(3, 10, labels = c("A", "B", "C")), 
                  c(phy.archipelago$tip.label))%>%as.data.frame.matrix

statuA <- rep("End", 10)
statuB <- c(rep("NoEnd", 8), rep("End", 2))
statuC <- c(rep("NoEnd", 2), rep("End", 8))

FamA <- gl(2, 5, labels = c("F1", "F2"))%>%as.vector
FamB <- gl(2, 5, labels = c("F1", "F3"))%>%as.vector
FamC <- gl(2, 5, labels = c("F1", "F4"))%>%as.vector


############################################################################
############################################################################
############################################################################
############################################################################

## data with colonization event, age of the archipelago ##

coloAge <- data.frame(Archipelago=c("A", "B", "C"), nb = c(2, 10, 5), age = c(10, 5, 15))


## community data to calculate the observed nnd with community matrix, traits, tree and data 
## for the statu and the family assigned to each species

df.archipelago <- data.frame("Statu" = c(statuA, statuB, statuC), 
                             "Family" = c(FamA, FamB, FamC))
rownames(df.archipelago) <- phy.archipelago$tip.label

data.community <- list("comm" = matarchi, 
                         "traits" = traits.archipelago, 
                         "phy" = phy.archipelago, 
                         "df" = df.archipelago)

## global data as a list with all information per archipelago including: (1) the tree, (2) the trait
## (3) the data with statu (End/NoEnd) and family, (4) the trait pool, (5) the phylo pool (6) the sigma and (7) the 
## bd rate for the families included in the pool

data.archip <- list("A" = list("phy" = sparchiA, "trait" = traitarchiA, 
                               "df" = data.frame(statuA, FamA), "phypool" = sppoolA, 
                               "traitpool"= traitpoolA, "rate" = list_rateA, "bd"= list_bdA), 
                    "B" = list("phy" = sparchiB, "trait" = traitarchiB, 
                               "df" = data.frame(statuB, FamB), "phypool" = sppoolB, 
                               "traitpool"= traitpoolB, "rate" = list_rateB, "bd"= list_bdB),
                    "C" = list("phy" = sparchiC, "trait" = traitarchiC, 
                               "df" = data.frame(statuC, FamC), "phypool" = sppoolC, 
                               "traitpool"= traitpoolC, "rate" = list_rateC, "bd"= list_bdC))


setwd("/Users/francoisrigal/Documents/POST_DOC_EXTINCTION/ARTICLES/Papier_BIRDS_ARCHIPELAGO/CodeCleanSubmitted")

save(coloAge, file="coloAge.Rdata")
save(data.community, file="data.community.Rdata")
save(data.archip, file="data.archip.Rdata")
save(sppoolGlobal, file="sppoolGlobal.Rdata")



