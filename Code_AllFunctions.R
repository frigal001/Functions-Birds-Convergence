
library(partitions)
library(mvMORPH)
library(phytools)
library(Rphylopars)
library(TreeSim)
library(picante)
library(adephylo)
library(ggplot2)


############################################################################
############################################################################
############################################################################
############################################################################

# function fct_phy_colo_spe_nnd_turn to calculate nnd between archipelago with the 
# reconstructed trees from the simulations

fct_phy_colo_spe_nnd_turn <- function(List_data, sppoolGlobal, nsim) {
  
  archipelago.names <- names(List_data)
  
  na <- 1:length(List_data)
  
  compare.l <- function(i, j){
    
    a1 <- List_data[[i]][[nsim]]
    a2 <- List_data[[j]][[nsim]]
    
    
    a1[[1]] <- a1[[1]][a1[[2]]$tip.label,]
    a2[[1]] <- a2[[1]][a2[[2]]$tip.label,]
    
    sp1 <- na.omit(c(rownames(a1[[1]][which(a1[[1]]$origin%in%c("colonization", "anagenesis")),]),
                     a1[[3]]$species.cladogenesis))
    
    
    sp2 <- na.omit(c(rownames(a2[[1]][which(a2[[1]]$origin%in%c("colonization", "anagenesis")),]),
                     a2[[3]]$species.cladogenesis))
    
    tree.between <- drop.tip(sppoolGlobal, which(!sppoolGlobal$tip.label%in%c(sp1, sp2)))
    
    Ast <- rownames(a1[[1]])
    Bst <- rownames(a2[[1]])
    
    all.sp <- (unique(c(sp1, sp2)))
    
    sp.com.a1 <- intersect(Ast, a2[[3]]$species.cladogenesis)
    sp.com.a2  <- intersect(Bst, a1[[3]]$species.cladogenesis)
    sp.in.com <- c(sp.com.a1, sp.com.a2)
    
    species.to.remove <- na.omit(c(a1[[3]]$species.cladogenesis, a2[[3]]$species.cladogenesis))
    species.to.remove  <- species.to.remove[which(!species.to.remove%in%sp.in.com)]
    
    #tree.between
    
    list_pair <- list(a1, a2)
    
    for (k in 1:2){
      
      if(is.na(list_pair[[k]][[3]]$placement.clades[1])){
        
        tree.between <- tree.between
        
      } else {
        
        n <- length(list_pair[[k]][[3]]$placement.clades)
        
        for (w in 1:n){
          tr <- list_pair[[k]][[3]]$tr.cladogenesis[[w]]
          placement.clades <- list_pair[[k]][[3]]$placement.clades[[w]]
          species.cladogenesis <- list_pair[[k]][[3]]$species.cladogenesis[w]
          where <- which(tree.between$tip.label == species.cladogenesis)
          
          nodes<-sapply(tree.between$tip.label,function(x,y) which(y==x),y=tree.between$tip.label)
          edge.lengths<-setNames(tree.between$edge.length[sapply(nodes,
                                                                 function(x,y) which(y==x),y=tree.between$edge[,2])],names(nodes))
          
          age.anc <- edge.lengths[species.cladogenesis]
          
          if(placement.clades >= age.anc) {
            placement <- age.anc
          } else {placement <- placement.clades}
          
          
          
          tree.between <- bind.tree(tree.between, tr, where = where, position = placement)
          
        }
      }
      
    }
    
    tree.between <- drop.tip(tree.between, which(tree.between$tip.label%in%species.to.remove))
    
    Ast <- rownames(a1[[1]])
    Bst <- rownames(a2[[1]])
    compmat <- cophenetic(tree.between)[Ast, Bst]
    Ann <- apply(as.matrix(compmat),1,min)
    Bnn <- apply(as.matrix(compmat),2,min)
    #Dnn <- mean(c(Ann, Bnn))
    turn <- min(c(mean(Ann),mean(Bnn)))
    turn
    
    return(turn)
  }
  
  nnd <- pairwise.table(compare.l, na, "none")
  
  dt <- rbind(rep(0, length(na)), cbind(nnd, rep(0, length(na)-1)))
  dts <- as.matrix(as.dist(dt))
  diag(dts) <- NA
  rownames(dts) <- (archipelago.names)
  colnames(dts) <- (archipelago.names)
  return(as.dist(dts))
  
}

# function fct_nnd_turnover to calculate nnd between archipelago with the 
# reconstructed traits from the simulations

fct_nnd_turnover <- function(mat, dist){
  
  sp.list <- list()
  for(i in rownames(mat)) {sp.list[[i]] <- colnames(mat)[mat[i,]==1]}
  
  na = 1:dim(mat)[1]
  compare.l <- function(i, j) {
    Asp     <- sp.list[[i]]
    Bsp     <- sp.list[[j]]
    compmat <- dist[Asp,Bsp]
    Ann     <- apply(as.matrix(compmat),1,min)
    Bnn     <- apply(as.matrix(compmat),2,min)
    Dnn     <- mean(c(Ann, Bnn))
    turn    <- min(c(mean(Ann),mean(Bnn)))
    return(turn)
  }
  
  h <- dim(mat)[1]
  turn.nnd <- pairwise.table(compare.l, na, "none")
  dt <- rbind(rep(0, h), cbind(turn.nnd, rep(0, h-1)))
  dts <- as.matrix(as.dist(dt))
  diag(dts) <- NA
  rownames(dts) <- names(sp.list)
  colnames(dts) <- names(sp.list)
  
  return(as.dist(dts))
  
}

# function to get the species age

fct_age <- function(phy){
  dage <- as.matrix(distTips(phy))
  diag(dage) <- NA
  age <- apply(dage, 2, function(x) min(x, na.rm = T)/2)
  species.ages <- data.frame(mrca.age=age, species = phy$tip.label)
  row.names(species.ages)<-phy$tip.label
  return(species.ages)
}




############################################################################
############################################################################
############################################################################
############################################################################

# the three core functions to generate random traits and trees for each archipelago based on random
# colonization and speciation while keeping the number of species and number of endemic species constant. Number of 
# colonization events was extracetd from the literature. Three functions are available but could be 
# probably merged to simplify the code. The functions need as inputs the number of species that will generate 
# cladogenetic events, number of species that will stay unchanged (native non-endemic) and the number of species 
# that will give birth to anagenetic species. 


fct_cladogenesis <- function(data.archip, n.cladogenesis, nb.do.clado, age, n.native, species.age, names){
  
  arc <- data.archip[[names]]

  species.cladogenesis <- n.cladogenesis
  age.clado <- NA
  for (k in 1:length(n.cladogenesis)){
    age.clado[k] <- species.age[which(rownames(species.age)%in%species.cladogenesis[k]),]$mrca.age
    }
    age.clado <- ifelse(age.clado > age, age, age.clado)

    # Prepare data to generate clado: the sigma and bd for the family the species belong to + age
    data_species_clado <- species.age[species.cladogenesis,]
    data_species_clado$family.clado <- arc$traitpool[species.cladogenesis,]$Fam
    data_species_clado$age <- age
    
    sigma.clado <- arc$rate[data_species_clado$family.clado]
    bd.clado <- arc$bd[data_species_clado$family.clado]
    
    cladogenesis <- list()
    placement.clades <- list()
  
    for (k in 1:length(nb.do.clado))
    {
     tr.clade <- sim.bd.taxa.age(nb.do.clado[k], 1, bd.clado[[k]][1], bd.clado[[k]][2], 1, age = age.clado[k], mrca = F)[[1]]
     tr.clade$root.edge <-NULL
     tr.clade$tip.label <- paste(names, "clado", k, 1:length(tr.clade$tip.label), sep = "_")
     placement.clades[[k]] <- max(node.age(tr.clade)$age)
     cladogenesis[[k]] <- tr.clade
    }
  
    sim.trait.clado <- list()

    for (k in 1:length(species.cladogenesis))
    {
      sim.trait.clado[[k]] <-mvSIM(cladogenesis[[k]], nsim = 1, model = "BM1", param = 
                                   list(theta = arc$traitpool[species.cladogenesis[k],][-ncol(arc$traitpool)]
                                        , sigma = sigma.clado[[k]]$pars$phylocov,
                                        ntraits=ncol(arc$traitpool)-1))
    }
  
  tree.colo <- drop.tip(arc$phypool, which(!arc$phypool$tip.label%in%c(n.native, species.cladogenesis)))
  plot(tree.colo)
  axisPhylo()
  
  for (i in 1:length(species.cladogenesis)){
    tr <- cladogenesis[[i]]
    where <- which(tree.colo$tip.label==species.cladogenesis[i])
    placement <- placement.clades[[i]]
    tree.colo <- bind.tree(tree.colo, tr, where = where, position = placement)
  }
  
  
  tree.colo <- drop.tip(tree.colo, species.cladogenesis)
  
  
  trait.clado <- do.call(rbind, sim.trait.clado)
  trait.native <-as.matrix(arc$traitpool[n.native,][-ncol(arc$traitpool)])
  rd.trait <- data.frame(rbind(trait.native, trait.clado))
  
  dim.nat <- ifelse(length(trait.native) == 0, 0, dim(trait.native)[1])
  dim.clado <- ifelse(length(trait.clado) == 0, 0, dim(trait.clado)[1])
  
  rd.trait$origin <- c(rep("colonization", dim.nat), rep("cladogenesis", dim.clado))
  
  tree.info <- list(cladogenesis, placement.clades, species.cladogenesis)
  names(tree.info) <- c("tr.cladogenesis", "placement.clades", "species.cladogenesis")
  
  final <- data.frame(rd.trait)

  return(list(final, tree.colo, tree.info))
  
}

fct_anagenesis <- function(data.archip, n.anagenesis, age, n.native, species.age, names){
  
  arc <- data.archip[[names]]
  
    species.anagenesis <- n.anagenesis
  age.ana <- NA
  for (k in 1:length(n.anagenesis)){
    age.ana[k] <- species.age[which(rownames(species.age)%in%n.anagenesis[k]),]$mrca.age
  }
  age.ana <- ifelse(age.ana > age, age, age.ana)
  
  # Prepare data to generate ana: the sigma and bd for the family the species belong to + age
  data_species_ana <- species.age[species.anagenesis,]
  data_species_ana$family.ana <- arc$traitpool[species.anagenesis,]$Fam
  data_species_ana$age <- age
  
  sigma.ana <- arc$rate[data_species_ana$family.ana]
  bd.ana <- arc$bd[data_species_ana$family.ana]
  
  anagenesis <- list()
  placement.clades.ana <- list()
  
  for (k in 1:length(n.anagenesis))
  {
    tr.anag <- sim.bd.taxa.age(2, 1, bd.ana[[k]][1], bd.ana[[k]][2], 1, age = age.ana[k], mrca = F)[[1]]
    tr.anag$root.edge <-NULL
    tr.anag$tip.label <- paste(names, "Ana", k, 1:length(tr.anag$tip.label), sep = "_")
    placement.clades.ana[[k]] <- max(dist.nodes(tr.anag))/2
    anagenesis[[k]] <- tr.anag
  }

  sim.trait.ana <- list()
  
  for (k in 1:length(n.anagenesis))
  {
    sim.trait.ana[[k]] <-mvSIM(anagenesis[[k]], nsim = 1, model = "BM1", param = 
                                   list(theta = arc$traitpool[species.anagenesis[k],][-ncol(arc$traitpool)]
                                        , sigma = sigma.ana[[k]]$pars$phylocov,
                                        ntraits=ncol(arc$traitpool)-1))[1,]
    }
  
  tree.colo <- drop.tip(arc$phypool, which(!arc$phypool$tip.label%in%c(n.native, species.anagenesis)))
  
  trait.ana <- do.call(rbind, sim.trait.ana)
  rownames(trait.ana) <- species.anagenesis
  trait.native <-as.matrix(arc$traitpool[n.native,][-ncol(arc$traitpool)])
  rd.trait <- data.frame(rbind(trait.native, trait.ana))
  
  dim.nat <- ifelse(length(trait.native) == 0, 0, dim(trait.native)[1])
  dim.ana <- ifelse(length(trait.ana) == 0, 0, dim(trait.ana)[1])
  
  rd.trait$origin <- c(rep("colonization", dim.nat), rep("anagenesis", dim.ana))
  
  tree.info <- list(NA, NA, NA)
  names(tree.info) <- c("clado.BM.shape", "placement.clades", "species.cladogenesis")
  
  final <- data.frame(rd.trait)
  
  return(list(final, tree.colo, tree.info))
}

fct_ana_clado <- function(data.archip, n.anagenesis, n.cladogenesis, nb.do.clado, age, n.native, species.age , names){
  
  arc <- data.archip[[names]]

  species.cladogenesis <- n.cladogenesis
  age.clado <- NA
  for (k in 1:length(n.cladogenesis)){
    
    age.clado[k] <- species.age[which(rownames(species.age)%in%species.cladogenesis[k]),]$mrca.age
  }
  age.clado <- ifelse(age.clado > age, age, age.clado)
  

  
  
  # Prepare data to generate clado: the sigma and bd for the family the species belong to + age
  data_species_clado <- species.age[species.cladogenesis,]
  data_species_clado$family.clado <- arc$traitpool[species.cladogenesis,]$Fam
  data_species_clado$age <- age
  
  sigma.clado <- arc$rate[data_species_clado$family.clado]
  bd.clado <- arc$bd[data_species_clado$family.clado]
  
  cladogenesis <- list()
  placement.clades <- list()
  
  for (k in 1:length(nb.do.clado))
  {
    tr.clade <- sim.bd.taxa.age(nb.do.clado[k], 1, bd.clado[[k]][1], bd.clado[[k]][2], 1, age = age.clado[k], mrca = F)[[1]]
    tr.clade$root.edge <-NULL
    tr.clade$tip.label <- paste(names, "clado", k, 1:length(tr.clade$tip.label), sep = "_")
    placement.clades[[k]] <- max(node.age(tr.clade)$age)
    cladogenesis[[k]] <- tr.clade
  }
  
  sim.trait.clado <- list()
  
  for (k in 1:length(species.cladogenesis))
  {
    sim.trait.clado[[k]] <-mvSIM(cladogenesis[[k]], nsim = 1, model = "BM1", param = 
                                   list(theta = arc$traitpool[species.cladogenesis[k],][-ncol(arc$traitpool)]
                                        , sigma = sigma.clado[[k]]$pars$phylocov,
                                        ntraits=ncol(arc$traitpool)-1))
  }
  
  

  species.anagenesis <- n.anagenesis
  age.ana <- NA
  for (k in 1:length(n.anagenesis)){

    age.ana[k] <- species.age[which(rownames(species.age)%in%n.anagenesis[k]),]$mrca.age
  }
  age.ana <- ifelse(age.ana > age, age, age.ana)
  
  # Prepare data to generate ana: the sigma and bd for the family the species belong to + age
  data_species_ana <- species.age[species.anagenesis,]
  data_species_ana$family.ana <- arc$traitpool[species.anagenesis,]$Fam
  data_species_ana$age <- age
  
  sigma.ana <- arc$rate[data_species_ana$family.ana]
  bd.ana <- arc$bd[data_species_ana$family.ana]
  
  anagenesis <- list()
  placement.clades.ana <- list()
  
  for (k in 1:length(n.anagenesis))
  {
    tr.anag <- sim.bd.taxa.age(2, 1, bd.ana[[k]][1], bd.ana[[k]][2], 1, age = age.ana[k], mrca = F)[[1]]
    tr.anag$root.edge <-NULL
    tr.anag$tip.label <- paste(names, "Ana", k, 1:length(tr.anag$tip.label), sep = "_")
    placement.clades.ana[[k]] <- max(dist.nodes(tr.anag))/2
    anagenesis[[k]] <- tr.anag
  }
  
  sim.trait.ana <- list()
  
  for (k in 1:length(n.anagenesis))
  {
    sim.trait.ana[[k]] <-mvSIM(anagenesis[[k]], nsim = 1, model = "BM1", param = 
                                 list(theta = arc$traitpool[species.anagenesis[k],][-ncol(arc$traitpool)]
                                      , sigma = sigma.ana[[k]]$pars$phylocov,
                                      ntraits=ncol(arc$traitpool)-1))[1,]
  }
  
  
  tree.colo <- drop.tip(arc$phypool, which(!arc$phypool$tip.label%in%c(n.native, species.cladogenesis, species.anagenesis)))
  
  for (i in 1:length(species.cladogenesis)){
    
    tr <- cladogenesis[[i]]
    where <- which(tree.colo$tip.label==species.cladogenesis[i])
    placement <- placement.clades[[i]]
    tree.colo <- bind.tree(tree.colo, tr, where = where, position = placement)
  }
  tree.colo <- drop.tip(tree.colo, species.cladogenesis)
  
  trait.ana <- do.call(rbind, sim.trait.ana)
  rownames(trait.ana) <- species.anagenesis
  trait.clado <- do.call(rbind, sim.trait.clado)
  
  trait.native <- as.matrix(arc$traitpool[n.native,][-ncol(arc$traitpool)])
  rd.trait <- data.frame(rbind(trait.native, trait.ana, trait.clado))
  
  dim.nat <- ifelse(length(trait.native) == 0, 0, dim(trait.native)[1])
  dim.ana <- ifelse(length(trait.ana) == 0, 0, dim(trait.ana)[1])
  dim.clado <- ifelse(length(trait.clado) == 0, 0, dim(trait.clado)[1])
  
  rd.trait$origin <- c(rep("colonization", dim.nat), 
                      rep("anagenesis", dim.ana),
                      rep("cladogenesis", dim.clado))
  
  tree.info <- list(cladogenesis, placement.clades, species.cladogenesis)
  names(tree.info) <- c("tr.cladogenesis", "placement.clades", "species.cladogenesis")
  
  final <- data.frame(rd.trait)
  
  return(list(final, tree.colo, tree.info))
  
}

############################################################################
############################################################################
############################################################################
############################################################################

# function prepare_data_simulation to prepare the new traits and trees of each archipelago that will 
# be used to calcualte the simulated nnd

prepare_data_simulation <- function(data.archip, coloAge, names, nsim, verbose = T) {
  
  age <- coloAge[coloAge$Archipelago == names,]$age
  nbco <- coloAge[coloAge$Archipelago == names,]$nb
  arc <- data.archip[[names]]
  
  phy_pool <- arc$phypool
  species.age <- fct_age(phy_pool)
  sp.region <- phy_pool$tip.label
  
  data.output <- list()
  
  spAna <- list()
  spClado <- list()
  
  for (j in 1:nsim){
    
      N <- length(arc$phy$tip.label);N
      End <- length(which(arc$df[,1]=="End"))
      sp.colonisation.native <-  N - End
      sp.colonisation.non.native <-  nbco - sp.colonisation.native
      ps <- restrictedparts(End, sp.colonisation.non.native, include.zero=FALSE, decreasing=F)
      distri.end <- ps[,sample(dim(ps)[2], 1)]
      
      
      sp.do.ana <- length(which(distri.end == 1))
      sp.do.clado <- length(which(distri.end > 1))
      nb.do.clado <- distri.end[distri.end>1]
      colonizers <- sample(sp.region, nbco)
      
      n.native <- sample(colonizers, sp.colonisation.native)
      n.anagenesis <- sample(colonizers[which(!colonizers%in%n.native)], sp.do.ana)
      n.cladogenesis <- colonizers[which(!colonizers%in%c(n.anagenesis, n.native))]
      
      spAna[[j]] <- n.anagenesis
      spClado[[j]] <- n.cladogenesis
      
      if (length(n.anagenesis) == 0 & length(n.cladogenesis) > 0) {
        outputs <- fct_cladogenesis(data.archip, n.cladogenesis, nb.do.clado, age, n.native, species.age, names)
        }
      
      if (length(n.anagenesis) > 0 & length(n.cladogenesis) == 0) {
        outputs <- fct_anagenesis(data.archip, n.anagenesis, age, n.native, species.age, names)
        }
      
      if (length(n.anagenesis) > 0 & length(n.cladogenesis) > 0) {
        outputs <- fct_ana_clado(data.archip, n.anagenesis, n.cladogenesis, nb.do.clado, age, n.native, species.age, names)
      }
      
    
      
      data.output[[j]] <- outputs
      
      if(verbose) print(j)
  }
  

  return(data.output)
  
  
}

# function Simulation_null_NND to calculate the simulated nnd

Simulation_null_NND <- function(List_data, sppoolGlobal){
  
  nsim <- length(List_data[[1]])
  
  names.archi <- names(List_data)
  
  
  List_data_end <- List_data
  
  for (i in 1:length(List_data_end))
  {
    for (j in 1:nsim)
    {
      l1 <- List_data_end[[i]][[j]][[1]]
      l1 <- l1[which(!l1$origin%in%("colonization")),]
      tr1 <- List_data_end[[i]][[j]][[2]]
      tr1 <- drop.tip(tr1, which(!tr1$tip.label%in%rownames(l1)))
      List_data_end[[i]][[j]][[1]] <- l1
      List_data_end[[i]][[j]][[2]] <- tr1
    }
  } 
  
  
  
  nnd.phy <- list()
  nnd.phy.end <- list()
  nnd.phy.nat <- list()
  
  nnd.trait <- list()
  nnd.trait.end<- list()
  nnd.trait.nat <- list()
  
  for (j in 1:nsim) {
    null.data <- list()
    for (i in 1:length(List_data)) {
      random.pool <- List_data[[i]][[j]][[1]]
      tree.pool <- List_data[[i]][[j]][[2]]
      random.pool$archipelago <- rep(names.archi[i], dim(random.pool)[1])
      random.pool$sp <- rownames(random.pool)
      null.data[[i]] <- random.pool
    }
    
    names(null.data) <- names.archi
    
    df.null.data <- do.call(rbind, null.data)
    df.null.data$full.sp <- rownames(df.null.data)
    
    df.null.colo <- df.null.data[df.null.data$origin == "colonization",]
    null.com.phy <- as.data.frame.matrix(table(df.null.colo$archipelago, df.null.colo$sp))
    tree.phy <- drop.tip(sppoolGlobal, which(!sppoolGlobal$tip.label%in%colnames(null.com.phy)))
    
    null.com <- as.data.frame.matrix(table(df.null.data$archipelago, df.null.data$full.sp))
    sp.nat <- unique(df.null.data[df.null.data$origin == "colonization",]$full.sp)
    sp.end <- unique(df.null.data[df.null.data$origin%in%c("anagenesis", "cladogenesis"),]$full.sp)
    
    null.com.end <- null.com[,sp.end]
    null.com.nat <- null.com[,sp.nat]
    
    null.com.end <- null.com.end[apply(null.com.end, 1, sum)>0,]
    null.com.end <- null.com.end[,apply(null.com.end, 2, sum)>0]
    
    null.com.nat <- null.com.nat[apply(null.com.nat, 1, sum)>0,]
    null.com.nat <- null.com.nat[,apply(null.com.nat, 2, sum)>0]
    
    names.end <- rownames(null.com.end)
    names.remove.end <- names(List_data_end)
    names.remove.end <- which(!names.remove.end%in%names.end)
    List_data_end[names.remove.end] <- NULL
    
    nnd.phy[[j]] <- fct_phy_colo_spe_nnd_turn(List_data, sppoolGlobal, j)
    nnd.phy.end[[j]] <- fct_phy_colo_spe_nnd_turn(List_data_end, sppoolGlobal, j)
    nnd.phy.nat[[j]] <- fct_nnd_turnover(null.com.phy, cophenetic(tree.phy))
    
    traits.names <- colnames(List_data[[1]][[j]][[1]])[-ncol(List_data[[1]][[j]][[1]])]
    
    null.trait <- data.frame(df.null.data[colnames(null.com),])[,traits.names]
    null.trait.end <- null.trait[colnames(null.com.end),]
    null.trait.nat <- null.trait[colnames(null.com.nat),]
    
    nnd.trait[[j]] <- fct_nnd_turnover(null.com, as.matrix(dist(null.trait)))
    nnd.trait.end[[j]] <- fct_nnd_turnover(null.com.end, as.matrix(dist(null.trait.end)))
    nnd.trait.nat[[j]] <- fct_nnd_turnover(null.com.nat, as.matrix(dist(null.trait.nat)))
    
    print(j) 
    
  }
  
  list_nnd_all <- list(nnd.trait, nnd.phy)
  names(list_nnd_all) <- c("nnd.trait", "nnd.phy")
  
  list_nnd_end <- list(nnd.trait.end, nnd.phy.end)
  names(list_nnd_end) <- c("nnd.trait.end", "nnd.phy.end")
  
  list_nnd_nat <- list(nnd.trait.nat, nnd.phy.nat)
  names(list_nnd_nat) <- c("nnd.trait.nat", "nnd.phy.nat")
  
  list_diss_null <- list(list_nnd_all, list_nnd_end, list_nnd_nat)
  
  names(list_diss_null) <- c("list_nnd_all", "list_nnd_end", "list_nnd_nat")
  
  return(list_diss_null)
}

# function observed_metrics_NND to  calculate the observed nnd

observed_metrics_NND <- function(list.archi){
  
  comm.end <- list.archi$comm[,list.archi$df$Statu == "End"]
  comm.nat <- list.archi$comm[,list.archi$df$Statu == "NoEnd"]
  
  comm.end <- comm.end[apply(comm.end, 1, sum)>0,]
  comm.nat <- comm.nat[apply(comm.nat, 1, sum)>0,]
  
  nnd.trait <- fct_nnd_turnover(list.archi$comm, as.matrix(dist(scale(list.archi$traits))))
  nnd.phy <- fct_nnd_turnover(list.archi$comm, cophenetic(list.archi$phy))
  
  nnd.trait.end <- fct_nnd_turnover(comm.end, as.matrix(dist(scale(list.archi$traits))))
  nnd.phy.end <- fct_nnd_turnover(comm.end, cophenetic(list.archi$phy))
  
  nnd.trait.nat <- fct_nnd_turnover(comm.nat, as.matrix(dist(scale(list.archi$traits))))
  nnd.phy.nat <- fct_nnd_turnover(comm.nat, cophenetic(list.archi$phy))
  

  list_nnd_all <- list(nnd.trait, nnd.phy)
  names(list_nnd_all) <- c("nnd.trait", "nnd.phy")
  
  list_nnd_end <- list(nnd.trait.end, nnd.phy.end)
  names(list_nnd_end) <- c("nnd.trait.end", "nnd.phy.end")
  
  list_nnd_nat <- list(nnd.trait.nat, nnd.phy.nat)
  names(list_nnd_nat) <- c("nnd.trait.nat", "nnd.phy.nat")
  
  list_diss_null <- list(list_nnd_all, list_nnd_end, list_nnd_nat)
  
  names(list_diss_null) <- c("list_nnd_all", "list_nnd_end", "list_nnd_nat")
  
  return(list_diss_null)
}

# function observed_metrics_NND to calculate the ses and pvalues

ses_function_complete <- function(list_diss_null, diss_obs){
  
  list_ses <- list()
  list_pv <- list()
  list_names <- list()
  list_vector_null <- list()
  list_obs_mean <- list()
  N <- length(list_diss_null)
  
  for (j in 1:N) {
    
    list_ses_part <- NA
    list_pv_part <- NA
    list_obs <- NA
    list_null <- list()
    
    M <- length(list_diss_null[[j]])
    
    for (i in 1:M)
    {
      
      mean.rd <- NA
      
      for (k in 1:length(list_diss_null[[j]][[i]]))
      {
        mean.rd[k] <- mean(list_diss_null[[j]][[i]][[k]])
      }
      
      list_ses_part[i] <- (mean(diss_obs[[j]][[i]]) - mean(mean.rd))/sd(mean.rd)
      list_pv_part[i] <-   sum(mean.rd <= mean(diss_obs[[j]][[i]]))/length(mean.rd)
      list_obs[i] <- mean(diss_obs[[j]][[i]]) 
      list_null[[i]] <- mean.rd
    }
    list_ses[[j]] <- list_ses_part
    list_pv[[j]] <- list_pv_part
    list_names[[j]] <- names(diss_obs[[j]])
    list_vector_null[[j]] <- do.call(rbind, list_null)
    list_obs_mean[[j]] <-  list_obs
  }
  
  df <- data.frame(metrics = unlist(list_names), obs = round(unlist(list_obs_mean), 3),
                   ses = round(unlist(list_ses),3), 
                   pv = round(unlist(list_pv), 4))
  
  null <- do.call(rbind, list_vector_null)
  rownames(df) <- df$metrics
  rownames(null) <- rownames(df)
  res <- list(df, null)
  names(res) <- c("df", "null")
  
  return(res)
}

############################################################################
############################################################################
############################################################################
############################################################################

# test with artifical data #
# see the code Code_CreateArtificialData #

load("coloAge.Rdata")
load("data.community.Rdata")
load("data.archip.Rdata")
load("sppoolGlobal.Rdata")

# step 1: we create a list for each archipelago with the N simulated data (traits and trees)

archipelago.names <- c("A", "B", "C")
List_data <- list()
for (i in 1:length(archipelago.names)){
  List_data[[i]] <- prepare_data_simulation(data.archip, coloAge, names = archipelago.names[i], nsim=100, verbose = T)
}
names(List_data) <- archipelago.names

# step 2: we simulate nnd for trait and phylo for all species, endemic and native non-endemic separately

sim.data <- Simulation_null_NND(List_data, sppoolGlobal)


#step3: we calculate the observed nnd for trait and phylo for all species, endemic and 
# native non-endemic separately 

diss_obs <- observed_metrics_NND(data.community)

# step 4: we compare the observed and simulated values to calculate the SES and the associated pvalues

res_ses <- ses_function_complete(sim.data, diss_obs)

res_ses$df # the results

##### violin plots #####

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






