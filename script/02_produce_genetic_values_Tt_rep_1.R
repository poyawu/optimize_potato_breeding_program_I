rm(list = ls()); gc()

library(AlphaSimR)
library(tibble)
library(dplyr)
source("../script/main_function.R")
############## procduce the genetic values for the target trait ##############
### take repetition = 1 as an example
r <- 1 
cat("rep:",r, "\n")

##### load the required files 
## crossplan 
load("random_crossplan_list_1000rep.RData")

## 2000 QTL effects 
load("01_2000QTL_A_D_effects.RData")

## genotypic information based on alphasimR format
load("01_100_parental_clones_AlphasimR_format.RData")

## obtain genetic map
genetic.map <- SP$.__enclos_env__$.__active__$genMap()

## the size of initial population 
n.seedling <- 300*3000


## the 300 crossplan at the rth repetition
crossplan <- crossplan.list[[r]]

## randomly determine 2000 QTL position
QTL.effect <- QTL.position.fun(genMap = genetic.map, effect.a = effect.a, effect.d = effect.d)
#### In the "QTL effect" ####
## 1st column: chromosome
## 2nd column: the location within a chromosome 
## 3st column: genetic position (morgen)
## 4th column: the location across all sites/markers
## 5th column: additive effects
## 6th column: dominance effects

## determine the number of progenies per cross evenly
n.progeny <- ceiling(n.seedling/nrow(crossplan))
n.progeny.vector <- rep(n.progeny, nrow(crossplan))
while(TRUE){
  if(sum(n.progeny.vector) > n.seedling){
    n.progeny.vector[which.max(n.progeny.vector)] <- max(n.progeny.vector) - 1
  }else if(sum(n.progeny.vector) == n.seedling)
    break
}


## produce genetic values of the target trait for an initital population within the rth repetition (the progenies of 300 crosses)
F1.infor <- do.call(dplyr::bind_rows, 
                    lapply(1:nrow(crossplan), function(i){
                      cat(i, "..")
                      plan <- matrix(crossplan[i,],1,2)
                      
                      ## new cross (F1)
                      F1.pop.tmp <- makeCross(P.pop, crossPlan = plan, nProgeny = n.progeny.vector[i], simParam = SP)
                      
                      ## obtain F1 genotypic information 
                      F1.sites.tmp <- pullSegSiteGeno(F1.pop.tmp, simParam = SP)
                      
                      ## extract QTL genotypic information in F1  
                      F1.QTL.tmp <- F1.sites.tmp[,QTL.effect$index]
                      
                      ## additive matrix 
                      matrix.A.tmp <- F1.QTL.tmp
                      
                      ## dominance matrix 
                      matrix.D.tmp <- apply(matrix.A.tmp, 2, function(x){
                        ifelse(x == 0 | x == 4, 0, 1)
                      })
                      
                      ga <- matrix.A.tmp %*% QTL.effect$effect.a 
                      rm(matrix.A.tmp)
                      
                      gd <- matrix.D.tmp %*% QTL.effect$effect.d
                      rm(matrix.D.tmp)
                      
                      T2.gv.dom.1 <- ga + gd

                      F1.infor.tmp <- data.frame(rep = rep(r, n.progeny.vector[i]), 
                                                 plan = rep(i, n.progeny.vector[i]), 
                                                 index = F1.pop.tmp@id, 
                                                 mother = F1.pop.tmp@mother,
                                                 father = F1.pop.tmp@father,
                                                 T2.gv.dom.1)
                      F1.infor.tmp <- tibble::as_tibble(F1.infor.tmp)
                      return(F1.infor.tmp)
                      rm(F1.infor.tmp); rm(F1.QTL.tmp); rm(F1.pop.tmp); rm(F1.sites.tmp); rm(plan); gc()
                    }))

save(F1.infor, file = paste("02_target_trait_genetic_values_rep_", r, ".RData", sep = ""))
