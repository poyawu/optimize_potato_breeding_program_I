rm(list = ls());gc()

library(AlphaSimR)
######## 01: prepare required datasets ########
##### I. create AlphaSimR format for 100 tetraploid potato clones 
## load haplotype information of 100 clones, genetic map information and centromere for 12 chromosome
load("100clones_haplotype.RData") ## haplotype.all
load("geneticmap_centromere.RData") ## geno.map, centromere

##　create MapPop format for 100 clones
founderPop <- newMapPop(genMap = geno.map, haplotypes = haplotype.all, inbred = F, ploidy = 4)
show(founderPop)
founderPop@centromere <- centromere

## create  "Simulation parameters" 
SP <- SimParam$new(founderPop)

## create "Pop" class for 100 clones
P.pop <- newPop(founderPop, simParam = SP)
show(P.pop)

save(P.pop, SP, file = "01_100_parental_clones_AlphasimR_format.RData")

##### II. simulate 2000 QTL additive and dominance effects 
shape.a <- 2
rate.a <- 5 
shape1.r <- 2
shape2.r <- 2

nQTL <- 2000

## produce additve effects 
effect.a <- rgamma(nQTL, shape = shape.a, rate = rate.a)

## produce ratio of d/a
ratio <- rbeta(nQTL, shape1 = shape1.r, shape2 = shape2.r)
effect.d <- effect.a*ratio

save(effect.a, effect.d, file = "01_2000QTL_A_D_effects.RData")

