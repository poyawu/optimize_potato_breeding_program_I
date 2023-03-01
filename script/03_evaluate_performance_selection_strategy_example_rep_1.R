rm(list =ls()); gc()

## load required functions 
source("../script/main_function.R")

######################################################################################################
######### Parameter settings in the simulation study: considering a specific scenario ################

T1.cost <- 1.4 ## phenotyping cost for Ta: auxiliary trait
T2.cost <- 25 ## phenotyping cost for Tt: target trait
geno.cost <- 25 ## genotyping cost 

GS.stage <- "GS-SL" ## selection strategy (could be PS, GS-SL, GS-SH, GS-A, GS-SL:SH, GS-SH:A, and GS-SL:SH:A)

## diﬀerent weights of genomic selection (GS) relative to phenotypic selection
alpha.1 <- 0.9  ## if GS applies at seedling stage
alpha.2 <- NA  ## if GS applies at single hills
alpha.3 <- NA  ## if GS applies at A clone

cor.T1.T2 <- -0.15 ## correlation between two traits 
pa <- 0.7 ## prediction accuracy of GS model for the target trait

h2.T1 <- 0.6 ## heritability for auxiliary trait
ratio.var.T2 <- c("vg"=1,"vge"=1,"ve"=0.5) ## variance components for the target trait

l1 <- 1 ## number of locations at seedling stage 
l2 <- 1 ## number of locations at single hills 
l3 <- 1 ## number of locations at A clone
l4 <- 2 ## number of locations at B clone
l5 <- 3 ## number of locations at A clone
l6 <- 4 ## number of locations at single hills

SL.to.SH <- 1/3 ## p1: selected proportion from seedlings to single hills 
SH.to.A <- 0.1 ## p2: selected proportion from single hills to A clone 
A.to.B <- 0.15 ## p3: selected proportion from A clone to B clone 
B.to.C <- 0.2  ## p4: selected proportion from B clone to C clone 
C.to.D <- 0.2 ## p5: selected proportion from C clone to D clone 


n.seedling.PS <- 300000  ## number of progenies at seedling stage in a standard PS 
######################################################################################################
######################################################################################################


## obtain the scenario information under a fixed budget
scenario.list <- produce.scenario.fun(T1.cost = T1.cost, T2.cost = T2.cost, 
                                      geno.cost = geno.cost, GS.stage = GS.stage,
                                      l1 = l1, l2 = l2, l3 = l3, l4 = l4, l5 = l5, l6 = l6,
                                      alpha.1 = alpha.1, alpha.2 = alpha.2, alpha.3 = alpha.3,
                                      n.seedling.PS = n.seedling.PS,
                                      SL.to.SH = SL.to.SH, SH.to.A = SH.to.A, 
                                      A.to.B = A.to.B, B.to.C = B.to.C, C.to.D = C.to.D)

## obtain the exact information at each stage under a fixed budget by adjusting the number of clones at each stage
if(GS.stage == "PS"){
  scenario <- scenario.list$PS 
}else{
  scenario <- scenario.list$GS 
}

## obtain the exact nummber of clones at each stage without the information of genotyping
scenario.input <- scenario[scenario$stage!="genotyping",] 

n.seedling <- scenario$n[scenario$stage == "seedling"]

## at repetition = 1, load the genetic values of initial population
r <- 1
load(paste("02_target_trait_genetic_values_rep_", r, ".RData", sep = ""))


## fix the size initial population as 300,000 under the basis of standard PS scheme 
F1.infor$cross <- paste(F1.infor$mother, F1.infor$father, sep =":")

unique.cross <- unique(F1.infor$cross)
index <- NULL
for(i in 1:length(unique.cross)){
  index <- c(index, sort(sample(which(F1.infor$cross == unique.cross[i]), 1000)))
}
F1.infor <- F1.infor[index, ]
F1.infor$pos <- 1:nrow(F1.infor)

## trait name: 
T2.type <- "T2.gv.dom.1"

## Randomly sample the reduced "n.seedling" from the initial population with an equal sample size for each cross population
n.progeny <- ceiling(n.seedling/length(unique.cross))
n.progeny.vector <- rep(n.progeny, length(unique.cross))
while(TRUE){
  if(sum(n.progeny.vector) > n.seedling){
    n.progeny.vector[which.max(n.progeny.vector)] <- max(n.progeny.vector) - 1
  }else if(sum(n.progeny.vector) == n.seedling)
    break
}
uni.cross <- unique(F1.infor$cross)
select <- NULL
for(i in 1:length(uni.cross)){
  select <- c(select, sort(sample(F1.infor$pos[F1.infor$cross == uni.cross[i]], n.progeny.vector[i])))
}

### evalulate the performance for a specific scenario
output.list <- evaluate.scenario.grid.function(T1.cost = T1.cost, T2.cost = T2.cost,
                                               geno.cost = geno.cost,
                                               F1.infor = F1.infor, GS.stage = GS.stage,
                                               scenario = scenario.input, 
                                               alpha.1 = alpha.1, alpha.2 = alpha.2, alpha.3 = alpha.3,
                                               select = select, 
                                               h2.T1 = h2.T1, ratio.var.T2 = ratio.var.T2,
                                               cor.T1.T2 = cor.T1.T2,
                                               pa = pa, r = r, T2.type = T2.type, print = F)

### Result: the performance for all stages
output.all.stage <- output.list$dat.output
output.all.stage
#@ regarding the output (data.frame)
# avg/var/median/median/max/min: the mean/variance/median/maximum/minimum of genetic values for the target trait at such stage 
# avg.rm.seedling: the mean genetic values for the target trait at such stage minus the mean genetic values for the target trait at seedling stage 
# delta.G.per.stage: the current mean genetic values for the target trait minus the previous mean genetic values for the target trait 

### Result: the performance at D clone stage only
output.D.stage <- output.list$dat.output.D
output.D.stage
