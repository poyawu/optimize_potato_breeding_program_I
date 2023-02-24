rm(list =ls())

#### VC 1
setwd("~/HPC/Archive/PotatoTools/PhDPoYa/result/202206end_simulation_1000rep/fixed_ratio/")
library(parallel)
require(doParallel)

source("../../../script/202206end_simulation_1000rep/main_function_202206.R")
load("../02_00_crossplan_re_order.RData")


GS.stage.seq <- c("GS:seedling", "GS:SH", "GS:A", "GS:seedling:SH", "GS:SH:A","GS:seedling:SH:A", "PS") # s



T1.cost <- 1.4
T2.cost <- 25

geno.cost <- 25

cor.seq <- c(0.3, 0.15, -0.15) # c

h2.T1 <- 0.6
ratio.var.T2 = c("vg"=1,"vge"=1,"ve"=0.5)



SL.to.SH <- 1/3
SH.to.A <- 0.1
A.to.B <- 0.15
B.to.C <- 0.2
C.to.D <- 0.2

n.seedling.PS <- 300000
pop <- "multiple"
T2.type <- "T2.gv.dom.1"


alpha.list <- list()
for(i in 1:length(GS.stage.seq)){
  GS <- GS.stage.seq[i]
  if(GS == "GS:seedling"){
    alpha.com <- tidyr::crossing(alpha.1 = seq(0.9, 0.4, -0.1), alpha.2 = NA, alpha.3 = NA)
  }
  if(GS == "GS:SH" ){
    alpha.com <- tidyr::crossing(alpha.1 = NA, alpha.2 = seq(0.9, 0.2, -0.1), alpha.3 = NA)
  }
  if(GS == "GS:A"){
    alpha.com <- tidyr::crossing(alpha.1 = NA, alpha.2 = NA, alpha.3 = seq(0.9, 0.2, -0.1))
  }
  if(GS == "GS:seedling:SH"){
    alpha.com <- tidyr::crossing(alpha.1 = seq(0.9, 0.4, -0.1), alpha.2 = seq(0.9, 0.2, -0.1), alpha.3 = NA)
  }
  if(GS == "GS:SH:A"){
    alpha.com <- tidyr::crossing(alpha.1 = NA, alpha.2 = seq(0.9, 0.2, -0.1), alpha.3 = seq(0.9, 0.2, -0.1))
  }
  if(GS == "GS:seedling:SH:A"){
    alpha.com <- tidyr::crossing(alpha.1 = seq(0.9, 0.4, -0.1), alpha.2 = seq(0.9, 0.2, -0.1), alpha.3 = seq(0.9, 0.2, -0.1))
  }
  if(GS == "PS"){
    alpha.com <- tidyr::crossing(alpha.1 = NA, alpha.2 = NA, alpha.3 = NA)
  }
  alpha.com <- as.data.frame(alpha.com)
  alpha.list[[GS]] <- alpha.com 
  
}
rm(i); rm(GS)


#rep.seq  <- 1:100



b <- proc.time()

cl <- makeCluster(20, type="FORK")

clusterExport(cl = cl, varlist = c("evaluate.scenario.grid.function", "produce.scenario.fun",
                                   "T2.type","ratio.var.T2",
                                   "cor.seq","h2.T1", "crossplan.list", 
                                   "alpha.list", "pop", "GS.stage.seq" ))



output.tmp <- parallel::parLapply(cl = cl, 1:1000, function(r){
  if(pop != "single"){
    crossplan <- crossplan.list[[r]]
  }
  
  load(paste("../", "02_02_T2gv_rep_",r,".RData", sep = ""))
  
  
  dat.summary <- NULL
  F1.infor$cross <- paste(F1.infor$mother, F1.infor$father, sep =":")
  
  unique.cross <- unique(F1.infor$cross)
  index <- NULL
  for(i in 1:length(unique.cross)){
    index <- c(index, sort(sample(which(F1.infor$cross == unique.cross[i]), 1000)))
  }
  F1.infor <- F1.infor[index, ]
  rm(unique.cross); rm(i); rm(index)
  
  F1.infor$pos <- 1:nrow(F1.infor)
  
  ##############################################################
  for(s in 1:length(GS.stage.seq)){
    GS.stage <- GS.stage.seq[s]
    alpha.com <- alpha.list[[GS.stage]]
    
    if(GS.stage == "PS"){
      pa.seq <- NA
    }else{
      pa.seq <- c(0.3, 0.5, 0.7)
    }
    
    for(a in 1:nrow(alpha.com)){
      alpha.tmp <- alpha.com[a, ]
      alpha.1 <- alpha.tmp[,"alpha.1"]
      alpha.2 <- alpha.tmp[,"alpha.2"]
      alpha.3 <- alpha.tmp[,"alpha.3"]
      
      scenario.list <- produce.scenario.fun(T1.cost = T1.cost, T2.cost = T2.cost, 
                                            geno.cost = geno.cost, GS.stage = GS.stage,
                                            alpha.1 = alpha.1, alpha.2 = alpha.2, alpha.3 = alpha.3,
                                            n.seedling.PS = n.seedling.PS,
                                            SL.to.SH = SL.to.SH, SH.to.A = SH.to.A, 
                                            A.to.B = A.to.B, B.to.C = B.to.C, C.to.D = C.to.D)
      if(GS.stage == "PS"){
        scenario <- scenario.list$PS 
      }else{
        scenario <- scenario.list$GS
      }
      
      scenario.input <- scenario[scenario$stage!="GS",]
      n.seedling <- scenario$n[scenario$stage == "seedling"]
      n.SH <- scenario$n[scenario$stage == "SH"]
      n.A <- scenario$n[scenario$stage == "A"]
      n.B <- scenario$n[scenario$stage == "B"]
      n.C <- scenario$n[scenario$stage == "C"]
      n.D <- scenario$n[scenario$stage == "D"]
      
      ## determine how many progenies should be selected in each combination
      if(pop == "multiple"){
        n.progeny <- ceiling(n.seedling/nrow(crossplan))
        n.progeny.vector <- rep(n.progeny, nrow(crossplan))
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
      }else{ ## pop = single, that is only one crossing 
        n.progeny <- n.seedling
        select <- sort(sample(1:nrow(F1.infor), n.progeny))
      }
      
      dat.tmp <- do.call(rbind, lapply(cor.seq, function(cor.T1.T2){
        tmp <- NULL
        for(p in 1:length(pa.seq)){
          pa <- pa.seq[p]
          tmp.list <- evaluate.scenario.grid.function(T1.cost = T1.cost, T2.cost = T2.cost,
                                                      geno.cost = geno.cost,
                                                      F1.infor = F1.infor, GS = GS.stage,
                                                      scenario = scenario.input, alpha.1 = alpha.1,
                                                      alpha.2 = alpha.2, alpha.3 = alpha.3,
                                                      select = select,
                                                      h2.T1 = h2.T1, ratio.var.T2 = ratio.var.T2,
                                                      cor.T1.T2 = cor.T1.T2,
                                                      pa = pa, r = r,
                                                      pop = pop, T2.type = T2.type, print = F)
          tmp <- rbind(tmp, tmp.list$dat.output)
        }
        return(tmp)
      }))
      
      dat.summary <- rbind(dat.summary, dat.tmp)
    } # end a (alpha)
  }# end s (GS)
  save(dat.summary, file = paste("03_01_T1cost_",T1.cost, "_genocost_",geno.cost,"summary_dom1_VC1_rep_", r, ".RData", sep = ""))

})

print(proc.time()-b)
stopCluster(cl)

