rm(list = ls());gc()

setwd("~/HPC/Archive/PotatoTools/PhDPoYa/result/202206end_simulation_1000rep/")

library(AlphaSimR)
library(tibble)
library(dplyr)
#library(parallel)
#library(doParallel)

source("~/HPC/Archive/PotatoTools/PhDPoYa/script/202206end_simulation_1000rep/main_function_202206.R")

## load parental information
load("01_produce_parent_re_order.RData")

## load cross plan 
load("02_00_crossplan_re_order.RData")

## load QTL effects 
load("02_01_fixed_QTLeffects.RData")


n.seedling <- 300*3000

rep <- length(crossplan.list)
rep

rep.seq <- 1:10
rep.seq

r <- 1 

# n.core <- detectCores()
# cat(n.core, "\n")

# if(n.core > 4){
#   n.core <- round(n.core * 5 / 6)
# }else{
#   n.core <- 1
# }
# cat(n.core, "\n")
# n.core <- 10
# 
# cl <- makeCluster(n.core, type="FORK")
# clusterExport(cl = cl, 
#               varlist = c("crossplan.list", "dat.infor", "geno.map", 
#                           "effect.a", "effect.d", "n.seedling",
#                           "P.pop", "SP"))
# clusterEvalQ(cl = cl, {library(AlphaSimR)
#   library(tibble)
#   library(dplyr)})


out <- lapply(rep.seq, function(r){
  b <- proc.time()
  crossplan <- crossplan.list[[r]]
  cat("rep:",r, "\n")
  
  QTL.effect <- QTL.effect.produce.function(dat.infor = dat.infor, nQTL = 2000, 
                                            genMap = geno.map, 
                                            effect.a = effect.a, effect.d = effect.d, plot = F)
  
  n.progeny <- ceiling(n.seedling/nrow(crossplan))
  n.progeny.vector <- rep(n.progeny, nrow(crossplan))
  while(TRUE){
    if(sum(n.progeny.vector) > n.seedling){
      n.progeny.vector[which.max(n.progeny.vector)] <- max(n.progeny.vector) - 1
    }else if(sum(n.progeny.vector) == n.seedling)
      break
  }
  
  #a <- proc.time()
  F1.infor <- do.call(bind_rows, 
                      lapply(1:nrow(crossplan), function(i){
                        cat(i, "..")
                        plan <- matrix(crossplan[i,],1,2)
                        F1.pop.tmp <- makeCross(P.pop, crossPlan = plan, nProgeny = n.progeny.vector[i], simParam = SP)
                        F1.sites.tmp <- pullSegSiteGeno(F1.pop.tmp, simParam = SP)
                        #F1.sites.tmp <- as_tibble(F1.sites.tmp)
                        #save(F1.pop.tmp, file = paste("../../result/20211229_simulation/all_geno_infor/02_02_all_sites_infor_rep_",r,"_plan_",i,".RData", sep = ""))
                        
                        F1.QTL.tmp <- F1.sites.tmp[,QTL.effect$index]
                        
                        ## to speed up the process, the genetic effects are simulated directly. The QTL matrix is not kept. 
                        matrix.A.tmp <- F1.QTL.tmp
                        
                        matrix.D.tmp <- apply(matrix.A.tmp, 2, function(x){
                          ifelse(x == 0 | x == 4, 0, 1)
                        })
                        
                        ga <- matrix.A.tmp %*% QTL.effect$effect.a 
                        rm(matrix.A.tmp)
                        
                        gd <- matrix.D.tmp %*% QTL.effect$effect.d
                        rm(matrix.D.tmp)
                        
                        T2.gv.dom.1 <- ga + gd
                        T2.gv.dom.2 <- ga + gd*2
                        
                        F1.infor.tmp <- data.frame(rep = rep(r, n.progeny.vector[i]), 
                                                   plan = rep(i, n.progeny.vector[i]), 
                                                   index = F1.pop.tmp@id, 
                                                   mother = F1.pop.tmp@mother,
                                                   father = F1.pop.tmp@father,
                                                   T2.gv.dom.1, T2.gv.dom.2)
                        
                        F1.infor.tmp <- as_tibble(F1.infor.tmp)
                        #F1.infor <- bind_rows(F1.infor, F1.infor.tmp)
                        return(F1.infor.tmp)
                        rm(F1.infor.tmp); rm(F1.QTL.tmp); rm(F1.pop.tmp); rm(F1.sites.tmp); rm(plan); gc()
                      }))
  cat("\n")
  #print(proc.time()-a)
  
  save(F1.infor,
       file = paste("02_02_T2gv_rep_",r,".RData", sep = ""))
  print(proc.time()-b)
  gc()
})


# 
# for(r in rep.seq){
#   crossplan <- crossplan.list[[r]]
#   cat("rep:",r, "\n")
#   
#   QTL.effect <- QTL.effect.produce.function(dat.infor = dat.infor, nQTL = 2000, 
#                                             genMap = geno.map, 
#                                             effect.a = effect.a, effect.d = effect.d, plot = F)
#   
#   n.progeny <- ceiling(n.seedling/nrow(crossplan))
#   n.progeny.vector <- rep(n.progeny, nrow(crossplan))
#   while(TRUE){
#     if(sum(n.progeny.vector) > n.seedling){
#       n.progeny.vector[which.max(n.progeny.vector)] <- max(n.progeny.vector) - 1
#     }else if(sum(n.progeny.vector) == n.seedling)
#       break
#   }
#   
#   # #F1.sites <- NULL
#   # #F1.QTL.matrix.A <- NULL
#   # F1.infor <- NULL
#   # a <- proc.time()
#   # for(i in 1:nrow(crossplan)){
#   #   cat(i," ..")
#   #   plan <- matrix(crossplan[i,],1,2)
#   #   F1.pop.tmp <- makeCross(P.pop, crossPlan = plan, nProgeny = n.progeny.vector[i], simParam = SP)
#   #   F1.sites.tmp <- pullSegSiteGeno(F1.pop.tmp, simParam = SP)
#   #   #F1.sites.tmp <- as_tibble(F1.sites.tmp)
#   #   
#   #   #save(F1.pop.tmp, file = paste("../../result/20211229_simulation/all_geno_infor/02_02_all_sites_infor_rep_",r,"_plan_",i,".RData", sep = ""))
#   #   
#   #   F1.QTL.tmp <- F1.sites.tmp[,QTL.effect$index]
#   #   
#   #   ## to speed up the process, the genetic effects are simulated directly. The QTL matrix is not kept. 
#   #   
#   #   matrix.A.tmp <- F1.QTL.tmp
#   #   
#   #   matrix.D.tmp <- apply(matrix.A.tmp, 2, function(x){
#   #     ifelse(x == 0 | x == 4, 0, 1)
#   #   })
#   #   
#   #   ga <- matrix.A.tmp %*% QTL.effect$effect.a 
#   #   rm(matrix.A.tmp)
#   #   
#   #   gd <- matrix.D.tmp %*% QTL.effect$effect.d
#   #   rm(matrix.D.tmp)
#   #   
#   #   T2.gv.dom.1 <- ga + gd
#   #   T2.gv.dom.2 <- ga + gd*2
#   #   
#   #   # if(i == 1){
#   #   #   #F1.sites <- F1.sites.tmp
#   #   #   F1.QTL.matrix.A <- F1.QTL.tmp
#   #   #   F1.QTL.matrix.A <- as_tibble(F1.QTL.matrix.A) 
#   #   # }else{
#   #   #   #F1.sites <- bind_rows(F1.sites, F1.sites.tmp)
#   #   #   #F1.QTL.matrix.A <- rbind(F1.QTL.matrix.A, F1.QTL.tmp)
#   #   #   F1.QTL.tmp <- as_tibble(F1.QTL.tmp)
#   #   #   F1.QTL.matrix.A <- bind_rows(F1.QTL.matrix.A, F1.QTL.tmp)
#   #   # }
#   #   
#   #   F1.infor.tmp <- data.frame(rep = rep(r, n.progeny.vector[i]), 
#   #                              plan = rep(i, n.progeny.vector[i]), 
#   #                              index = F1.pop.tmp@id, 
#   #                              mother = F1.pop.tmp@mother,
#   #                              father = F1.pop.tmp@father,
#   #                              T2.gv.dom.1, T2.gv.dom.2)
#   #   
#   #   F1.infor.tmp <- as_tibble(F1.infor.tmp)
#   #   F1.infor <- bind_rows(F1.infor, F1.infor.tmp)
#   #   rm(F1.infor.tmp); rm(F1.QTL.tmp); rm(F1.pop.tmp); rm(F1.sites.tmp); rm(plan); gc()
#   # }
#   # cat("\n")
#   # print(proc.time()-a)
#   # 
#   a <- proc.time()
#   F1.infor <- do.call(bind_rows, 
#                       lapply(1:nrow(crossplan), function(i){
#                         cat(i, "..")
#                         plan <- matrix(crossplan[i,],1,2)
#                         F1.pop.tmp <- makeCross(P.pop, crossPlan = plan, nProgeny = n.progeny.vector[i], simParam = SP)
#                         F1.sites.tmp <- pullSegSiteGeno(F1.pop.tmp, simParam = SP)
#                         #F1.sites.tmp <- as_tibble(F1.sites.tmp)
#                         #save(F1.pop.tmp, file = paste("../../result/20211229_simulation/all_geno_infor/02_02_all_sites_infor_rep_",r,"_plan_",i,".RData", sep = ""))
#                         
#                         F1.QTL.tmp <- F1.sites.tmp[,QTL.effect$index]
#                         
#                         ## to speed up the process, the genetic effects are simulated directly. The QTL matrix is not kept. 
#                         matrix.A.tmp <- F1.QTL.tmp
#                         
#                         matrix.D.tmp <- apply(matrix.A.tmp, 2, function(x){
#                           ifelse(x == 0 | x == 4, 0, 1)
#                         })
#                         
#                         ga <- matrix.A.tmp %*% QTL.effect$effect.a 
#                         rm(matrix.A.tmp)
#                         
#                         gd <- matrix.D.tmp %*% QTL.effect$effect.d
#                         rm(matrix.D.tmp)
#                         
#                         T2.gv.dom.1 <- ga + gd
#                         T2.gv.dom.2 <- ga + gd*2
#                         
#                         F1.infor.tmp <- data.frame(rep = rep(r, n.progeny.vector[i]), 
#                                                    plan = rep(i, n.progeny.vector[i]), 
#                                                    index = F1.pop.tmp@id, 
#                                                    mother = F1.pop.tmp@mother,
#                                                    father = F1.pop.tmp@father,
#                                                    T2.gv.dom.1, T2.gv.dom.2)
#                         
#                         F1.infor.tmp <- as_tibble(F1.infor.tmp)
#                         #F1.infor <- bind_rows(F1.infor, F1.infor.tmp)
#                         return(F1.infor.tmp)
#                         rm(F1.infor.tmp); rm(F1.QTL.tmp); rm(F1.pop.tmp); rm(F1.sites.tmp); rm(plan); gc()
#                       }))
#   cat("\n")
#   print(proc.time()-a)
#   
#   save(F1.infor,
#        file = paste("02_02_T2gv_rep_",r,".RData", sep = ""))
#   
#   # rm(F1.QTL.matrix.A); rm(F1.QTL.matrix.D);
#   gc()
# }
