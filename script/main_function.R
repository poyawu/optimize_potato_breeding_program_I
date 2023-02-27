####### determine the position for each QTL ##########
QTL.position.function <- function(nQTL, genMap, effect.a, effect.d){
  
  if(length(effect.a) != length(effect.d)){
    stop("Please check the additive and dominance effects. Their numbers are different.")
  }
  nQTL <- length(effect.a)
  ## calculate # of site per chromosome for the selected SNP
  num.QTL.chr <- sapply(genMap, length)
  
  ## determine how many sites are used as QTL for each chromosome
  select.num.QTL.chr <- round(num.QTL.chr/sum(num.QTL.chr)*nQTL, 0)
  nchr <- length(num.QTL.chr)
  while(sum(select.num.QTL.chr) != nQTL){
    if(sum(select.num.QTL.chr) > nQTL){
      select.num.QTL.chr[which.max(select.num.QTL.chr)] <- max(select.num.QTL.chr)-1
    }else if(sum(select.num.QTL.chr) < nQTL){
      select.num.QTL.chr[which.min(select.num.QTL.chr)] <- min(select.num.QTL.chr)+1
    }
  }
  ## randamly determine which positions are QTL based on genetic map 
  pos <- sapply(1:length(genMap), function(i){
    map.round <- round(genMap[[i]], 3)
    freq <- as.data.frame(table(round(genMap[[i]], 3)), stringsAsFactors = F)
    index <- sort(sample(1:nrow(freq), select.num.QTL.chr[i], replace = F))
    index.map <- freq$Var1[index]
    pos.tmp <- sapply(index.map, function(x) sample(which(map.round == as.numeric(x)), 1))
    names(pos.tmp) <- NULL
    return(pos.tmp)
  })
  
  genMap.dataframe <- lapply(1:nchr, function(x){
    data.frame(chr = rep(x, length(genMap[[x]])),
               site = 1:length(genMap[[x]]), 
               pos = genMap[[x]])
  })
  
  ## combine all infor across 12 chromosome into data.frame
  genMap.dataframe <- do.call(rbind, genMap.dataframe)
  genMap.dataframe$index <- 1:nrow(genMap.dataframe)
  
  QTL.summary <- lapply(1:nchr, function(x){
    sub <- genMap.dataframe[genMap.dataframe$chr == x, ]
    sub[match(pos[[x]], sub$site),]
  })
  
  QTL.summary <- do.call(rbind, QTL.summary)
  
  QTL.summary$effect.a <- effect.a
  QTL.summary$effect.d <- effect.d
  
  
  return(QTL.summary)
}

#####################################################################################################################
############### produce the information (# of clones, # of rep, # of loc...) for the specific scenario ##############
#####################################################################################################################
#### strategies (GS.stage): 
# "none": only PS
# "seedling": apply in seedling
# "SH": apply in SH
# "A": apply in A
# "seedling.SH.A": apply in seedling, SH, and A using the same EGV (same scenario as "seedling")
# "seedling.SH": apply in seedling, SH, and A using the same EGV (same scenario as "seedling")
# "SH.A": apply in SH, and A using the same GEBV (same scenario as "SH")

produce.scenario.fun <- function(T1.cost, T2.cost, geno.cost,
                                 GS.stage, 
                                 n.seedling.PS = 300000,
                                 l1, l2, l3, l4, l5, l6, 
                                 alpha.1 = alpha.1, alpha.2 = alpha.2, alpha.3 = alpha.3, 
                                 SL.to.SH = SL.to.SH, SH.to.A = SH.to.A, 
                                 A.to.B = A.to.B, B.to.C = B.to.C, C.to.D = C.to.D){
  
  ## set parameters for PS
  SL.to.SH.PS <- 1/3
  SH.to.A.PS <- 0.1
  A.to.B.PS <- 0.15
  B.to.C.PS <- 0.2
  C.to.D.PS <- 0.2
  
  ## number of clones in each stage (for PS) 
  n.SH.PS <- n.seedling.PS * SL.to.SH.PS
  n.A.PS <- n.SH.PS * SH.to.A.PS 
  n.B.PS <- n.A.PS * A.to.B.PS
  n.C.PS <- n.B.PS * B.to.C.PS
  n.D.PS <- n.C.PS * C.to.D.PS
  
  alpha <- switch(GS.stage, 
                  "GS-SL" = alpha.1,
                  "GS-SL:SH:A" = alpha.1,
                  "GS-SL:SH" = alpha.1,
                  "GS-SH" = alpha.2,
                  "GS-SH:A" = alpha.2,
                  "GS-A" = alpha.3)
  
  scenario.PS <- data.frame(stage = c("seedling", "SH", LETTERS[1:4]),
                            n = c(n.seedling.PS, n.SH.PS, n.A.PS, n.B.PS, n.C.PS, n.D.PS),
                            select.cost = c(T1.cost, T1.cost, T1.cost, T2.cost, T2.cost, T2.cost),
                            n.rep = rep(1, 6), 
                            n.loc = c(l1,l2,l3,l4,l5,l6),
                            stringsAsFactors = F)
  
  scenario.PS$cost <- scenario.PS$n * scenario.PS$select.cost * scenario.PS$n.rep * scenario.PS$n.loc 
  
  
  if(grepl("GS-SL", GS.stage)){
    scenario.GS <- rbind(scenario.PS[1, ], rep(NA, ncol(scenario.PS)), scenario.PS[2:6,])
    rownames(scenario.GS) <- NULL
    scenario.GS$stage[is.na(scenario.GS$stage)] <- "GS"
    n.geno <- scenario.PS$n[scenario.PS$stage == "seedling"]*alpha
    
    scenario.GS[scenario.GS$stage == "GS",2:5] <- c(n.geno, geno.cost, 1, 1)
    scenario.GS$cost[scenario.GS$stage == "GS"] <- scenario.GS$n[scenario.GS$stage == "GS"]*scenario.GS$select.cost[scenario.GS$stage == "GS"]
    alpha.2 <- scenario.GS$n[which(scenario.GS$stage == "GS")+1]/scenario.GS$n[scenario.GS$stage == "GS"]
    
    adj.ratio <- sum(scenario.PS$cost)/sum(scenario.GS$cost)
    
    scenario <- scenario.GS
    scenario$n[scenario$stage == "seedling"] <- scenario$n[scenario$stage == "seedling"] * adj.ratio 
    scenario$n[scenario$stage == "GS"] <-  scenario$n[scenario$stage == "seedling"] * alpha
    scenario$n[scenario$stage == "SH"] <- scenario$n[scenario$stage == "GS"] * alpha.2
    scenario$n[scenario$stage == "A"] <- scenario$n[scenario$stage == "SH"] * SH.to.A
    scenario$n[scenario$stage == "B"] <- scenario$n[scenario$stage == "A"] * A.to.B
    scenario$n[scenario$stage == "C"] <- scenario$n[scenario$stage == "B"] * B.to.C
    scenario$n[scenario$stage == "D"] <- scenario$n[scenario$stage == "C"] * C.to.D
    scenario$n <- round(scenario$n, 0)
    
    scenario$cost <- scenario$n*scenario$select.cost*scenario$n.rep*scenario$n.loc
  }
  
  if(grepl("GS-SH", GS.stage)){
    scenario.GS <- rbind(scenario.PS[1:2, ], rep(NA, ncol(scenario.PS)), scenario.PS[3:6,])
    n.geno <- scenario.PS$n[scenario.PS$stage == "SH"]*alpha
    rownames(scenario.GS) <- NULL
    
    scenario.GS$stage[is.na(scenario.GS$stage)] <- "GS"
    scenario.GS[scenario.GS$stage == "GS",2:5] <- c(n.geno, geno.cost, 1, 1)
    scenario.GS$cost[scenario.GS$stage == "GS"] <- scenario.GS$n[scenario.GS$stage == "GS"]*scenario.GS$select.cost[scenario.GS$stage == "GS"]
    alpha.2 <- scenario.GS$n[which(scenario.GS$stage == "GS")+1]/scenario.GS$n[scenario.GS$stage == "GS"]
    
    adj.ratio <- sum(scenario.PS$cost)/sum(scenario.GS$cost)
    
    scenario <- scenario.GS
    scenario$n[scenario$stage == "seedling"] <- scenario$n[scenario$stage == "seedling"] * adj.ratio 
    scenario$n[scenario$stage == "SH"] <- scenario$n[scenario$stage == "seedling"] * SL.to.SH
    scenario$n[scenario$stage == "GS"] <-  scenario$n[scenario$stage == "SH"] * alpha
    scenario$n[scenario$stage == "A"] <- scenario$n[scenario$stage == "GS"] * alpha.2
    scenario$n[scenario$stage == "B"] <- scenario$n[scenario$stage == "A"] * A.to.B
    scenario$n[scenario$stage == "C"] <- scenario$n[scenario$stage == "B"] * B.to.C
    scenario$n[scenario$stage == "D"] <- scenario$n[scenario$stage == "C"] * C.to.D
    scenario$n <- round(scenario$n, 0)
    
    scenario$cost <- scenario$n*scenario$select.cost*scenario$n.rep*scenario$n.loc
  }
  
  if(GS.stage == "GS-A"){
    scenario.GS <- rbind(scenario.PS[1:3, ], rep(NA, ncol(scenario.PS)), scenario.PS[4:6,])
    n.geno <- scenario.PS$n[scenario.PS$stage == "A"]*alpha
    rownames(scenario.GS) <- NULL
    
    scenario.GS$stage[is.na(scenario.GS$stage)] <- "GS"
    scenario.GS[scenario.GS$stage == "GS",2:5] <- c(n.geno, geno.cost, 1, 1)
    scenario.GS$cost[scenario.GS$stage == "GS"] <- scenario.GS$n[scenario.GS$stage == "GS"]*scenario.GS$select.cost[scenario.GS$stage == "GS"]
    alpha.2 <- scenario.GS$n[which(scenario.GS$stage == "GS")+1]/scenario.GS$n[scenario.GS$stage == "GS"]
    
    adj.ratio <- sum(scenario.PS$cost)/sum(scenario.GS$cost)
    
    scenario <- scenario.GS
    scenario$n[scenario$stage == "seedling"] <- scenario$n[scenario$stage == "seedling"] * adj.ratio 
    scenario$n[scenario$stage == "SH"] <- scenario$n[scenario$stage == "seedling"] * SL.to.SH
    scenario$n[scenario$stage == "A"] <-  scenario$n[scenario$stage == "SH"] * SH.to.A
    scenario$n[scenario$stage == "GS"] <- scenario$n[scenario$stage == "A"] * alpha
    scenario$n[scenario$stage == "B"] <- scenario$n[scenario$stage == "GS"] * alpha.2
    scenario$n[scenario$stage == "C"] <- scenario$n[scenario$stage == "B"] * B.to.C
    scenario$n[scenario$stage == "D"] <- scenario$n[scenario$stage == "C"] * C.to.D
    scenario$n <- round(scenario$n, 0) 
    
    scenario$cost <- scenario$n*scenario$select.cost*scenario$n.rep*scenario$n.loc
  }
  if(GS.stage!="PS"){
    scenario$stage[scenario$stage == "GS"] <- "genotyping"
    scenario.GS$stage[scenario.GS$stage == "GS"] <- "genotyping"
    return(list(PS = scenario.PS, 
                GS.origin = scenario.GS, 
                GS = scenario,
                alpha.2 = alpha.2))
  }else{
    return(list(PS = scenario.PS))
  }
}

##################################################################################
###### evaluate the performance for the specific scenario
evaluate.scenario.grid.function <- function(T1.cost = T1.cost, T2.cost = T2.cost, geno.cost = geno.cost, 
                                            F1.infor = F1.infor, GS.stage, 
                                            scenario, alpha.1, alpha.2, alpha.3,
                                            select = select,
                                            h2.T1 = h2.T1, ratio.var.T2 = c("vg"=1,"vge"=1,"ve"=0.5), 
                                            cor.T1.T2 = cor.T1.T2, 
                                            pa = pa, r, T2.type = T2.type, print = T){
  
  library(tibble)
  library(dplyr)
  
  if(GS.stage=="PS"){
    B <- sum(scenario$n * scenario$select.cost* scenario$n.loc)
  }else{
    alpha <- switch(GS.stage, 
                    "GS-SL" = alpha.1,
                    "GS-SL:SH:A" = alpha.1,
                    "GS-SL:SH" = alpha.1,
                    "GS-SH" = alpha.2,
                    "GS-SH:A" = alpha.2,
                    "GS-A" = alpha.3)
  }
  
  a1 <- as.character(round(alpha.1,2))
  a2 <- as.character(round(alpha.2,2))
  a3 <- as.character(round(alpha.3,2))
  
  scenario.name <- paste("T1cost_", T1.cost,"_T2cost_", T2.cost, "_genocost_", geno.cost, 
                         "_GS_", GS.stage, "_corT1T2_", cor.T1.T2, "_pa_", pa, 
                         "_a1_", a1, "_a2_", a2, "_a3_", a3, "_loc_", paste(scenario$n.loc, collapse = ""), sep = "")
  
  if(print){
    cat(scenario.name, "\n")
  }
  
  n.seedling <- scenario$n[scenario$stage == "seedling"]
  n.SH <- scenario$n[scenario$stage == "SH"]
  n.A <- scenario$n[scenario$stage == "A"]
  n.B <- scenario$n[scenario$stage == "B"]
  n.C <- scenario$n[scenario$stage == "C"]
  n.D <- scenario$n[scenario$stage == "D"]
  
  ######## variance for target trait 

  T2.gv <- F1.infor %>% pull(T2.type) # T2.gv <- F1.infor[, T2.type]
  vg.T2 <- var(T2.gv)
  
  vge.T2 <- vg.T2 * ratio.var.T2[names(ratio.var.T2) == "vge"]/ratio.var.T2[names(ratio.var.T2) == "vg"]
  ve.T2 <- vg.T2 * ratio.var.T2[names(ratio.var.T2) == "ve"]/ratio.var.T2[names(ratio.var.T2) == "vg"]
  
  ve.fake.seq <- ifelse(scenario$n.loc > 1, vge.T2/scenario$n.loc, 0) + ve.T2/(scenario$n.loc * scenario$n.rep)
  names(ve.fake.seq) <- as.character(scenario$stage)
  
  ######## obtain the genetic values, genetic variance and error variance of Ta
  ## GT1 = GT2 + e, e~N(0, var.T1.T2), where var.T1.T2 is determined by the degree of cor(GT1, GT2) 
  avg.T2 <- mean(T2.gv)
  #var.T1.T2 <- (1/(nrow(T2.gv)-2))*(1 - cor.T1.T2^2)/cor.T1.T2^2*sum((T2.gv - avg.T2)^2)
  var.T1.T2 <- (1/(length(T2.gv)-2))*(1 - cor.T1.T2^2)/cor.T1.T2^2*sum((T2.gv - avg.T2)^2)
  if(cor.T1.T2 > 0){
    T1.gv <- T2.gv + rnorm(nrow(F1.infor), mean = 0, sd = sqrt(var.T1.T2))
  }else{
    T1.gv <- T2.gv *(-1)+ rnorm(nrow(F1.infor), mean = 0, sd = sqrt(var.T1.T2))
  }
  
  ## variance for Ta 
  vg.T1 <- var(T1.gv)
  ve.T1 <- (1- h2.T1)/h2.T1*vg.T1
  
  F1.infor$T1.gv <- T1.gv
  
  rm(T1.gv); rm(T2.gv)
  
  ## select required progenies in this scenario (== n.seedling)
  F1.infor.tmp <- F1.infor
  F1.infor.tmp <- F1.infor.tmp[select, ]
  
  
  ### select GT2 and generate GT1 based on their correlations
  T2.gv <- F1.infor.tmp %>% pull(T2.type)#T2.gv <- F1.infor.tmp[,T2.type]
  
  if(nrow(F1.infor.tmp)!=n.seedling){
    stop("check T2.gv")
  }
  
  
  ###### 1. Seedling stage
  dat.tmp <- data.frame(F1.infor.tmp, T2.gv)
  rm(F1.infor.tmp)
  
  ## generate phenotypic values: T1
  dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))

  dat.F1 <- dat.tmp; rm(dat.tmp)
  dat.F1 <- dat.F1[order(dat.F1$T1.pv, decreasing = T), ] ## ordering by pheno.T1
  
  dat.F1$cross <- factor(dat.F1$cross)
  
  if(GS.stage == "PS"){
    ##### 2. F1 -> SH: based on pheno.T1 ####
    dat.tmp <- dat.F1[1:n.SH, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))

    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] ## ordering by T1 pV
    
    #### 3. SH -> A: based on T1
    dat.tmp <- dat.SH[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ] ## order by T1 pv
    
    ##### 4. A -> B: based on T1
    dat.tmp <- dat.A[1:n.B, ]
    dim(dat.tmp)
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))
    
    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] ## order by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]

    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
  }#end PS
  
  if(GS.stage == "GS-SL"){
    n.GS <- round(n.seedling * alpha.1, 0)
    
    ## 2.1 F1 ->  middle F1: based on T1 (for genotyped)
    dat.tmp <- dat.F1[1:n.GS, ]

    ## generate the GEBV
    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv))

    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by GEBV
    
    ## 2.2 middle F1 -> SH (based on GS)
    dat.tmp <- dat.tmp[1:n.SH, ]
    
    ## generate phenotype: T2
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] # order by T1 pv
    
    #### 3. SH -> A: based on T1
    dat.tmp <- dat.SH[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ]  ## order by T1 pv
    
    ##### 4. A -> B: based on T1
    dat.tmp <- dat.A[1:n.B, ]
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))
    
    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]
    dim(dat.tmp)
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]

    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
    
  }# end if seedling
  
  if(GS.stage == "GS-SH"){
    n.GS <- round(n.SH * alpha.2, 0)
    
    ##### 2. F1 -> SH: based on T1 ####
    dat.tmp <- dat.F1[1:n.SH, ]
 
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))

    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] ## ordering by T1 pV
    
    ## 3.1 SH -> middle SH (for genotyped)
    dat.tmp <- dat.SH[1:n.GS, ]

    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv))
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by gebv
    
    ## 3.2 middle SH -> A (based on GS)
    dat.tmp <- dat.tmp[1:n.A, ]

    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ] # order by T1 pv
    
    ##### 4. A -> B: based on T1
    dat.tmp <- dat.A[1:n.B, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))

    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))

    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] ## order by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]

    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
    
  }# end if SH
  
  if(GS.stage == "GS-A"){
    n.GS <- round(n.A*alpha.3, 0)
    
    ##### 2. F1 -> SH: based on T1 ####
    dat.tmp <- dat.F1[1:n.SH, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] ## ordering by T1 pV
    
    #### 3. SH -> A: based on T1
    dat.tmp <- dat.SH[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ] ## order by T1 pv
    
    ## 4.1 A -> middle A (for genotyped)
    dat.tmp <- dat.A[1:n.GS, ]
    #dim(dat.tmp)
    
    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv)); cor(dat.tmp$T2.ebv, dat.tmp$T2.gv)
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] ## order by gebv
    
    
    ## 4.2 middle A -> B (based on GS)
    dat.tmp <- dat.tmp[1:n.B, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))
    
    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]
    dim(dat.tmp)
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ]
        ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]

    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
    
  }# end if A
  
  if(GS.stage == "GS-SL:SH:A"){
    n.GS <- round(n.seedling*alpha.1, 0)
    
    ## 2.1 F1 ->  middle F1: based on pheno T1 
    dat.tmp <- dat.F1[1:n.GS, ] # n.GS = # of individual are gentyped
    #dim(dat.tmp)
    
    ## generate the GEBV
    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv))
    #cor(dat.tmp$T2.ebv, dat.tmp$T2.gv)
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by GEBV
    
    ## 2.2 middle F1 -> SH (based on GS, GEBV)
    dat.tmp <- dat.tmp[1:n.SH, ]
    
    ## generate phenotype: T1 
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1)); cor(dat.tmp$T1.pv, dat.tmp$T1.gv)^2
    
    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] # order by T1 pv
    
    
    #### 3.1 SH -> middle SH: based on pheno T1
    n.middle <- round(n.SH * alpha.2, 0)
    
    dat.tmp <- dat.SH[1:n.middle, ]; rm(n.middle)
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by gebv
    
    #### 3.2 middle SH -> A (based on GS)
    dat.tmp <- dat.tmp[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ]  ## order by T1 pv
    
    
    #### 4.1 A -> middle A (for genotyped)
    n.middle <- round(n.A * alpha.3, 0)
    
    dat.tmp <- dat.A[1:n.middle, ]; rm(n.middle)
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] ## order by gebv
    
    ## 4.2 middle A -> B (based on GS, GEBV)
    dat.tmp <- dat.tmp[1:n.B, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))
    
    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]
    dim(dat.tmp)
    
    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
    
  }# end if seedling.SH.A
  
  if(GS.stage == "GS-SL:SH"){
    n.GS <- round(n.seedling*alpha.1, 0)
    
    ## 2.1 F1 ->  middle F1: based on pheno T1 
    dat.tmp <- dat.F1[1:n.GS, ] # n.GS = # of individual are gentyped
    #dim(dat.tmp)
    
    ## generate the GEBV
    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv))
    #cor(dat.tmp$T2.ebv, dat.tmp$T2.gv)
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by GEBV
    
    ## 2.2 middle F1 -> SH (based on GS, GEBV)
    dat.tmp <- dat.tmp[1:n.SH, ]
    
    ## generate phenotype: T1
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] # order by T1 pv
    
    #### 3.1 SH -> middle SH: based on pheno T1
    n.middle <- round(n.SH * alpha.2, 0)
    
    dat.tmp <- dat.SH[1:n.middle, ]; rm(n.middle)
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by gebv
    
    #### 3.2 middle SH -> A (based on GS)
    dat.tmp <- dat.tmp[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ]  ## order by T1 pv
    
    ##### 4. A -> B: based on T1
    dat.tmp <- dat.A[1:n.B, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))

    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]

    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]
    
    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
  }# end if seedling.SH
  
  if(GS.stage == "GS-SH:A"){
    n.GS <- round(n.SH * alpha.2, 0)
    
    ##### 2. F1 -> SH: based on T1 ####
    dat.tmp <- dat.F1[1:n.SH, ]

    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))

    dat.SH <- dat.tmp; rm(dat.tmp)
    dat.SH <- dat.SH[order(dat.SH$T1.pv, decreasing = T), ] ## ordering by T1 pV
    
    ## 3.1 SH -> middle SH (for genotyped)
    dat.tmp <- dat.SH[1:n.GS, ]

    var.ebv.bv <- (1/(nrow(dat.tmp)-2))*(1 - pa^2)/pa^2*sum((dat.tmp$T2.gv - mean(dat.tmp$T2.gv))^2)
    dat.tmp$T2.ebv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(var.ebv.bv))
    
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] # order by gebv
    
    #### 3.2 middle SH -> A (based on GS)
    dat.tmp <- dat.tmp[1:n.A, ]
    
    ## generate the phenotypic values T1:
    dat.tmp$T1.pv <- dat.tmp$T1.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.T1))
    
    dat.A <- dat.tmp; rm(dat.tmp)
    dat.A <- dat.A[order(dat.A$T1.pv, decreasing = T), ]  ## order by T1 pv
    
    #### 4.1 A -> middle A (for genotyped)
    n.middle <- round(n.A * alpha.3, 0)
    
    dat.tmp <- dat.A[1:n.middle, ]; rm(n.middle)
    dat.tmp <- dat.tmp[order(dat.tmp$T2.ebv, decreasing = T),] ## order by gebv
    
    ## 4.2 middle A -> B (based on GS, GEBV)
    dat.tmp <- dat.tmp[1:n.B, ]
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "B"]))
    
    dat.B <- dat.tmp; rm(dat.tmp)
    dat.B <- dat.B[order(dat.B$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    #### 5. B -> C : based on T2
    dat.tmp <- dat.B[1:n.C, ]
    dim(dat.tmp)
    
    ## generate phenotypes: T2 and no longer update T1
    dat.tmp$T2.pv <- dat.tmp$T2.gv + rnorm(nrow(dat.tmp), 0, sqrt(ve.fake.seq[names(ve.fake.seq) == "C"]))
    
    dat.C <- dat.tmp; rm(dat.tmp)
    dat.C <- dat.C[order(dat.C$T2.pv, decreasing = T), ] # ordering by T2 pv
    
    ## 6. C -> D : based on T2
    dat.tmp <- dat.C[1:n.D, ]
    
    dat.D <- dat.tmp; rm(dat.tmp)
    dat.D <- dat.D[order(dat.D$T2.pv, decreasing = T), ] ## order by T2 pv
    
  }# end if SH.A
  
  
  ### calculate mean, sd, max, min...
  calculate.output.fun <-  function(gv.T2.seedling, gv.T2, stage){
    out <- data.frame(avg = mean(gv.T2),
                      avg.rm.seedling = mean(gv.T2) - mean(gv.T2.seedling), 
                      var = var(gv.T2),
                      median = median(gv.T2),
                      max = max(gv.T2), 
                      min = min(gv.T2))
    return(out)
  }
  
  infor <- scenario

  n.stage <- nrow(infor)
  
  select.method <- GS.stage
  
  out <- NULL 
  for(j in 1:n.stage){
    if(infor$stage[j] == "seedling"){
      stage <- "seedling"
      gv.T2 <- dat.F1$T2.gv
    }else{
      stage <- infor$stage[j] 
      gv.T2 <- get(paste("dat.", stage, sep = ""))[,"T2.gv"]
    }
    out <- rbind(out, calculate.output.fun(gv.T2.seedling = dat.F1$T2.gv, gv.T2 = gv.T2, stage = stage))
  }
  
  out$delta.G.per.stage <- c(0, diff(out$avg))
  
  dat.output <- data.frame(infor, out,
                           scenario = rep(scenario.name, n.stage),
                           T1cost = rep(T1.cost, n.stage),
                           T2cost = rep(T2.cost, n.stage),
                           genocost = rep(geno.cost, n.stage),
                           n.seedling = rep(n.seedling, n.stage),
                           select = rep(select.method, n.stage),
                           cor = rep(cor.T1.T2, n.stage),
                           pa = rep(pa, n.stage),
                           alpha.1 = rep(a1, n.stage),
                           alpha.2 = rep(a2, n.stage),
                           alpha.3 = rep(a3, n.stage),
                           rep = rep(r, n.stage), stringsAsFactors = FALSE)
  

  dat.output.D <- dat.output[dat.output$stage == "D",]
  
  return(list(dat.output = dat.output,
              dat.output.D = dat.output.D))
  
}