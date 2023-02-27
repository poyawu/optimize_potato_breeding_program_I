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