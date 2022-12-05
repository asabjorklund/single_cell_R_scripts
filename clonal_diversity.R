## Modified from scRepertioir function:
#https://raw.githubusercontent.com/ncborcherding/scRepertoire/bc9b6840e78bed4951df715ca9abdaab825d4934/R/diversity.R

diversity_est = function(x, clone.col = "CLONE"){
  # run the different diversity estimates in the vegan package for one group.
  require(vegan)
  t = as.data.frame(table(x[,clone.col]))
  w <- vegan::diversity(t[,"Freq"], index = "shannon")
  x <- vegan::diversity(t[,"Freq"], index = "invsimpson")
  y <- vegan::estimateR(t[,"Freq"])[2] #Chao
  z <- vegan::estimateR(t[,"Freq"])[4] #ACE
  z2 <- vegan::diversity(t[,"Freq"], index = "shannon")/log(length(t[,"Freq"]))
  out <- c(w,x,y,z, z2)
  names(out) = c("Shannon","invSimpson","Chao","ACE","invPielou")
  return(out)
}

calculate_div = function(x, clone.col = "CLONE", n.sample = NULL, subsample.min = FALSE,
                         n.iter = 100, coldef = NULL){
  # run diversity_est on a list of samples.
  # subsample if necessary, specify size of sample with "subsample"
  # if subsample.min = TRUE - use the size of the smallest group.
  # use column "clone.col" to define the clonal info.
  if (is.null(n.sample) && subsample.min == FALSE){
    div = data.frame(Reduce(rbind,lapply(x,diversity_est)))
    div$group = names(x)
    m = reshape2::melt(div,id.vars = "group")
    p = ggplot(m, aes(x=group, y=value, fill = group)) + facet_wrap(~variable, scales = "free") + geom_bar(stat = "identity") + ggtitle("No subsampling") + theme_classic() + RotatedAxis()
    print(p)
  }else{
    cut = n.sample
    if (is.null(cut)){
      cut = min(sapply(x,nrow))
    }
    iterations = list()
    set.seed(1)
    for (i in 1:n.iter){
      x.sub = lapply(x, function(y) y[sample(1:nrow(y),min(cut,nrow(y))),])
      div = data.frame(Reduce(rbind,lapply(x.sub,diversity_est, clone.col = clone.col)))
      div$group = names(x)
      iterations[[i]] = div
    }
    div = Reduce(rbind,iterations)
    m = reshape2::melt(div,id.vars = "group")
    p = ggplot(m, aes(x=group, y=value, fill = group)) + facet_wrap(~variable, scales = "free") + geom_boxplot()  + ggtitle(sprintf("Subsampling %d cells, %d times",cut,n.iter))  + theme_classic() + RotatedAxis()
    if (!is.null(coldef)){ p = p + scale_fill_manual(values = coldef)}
    print(p)
  }
  invisible(div)
}