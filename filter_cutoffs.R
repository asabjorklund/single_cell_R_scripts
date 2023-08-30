# filter cells in seurat object
# cells grouped by batch.key
# qc stat in meta.name meta.data slot.
# define cutoff by median +/- mad.cut* mad 
# upper TRUE/FALSE to cut upper/lower tail.

filter_cutoffs = function(sobj, meta.name, batch.key = "orig.ident", 
                          upper = TRUE, mad.cut = 2, 
                          plot = T, with.points = FALSE){
  d = sobj@meta.data[,meta.name]
  names(d) = colnames(sobj)
  per.batch = split(d, sobj@meta.data[,batch.key])
  stats = lapply(per.batch, function(x){
    m  = median(x)
    s = mad(x)
    if (upper){
      cut = m + mad.cut*s
      cells = names(x)[x < cut]
    }else { 
      cut = m - mad.cut*s
      cells = names(x)[x > cut]
    }  
    return(list(cells = cells, cutoff=cut))
  })
  keep.cells = unlist(lapply(stats, function(x) x$cells)  )
  cutoffs = unlist(lapply(stats, function(x) x$cutoff)  )
  print(cutoffs)
  
  if (plot){
    meta = data.frame(batch = sobj@meta.data[,batch.key], data = d, cutoff = cutoffs[sobj@meta.data[,batch.key]])
    p1 <- ggplot(meta, aes(x=batch, y=data, fill = batch)) + 
      geom_violin(scale = "width", trim=TRUE, adjust = 1) + 
      geom_errorbar(width=0.8, aes(ymax=cutoff, ymin=cutoff)) +
      theme_classic() + 
      RotatedAxis() + NoLegend() + 
      ggtitle(meta.name)
    if (with.points){ 
      p1 = p1 + geom_jitter(height = 0, size = .1) + theme_classic()
    }
    print(p1)
  }
  
  return(keep.cells)
}
