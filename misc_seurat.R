run_seurat = function(tmp, nd=30, res=0.8, regress = NULL, scale = T, sct = F, reduction = "pca"){
    if (scale & !sct){
      tmp = ScaleData(tmp, vars.to.regress = regress, verbose = F)
    }
    if (sct){
      tmp = SCTransform(tmp, vars.to.regress = regress, verbose = F)
    }
    if (reduction == "pca"){
      tmp <- RunPCA(tmp, verbose = F)
    }
    tmp <- FindNeighbors(tmp, dims = 1:nd, verbose = F, reduction = reduction)
    tmp <- FindClusters(tmp, resolution = res, verbose = F)
    tmp <- RunUMAP(tmp, dims = 1:nd, verbose = F, reduction = reduction)
    return(tmp)
}


plot_seurat = function(x, reduction = "umap") {
p1 = DimPlot(x, group.by = "seurat_clusters", label = T, reduction = reduction) + NoAxes()
p2 = DimPlot(x, group.by = "Timepoint", reduction = reduction) + NoAxes()
p3 = DimPlot(x, group.by = "Tomato_sort", reduction = reduction) + NoAxes()
p4 = DimPlot(x, group.by = "Phase", reduction = reduction) + NoAxes()
p5 = FeaturePlot(x, features = "nFeature_RNA", reduction = reduction) + NoAxes()
p6 = FeaturePlot(x, features = "percent_mito", reduction = reduction) + NoAxes()
grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3)

VlnPlot(x, features = c("nFeature_RNA", "percent_mito", "S.Score","G2M.Score"), pt.size = 0, ncol = 2)
}



plot_markers = function(tmp, reduction = "umap") {
  tmp@active.assay = "RNA"
  plots = list()
  for (i in 1:length(marker.genes)){
    g = marker.genes[i]
    plots[[g]] = FeaturePlot(tmp, features = g, order = T, reduction = reduction) + NoAxes() + NoLegend() + ggtitle(g) + theme(plot.title = element_text(size = 10)) + cc
  }
  grid.arrange(grobs=plots, ncol=6)
}

run.deg = function(tmp, filename, topN = 10, method = "wilcox", latent.vars = NULL, 
                   force = FALSE, subsampleN = "min") {
  # run FindMarkers for all vs rest in ActiveIdent
  # subsample with:
  # None if subsampeN = Inf,
  # smallest group if subsampleN = "min"
  # subsampleN  if N = integer.
  # return full matrix and topN genes per group.
  DefaultAssay(tmp) = "RNA"
  if (file.exists(filename) & !force ) {
    markers = read.csv(filename)
  }else{
    if (is.integer(subsampleN)){
      tmp = tmp[,WhichCells(tmp,downsample = subsampleN)]
    }else if (subsampleN == "min"){
      n = min(table(tmp@active.ident))
      tmp = tmp[,WhichCells(tmp,downsample = n)]
    }
    
    markers = FindAllMarkers(tmp,assay = "RNA", only.pos = T, test.use = method, latent.vars = latent.vars)
    write.csv(markers, file=filename)
  }
  per.cl = split(markers$gene, markers$cluster)
  topM = lapply(per.cl, function(x) x[1:topN])
  topM = unique(unlist(topM))
  return(list(m=markers,top=topM))
}


run.deg.pairwise = function(tmp, file_prefix , id.col = "orig.ident", topN = 10, pval.cut = 0.01, 
                            method = "wilcox", latent.vars = NULL, force = FALSE) {
  # run pairwise DEG analysis for all possible pairs in a data column
  DefaultAssay(tmp) = "RNA"
  tmp = SetIdent(tmp, value = id.col) 
  idents = sort(unique(tmp@meta.data[,id.col]))
  cat("Running DEG analysis for: ", idents, "\n")
  nI = length(idents)
  output = list()
  for (i in 1:(nI-1)){
    id1 = idents[i]
    for (j in (i+1):nI){
      id2 = idents[j]
      filename = paste(c(file_prefix,id1,id2,"csv"), collapse = "." )
      if (file.exists(filename) & !force ) {
        mm = read.csv(filename, row.names = 1)
      }else{
        cat("DE", id1, "vs", id2, "\n")
        mm = FindMarkers(tmp, id1,id2, assay = "RNA",  test.use = method, latent.vars = latent.vars)
        write.csv(mm, file=filename)
      }
      output[[paste0(id1,"_vs_",id2)]] = list(
        m = mm,
        up = rownames(mm[mm$avg_log2FC>0,])[1:topN],
        down = rownames(mm[mm$avg_log2FC<0,])[1:topN],
        up.sign = rownames(mm[mm$avg_log2FC>0 & mm$p_val_adj < pval.cut,]),
        down.sign = rownames(mm[mm$avg_log2FC<0 & mm$p_val_adj < pval.cut,])
      )
    }
  }
  return(output)  
}







