# Boxplot with top expressed genes from a Seurat object. 

top_expressed_genes = function(sobject, nPlot=20, title="Top expressed genes", assay = "RNA"){
  C = sobject@assays[[assay]]@counts
  C@x = C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[nPlot:1]
  boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "Fraction counts per cell", 
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main=title)
  invisible(rownames(C)[most_expressed])
}
